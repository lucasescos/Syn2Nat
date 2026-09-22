#!/usr/bin/env python3
"""
Scaled ESMC Representation Extraction Pipeline.
Extracts both:
1. Sequence-level representations ([2560] float32 vectors in Parquet)
2. Per-residue representations ([L, 2560] float16 matrices in compressed NPZ)
with multi-threaded batching and fault-tolerant checkpointing.
"""

import argparse
import concurrent.futures
import json
import os
import sys
import time
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import socket
import threading
import numpy as np
import pandas as pd
import torch

# Prevent any socket from hanging indefinitely on network drops
socket.setdefaulttimeout(35.0)

ROOT_DIR = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT_DIR))

from esm.sdk import esmc_client
from esm.sdk.api import ESMProtein, LogitsConfig
from bindpred.esmc_client import get_api_key, get_api_keys


def parse_args():
    parser = argparse.ArgumentParser(description="Extract ESMC sequence and per-residue representations.")
    parser.add_argument(
        "--dataset",
        choices=["microproteins", "human", "hmpa", "test"],
        default="microproteins",
        help="Dataset to process (microproteins: 7,264; human: 20,652; hmpa: 617,463; test: 5 samples)",
    )
    parser.add_argument(
        "--model",
        default="esmc-6b-2024-12",
        help="ESMC model variant (default: esmc-6b-2024-12)",
    )
    parser.add_argument(
        "--layer",
        type=int,
        default=79,
        help="Layer to extract (79 = penultimate layer for 6B)",
    )
    parser.add_argument(
        "--batch-size",
        type=int,
        default=200,
        help="Number of proteins per batch file (default: 200)",
    )
    parser.add_argument(
        "--workers",
        type=int,
        default=6,
        help="Concurrent API workers (default: 6)",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=ROOT_DIR / "bindpred" / "representations",
        help="Directory to store extracted representations",
    )
    parser.add_argument(
        "--limit",
        type=int,
        default=None,
        help="Optional limit on number of sequences to process",
    )
    return parser.parse_args()


class RepresentationExtractor:
    def __init__(self, model_name: str = "esmc-6b-2024-12", layer: int = 79, workers: int = 6):
        self.model_name = model_name
        self.layer = layer
        self.workers = workers
        self.tokens = get_api_keys()
        if not self.tokens:
            raise ValueError("No Biohub API keys found in environment or .env files.")
        self.exhausted_tokens = set()
        self._thread_counter = 0
        self._lock = threading.Lock()
        self._local = threading.local()
        print(f"[*] Initialized RepresentationExtractor with {len(self.tokens)} configured API key(s).")

    def _get_active_tokens(self) -> List[str]:
        return [t for t in self.tokens if t not in self.exhausted_tokens]

    def _mark_token_exhausted(self, token: str):
        with self._lock:
            if token not in self.exhausted_tokens:
                self.exhausted_tokens.add(token)
                active = self._get_active_tokens()
                masked = f"...{token[-6:]}" if len(token) >= 6 else "***"
                print(
                    f"\n[KEY EXHAUSTED] Biohub API Key ({masked}) hit daily credit limit. Blacklisted from pool. ({len(active)} active keys remaining)",
                    file=sys.stderr,
                )
                if not active:
                    print("\n[CRITICAL] All configured Biohub API keys have exceeded daily credit quotas!", file=sys.stderr)

    def _get_client_and_token(self) -> Tuple[Any, str]:
        with self._lock:
            active = self._get_active_tokens()
            if not active:
                raise RuntimeError("All Biohub API keys exhausted daily credit limit.")
            cur_token = getattr(self._local, "client_token", None)
            if cur_token not in active:
                self._thread_counter += 1
                assigned_token = active[self._thread_counter % len(active)]
                self._local.client = esmc_client(
                    model=self.model_name,
                    url="https://biohub.ai",
                    token=assigned_token,
                    request_timeout=25.0,
                )
                self._local.client_token = assigned_token

        return self._local.client, self._local.client_token

    def _rotate_thread_client(self, failed_token: Optional[str] = None):
        with self._lock:
            active = self._get_active_tokens()
            if not active:
                raise RuntimeError("All Biohub API keys exhausted daily credit limit.")
            other_active = [t for t in active if t != failed_token]
            next_token = other_active[0] if other_active else active[0]

        self._local.client = esmc_client(
            model=self.model_name,
            url="https://biohub.ai",
            token=next_token,
            request_timeout=25.0,
        )
        self._local.client_token = next_token

    def extract_single(self, target: Dict[str, Any], max_retries: int = 4) -> Optional[Dict[str, Any]]:
        """Extract both per-residue and sequence-level vectors for a single target."""
        target_id = target["id"]
        full_seq = "".join(target["sequence"].split()).upper()

        if not full_seq:
            return None

        orig_len = len(full_seq)
        seq = full_seq[:2048]  # ESMC context window cap

        for attempt in range(1, max_retries + 1):
            token = None
            try:
                client, token = self._get_client_and_token()
                protein = ESMProtein(sequence=seq)
                tensor_input = client.encode(protein)
                config = LogitsConfig(
                    sequence=True,
                    return_hidden_states=True,
                    ith_hidden_layer=self.layer,
                )
                output = client.logits(tensor_input, config)

                if hasattr(output, "error_code"):
                    raise RuntimeError(f"ESMProteinError {output.error_code}: {output.error_msg}")

                # Raw shape: [1, 1, L + 2, 2560]
                raw_hs = output.hidden_states.float().squeeze(0).squeeze(0)  # [L + 2, 2560]
                # Slice off <BOS> (index 0) and <EOS> (index -1)
                residue_tensor = raw_hs[1:-1]  # [L, 2560]

                # Mean pool along sequence length L to get sequence-level vector
                seq_vector = residue_tensor.mean(dim=0).detach().cpu().numpy()  # [2560] float32

                # Convert per-residue matrix to float16 to save 50% disk space
                residue_matrix = residue_tensor.detach().cpu().numpy().astype(np.float16)  # [L, 2560] float16

                return {
                    "id": target_id,
                    "gene": target.get("gene", ""),
                    "length": len(seq),
                    "orig_length": orig_len,
                    "seq_vector": seq_vector,
                    "residue_matrix": residue_matrix,
                }
            except Exception as e:
                err_str = str(e)
                # Reset cached client on error to ensure a clean new connection on retry
                if hasattr(self._local, "client"):
                    del self._local.client
                if hasattr(self._local, "client_token"):
                    del self._local.client_token

                # Check if key exceeded daily credit limit
                if any(x in err_str.lower() for x in ["daily credit limit", "credit limit", "exceeded your daily credit"]):
                    if token:
                        self._mark_token_exhausted(token)
                    if not self._get_active_tokens():
                        raise RuntimeError("All Biohub API keys have exhausted their daily credit quota.")
                    time.sleep(1.0)
                    continue

                # Check if transient rate limited or general 429/quota error, rotate token
                if "429" in err_str or "quota" in err_str.lower() or "rate" in err_str.lower() or "403" in err_str:
                    self._rotate_thread_client(failed_token=token)
                    time.sleep(1.0 + (attempt * 0.5))
                    continue

                if attempt == max_retries:
                    print(f"\n[ERROR] Failed target {target_id} after {max_retries} attempts: {e}", file=sys.stderr)
                    return None
                time.sleep(attempt * 2)

        return None




def load_dataset(dataset_type: str, limit: Optional[int] = None) -> List[Dict[str, Any]]:
    targets = []
    if dataset_type in ("microproteins", "test"):
        path = ROOT_DIR / "data" / "parsed" / "orfs.parquet"
        df = pd.read_parquet(path)
        for _, row in df.iterrows():
            targets.append({
                "id": str(row["orf_id"]),
                "gene": str(row.get("gene", "")),
                "sequence": str(row["sequence"]),
            })
    elif dataset_type == "human":
        path = ROOT_DIR / "data" / "targets" / "human_canonical_proteome.parquet"
        df = pd.read_parquet(path)
        for _, row in df.iterrows():
            targets.append({
                "id": str(row["uniprot_id"]),
                "gene": str(row.get("gene_symbol", "")),
                "sequence": str(row["sequence"]),
            })
    elif dataset_type == "hmpa":
        path = ROOT_DIR / "data" / "curated" / "hmpa_617k_microproteins.parquet"
        df = pd.read_parquet(path)
        for _, row in df.iterrows():
            targets.append({
                "id": str(row["id"]),
                "gene": str(row.get("gene", "")),
                "sequence": str(row["sequence"]),
            })

    if dataset_type == "test":
        targets = targets[:5]

    return targets


def get_completed_ids(dataset_dir: Path) -> set:
    """Read all completed IDs from existing batch metadata files."""
    completed = set()
    batches_dir = dataset_dir / "batches"
    if not batches_dir.exists():
        return completed

    for parquet_file in batches_dir.glob("batch_*.parquet"):
        try:
            df = pd.read_parquet(parquet_file, columns=["id"])
            completed.update(df["id"].tolist())
        except Exception:
            pass
    return completed


def main():
    args = parse_args()
    print("=" * 80)
    print("Scaled ESMC Representation Extraction Pipeline")
    print(f"Model:      {args.model} (Layer {args.layer})")
    print(f"Dataset:    {args.dataset}")
    print(f"Workers:    {args.workers} concurrent threads")
    print(f"Batch Size: {args.batch_size}")
    print("=" * 80)

    if args.dataset in ("microproteins", "test"):
        dataset_name = "microproteins"
    elif args.dataset == "hmpa":
        dataset_name = "hmpa"
    else:
        dataset_name = "human"
    dataset_dir = args.output_dir / dataset_name
    batches_dir = dataset_dir / "batches"
    batches_dir.mkdir(parents=True, exist_ok=True)

    # 1. Load targets
    all_targets = load_dataset(args.dataset, limit=args.limit)
    print(f"Loaded {len(all_targets):,} total targets for {args.dataset}.")

    # 2. Check completed targets
    completed_ids = get_completed_ids(dataset_dir)
    print(f"Already completed targets found: {len(completed_ids):,}")

    remaining_targets = [t for t in all_targets if t["id"] not in completed_ids]
    if args.limit is not None:
        remaining_targets = remaining_targets[:args.limit]
    print(f"Remaining targets to extract:   {len(remaining_targets):,}")

    if not remaining_targets:
        print("\nAll targets already extracted! Nothing to do.")
        return

    extractor = RepresentationExtractor(
        model_name=args.model, layer=args.layer, workers=args.workers
    )

    # Determine next batch index
    existing_batches = list(batches_dir.glob("batch_*.parquet"))
    next_batch_idx = len(existing_batches)

    total_to_process = len(remaining_targets)
    processed_count = 0
    start_time = time.time()

    # Process in batches
    for b_start in range(0, total_to_process, args.batch_size):
        b_end = min(b_start + args.batch_size, total_to_process)
        batch_targets = remaining_targets[b_start:b_end]
        batch_idx = next_batch_idx + (b_start // args.batch_size)

        print(f"\n>>> Starting Batch {batch_idx:04d} ({len(batch_targets)} targets, indices {b_start}..{b_end - 1})")
        b_t0 = time.time()

        batch_results = []
        all_exhausted = False
        with concurrent.futures.ThreadPoolExecutor(max_workers=args.workers) as executor:
            futures = {
                executor.submit(extractor.extract_single, t): t["id"]
                for t in batch_targets
            }
            for future in concurrent.futures.as_completed(futures):
                try:
                    res = future.result()
                    if res is not None:
                        batch_results.append(res)
                except RuntimeError as re:
                    if "exhausted" in str(re).lower():
                        print(f"\n[STOP] Halting batch: {re}", file=sys.stderr)
                        all_exhausted = True
                        for f in futures:
                            f.cancel()
                        break
                    print(f"\n[ERROR] Worker exception: {re}", file=sys.stderr)
                except Exception as ex:
                    print(f"\n[ERROR] Worker exception: {ex}", file=sys.stderr)

                processed_count += 1
                if processed_count % 10 == 0 or processed_count == total_to_process:
                    elapsed = time.time() - start_time
                    rate = processed_count / (elapsed + 1e-6)
                    remaining = (total_to_process - processed_count) / (rate + 1e-6)
                    print(
                        f"   Progress: {processed_count}/{total_to_process} ({processed_count/total_to_process*100:.1f}%) | "
                        f"Rate: {rate:.1f} seq/s | ETA: {remaining/60:.1f} min",
                        end="\r",
                        flush=True,
                    )

        if not batch_results:
            if all_exhausted:
                print("\n[STOP] All configured API keys are exhausted for today.", file=sys.stderr)
                sys.exit(3)
            print(f"\n[WARN] No successful results for Batch {batch_idx:04d}")
            continue

        # Save Batch
        # 1. Per-residue matrices -> compressed NPZ
        npz_dict = {res["id"]: res["residue_matrix"] for res in batch_results}
        npz_file = batches_dir / f"batch_{batch_idx:04d}_residues.npz"
        np.savez_compressed(npz_file, **npz_dict)

        # 2. Sequence-level vectors -> Parquet table
        parquet_rows = [
            {
                "id": res["id"],
                "gene": res["gene"],
                "length": res["length"],
                "orig_length": res.get("orig_length", res["length"]),
                "vector": res["seq_vector"].tolist(),
            }
            for res in batch_results
        ]
        df_batch = pd.DataFrame(parquet_rows)
        parquet_file = batches_dir / f"batch_{batch_idx:04d}_sequence.parquet"
        df_batch.to_parquet(parquet_file, index=False)

        batch_time = time.time() - b_t0
        print(f"\n   [SAVED] Batch {batch_idx:04d}: {len(batch_results)} targets in {batch_time:.1f}s -> {parquet_file.name}, {npz_file.name}")

        # Live status tracking
        try:
            total_done = len(completed_ids) + processed_count
            status_data = {
                "dataset": args.dataset,
                "total_catalog": len(all_targets),
                "total_completed": total_done,
                "percent_complete": round((total_done / len(all_targets)) * 100, 2),
                "last_batch": batch_idx,
                "last_batch_size": len(batch_results),
                "rate_seq_per_sec": round(processed_count / max(time.time() - start_time, 0.001), 2),
                "updated_at": time.time(),
            }
            with open(dataset_dir / "status.json", "w") as sf:
                json.dump(status_data, sf, indent=2)
        except Exception:
            pass

        if all_exhausted:
            print("\n[STOP] All configured API keys are exhausted for today. Saved current batch before exit.", file=sys.stderr)
            sys.exit(3)

    total_time = time.time() - start_time
    print("\n" + "=" * 80)
    print(f"[COMPLETE] Processed {processed_count} targets in {total_time/60:.1f} minutes.")
    print(f"Outputs written to: {dataset_dir}")
    print("=" * 80)


if __name__ == "__main__":
    main()
