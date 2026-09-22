#!/usr/bin/env python3
"""
High-Throughput Biohub ESMFold2 Cofolding Pipeline for Microprotein-Partner Complexes.

Folds multi-chain protein complexes (Chain A: Microprotein, Chain B: Interacting Partner)
using Biohub ESMFold2 `fold_all_atom` API.
Outputs PDB coordinate files, PAE matrices, and structured metrics (ipTM, pTM, per-chain pLDDT).
Includes persistent checkpointing, length-limit handling (768 aa online limit), and exponential backoff.
"""

import argparse
import json
import os
import sys
import time
from pathlib import Path
import pandas as pd
import numpy as np

# Ensure esm SDK is available
try:
    from esm.sdk.forge import SequenceStructureForgeInferenceClient, FoldingConfig
    from esm.utils.structure.input_builder import StructurePredictionInput, ProteinInput
    from esm.sdk.api import ESMProteinError
except ImportError as e:
    print(f"[ERROR] Failed to import ESM SDK: {e}", file=sys.stderr)
    sys.exit(1)


def load_api_keys() -> list[str]:
    """Retrieve all configured ESM API keys (primary, alt, and numbered) from environment or .env files."""
    keys = []
    for var, val in os.environ.items():
        if var.startswith("ESM_API_KEY") and val.strip() and val.strip() not in keys:
            keys.append(val.strip())

    for p in [Path.cwd() / ".env", Path.home() / ".env"]:
        if p.is_file():
            try:
                with open(p, "r", encoding="utf-8") as f:
                    for line in f:
                        line = line.strip()
                        if line.startswith("#") or "=" not in line:
                            continue
                        var, val = line.split("=", 1)
                        var = var.strip()
                        val = val.strip().strip("\"'")
                        if var.startswith("ESM_API_KEY") and val and val not in keys:
                            keys.append(val)
            except Exception:
                pass
    return keys


def load_api_key() -> str:
    """Retrieve primary or alternative ESM API key."""
    keys = load_api_keys()
    if keys:
        return keys[0]
    print("[ERROR] Neither ESM_API_KEY nor ESM_API_KEY_ALT found in environment, .env, or ~/.env", file=sys.stderr)
    sys.exit(1)


class RotatingESMClient:
    """Multi-key client with automatic rotation when rate limits or quotas are hit."""
    def __init__(self, model: str, tokens: list[str]):
        self.model = model
        self.tokens = tokens
        self.current_idx = 0
        self.clients = [
            SequenceStructureForgeInferenceClient(model=model, token=t)
            for t in tokens
        ]

    @property
    def current_client(self):
        return self.clients[self.current_idx]

    def rotate(self):
        if len(self.clients) > 1:
            self.current_idx = (self.current_idx + 1) % len(self.clients)
            print(f"\n    [KEY ROTATION] Biohub API quota/rate limit reached. Switched to key #{self.current_idx + 1}/{len(self.clients)}", file=sys.stderr)

    def fold_all_atom(self, *args, **kwargs):
        return self.current_client.fold_all_atom(*args, **kwargs)


def load_checkpoint(state_file: Path) -> dict:
    """Load previously completed or skipped pairs from the JSONL checkpoint file."""
    completed = {}
    if not state_file.is_file():
        return completed
    with open(state_file, "r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            try:
                rec = json.loads(line)
                if "pair_id" in rec:
                    completed[rec["pair_id"]] = rec
            except Exception:
                pass
    return completed


def save_checkpoint_record(state_file: Path, record: dict):
    """Append a single record to the state ledger."""
    with open(state_file, "a", encoding="utf-8") as f:
        f.write(json.dumps(record) + "\n")


def fold_complex_with_retry(client, spi, config, max_retries=4, base_backoff=3.0):
    """Execute fold_all_atom with exponential backoff and automatic key rotation on rate limits."""
    for attempt in range(1, max_retries + 1):
        try:
            res = client.fold_all_atom(spi, config=config)
            if isinstance(res, ESMProteinError):
                if res.error_code == 422:
                    raise ValueError(f"ESMProteinError [422]: {res.error_msg}")
                if res.error_code in (429, 403) or "quota" in str(res.error_msg).lower() or "rate" in str(res.error_msg).lower():
                    if hasattr(client, "rotate"):
                        client.rotate()
                        time.sleep(1.0)
                        continue
                raise RuntimeError(f"ESMProteinError [{res.error_code}]: {res.error_msg}")
            return res
        except ValueError as e:
            # Fatal validation error (e.g. sequence length), do not retry
            raise e
        except Exception as e:
            err_str = str(e)
            if "422" in err_str or "exceeds maximum allowed" in err_str:
                raise ValueError(err_str)
            if "429" in err_str or "quota" in err_str.lower() or "rate limit" in err_str.lower() or "403" in err_str:
                if hasattr(client, "rotate"):
                    client.rotate()
                    time.sleep(1.0)
                    continue
            if attempt == max_retries:
                raise e
            wait_time = base_backoff * (2 ** (attempt - 1))
            print(f"    [WARN] Attempt {attempt} failed: {err_str[:120]}... Backing off {wait_time:.1f}s")
            time.sleep(wait_time)



def sync_summary_tables(candidates_df: pd.DataFrame, state_file: Path, out_parquet: Path, out_tsv: Path):
    """Merge checkpoint state with full metadata and update summary tables."""
    completed_records = load_checkpoint(state_file)
    if not completed_records:
        return

    records_list = list(completed_records.values())
    state_df = pd.DataFrame(records_list)

    cols_from_state = [
        "pair_id", "status", "ptm", "iptm", "mean_plddt",
        "microprotein_plddt", "partner_plddt", "binding_classification",
        "pdb_path", "pae_path", "model_used", "elapsed_seconds", "folded_at_utc"
    ]
    cols_to_merge = [c for c in cols_from_state if c in state_df.columns]
    state_sub = state_df[cols_to_merge].drop_duplicates(subset=["pair_id"], keep="last")

    merged = pd.merge(candidates_df, state_sub, on="pair_id", how="left")
    merged.to_parquet(out_parquet, index=False)
    merged.to_csv(out_tsv, sep="\t", index=False)
    print(f"[SYNC] Summary tables updated: {len(state_sub)} completed records stored.")


def main():
    parser = argparse.ArgumentParser(description="Run ESMFold2 cofolding campaign for microprotein-partner complexes")
    parser.add_argument("--input-file", default=None, help="Custom input parquet/tsv file of candidate pairs")
    parser.add_argument("--output-prefix", default=None, help="Custom output prefix for summary parquet/tsv")
    parser.add_argument("--phase", choices=["1", "2", "all"], default="all",
                        help="Phase 1 (partner_rank=1, n=50), Phase 2 (partner_rank>1, n=450), or all (n=500)")
    parser.add_argument("--limit", type=int, default=None, help="Optional limit on number of complexes to process")
    parser.add_argument("--model", default="esmfold2-fast-2026-05", choices=["esmfold2-fast-2026-05", "esmfold2-2026-05"],
                        help="ESMFold2 model version")
    parser.add_argument("--num-loops", type=int, default=20, help="Number of trunk loops (default: 20)")
    parser.add_argument("--num-steps", type=int, default=100, help="Diffusion sampling steps (default: 100)")
    parser.add_argument("--delay", type=float, default=0.5, help="Polite delay between API calls in seconds (default: 0.5)")
    parser.add_argument("--max-length", type=int, default=768, help="Max complex length accepted by API (default: 768)")
    parser.add_argument("--include-pae", action="store_true", default=True, help="Include PAE matrix")

    args = parser.parse_args()
    api_keys = load_api_keys()
    if not api_keys:
        print("[ERROR] Neither ESM_API_KEY nor ESM_API_KEY_ALT found in environment, .env, or ~/.env", file=sys.stderr)
        sys.exit(1)

    base_dir = Path.cwd()
    if args.input_file:
        in_path = Path(args.input_file)
        if in_path.suffix == ".parquet":
            df = pd.read_parquet(in_path)
        else:
            df = pd.read_csv(in_path, sep="\t")
        subset = df.copy()
        phase_desc = f"Custom Input File ({in_path.name})"
        out_name = args.output_prefix if args.output_prefix else in_path.stem.replace("_candidates_", "_cofolding_summary_")
        out_parquet = base_dir / "reports" / f"{out_name}.parquet"
        out_tsv = base_dir / "reports" / f"{out_name}.tsv"
    else:
        data_path = base_dir / "reports" / "top50_microprotein_binding_partners_top10.parquet"
        if not data_path.is_file():
            data_path = base_dir / "reports" / "top50_microprotein_binding_partners_top10.tsv"
            df = pd.read_csv(data_path, sep="\t")
        else:
            df = pd.read_parquet(data_path)

        out_parquet = base_dir / "reports" / "esmfold2_cofolding_summary.parquet"
        out_tsv = base_dir / "reports" / "esmfold2_cofolding_summary.tsv"

        # Filter by phase
        if args.phase == "1":
            subset = df[df["partner_rank"] == 1].copy()
            phase_desc = "Phase 1 (Top 1 Cognate Partner per microprotein, n=50)"
        elif args.phase == "2":
            subset = df[df["partner_rank"] > 1].copy()
            phase_desc = "Phase 2 (Partners #2-#10, n=450)"
        else:
            subset = df.copy()
            phase_desc = "Full Portfolio (All 500 Candidate Complexes)"

    if args.limit:
        subset = subset.head(args.limit)

    print(f"=== ESMFold2 Cofolding Campaign ===")
    print(f"Selected Scope:  {phase_desc}")
    print(f"Target Count:    {len(subset)} complexes")
    print(f"Model:           {args.model} (loops={args.num_loops}, steps={args.num_steps})")
    print(f"Max Seq Length:  {args.max_length} aa (Biohub Platform Limit)")

    # Setup directories
    struct_dir = base_dir / "structures"
    pdb_dir = struct_dir / "pdbs"
    pae_dir = struct_dir / "pae"
    pdb_dir.mkdir(parents=True, exist_ok=True)
    pae_dir.mkdir(parents=True, exist_ok=True)

    state_file = struct_dir / "cofolding_state.jsonl"
    out_parquet = base_dir / "reports" / "esmfold2_cofolding_summary.parquet"
    out_tsv = base_dir / "reports" / "esmfold2_cofolding_summary.tsv"

    completed_map = load_checkpoint(state_file)
    print(f"Already in Ledger: {len(completed_map)} complexes")

    client = RotatingESMClient(model=args.model, tokens=api_keys)
    folding_config = FoldingConfig(
        num_loops=args.num_loops,
        num_sampling_steps=args.num_steps,
        include_pae=args.include_pae,
        include_pair_chains_iptm=True
    )

    t_start = time.time()
    n_processed = 0
    n_skipped = 0
    n_failed = 0

    for i, (_, row) in enumerate(subset.iterrows(), start=1):
        pair_id = row["pair_id"]
        total_len = row["esmfold_total_length"]
        pdb_target = pdb_dir / f"{pair_id}.pdb"
        pae_target = pae_dir / f"{pair_id}_pae.json"

        # Check if already processed (successful PDB or known skip)
        if pair_id in completed_map:
            prev = completed_map[pair_id]
            if prev.get("status") == "success" and pdb_target.is_file() and pdb_target.stat().st_size > 0:
                n_skipped += 1
                continue
            if prev.get("status") == "exceeds_api_limit_768":
                n_skipped += 1
                continue

        # Check API length limit
        if total_len > args.max_length:
            print(f"\n[{i}/{len(subset)}] [SKIP] {pair_id} ({row['microprotein_gene']} + {row['partner_gene_symbol']})")
            print(f"  Length: {row['microprotein_length']} aa + {row['partner_length']} aa = {total_len} aa exceeds API limit ({args.max_length} aa).")
            record = {
                "pair_id": pair_id,
                "microprotein_gene": row["microprotein_gene"],
                "partner_gene_symbol": row["partner_gene_symbol"],
                "status": "exceeds_api_limit_768",
                "ptm": None,
                "iptm": None,
                "mean_plddt": None,
                "microprotein_plddt": None,
                "partner_plddt": None,
                "binding_classification": "Exceeds API Context Limit (>768 aa)",
                "pdb_path": None,
                "pae_path": None,
                "model_used": args.model,
                "elapsed_seconds": 0.0,
                "folded_at_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())
            }
            save_checkpoint_record(state_file, record)
            completed_map[pair_id] = record
            n_skipped += 1
            continue

        print(f"\n[{i}/{len(subset)}] Folding {pair_id} ({row['microprotein_gene']} + {row['partner_gene_symbol']}) ...")
        print(f"  Length: {row['microprotein_length']} aa + {row['partner_length']} aa = {total_len} aa | Tier: {row['esmfold_feasibility']}")

        spi = StructurePredictionInput(sequences=[
            ProteinInput(id="A", sequence=row["microprotein_sequence"]),
            ProteinInput(id="B", sequence=row["partner_sequence"])
        ])

        try:
            t0 = time.time()
            res = fold_complex_with_retry(client, spi, folding_config)
            t_elapsed = time.time() - t0

            # Extract metrics
            ptm_val = float(res.ptm) if res.ptm is not None else None
            iptm_val = float(res.iptm) if res.iptm is not None else None

            # pLDDT extraction
            plddt_raw = res.plddt
            if hasattr(plddt_raw, "tolist"):
                plddt_list = plddt_raw.tolist()
            elif isinstance(plddt_raw, list):
                plddt_list = plddt_raw
            else:
                plddt_list = []

            valid_plddt = [x for x in plddt_list if x is not None and not np.isnan(x)]
            mean_plddt = (sum(valid_plddt) / len(valid_plddt)) if valid_plddt else None

            # Chain specific pLDDT
            l_a = row["microprotein_length"]
            plddt_a_list = [x for x in plddt_list[:l_a] if x is not None and not np.isnan(x)]
            plddt_b_list = [x for x in plddt_list[l_a:] if x is not None and not np.isnan(x)]
            plddt_a = (sum(plddt_a_list) / len(plddt_a_list)) if plddt_a_list else None
            plddt_b = (sum(plddt_b_list) / len(plddt_b_list)) if plddt_b_list else None

            # Classification
            if iptm_val is not None:
                if iptm_val >= 0.75:
                    b_class = "High Confidence (ipTM >= 0.75)"
                elif iptm_val >= 0.60:
                    b_class = "Plausible Interface (0.60 <= ipTM < 0.75)"
                elif iptm_val >= 0.45:
                    b_class = "Weak / Transient (0.45 <= ipTM < 0.60)"
                else:
                    b_class = "Non-Interacting / Unbound (ipTM < 0.45)"
            else:
                b_class = "Undetermined"

            # Export PDB
            try:
                pc = res.complex.to_protein_complex()
                pdb_str = pc.to_pdb_string()
                with open(pdb_target, "w", encoding="utf-8") as f_pdb:
                    f_pdb.write(pdb_str)
            except Exception as e_pdb:
                print(f"    [WARN] PDB export failed, fallback to mmCIF: {e_pdb}")
                try:
                    cif_target = pdb_dir / f"{pair_id}.cif"
                    with open(cif_target, "w", encoding="utf-8") as f_cif:
                        f_cif.write(res.complex.to_mmcif())
                    pdb_target = cif_target
                except Exception as e_cif:
                    print(f"    [ERROR] Coordinate export failed completely: {e_cif}")

            # Export PAE
            if args.include_pae and hasattr(res, "pae") and res.pae is not None:
                pae_raw = res.pae
                pae_matrix = pae_raw.tolist() if hasattr(pae_raw, "tolist") else pae_raw
                with open(pae_target, "w", encoding="utf-8") as f_pae:
                    json.dump(pae_matrix, f_pae)

            record = {
                "pair_id": pair_id,
                "microprotein_gene": row["microprotein_gene"],
                "partner_gene_symbol": row["partner_gene_symbol"],
                "status": "success",
                "ptm": round(ptm_val, 4) if ptm_val else None,
                "iptm": round(iptm_val, 4) if iptm_val else None,
                "mean_plddt": round(mean_plddt, 2) if mean_plddt else None,
                "microprotein_plddt": round(plddt_a, 2) if plddt_a else None,
                "partner_plddt": round(plddt_b, 2) if plddt_b else None,
                "binding_classification": b_class,
                "pdb_path": str(pdb_target),
                "pae_path": str(pae_target) if pae_target.exists() else None,
                "model_used": args.model,
                "elapsed_seconds": round(t_elapsed, 2),
                "folded_at_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())
            }

            save_checkpoint_record(state_file, record)
            completed_map[pair_id] = record
            n_processed += 1

            plddt_a_str = f"{plddt_a:.1f}" if plddt_a else "N/A"
            plddt_b_str = f"{plddt_b:.1f}" if plddt_b else "N/A"
            iptm_str = f"{iptm_val:.3f}" if iptm_val else "N/A"
            ptm_str = f"{ptm_val:.3f}" if ptm_val else "N/A"
            print(f"  [OK] Done in {t_elapsed:.1f}s | ipTM: {iptm_str} | pTM: {ptm_str} | pLDDT: uProt={plddt_a_str}, Partner={plddt_b_str} | {b_class}")

            # Intermediate sync every 5 predictions
            if n_processed % 5 == 0:
                sync_summary_tables(df, state_file, out_parquet, out_tsv)

            if args.delay > 0:
                time.sleep(args.delay)

        except Exception as err:
            print(f"  [FAIL] Failed to fold {pair_id}: {err}")
            fail_record = {
                "pair_id": pair_id,
                "microprotein_gene": row["microprotein_gene"],
                "partner_gene_symbol": row["partner_gene_symbol"],
                "status": "error",
                "error_message": str(err),
                "folded_at_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())
            }
            save_checkpoint_record(state_file, fail_record)
            n_failed += 1

    # Final synchronization
    sync_summary_tables(df, state_file, out_parquet, out_tsv)
    total_time = time.time() - t_start
    print(f"\n=== Execution Complete ===")
    print(f"Total Elapsed:   {total_time/60:.2f} minutes")
    print(f"Processed:       {n_processed}")
    print(f"Skipped (Prior): {n_skipped}")
    print(f"Failed:          {n_failed}")


if __name__ == "__main__":
    main()
