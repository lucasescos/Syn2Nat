#!/usr/bin/env python3
"""
Step 2: High-Throughput AlphaGenome AVI Batch Query for 617k HMPA Microproteins.

Queries dense AlphaGenome Variant Impact (AVI) purifying selection scores across
microprotein coding intervals using Google DeepMind's AlphaGenome Atlas API.

Features:
- Checkpointed Parquet batching (e.g. 500 ORFs per batch).
- Fast resumption: automatically detects already completed ORFs and resumes where left off.
- Concurrent multi-threaded queries with ThreadPoolExecutor.
- Automatic exponential backoff on gRPC rate/quota limits (RESOURCE_EXHAUSTED).
- Real-time JSON telemetry status monitor (`hmpa_avi_status.json`).
"""

import argparse
import datetime
import json
import os
import sys
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import dotenv
import numpy as np
import pandas as pd
from tqdm import tqdm

from alphagenome.atlas import atlas
from alphagenome.data import genome

BASE_DIR = Path(__file__).resolve().parent.parent
DATA_DIR = BASE_DIR / "data"
CURATED_FILE = DATA_DIR / "curated" / "hmpa_microprotein_intervals_all_617k.parquet"
SCREENING_DIR = DATA_DIR / "screening"
BATCHES_DIR = SCREENING_DIR / "hmpa_avi_batches"
STATUS_FILE = SCREENING_DIR / "hmpa_avi_status.json"


def phred_from_quantiles(quantiles: np.ndarray) -> np.ndarray:
    """Calculate calibrated Phred score: -10 * log10(1 - quantile)."""
    return -10.0 * np.log10(np.clip(1.0 - quantiles, 1e-7, 1.0))


def query_single_microprotein(client, row: dict) -> dict:
    """
    Query dense AVI scores across a single microprotein genomic interval.
    Returns summary metrics dict or error info.
    """
    chrom = str(row["chrom"])
    st_0b = int(row["start_0b"])
    en_0b = int(row["end_0b"])
    smorf_id = str(row["smORF_id"])

    interval_obj = genome.Interval(chrom, st_0b, en_0b)

    max_retries = 6
    res = None
    for attempt in range(max_retries):
        try:
            res = client.query_interval(
                interval_obj,
                requested_scorers=["AVI_SCORE"],
                progress_bar=False,
            )
            break
        except Exception as e:
            err_str = str(e)
            if "RESOURCE_EXHAUSTED" in err_str or "Quota exceeded" in err_str:
                wait_time = 15 * (2**attempt)
                tqdm.write(f"[{smorf_id}] Quota limited. Pausing {wait_time}s before retry ({attempt+1}/{max_retries})...")
                time.sleep(wait_time)
            elif "DEADLINE_EXCEEDED" in err_str or "UNAVAILABLE" in err_str:
                wait_time = 5 * (attempt + 1)
                tqdm.write(f"[{smorf_id}] Transient connection error. Retrying in {wait_time}s...")
                time.sleep(wait_time)
            else:
                tqdm.write(f"[{smorf_id}] Query error: {err_str}")
                return {
                    "smORF_id": smorf_id,
                    "error": err_str,
                    "status": "error",
                }

    if res is None or "AVI_SCORE" not in res:
        return {
            "smORF_id": smorf_id,
            "error": f"Failed after {max_retries} attempts",
            "status": "failed",
        }

    adata = res["AVI_SCORE"]
    quantiles = adata.layers["quantiles"].flatten()
    if len(quantiles) == 0:
        return {
            "smORF_id": smorf_id,
            "error": "No variants returned in interval",
            "status": "empty",
        }

    phreds = phred_from_quantiles(quantiles)

    # Locate peak deleterious hotspot variant
    variants = adata.obs["variant"].tolist() if "variant" in adata.obs else []
    max_idx = int(np.argmax(phreds))
    if max_idx < len(variants):
        v = variants[max_idx]
        if hasattr(v, "chromosome"):
            peak_variant = f"{v.chromosome}:{v.position}:{v.reference_bases}>{v.alternate_bases}"
        else:
            peak_variant = str(v)
    else:
        peak_variant = None

    # Construct clickable AlphaGenome Atlas Deep-Link
    min_st_1b = st_0b + 1
    max_en_1b = en_0b
    atlas_url = (
        f"https://deepmind.google.com/science/alphagenome/atlas"
        f"?q={chrom}:{min_st_1b}-{max_en_1b}&m=locus&lItems=avi,section:RNA_SEQ,section:DNASE"
    )

    return {
        "smORF_id": smorf_id,
        "chrom": chrom,
        "start_1b": int(row["start_1b"]),
        "end_1b": int(row["end_1b"]),
        "span_bp": int(row["span_bp"]),
        "length_aa": int(row["length_aa"]),
        "plddt": float(row["plddt"]),
        "nearest_gene": str(row["nearest_gene"]),
        "n_variants_scored": len(phreds),
        "median_avi_phred": round(float(np.median(phreds)), 3),
        "mean_avi_phred": round(float(np.mean(phreds)), 3),
        "max_avi_phred": round(float(phreds[max_idx]), 3),
        "p90_avi_phred": round(float(np.percentile(phreds, 90)), 3),
        "peak_variant": peak_variant,
        "atlas_url": atlas_url,
        "status": "success",
    }


def get_completed_ids() -> set:
    """Scan existing batch files and return set of completed smORF_ids."""
    if not BATCHES_DIR.exists():
        return set()

    completed = set()
    for batch_path in BATCHES_DIR.glob("batch_*.parquet"):
        try:
            df = pd.read_parquet(batch_path, columns=["smORF_id", "status"])
            succ = df[df["status"] == "success"]["smORF_id"].tolist()
            completed.update(succ)
        except Exception:
            continue
    return completed


def update_status(total_catalog: int, completed_count: int, start_time: float, active_workers: int):
    """Write live telemetry JSON status."""
    elapsed = max(time.time() - start_time, 1.0)
    rate = completed_count / elapsed if completed_count > 0 else 0.0
    remaining = max(total_catalog - completed_count, 0)
    eta_hours = (remaining / rate / 3600) if rate > 0 else None

    status_data = {
        "total_catalog": total_catalog,
        "completed_count": completed_count,
        "pct_completed": round((completed_count / total_catalog) * 100, 2) if total_catalog > 0 else 0.0,
        "rate_orfs_per_sec": round(rate, 2),
        "eta_hours": round(eta_hours, 1) if eta_hours else None,
        "active_workers": active_workers,
        "last_updated": datetime.datetime.now().isoformat(),
    }

    try:
        with open(STATUS_FILE, "w", encoding="utf-8") as f:
            json.dump(status_data, f, indent=2)
    except Exception:
        pass


def main():
    parser = argparse.ArgumentParser(description="HMPA AlphaGenome AVI High-Throughput Batch Screener.")
    parser.add_argument("--workers", type=int, default=6, help="Concurrent worker threads (default: 6)")
    parser.add_argument("--batch-size", type=int, default=500, help="Checkpointed batch size (default: 500)")
    parser.add_argument("--max-orfs", type=int, default=None, help="Stop after scoring N orfs (optional pilot)")
    args = parser.parse_args()

    # Load environment
    dotenv.load_dotenv(BASE_DIR / ".env")
    api_key = os.environ.get("ALPHAGENOME_API_KEY")
    if not api_key:
        print("ERROR: ALPHAGENOME_API_KEY is missing from environment or .env!", file=sys.stderr)
        sys.exit(1)

    BATCHES_DIR.mkdir(parents=True, exist_ok=True)
    SCREENING_DIR.mkdir(parents=True, exist_ok=True)

    print("=" * 80)
    print("HMPA AlphaGenome AVI High-Throughput Screening Daemon")
    print(f"Workers: {args.workers} | Batch Size: {args.batch_size} | Max ORFs: {args.max_orfs or 'All 617k'}")
    print("=" * 80)

    print(f"Loading curated intervals: {CURATED_FILE}...")
    df_catalog = pd.read_parquet(CURATED_FILE)
    total_in_file = len(df_catalog)
    print(f"Loaded {total_in_file:,} candidate microprotein intervals.")

    # Detect previously completed ORFs
    completed_ids = get_completed_ids()
    print(f"Detected {len(completed_ids):,} already completed smORFs in {BATCHES_DIR}.")

    # Filter for uncompleted ORFs
    uncompleted_df = df_catalog[~df_catalog["smORF_id"].isin(completed_ids)].reset_index(drop=True)
    print(f"Remaining candidates to score: {len(uncompleted_df):,}")

    if len(uncompleted_df) == 0:
        print("All candidates are already completed! Nothing to query.")
        return

    if args.max_orfs:
        uncompleted_df = uncompleted_df.head(args.max_orfs)
        print(f"Subsetting to first {len(uncompleted_df):,} candidates for this session.")

    client = atlas.create(api_key)
    print("AlphaGenome Atlas client initialized successfully.")

    start_time = time.time()
    records_buffer = []
    total_completed = len(completed_ids)
    batch_counter = len(list(BATCHES_DIR.glob("batch_*.parquet")))

    rows_to_process = uncompleted_df.to_dict(orient="records")

    with ThreadPoolExecutor(max_workers=args.workers) as executor:
        # Submit in sliding chunks to keep memory compact
        chunk_size = args.batch_size * 2
        pbar = tqdm(total=len(rows_to_process), desc="Scoring smORFs", unit="orf")

        for i in range(0, len(rows_to_process), chunk_size):
            chunk = rows_to_process[i : i + chunk_size]
            future_to_row = {
                executor.submit(query_single_microprotein, client, row): row
                for row in chunk
            }

            for future in as_completed(future_to_row):
                row_res = future.result()
                if row_res:
                    records_buffer.append(row_res)
                    if row_res.get("status") == "success":
                        total_completed += 1

                pbar.update(1)

                # Flush batch to Parquet
                if len(records_buffer) >= args.batch_size:
                    batch_counter += 1
                    batch_df = pd.DataFrame(records_buffer)
                    out_path = BATCHES_DIR / f"batch_{batch_counter:06d}.parquet"
                    batch_df.to_parquet(out_path, index=False, compression="snappy")
                    records_buffer = []
                    update_status(total_in_file, total_completed, start_time, args.workers)

        # Flush any remaining in buffer
        if records_buffer:
            batch_counter += 1
            batch_df = pd.DataFrame(records_buffer)
            out_path = BATCHES_DIR / f"batch_{batch_counter:06d}.parquet"
            batch_df.to_parquet(out_path, index=False, compression="snappy")
            records_buffer = []
            update_status(total_in_file, total_completed, start_time, args.workers)

        pbar.close()

    elapsed = time.time() - start_time
    print(f"\nDone! Scored in {elapsed:.1f}s. Total completed catalog: {total_completed:,} / {total_in_file:,}.")


if __name__ == "__main__":
    main()
