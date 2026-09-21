#!/usr/bin/env python3
"""
Step 2: Direct AlphaGenome API Batch Query for Microproteins.
Queries the AlphaGenome Atlas dense AVI scores across microprotein coding intervals.
Computes median, mean, max, and percentile Phred scores, with automatic batch checkpointing.
Treats all microproteins as equals.
"""

import argparse
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
CURATED_FILE = DATA_DIR / "curated" / "microprotein_intervals_all_7264.parquet"
SCREENING_DIR = DATA_DIR / "screening"
BATCHES_DIR = SCREENING_DIR / "avi_batches"
CONSOLIDATED_FILE = SCREENING_DIR / "microprotein_avi_scores.parquet"


def score_single_orf(client, row):
    """
    Query dense AVI scores across all exon intervals of a single ORF.
    Returns summary statistics dictionary or None on error.
    """
    chrom = row["chrom"]
    intervals_0b = row["intervals_0based"].split(";")
    
    all_quantiles = []
    all_variants = []

    for iv_str in intervals_0b:
        _, coords = iv_str.split(":")
        st, en = map(int, coords.split("-"))
        interval_obj = genome.Interval(chrom, st, en)

        # Query dense precomputed AVI score with retry on quota exhaustion
        max_retries = 5
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
                    wait_time = 30 * (attempt + 1)
                    print(f"\n[Quota limit reached] Pausing {wait_time}s before retry (attempt {attempt+1}/{max_retries})...")
                    time.sleep(wait_time)
                else:
                    raise e

        if res is None:
            raise RuntimeError(f"Failed to query interval {iv_str} after {max_retries} retries.")

        adata = res["AVI_SCORE"]
        q = adata.layers["quantiles"]
        all_quantiles.append(q.flatten())
        if "variant" in adata.obs:
            all_variants.extend(adata.obs["variant"].tolist())

    if not all_quantiles:
        return None

    quantiles = np.concatenate(all_quantiles)
    # Phred calibration: Phred = -10 * log10(1 - quantile)
    phreds = -10.0 * np.log10(np.clip(1.0 - quantiles, 1e-7, 1.0))

    # Identify peak variant
    max_idx = int(np.argmax(phreds))
    if max_idx < len(all_variants):
        v = all_variants[max_idx]
        if hasattr(v, "chromosome"):
            peak_variant = f"{v.chromosome}:{v.position}:{v.reference_bases}>{v.alternate_bases}"
        else:
            peak_variant = str(v)
    else:
        peak_variant = None
    max_phred = float(phreds[max_idx])

    # Construct clickable AlphaGenome Atlas URL
    min_st_1b = row["min_start_0based"] + 1
    max_en_1b = row["max_end_0based"]
    atlas_url = (
        f"https://deepmind.google.com/science/alphagenome/atlas"
        f"?q={chrom}:{min_st_1b}-{max_en_1b}&m=locus&lItems=avi,section:RNA_SEQ,section:DNASE"
    )

    return {
        "orf_id": row["orf_id"],
        "gene": row["gene"],
        "transcript": row["transcript"],
        "orf_biotype": row["orf_biotype"],
        "final_tier": row["final_tier"],
        "is_peptidein": row["is_peptidein"],
        "strand": row["strand"],
        "chrom": chrom,
        "n_exons": row["n_exons"],
        "length_aa": row["length_aa"],
        "length_bp": row["length_bp"],
        "n_variants_scored": len(phreds),
        "median_avi_phred": round(float(np.median(phreds)), 3),
        "mean_avi_phred": round(float(np.mean(phreds)), 3),
        "max_avi_phred": round(max_phred, 3),
        "p90_avi_phred": round(float(np.percentile(phreds, 90)), 3),
        "pct_phred_ge_15": round(float(np.mean(phreds >= 15.0) * 100), 2),
        "pct_phred_ge_20": round(float(np.mean(phreds >= 20.0) * 100), 2),
        "pct_phred_ge_30": round(float(np.mean(phreds >= 30.0) * 100), 2),
        "peak_variant": peak_variant,
        "atlas_url": atlas_url,
    }


def main():
    parser = argparse.ArgumentParser(description="Query AlphaGenome AVI scores for microproteins.")
    parser.add_argument("--batch-size", type=int, default=50, help="Number of ORFs per checkpoint batch.")
    parser.add_argument("--limit", type=int, default=None, help="Optional max ORFs to query (for quick runs).")
    parser.add_argument("--workers", type=int, default=6, help="Parallel worker threads for gRPC queries.")
    parser.add_argument("--start-idx", type=int, default=0, help="Starting index in the catalog.")
    args = parser.parse_args()

    print("=" * 70)
    print("Step 2: Direct AlphaGenome API Batch Query for Microproteins")
    print("=" * 70)

    # Load credentials
    dotenv.load_dotenv(BASE_DIR / ".env")
    dotenv.load_dotenv(Path.home() / ".env")
    api_key = os.environ.get("ALPHAGENOME_API_KEY")
    if not api_key:
        print("ERROR: ALPHAGENOME_API_KEY is missing from .env!")
        sys.exit(1)

    print("Connecting to AlphaGenome Atlas gRPC API...")
    client = atlas.create(api_key)
    print("Authenticated successfully.")

    BATCHES_DIR.mkdir(parents=True, exist_ok=True)

    print(f"Loading standardized catalog: {CURATED_FILE}...")
    df_all = pd.read_parquet(CURATED_FILE)
    print(f"Total catalog ORFs: {len(df_all):,}")

    # Check already scored ORFs across all saved batches
    scored_ids = set()
    existing_batches = list(BATCHES_DIR.glob("batch_*.parquet"))
    if existing_batches:
        print(f"Found {len(existing_batches)} existing batch files in {BATCHES_DIR}. Checking completed ORFs...")
        for b in existing_batches:
            try:
                b_df = pd.read_parquet(b)
                scored_ids.update(b_df["orf_id"].tolist())
            except Exception as e:
                print(f"Warning reading {b.name}: {e}")
        print(f"Already scored: {len(scored_ids):,} ORFs.")

    # Filter to remaining uncompleted ORFs
    df_todo = df_all[~df_all["orf_id"].isin(scored_ids)].copy()
    if args.start_idx > 0:
        df_todo = df_todo.iloc[args.start_idx:]
    if args.limit:
        df_todo = df_todo.head(args.limit)

    print(f"ORFs to query in this execution run: {len(df_todo):,}")
    if len(df_todo) == 0:
        print("All requested ORFs are already scored! Proceeding to consolidation...")
        consolidate_results()
        return

    # Process in batches
    batch_records = []
    total_processed = 0
    t0 = time.time()

    def process_row(r):
        try:
            return score_single_orf(client, r)
        except Exception as err:
            return {"orf_id": r["orf_id"], "error": str(err)}

    # Batch iterator
    batch_indices = [
        int(f.stem.split("_")[1])
        for f in existing_batches
        if len(f.stem.split("_")) > 1 and f.stem.split("_")[1].isdigit()
    ]
    batch_num = max(batch_indices) if batch_indices else 0
    rows = [row for _, row in df_todo.iterrows()]

    for i in range(0, len(rows), args.batch_size):
        chunk = rows[i : i + args.batch_size]
        chunk_results = []
        batch_num += 1

        print(f"\n--- Querying Batch {batch_num} ({len(chunk)} ORFs, workers={args.workers}) ---")
        batch_t0 = time.time()

        with ThreadPoolExecutor(max_workers=args.workers) as executor:
            future_to_orf = {executor.submit(process_row, r): r["orf_id"] for r in chunk}
            for fut in as_completed(future_to_orf):
                res = fut.result()
                if res and "error" not in res:
                    chunk_results.append(res)
                elif res and "error" in res:
                    print(f"  [!] Error on {res['orf_id']}: {res['error']}")

        batch_elapsed = time.time() - batch_t0
        total_processed += len(chunk_results)

        # Save batch checkpoint
        if chunk_results:
            batch_df = pd.DataFrame(chunk_results)
            batch_path = BATCHES_DIR / f"batch_{batch_num:05d}_{int(time.time())}.parquet"
            batch_df.to_parquet(batch_path, index=False)
            print(f"Saved checkpoint: {batch_path.name} ({len(batch_df)} ORFs) in {batch_elapsed:.1f}s ({batch_elapsed/len(chunk):.2f}s/ORF)")

            # Quick preview of top candidate in this batch
            top_in_batch = batch_df.sort_values("median_avi_phred", ascending=False).iloc[0]
            print(f"  Batch Top Constraint: {top_in_batch['orf_id']} ({top_in_batch['gene']}) - Median Phred: {top_in_batch['median_avi_phred']}, Max Phred: {top_in_batch['max_avi_phred']}")

    total_elapsed = time.time() - t0
    print("\n" + "=" * 70)
    print(f"Finished querying run: {total_processed} ORFs in {total_elapsed:.1f} s ({total_elapsed/max(1, total_processed):.2f} s/ORF)")
    print("=" * 70)

    consolidate_results()


def consolidate_results():
    """Consolidate all batch files into a single master dataset."""
    all_batch_files = sorted(BATCHES_DIR.glob("batch_*.parquet"))
    if not all_batch_files:
        print("No batch files to consolidate.")
        return

    print(f"\nConsolidating {len(all_batch_files)} batches...")
    dfs = [pd.read_parquet(f) for f in all_batch_files]
    df_master = pd.concat(dfs, ignore_index=True).drop_duplicates(subset=["orf_id"])
    df_master = df_master.sort_values("median_avi_phred", ascending=False)

    df_master.to_parquet(CONSOLIDATED_FILE, index=False)
    tsv_file = SCREENING_DIR / "microprotein_avi_scores.tsv"
    df_master.to_csv(tsv_file, sep="\t", index=False)

    print(f"Successfully consolidated {len(df_master):,} scored microproteins:")
    print(f" - Master Parquet: {CONSOLIDATED_FILE}")
    print(f" - Master TSV:     {tsv_file}")

    # Display summary statistics
    print("\nScore Distribution (Phred scale):")
    print(df_master[["median_avi_phred", "mean_avi_phred", "max_avi_phred", "p90_avi_phred"]].describe())

    # Top 5 standout candidates
    print("\nTop 5 Standout Microproteins under Highest Negative Selection (Median AVI Phred):")
    top5 = df_master.head(5)[["orf_id", "gene", "final_tier", "orf_biotype", "length_aa", "median_avi_phred", "max_avi_phred", "atlas_url"]]
    print(top5.to_string(index=False))


if __name__ == "__main__":
    main()
