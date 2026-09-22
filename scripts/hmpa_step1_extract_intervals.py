#!/usr/bin/env python3
"""
Step 1: Rapid Coordinate Extraction & Standardization for 617,463 HMPA Microproteins.
Reads hmpa_structure_alphafold_scores.csv, validates genomic coordinates,
converts to 0-based half-open intervals for AlphaGenome API, and saves to curated Parquet.
"""

import sys
import time
from pathlib import Path
import numpy as np
import pandas as pd

BASE_DIR = Path(__file__).resolve().parent.parent
DATA_DIR = BASE_DIR / "data"
RAW_FILE = DATA_DIR / "raw" / "hmpa" / "hmpa_structure_alphafold_scores.csv"
CURATED_DIR = DATA_DIR / "curated"
OUTPUT_PARQUET = CURATED_DIR / "hmpa_microprotein_intervals_all_617k.parquet"


def main():
    print("=" * 70)
    print("HMPA Step 1: Standardizing Microprotein Genomic Intervals (617k)")
    print("=" * 70)

    CURATED_DIR.mkdir(parents=True, exist_ok=True)

    if not RAW_FILE.exists():
        print(f"ERROR: Raw HMPA file not found at {RAW_FILE}", file=sys.stderr)
        sys.exit(1)

    t0 = time.time()
    print(f"Reading raw HMPA table: {RAW_FILE}...")
    cols = [
        "smORF_id",
        "Chrome",
        "Start",
        "End",
        "length",
        "plddt",
        "nearest_gene",
        "Sequence",
    ]
    df = pd.read_csv(RAW_FILE, usecols=cols, low_memory=False)
    print(f"Loaded {len(df):,} rows in {time.time() - t0:.1f}s.")

    # Clean non-numeric or malformed rows
    df["Start_num"] = pd.to_numeric(df["Start"], errors="coerce")
    df["End_num"] = pd.to_numeric(df["End"], errors="coerce")
    df["length_num"] = pd.to_numeric(df["length"], errors="coerce")
    df["plddt_num"] = pd.to_numeric(df["plddt"], errors="coerce")

    valid_mask = (
        df["Start_num"].notna()
        & df["End_num"].notna()
        & df["Chrome"].str.startswith("chr", na=False)
        & (df["End_num"] > df["Start_num"])
    )

    invalid_count = len(df) - valid_mask.sum()
    print(f"Valid coordinates: {valid_mask.sum():,} | Filtered non-numeric/invalid: {invalid_count:,}")

    clean_df = df[valid_mask].copy()

    # Cast coordinates
    # Note: Genomic coordinates in HMPA are 1-based closed (e.g. Start 9943442, End 9943508).
    # For AlphaGenome API: genome.Interval requires 0-based half-open (start - 1, end).
    start_1b = clean_df["Start_num"].astype(np.int64)
    end_1b = clean_df["End_num"].astype(np.int64)
    start_0b = (start_1b - 1).astype(np.int64)
    end_0b = end_1b

    clean_df["chrom"] = clean_df["Chrome"].astype(str)
    clean_df["start_1b"] = start_1b
    clean_df["end_1b"] = end_1b
    clean_df["start_0b"] = start_0b
    clean_df["end_0b"] = end_0b
    clean_df["span_bp"] = end_0b - start_0b
    clean_df["length_aa"] = clean_df["length_num"].fillna(0).astype(int)
    clean_df["plddt"] = clean_df["plddt_num"].round(2)
    clean_df["nearest_gene"] = clean_df["nearest_gene"].fillna("").astype(str)

    # Format 0-based interval string: e.g. "chr1:9943441-9943508"
    clean_df["interval_0b"] = (
        clean_df["chrom"] + ":" + clean_df["start_0b"].astype(str) + "-" + clean_df["end_0b"].astype(str)
    )

    out_cols = [
        "smORF_id",
        "chrom",
        "start_1b",
        "end_1b",
        "start_0b",
        "end_0b",
        "span_bp",
        "length_aa",
        "plddt",
        "nearest_gene",
        "interval_0b",
        "Sequence",
    ]
    out_df = clean_df[out_cols].reset_index(drop=True)

    print(f"Writing curated parquet: {OUTPUT_PARQUET}...")
    out_df.to_parquet(OUTPUT_PARQUET, index=False, compression="snappy")
    size_mb = OUTPUT_PARQUET.stat().st_size / (1024 * 1024)
    print(f"Successfully saved {len(out_df):,} curated intervals ({size_mb:.2f} MB) in {time.time() - t0:.1f}s.")


if __name__ == "__main__":
    main()
