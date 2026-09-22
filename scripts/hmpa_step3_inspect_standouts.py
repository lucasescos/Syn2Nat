#!/usr/bin/env python3
"""
Step 3: Progressive Standout Inspector for HMPA AlphaGenome Screening.
Aggregates completed batches from data/screening/hmpa_avi_batches/,
calculates purifying selection distribution, and displays top ranked microproteins.
"""

import argparse
from pathlib import Path
import pandas as pd
import numpy as np

BASE_DIR = Path(__file__).resolve().parent.parent
DATA_DIR = BASE_DIR / "data"
BATCHES_DIR = DATA_DIR / "screening" / "hmpa_avi_batches"
TOTAL_CATALOG = 617262


def phred_to_percentile_str(phred: float) -> str:
    """Format Phred score into human-readable top percentile string."""
    if phred >= 40:
        return f"{phred:.2f} (Top <0.01%)"
    elif phred >= 30:
        return f"{phred:.2f} (Top <0.10%)"
    elif phred >= 20:
        return f"{phred:.2f} (Top <1.0%)"
    elif phred >= 15:
        return f"{phred:.2f} (Top <3.2%)"
    elif phred >= 10:
        return f"{phred:.2f} (Top <10%)"
    else:
        return f"{phred:.2f} (Bottom 90%)"


def main():
    parser = argparse.ArgumentParser(description="Inspect standout microproteins from HMPA screening.")
    parser.add_argument("--top-n", type=int, default=20, help="Number of top candidates to display (default: 20)")
    parser.add_argument("--min-plddt", type=float, default=None, help="Filter for minimum AlphaFold pLDDT (e.g. 70.0)")
    args = parser.parse_args()

    batch_files = sorted(BATCHES_DIR.glob("batch_*.parquet"))
    if not batch_files:
        print(f"No batch files found in {BATCHES_DIR}. Please start screening first!")
        return

    print(f"Reading {len(batch_files)} batch files from {BATCHES_DIR}...")
    dfs = [pd.read_parquet(f) for f in batch_files]
    df = pd.concat(dfs, ignore_index=True)
    df = df[df["status"] == "success"].drop_duplicates(subset=["smORF_id"]).reset_index(drop=True)

    pct_done = (len(df) / TOTAL_CATALOG) * 100

    print("=" * 85)
    print(f"HMPA AlphaGenome Screening Progress: {len(df):,} / {TOTAL_CATALOG:,} ({pct_done:.2f}%)")
    print("=" * 85)

    if len(df) == 0:
        print("No successful records found yet.")
        return

    # Constraint summary
    top_1pct = (df["median_avi_phred"] >= 20.0).sum()
    top_3pct = (df["median_avi_phred"] >= 15.0).sum()
    top_10pct = (df["median_avi_phred"] >= 10.0).sum()

    print(f"\nConstraint Distribution Across Scored Candidates ({len(df):,} total):")
    print(f" - Top 1% Constraint  (Median Phred >= 20): {top_1pct:,} ({top_1pct/len(df)*100:.1f}%)")
    print(f" - Top 3% Constraint  (Median Phred >= 15): {top_3pct:,} ({top_3pct/len(df)*100:.1f}%)")
    print(f" - Top 10% Constraint (Median Phred >= 10): {top_10pct:,} ({top_10pct/len(df)*100:.1f}%)")
    print(f" - Neutral / Weak     (Median Phred < 10):  {len(df) - top_10pct:,} ({(len(df)-top_10pct)/len(df)*100:.1f}%)")

    # Optional pLDDT filter
    filtered_df = df
    if args.min_plddt:
        filtered_df = df[df["plddt"] >= args.min_plddt]
        print(f"\nFiltered by pLDDT >= {args.min_plddt}: {len(filtered_df):,} candidates match.")

    # Sort by Median AVI Phred (Purifying Selection)
    df_sorted = filtered_df.sort_values("median_avi_phred", ascending=False).reset_index(drop=True)

    print(f"\n{'='*85}")
    print(f"TOP {min(args.top_n, len(df_sorted))} STANDOUT MICROPROTEINS UNDER HIGHEST PURIFYING SELECTION")
    print(f"{'='*85}")

    for idx, row in df_sorted.head(args.top_n).iterrows():
        med_str = phred_to_percentile_str(row["median_avi_phred"])
        max_str = phred_to_percentile_str(row["max_avi_phred"])
        plddt_badge = f"pLDDT: {row['plddt']:.1f}" if pd.notna(row['plddt']) else "pLDDT: N/A"
        gene_badge = row['nearest_gene'] if row['nearest_gene'] else "Intergenic"

        print(f"\n#{idx+1}: {row['smORF_id']} | Gene: {gene_badge} | {plddt_badge} | Length: {row['length_aa']} aa | Span: {row['span_bp']} bp")
        print(f"    Median AVI Phred: {med_str}")
        print(f"    Peak Hotspot:     {max_str} (Variant: {row['peak_variant']})")
        print(f"    Atlas Deep-Link:  {row['atlas_url']}")


if __name__ == "__main__":
    main()
