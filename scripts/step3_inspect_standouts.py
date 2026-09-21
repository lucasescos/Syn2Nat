#!/usr/bin/env python3
"""
Step 3: Inspect Standout Microproteins from AlphaGenome Screening.
Aggregates completed batches, identifies top outliers under negative selection,
and outputs formatted tables with clickable AlphaGenome Atlas URLs.
"""

import argparse
from pathlib import Path
import pandas as pd
import numpy as np

BASE_DIR = Path(__file__).resolve().parent.parent
DATA_DIR = BASE_DIR / "data"
BATCHES_DIR = DATA_DIR / "screening" / "avi_batches"
CURATED_ALL = DATA_DIR / "curated" / "microprotein_intervals_all_7264.parquet"

def phred_to_percentile_str(phred):
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
    parser = argparse.ArgumentParser(description="Inspect standout microproteins.")
    parser.add_argument("--top-n", type=int, default=20, help="Number of top candidates to display.")
    parser.add_argument("--min-phred", type=float, default=15.0, help="Minimum median AVI Phred threshold.")
    args = parser.parse_args()

    batch_files = sorted(BATCHES_DIR.glob("batch_*.parquet"))
    if not batch_files:
        print(f"No batch files found in {BATCHES_DIR}")
        return

    dfs = [pd.read_parquet(f) for f in batch_files]
    df = pd.concat(dfs, ignore_index=True).drop_duplicates(subset=["orf_id"])

    total_catalog = 7264
    pct_done = len(df) / total_catalog * 100

    print("=" * 80)
    print(f"AlphaGenome Microprotein Screening Progress: {len(df):,} / {total_catalog:,} ({pct_done:.1f}%)")
    print("=" * 80)

    # Constraint summary
    top_1pct = (df["median_avi_phred"] >= 20.0).sum()
    top_3pct = (df["median_avi_phred"] >= 15.0).sum()
    top_10pct = (df["median_avi_phred"] >= 10.0).sum()

    print(f"\nConstraint Distribution Across Scored Candidates:")
    print(f" - Top 1% Constraint  (Median Phred >= 20): {top_1pct} ({top_1pct/len(df)*100:.1f}%)")
    print(f" - Top 3% Constraint  (Median Phred >= 15): {top_3pct} ({top_3pct/len(df)*100:.1f}%)")
    print(f" - Top 10% Constraint (Median Phred >= 10): {top_10pct} ({top_10pct/len(df)*100:.1f}%)")
    print(f" - Neutral / Weak     (Median Phred < 10):  {len(df) - top_10pct} ({(len(df)-top_10pct)/len(df)*100:.1f}%)")

    # Ranking by Median AVI Phred (Negative Selection)
    df_sorted = df.sort_values("median_avi_phred", ascending=False).reset_index(drop=True)

    print(f"\n{'='*80}")
    print(f"TOP {min(args.top_n, len(df_sorted))} STANDOUT MICROPROTEINS UNDER HIGHEST NEGATIVE SELECTION")
    print(f"{'='*80}")

    cols_show = ["orf_id", "gene", "final_tier", "orf_biotype", "length_aa", "median_avi_phred", "max_avi_phred", "atlas_url"]
    
    for idx, row in df_sorted.head(args.top_n).iterrows():
        tier_badge = f"Tier {row['final_tier']}"
        med_str = phred_to_percentile_str(row['median_avi_phred'])
        max_str = phred_to_percentile_str(row['max_avi_phred'])
        print(f"\n#{idx+1}: {row['orf_id']} | Gene: {row['gene']} | {tier_badge} | Biotype: {row['orf_biotype']} | Length: {row['length_aa']} aa")
        print(f"    Median AVI Phred: {med_str}")
        print(f"    Peak Hotspot:     {max_str} (Variant: {row['peak_variant']})")
        print(f"    Atlas Deep-Link:  {row['atlas_url']}")

if __name__ == "__main__":
    main()
