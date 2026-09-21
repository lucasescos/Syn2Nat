#!/usr/bin/env python3
"""
Step 1: Rapid Coordinate Extraction for all 7,264 GENCODE non-canonical ORFs.
Parses hg38 exonic intervals from supp_table11_orbl.xlsx, merges with tier/biotype metadata,
and standardizes 0-based half-open intervals for AlphaGenome API querying.
"""

import os
import re
from pathlib import Path
import pandas as pd

BASE_DIR = Path(__file__).resolve().parent.parent
DATA_DIR = BASE_DIR / "data"
RAW_FILE = DATA_DIR / "raw" / "supp_table11_orbl.xlsx"
PARSED_ORFS = DATA_DIR / "parsed" / "orfs.parquet"
PARSED_TIERS = DATA_DIR / "parsed" / "orf_tiers.parquet"
CURATED_DIR = DATA_DIR / "curated"

def parse_intervals(interval_str):
    """
    Parse '+' separated intervals like 'chr1:847672-847806+chr1:850181-850321'
    Returns:
      chrom: str
      exons_1based: list of (start, end) tuples (1-based closed)
      exons_0based: list of (start-1, end) tuples (0-based half-open)
      total_bp: int
    """
    parts = interval_str.strip().split('+')
    exons_1based = []
    exons_0based = []
    chrom = None
    total_bp = 0
    
    for p in parts:
        c, coords = p.split(':')
        if chrom is None:
            chrom = c
        st, en = map(int, coords.split('-'))
        exons_1based.append((st, en))
        exons_0based.append((st - 1, en))
        total_bp += (en - st + 1)
        
    return chrom, exons_1based, exons_0based, total_bp

def main():
    print("=" * 70)
    print("Step 1: Extracting & Standardizing Microprotein Genomic Intervals")
    print("=" * 70)

    CURATED_DIR.mkdir(parents=True, exist_ok=True)

    print(f"Reading raw supplementary Table 11: {RAW_FILE}...")
    df_raw = pd.read_excel(RAW_FILE, skiprows=5)
    print(f"Loaded {len(df_raw)} raw rows from Table 11.")

    print(f"Reading parsed tier annotations: {PARSED_TIERS}...")
    df_tiers = pd.read_parquet(PARSED_TIERS)

    # Merge on orf_id == Ribo-Seq_ORF
    merged = pd.merge(
        df_raw,
        df_tiers[['orf_id', 'final_tier', 'is_peptidein', 'hla_detected', 'nonhla_detected', 'length']],
        left_on='Ribo-Seq_ORF',
        right_on='orf_id',
        how='inner'
    )
    print(f"Successfully matched {len(merged)} ORFs with tier annotations.")

    records = []
    for _, row in merged.iterrows():
        raw_iv = str(row['IntervalsWithStop'])
        chrom, exons_1b, exons_0b, total_bp = parse_intervals(raw_iv)
        
        # Format 0-based intervals string: e.g. "chr1:847671-847806;chr1:850180-850321"
        intervals_0based_str = ";".join([f"{chrom}:{st}-{en}" for st, en in exons_0b])
        
        records.append({
            'orf_id': row['orf_id'],
            'gene': str(row['GeneNameV42']),
            'transcript': str(row['TranscriptV42']),
            'orf_biotype': str(row['BiotypeV42']),
            'final_tier': str(row['final_tier']),
            'is_peptidein': bool(row['is_peptidein']),
            'strand': str(row['Strand']),
            'chrom': chrom,
            'n_exons': len(exons_0b),
            'length_aa': int(row['length']),
            'length_bp': total_bp,
            'intervals_1based': raw_iv,
            'intervals_0based': intervals_0based_str,
            'min_start_0based': exons_0b[0][0],
            'max_end_0based': exons_0b[-1][1],
            'span_bp': exons_0b[-1][1] - exons_0b[0][0]
        })

    df_out = pd.DataFrame(records)
    out_parquet = CURATED_DIR / "microprotein_intervals_all_7264.parquet"
    out_tsv = CURATED_DIR / "microprotein_intervals_all_7264.tsv"

    df_out.to_parquet(out_parquet, index=False)
    df_out.to_csv(out_tsv, sep="\t", index=False)

    print(f"\nExtracted and saved {len(df_out)} standardized microprotein intervals:")
    print(f" - Parquet: {out_parquet} ({out_parquet.stat().st_size / 1024:.1f} KB)")
    print(f" - TSV:     {out_tsv} ({out_tsv.stat().st_size / 1024:.1f} KB)")

    print("\nSummary by Final Tier:")
    print(df_out['final_tier'].value_counts())
    print("\nSummary by Biotype:")
    print(df_out['orf_biotype'].value_counts())

if __name__ == "__main__":
    main()
