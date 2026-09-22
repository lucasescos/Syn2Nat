#!/usr/bin/env python3
"""
Filter HMPA Microproteins by Genomic Autonomy (Remove Canonical CDS intORFs).

Parses canonical CDS intervals from GENCODE v46 (hg38), performs fast interval
intersection against all 617,262 HMPA smORFs, and partitions the catalog into:
- Autonomous smORFs (lncRNA-ORFs, 5' UTR uORFs, 3' UTR dORFs, intergenic)
- Canonical CDS-overlapping intORFs (hitchhiking confounder)

Saves curated autonomous candidates to Parquet.
"""

import gzip
import re
import sys
import time
from bisect import bisect_left, bisect_right
from pathlib import Path
import numpy as np
import pandas as pd

BASE_DIR = Path(__file__).resolve().parent.parent
DATA_DIR = BASE_DIR / "data"
GTF_FILE = DATA_DIR / "raw" / "gencode.v46.basic.annotation.gtf.gz"
HMPA_PARQUET = DATA_DIR / "curated" / "hmpa_microprotein_intervals_all_617k.parquet"
OUTPUT_AUTONOMOUS = DATA_DIR / "curated" / "hmpa_autonomous_intervals.parquet"
OUTPUT_INTORFS = DATA_DIR / "curated" / "hmpa_intorfs_intervals.parquet"


def parse_canonical_cds(gtf_path: Path):
    """
    Parse all CDS features from GENCODE GTF.
    Returns dict: {chrom: [(start, end, gene_name), ...]}
    """
    print(f"Parsing canonical CDS intervals from {gtf_path}...")
    cds_by_chrom = {}
    gene_re = re.compile(r'gene_name "([^"]+)"')

    count = 0
    with gzip.open(gtf_path, "rt", encoding="utf-8") as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9 or parts[2] != "CDS":
                continue

            chrom = parts[0]
            start = int(parts[3])
            end = int(parts[4])
            m = gene_re.search(parts[8])
            gene = m.group(1) if m else ""

            if chrom not in cds_by_chrom:
                cds_by_chrom[chrom] = []
            cds_by_chrom[chrom].append((start, end, gene))
            count += 1

    print(f"Loaded {count:,} raw CDS segments across {len(cds_by_chrom)} chromosomes.")

    # Merge overlapping CDS intervals per chromosome for ultra-fast lookup
    merged_cds = {}
    for chrom, ivs in cds_by_chrom.items():
        # Sort by start coordinate
        ivs.sort(key=lambda x: x[0])
        merged = []
        for st, en, gene in ivs:
            if not merged:
                merged.append([st, en, {gene}])
            else:
                prev_st, prev_en, prev_genes = merged[-1]
                if st <= prev_en + 1:  # Overlapping or adjacent
                    merged[-1][1] = max(prev_en, en)
                    merged[-1][2].add(gene)
                else:
                    merged.append([st, en, {gene}])
        
        # Split into numpy arrays for bisect search
        starts = np.array([m[0] for m in merged], dtype=np.int64)
        ends = np.array([m[1] for m in merged], dtype=np.int64)
        genes = [",".join(sorted(m[2])) for m in merged]
        merged_cds[chrom] = (starts, ends, genes)

    total_merged = sum(len(v[0]) for v in merged_cds.values())
    print(f"Collapsed into {total_merged:,} non-redundant canonical CDS intervals.")
    return merged_cds


def check_cds_overlap(merged_cds, chrom: str, st: int, en: int):
    """
    Check if interval [st, en] overlaps any merged CDS on chrom.
    Returns (is_intorf: bool, overlapping_host_genes: str)
    """
    if chrom not in merged_cds:
        return False, ""

    starts, ends, genes = merged_cds[chrom]
    # Binary search for any interval that could overlap [st, en]
    # An interval [c_st, c_en] overlaps [st, en] iff c_st <= en and c_en >= st
    # Find rightmost interval with c_st <= en
    idx_right = bisect_right(starts, en)
    if idx_right == 0:
        return False, ""

    # Check candidates up to idx_right
    overlapping_genes = set()
    is_overlapping = False
    # Only need to scan backwards while ends could be >= st
    for i in range(idx_right - 1, -1, -1):
        c_st = starts[i]
        c_en = ends[i]
        if c_en < st:
            # Since intervals are sorted by start, earlier ones might still have large ends,
            # but in merged non-overlapping intervals, c_en is monotonically increasing!
            break
        if c_st <= en and c_en >= st:
            is_overlapping = True
            overlapping_genes.add(genes[i])

    if is_overlapping:
        return True, ";".join(sorted(overlapping_genes))
    return False, ""


def main():
    print("=" * 80)
    print("HMPA Catalog De-confounding: Removing Canonical CDS intORFs")
    print("=" * 80)

    t0 = time.time()
    merged_cds = parse_canonical_cds(GTF_FILE)

    print(f"\nLoading curated HMPA intervals from {HMPA_PARQUET}...")
    df = pd.read_parquet(HMPA_PARQUET)
    print(f"Loaded {len(df):,} microproteins in {time.time() - t0:.1f}s.")

    print("\nScanning for CDS overlap across all 617,262 microproteins...")
    is_intorf_list = []
    host_genes_list = []

    chroms = df["chrom"].values
    starts = df["start_1b"].values
    ends = df["end_1b"].values

    for i in range(len(df)):
        is_ov, genes = check_cds_overlap(merged_cds, chroms[i], starts[i], ends[i])
        is_intorf_list.append(is_ov)
        host_genes_list.append(genes)

    df["is_intorf"] = is_intorf_list
    df["cds_host_genes"] = host_genes_list
    df["is_autonomous"] = ~df["is_intorf"]

    intorf_count = df["is_intorf"].sum()
    auto_count = df["is_autonomous"].sum()

    print("\n" + "=" * 80)
    print(f"Genomic Classification Results:")
    print(f" - Total Microproteins Analyzed : {len(df):,}")
    print(f" - Canonical intORFs (Filtered) : {intorf_count:,} ({intorf_count / len(df) * 100:.1f}%)")
    print(f" - Autonomous Microproteins     : {auto_count:,} ({auto_count / len(df) * 100:.1f}%)")
    print("=" * 80)

    # Save Autonomous dataset
    df_auto = df[df["is_autonomous"]].reset_index(drop=True)
    df_auto.to_parquet(OUTPUT_AUTONOMOUS, index=False, compression="snappy")
    print(f"\nSaved {len(df_auto):,} Autonomous microproteins to: {OUTPUT_AUTONOMOUS}")

    # Save intORFs dataset for records
    df_int = df[df["is_intorf"]].reset_index(drop=True)
    df_int.to_parquet(OUTPUT_INTORFS, index=False, compression="snappy")
    print(f"Saved {len(df_int):,} intORFs to: {OUTPUT_INTORFS}")

    print(f"\nCompleted in {time.time() - t0:.1f}s.")


if __name__ == "__main__":
    main()
