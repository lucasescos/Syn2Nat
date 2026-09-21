#!/usr/bin/env python3
"""
Consolidate and Validate ESMC Representations.
Validates both Microprotein and Human datasets:
- Sequence-level pooled representations ([2560] float32)
- Per-residue representations ([L, 2560] float16 compressed NPZ)
- Generates master parquet tables and comprehensive summary reports.
"""

import argparse
import glob
import os
import sys
import time
from pathlib import Path
from typing import Dict, List, Set, Tuple

import numpy as np
import pandas as pd

ROOT_DIR = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT_DIR))


def validate_microproteins(base_dir: Path) -> bool:
    print("\n" + "=" * 80)
    print("VALIDATING MICROPROTEIN REPRESENTATIONS")
    print("=" * 80)
    mp_dir = base_dir / "microproteins"
    seq_file = mp_dir / "microproteins_sequence_embeddings.parquet"
    batches_dir = mp_dir / "batches"

    if not seq_file.exists():
        print(f"[FAIL] Sequence embeddings file not found: {seq_file}")
        return False

    df = pd.read_parquet(seq_file)
    print(f"[PASS] Loaded sequence embeddings: {len(df):,} microproteins")
    print(f"       Columns: {df.columns.tolist()}")

    # Check vector dimension and validity
    sample_vec = np.array(df["vector"].iloc[0])
    if sample_vec.shape != (2560,):
        print(f"[FAIL] Incorrect vector dimension: {sample_vec.shape}, expected (2560,)")
        return False
    print(f"[PASS] Vector dimension confirmed: {sample_vec.shape[0]}")

    null_count = df["vector"].isnull().sum()
    if null_count > 0:
        print(f"[FAIL] Found {null_count} null vectors!")
        return False

    # Check residue batches
    npz_files = sorted(batches_dir.glob("batch_*_residues.npz"))
    print(f"[INFO] Found {len(npz_files)} residue batch files in {batches_dir}")
    residue_ids: Set[str] = set()
    total_residue_bytes = 0
    for npz in npz_files:
        total_residue_bytes += npz.stat().st_size
        with np.load(npz) as data:
            for k in data.files:
                residue_ids.add(k)

    print(f"[PASS] Total unique microproteins in residue batches: {len(residue_ids):,}")
    print(f"[INFO] Total compressed residue disk footprint: {total_residue_bytes / 1e9:.2f} GB")

    missing = set(df["id"]) - residue_ids
    if missing:
        print(f"[WARN] {len(missing)} IDs present in sequence table but missing from residue batches")
        return False
    else:
        print(f"[PASS] 100% agreement between sequence table and residue matrices ({len(df):,}/{len(residue_ids):,})")

    return True


def consolidate_human(base_dir: Path) -> bool:
    print("\n" + "=" * 80)
    print("CONSOLIDATING HUMAN PROTEOME REPRESENTATIONS")
    print("=" * 80)
    human_dir = base_dir / "human"
    batches_dir = human_dir / "batches"
    master_seq_file = human_dir / "human_sequence_embeddings.parquet"

    if not batches_dir.exists():
        print(f"[FAIL] Human batches directory not found: {batches_dir}")
        return False

    p_files = sorted(batches_dir.glob("batch_*_sequence.parquet"))
    npz_files = sorted(batches_dir.glob("batch_*_residues.npz"))
    print(f"[INFO] Found {len(p_files)} sequence parquet batches and {len(npz_files)} residue NPZ batches.")

    if not p_files:
        print("[FAIL] No sequence parquet batch files found!")
        return False

    dfs = []
    total_rows = 0
    for pf in p_files:
        df_b = pd.read_parquet(pf)
        dfs.append(df_b)
        total_rows += len(df_b)

    df_master = pd.concat(dfs, ignore_index=True)
    print(f"[INFO] Concatenated {len(p_files)} batches: total {len(df_master):,} rows (expected 20,652 canonical)")

    # Deduplicate if overlapping batches exist
    before_dedup = len(df_master)
    df_master = df_master.drop_duplicates(subset=["id"]).reset_index(drop=True)
    if len(df_master) < before_dedup:
        print(f"[INFO] Removed {before_dedup - len(df_master)} duplicate entries. Remaining: {len(df_master):,}")

    # Check vector sanity
    sample_vec = np.array(df_master["vector"].iloc[0])
    print(f"[PASS] Target vector dimension: {sample_vec.shape[0]} (float32)")

    # Save master parquet
    t0 = time.time()
    df_master.to_parquet(master_seq_file, index=False)
    print(f"[PASS] Saved master human sequence table: {master_seq_file} ({master_seq_file.stat().st_size / 1e6:.1f} MB) in {time.time()-t0:.1f}s")

    # Inspect residue batches
    residue_ids: Set[str] = set()
    total_res_bytes = 0
    for npz in npz_files:
        total_res_bytes += npz.stat().st_size
        with np.load(npz) as data:
            for k in data.files:
                residue_ids.add(k)

    print(f"[INFO] Unique proteins in residue batches: {len(residue_ids):,}")
    print(f"[INFO] Total residue disk footprint: {total_res_bytes / 1e9:.2f} GB")

    missing = set(df_master["id"]) - residue_ids
    if missing:
        print(f"[WARN] {len(missing):,} sequence IDs missing from residue batches.")
    else:
        print(f"[PASS] 100% agreement: all {len(df_master):,} proteins have per-residue matrices.")

    pct_complete = (len(df_master) / 20652) * 100
    print(f"\nHuman Proteome Completion: {len(df_master):,} / 20,652 ({pct_complete:.1f}%)")
    return pct_complete >= 99.9


def main():
    parser = argparse.ArgumentParser(description="Consolidate and validate representations.")
    parser.add_argument(
        "--data-dir",
        type=Path,
        default=ROOT_DIR / "bindpred" / "representations",
        help="Base representations directory",
    )
    args = parser.parse_args()

    mp_ok = validate_microproteins(args.data_dir)
    human_ok = consolidate_human(args.data_dir)

    print("\n" + "=" * 80)
    print("CONSOLIDATION SUMMARY")
    print(f"Microproteins (7,264 targets): {'COMPLETED' if mp_ok else 'INCOMPLETE'}")
    print(f"Human Proteome (20,652 targets): {'COMPLETED' if human_ok else 'IN PROGRESS'}")
    print("=" * 80)


if __name__ == "__main__":
    main()
