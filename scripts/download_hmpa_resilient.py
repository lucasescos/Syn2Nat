#!/usr/bin/env python3
"""
Resilient, Resume-Aware Downloader for HMPA Datasets (Zenodo Record 19547725).
Handles intermittent connection resets by repeatedly resuming with curl -C -
until the target byte size is fully reached.
"""

import os
import subprocess
import sys
import time
from pathlib import Path

DEST_DIR = Path(__file__).resolve().parent.parent / "data" / "raw" / "hmpa"
DEST_DIR.mkdir(parents=True, exist_ok=True)

TARGETS = [
    {
        "filename": "hmpa_structure_alphafold_scores.csv",
        "url": "https://zenodo.org/records/19547725/files/hmpa_structure_alphafold_scores.csv?download=1",
        "expected_bytes": 146198665,  # 139.43 MB
        "desc": "Master 617,462 microprotein sequences, coordinates, pLDDT, nearest gene",
    },
    {
        "filename": "hmpa_pathogenicity_esm2_scores.csv",
        "url": "https://zenodo.org/records/19547725/files/hmpa_pathogenicity_esm2_scores.csv?download=1",
        "expected_bytes": 14032997,  # 13.38 MB
        "desc": "ESM-2 mutational pathogenicity scores",
    },
    {
        "filename": "hmpa_structure_classification.csv",
        "url": "https://zenodo.org/records/19547725/files/hmpa_structure_classification.csv?download=1",
        "expected_bytes": 45183861,  # 43.09 MB
        "desc": "CATH structural fold classifications",
    },
    {
        "filename": "hmpa_structure_clusters_foldseek.tsv",
        "url": "https://zenodo.org/records/19547725/files/hmpa_structure_clusters_foldseek.tsv?download=1",
        "expected_bytes": 14589255,  # 13.91 MB
        "desc": "Foldseek structural cluster assignments",
    },
    {
        "filename": "hmpa_sequence_clusters_esm2_leiden.csv",
        "url": "https://zenodo.org/records/19547725/files/hmpa_sequence_clusters_esm2_leiden.csv?download=1",
        "expected_bytes": 10609440,  # 10.12 MB
        "desc": "ESM-2 sequence clusters",
    },
]


def download_resilient(target: dict):
    out_path = DEST_DIR / target["filename"]
    expected = target["expected_bytes"]
    print(f"\n=======================================================")
    print(f"Target: {target['filename']} ({target['desc']})")
    print(f"Expected Size: {expected / (1024*1024):.2f} MB")
    print(f"=======================================================")

    attempts = 0
    max_attempts = 100

    while attempts < max_attempts:
        current_size = out_path.stat().st_size if out_path.exists() else 0
        if current_size >= expected:
            print(f"[*] Complete: {target['filename']} ({current_size / (1024*1024):.2f} MB)")
            return True

        if current_size > 0:
            print(f"[*] Resuming from {current_size / (1024*1024):.2f} MB / {expected / (1024*1024):.2f} MB ({(current_size/expected)*100:.1f}%)")

        cmd = [
            "curl.exe",
            "-L",
            "-C", "-",
            "-A", "Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/120.0.0.0 Safari/537.36",
            "-e", "https://zenodo.org/records/19547725",
            "--connect-timeout", "15",
            "--max-time", "120",
            "-o", str(out_path),
            target["url"]
        ]

        try:
            res = subprocess.run(cmd, capture_output=False)
        except Exception as e:
            print(f"[!] Subprocess error: {e}")

        new_size = out_path.stat().st_size if out_path.exists() else 0
        if new_size >= expected:
            print(f"[SUCCESS] Finished {target['filename']} ({new_size / (1024*1024):.2f} MB)")
            return True

        attempts += 1
        time.sleep(2)

    print(f"[ERROR] Reached max attempts for {target['filename']}")
    return False


def main():
    print(f"Starting resilient download of HMPA Zenodo assets into {DEST_DIR}...")
    for t in TARGETS:
        download_resilient(t)
    print("\nAll downloads completed!")


if __name__ == "__main__":
    main()
