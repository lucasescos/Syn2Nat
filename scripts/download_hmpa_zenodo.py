#!/usr/bin/env python3
"""
Download curated HMPA 600k+ microprotein datasets from Zenodo (Record 19547725).
Uses curl with resume and browser User-Agent headers to ensure reliable downloading.
"""

import subprocess
import sys
import time
from pathlib import Path

DEST_DIR = Path(__file__).resolve().parent.parent / "data" / "raw" / "hmpa"
DEST_DIR.mkdir(parents=True, exist_ok=True)

ZENODO_FILES = [
    {
        "filename": "hmpa_structure_alphafold_scores.csv",
        "url": "https://zenodo.org/records/19547725/files/hmpa_structure_alphafold_scores.csv?download=1",
        "desc": "Master table with 617,462 microprotein sequences, coordinates, pLDDT, nearest gene",
        "approx_mb": 139.4,
    },
    {
        "filename": "hmpa_pathogenicity_esm2_scores.csv",
        "url": "https://zenodo.org/records/19547725/files/hmpa_pathogenicity_esm2_scores.csv?download=1",
        "desc": "ESM-2 mutational pathogenicity scores",
        "approx_mb": 13.4,
    },
    {
        "filename": "hmpa_structure_classification.csv",
        "url": "https://zenodo.org/records/19547725/files/hmpa_structure_classification.csv?download=1",
        "desc": "CATH fold classifications (Mainly alpha, beta, etc.)",
        "approx_mb": 43.1,
    },
    {
        "filename": "hmpa_structure_clusters_foldseek.tsv",
        "url": "https://zenodo.org/records/19547725/files/hmpa_structure_clusters_foldseek.tsv?download=1",
        "desc": "Foldseek structural cluster assignments",
        "approx_mb": 13.9,
    },
    {
        "filename": "hmpa_sequence_clusters_esm2_leiden.csv",
        "url": "https://zenodo.org/records/19547725/files/hmpa_sequence_clusters_esm2_leiden.csv?download=1",
        "desc": "ESM-2 sequence clusters (17 Leiden clusters)",
        "approx_mb": 10.1,
    },
]


def download_with_curl(item: dict):
    out_file = DEST_DIR / item["filename"]
    print(f"\n[DOWNLOAD] {item['filename']} (~{item['approx_mb']} MB)")
    print(f"  Description: {item['desc']}")
    print(f"  Destination: {out_file}")

    if out_file.exists() and out_file.stat().st_size > (item["approx_mb"] * 0.9 * 1024 * 1024):
        print(f"  [OK] Already downloaded and complete ({out_file.stat().st_size / (1024*1024):.2f} MB). Skipping.")
        return

    cmd = [
        "curl.exe",
        "-L",
        "-C", "-",
        "-A", "Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/120.0.0.0 Safari/537.36",
        "-e", "https://zenodo.org/records/19547725",
        "--retry", "5",
        "--retry-delay", "3",
        "-o", str(out_file),
        item["url"]
    ]

    t0 = time.time()
    res = subprocess.run(cmd)
    if res.returncode == 0:
        elapsed = time.time() - t0
        mb = out_file.stat().st_size / (1024 * 1024)
        print(f"  [SUCCESS] {item['filename']} finished: {mb:.2f} MB in {elapsed:.1f}s")
    else:
        print(f"  [ERROR] Failed to download {item['filename']}, returncode: {res.returncode}")


def main():
    print(f"=== Downloading HMPA Microprotein Atlas Data from Zenodo ===")
    print(f"Target Directory: {DEST_DIR}")
    for item in ZENODO_FILES:
        download_with_curl(item)

    print("\n=== Download Summary ===")
    for f in sorted(DEST_DIR.iterdir()):
        if not f.name.endswith(".part"):
            print(f"  {f.name:<45} : {f.stat().st_size / (1024*1024):8.2f} MB")


if __name__ == "__main__":
    main()
