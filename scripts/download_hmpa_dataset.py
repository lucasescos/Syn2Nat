#!/usr/bin/env python3
"""
Download core data files from the Human Microprotein Atlas (HMPA) (Communications Chemistry, 2026).
Downloads master sequences, coordinates, annotations, and classifications for 617,462 human microproteins.
"""

import os
import sys
import time
import urllib.request
from pathlib import Path

DEST_DIR = Path(__file__).resolve().parent.parent / "data" / "raw" / "hmpa"
DEST_DIR.mkdir(parents=True, exist_ok=True)

FILES = [
    {
        "name": "smORFs.pro.fa",
        "url": "http://www.cuilab.cn/microaf/static/mipro_info/smORFs.pro.fa",
        "expected_size_mb": 26.83,
        "desc": "FASTA format protein sequences for all 617,462 human microproteins",
    },
    {
        "name": "smORFs.merge.hg38.bed",
        "url": "http://www.cuilab.cn/microaf/static/mipro_info/smORFs.merge.hg38.bed",
        "expected_size_mb": 25.0,
        "desc": "Genomic coordinates in hg38 BED format",
    },
    {
        "name": "micro_complete_info.csv",
        "url": "http://www.cuilab.cn/microaf/static/mipro_info/micro_complete_info.csv",
        "expected_size_mb": 324.99,
        "desc": "Master table with 617,462 microproteins, sequences, AF2 pLDDT, properties",
    },
    {
        "name": "data_source.tsv",
        "url": "http://www.cuilab.cn/microaf/static/mipro_info/data_source.tsv",
        "expected_size_mb": 81.76,
        "desc": "Source ID cross-reference (SmProt, sORFs.org, Thomas et al.)",
    },
    {
        "name": "microprotein_subcellular_location.csv",
        "url": "http://www.cuilab.cn/microaf/static/subcellular_location/microprotein_subcellular_location.csv",
        "expected_size_mb": 52.78,
        "desc": "LocPro subcellular localization probabilities across 10 compartments",
    },
    {
        "name": "hmpa_pathogenicity_esm2_scores.csv",
        "url": "https://zenodo.org/api/records/19547725/files/hmpa_pathogenicity_esm2_scores.csv/content",
        "expected_size_mb": 13.38,
        "desc": "ESM-2 mutational pathogenicity scores (Zenodo)",
    },
    {
        "name": "hmpa_structure_classification.csv",
        "url": "https://zenodo.org/api/records/19547725/files/hmpa_structure_classification.csv/content",
        "expected_size_mb": 43.09,
        "desc": "CATH fold classifications (Mainly alpha, beta, etc.)",
    },
    {
        "name": "hmpa_structure_clusters_foldseek.tsv",
        "url": "https://zenodo.org/api/records/19547725/files/hmpa_structure_clusters_foldseek.tsv/content",
        "expected_size_mb": 13.91,
        "desc": "Foldseek structural cluster assignments",
    },
]


def download_file(item: dict):
    out_path = DEST_DIR / item["name"]
    if out_path.exists() and out_path.stat().st_size > 1024 * 1024:
        print(f"[SKIP] {item['name']} already exists ({out_path.stat().st_size / (1024*1024):.2f} MB).")
        return

    print(f"\n[DOWNLOADING] {item['name']} ({item['desc']})")
    print(f"  URL: {item['url']}")
    print(f"  Target: {out_path}")

    headers = {"User-Agent": "Mozilla/5.0 (Windows NT 10.0; Win64; x64)"}
    req = urllib.request.Request(item["url"], headers=headers)

    tmp_path = out_path.with_suffix(out_path.suffix + ".part")
    start_time = time.time()

    try:
        with urllib.request.urlopen(req, timeout=60) as resp, open(tmp_path, "wb") as f:
            total_bytes = int(resp.headers.get("Content-Length", 0))
            downloaded = 0
            last_print = 0

            while True:
                chunk = resp.read(1024 * 512)  # 512 KB chunks
                if not chunk:
                    break
                f.write(chunk)
                downloaded += len(chunk)

                if time.time() - last_print > 2 or downloaded == total_bytes:
                    last_print = time.time()
                    elapsed = last_print - start_time
                    speed = (downloaded / (1024 * 1024)) / max(elapsed, 0.001)
                    if total_bytes > 0:
                        pct = (downloaded / total_bytes) * 100
                        print(f"  Progress: {pct:5.1f}% ({downloaded / (1024*1024):.1f}/{total_bytes / (1024*1024):.1f} MB) at {speed:.2f} MB/s", end="\r")
                    else:
                        print(f"  Downloaded: {downloaded / (1024*1024):.1f} MB at {speed:.2f} MB/s", end="\r")

        print()
        tmp_path.replace(out_path)
        print(f"[SUCCESS] {item['name']} saved ({out_path.stat().st_size / (1024*1024):.2f} MB) in {time.time() - start_time:.1f}s")
    except Exception as e:
        print(f"\n[ERROR] Failed to download {item['name']}: {e}")
        if tmp_path.exists():
            tmp_path.unlink()


def main():
    print(f"Starting HMPA dataset download to {DEST_DIR}...")
    for item in FILES:
        download_file(item)
    print("\nAll downloads finished. Checking contents:")
    for f in sorted(DEST_DIR.iterdir()):
        print(f"  - {f.name}: {f.stat().st_size / (1024*1024):.2f} MB")


if __name__ == "__main__":
    main()
