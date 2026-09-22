#!/usr/bin/env python3
"""
Supervisor Daemon for Autonomous HMPA Microprotein AlphaGenome Screening.
Continuously runs cycles of `query_autonomous_avi.py` across the 93,929 genuine autonomous
microproteins (non-intORFs), prioritized by AlphaFold structure confidence (pLDDT).
"""

import argparse
import datetime
import json
import os
import subprocess
import sys
import time
from datetime import timedelta
from pathlib import Path
import pandas as pd

if sys.platform == "win32":
    try:
        sys.stdout.reconfigure(encoding="utf-8", errors="replace")
        sys.stderr.reconfigure(encoding="utf-8", errors="replace")
    except Exception:
        pass

ROOT_DIR = Path(__file__).resolve().parent.parent
SCREENING_DIR = ROOT_DIR / "data" / "screening"
BATCHES_DIR = SCREENING_DIR / "hmpa_autonomous_avi_batches"
STATUS_FILE = SCREENING_DIR / "hmpa_autonomous_avi_status.json"
TOTAL_CATALOG = 93929


def parse_args():
    parser = argparse.ArgumentParser(description="Supervisor Daemon for Autonomous HMPA Screening.")
    parser.add_argument("--chunk-size", type=int, default=500, help="Microproteins per cycle (default: 500)")
    parser.add_argument("--workers", type=int, default=3, help="Concurrent API workers (default: 3)")
    parser.add_argument("--batch-size", type=int, default=250, help="Checkpoint batch size (default: 250)")
    parser.add_argument("--max-retries", type=int, default=50, help="Max consecutive failures")
    return parser.parse_args()


def get_current_completed_count() -> int:
    """Fast count of completed targets from batch files."""
    if not BATCHES_DIR.exists():
        return 0

    count = 0
    for pf in BATCHES_DIR.glob("batch_*.parquet"):
        try:
            df = pd.read_parquet(pf, columns=["smORF_id", "status"])
            count += len(df[df["status"] == "success"])
        except Exception:
            pass
    return count


def format_duration(seconds: float) -> str:
    return str(timedelta(seconds=int(seconds)))


def main():
    args = parse_args()
    print("=" * 80)
    print("[*] Autonomous HMPA AlphaGenome Screening Supervisor Daemon")
    print(f"Total Target Catalog : {TOTAL_CATALOG:,} autonomous microproteins (non-intORFs)")
    print(f"Cycle Chunk Size     : {args.chunk_size} ORFs/cycle")
    print(f"Concurrency Workers  : {args.workers} threads (quota-balanced)")
    print(f"Checkpoint Batch Size: {args.batch_size} ORFs/batch")
    print(f"Batches Directory    : {BATCHES_DIR}")
    print("=" * 80)

    consecutive_failures = 0
    cycle_num = 1
    session_start_time = time.time()
    initial_completed = get_current_completed_count()
    print(f"[*] Initial completed autonomous microproteins: {initial_completed:,} ({(initial_completed/TOTAL_CATALOG)*100:.2f}%)\n")

    while True:
        completed = get_current_completed_count()
        remaining = TOTAL_CATALOG - completed

        if remaining <= 0:
            print("\n" + "=" * 80)
            print("[SUCCESS] ALL 93,929 AUTONOMOUS MICROPROTEINS HAVE BEEN FULLY SCORED!")
            print("=" * 80)
            break

        elapsed_sec = time.time() - session_start_time
        newly_done = max(completed - initial_completed, 0)
        rate = newly_done / max(elapsed_sec, 0.001)
        eta_str = format_duration(remaining / max(rate, 0.01)) if rate > 0.01 else "calculating..."

        now_str = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        print(f"[{now_str}] Cycle {cycle_num} | Done: {completed:,}/{TOTAL_CATALOG:,} ({(completed/TOTAL_CATALOG)*100:.2f}%) | Rate: {rate:.2f} orf/s | ETA: {eta_str}")

        chunk = min(args.chunk_size, remaining)
        cmd = [
            "uv",
            "run",
            "--with",
            "alphagenome,python-dotenv,pandas,pyarrow,tqdm,numpy",
            "python",
            str(ROOT_DIR / "scripts" / "query_autonomous_avi.py"),
            "--max-orfs", str(chunk),
            "--batch-size", str(args.batch_size),
            "--workers", str(args.workers),
        ]

        t0 = time.time()
        try:
            res = subprocess.run(cmd, cwd=str(ROOT_DIR))
            cycle_duration = time.time() - t0

            if res.returncode == 0:
                consecutive_failures = 0
                cycle_num += 1
                time.sleep(2)  # Healthy pause between cycles
            else:
                consecutive_failures += 1
                backoff = min(20 * (2 ** (consecutive_failures - 1)), 300)
                print(f"[!] Cycle {cycle_num} exited with code {res.returncode}. Retry {consecutive_failures}/{args.max_retries} in {backoff}s...")
                time.sleep(backoff)

        except KeyboardInterrupt:
            print("\n[!] Supervisor interrupted by user. Graceful shutdown.")
            break
        except Exception as e:
            consecutive_failures += 1
            backoff = min(20 * (2 ** (consecutive_failures - 1)), 300)
            print(f"[!] Cycle {cycle_num} failed with error: {e}. Retry in {backoff}s...")
            time.sleep(backoff)


if __name__ == "__main__":
    main()
