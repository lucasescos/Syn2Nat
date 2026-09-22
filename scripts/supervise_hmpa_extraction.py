#!/usr/bin/env python3
"""
Self-Healing Supervisor Daemon for Large-Scale Microprotein Representation Extraction.
Continuously runs extraction cycles of `bindpred/extract_representations.py` until
all 617,463 human microproteins from the HMPA dataset are fully extracted.

Features:
- Subprocess cycle management (refreshes memory & OS resources every N sequences).
- Auto-restart with exponential backoff upon transient network/API drops.
- Heartbeat logging (completed count, percentage, throughput, ETA).
- Graceful termination on Ctrl+C.
"""

import argparse
import json
import os
import subprocess
import sys
import time
from datetime import datetime, timedelta
from pathlib import Path

if sys.platform == "win32":
    try:
        sys.stdout.reconfigure(encoding="utf-8", errors="replace")
        sys.stderr.reconfigure(encoding="utf-8", errors="replace")
    except Exception:
        pass

ROOT_DIR = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT_DIR))
HMPA_DIR = ROOT_DIR / "bindpred" / "representations" / "hmpa"
BATCHES_DIR = HMPA_DIR / "batches"
STATUS_FILE = HMPA_DIR / "status.json"
TOTAL_CATALOG = 617463


def parse_args():
    parser = argparse.ArgumentParser(description="Supervisor Daemon for HMPA Representation Extraction.")
    parser.add_argument("--chunk-size", type=int, default=500, help="Number of microproteins per cycle (default: 500)")
    parser.add_argument("--workers", type=int, default=16, help="Concurrent API workers (default: 16)")
    parser.add_argument("--batch-size", type=int, default=100, help="Batch size per checkpoint file (default: 100)")
    parser.add_argument("--max-retries", type=int, default=50, help="Max consecutive failure retries before pausing")
    return parser.parse_args()


def get_current_completed_count() -> int:
    """Fast count of completed targets from status.json or batch files."""
    if STATUS_FILE.exists():
        try:
            with open(STATUS_FILE, "r") as f:
                data = json.load(f)
                return data.get("total_completed", 0)
        except Exception:
            pass

    if not BATCHES_DIR.exists():
        return 0

    import pandas as pd
    count = 0
    for pf in BATCHES_DIR.glob("batch_*_sequence.parquet"):
        try:
            df = pd.read_parquet(pf, columns=["id"])
            count += len(df)
        except Exception:
            pass
    return count


def format_duration(seconds: float) -> str:
    return str(timedelta(seconds=int(seconds)))


def main():
    args = parse_args()
    from bindpred.esmc_client import get_api_keys
    configured_keys = get_api_keys()

    print("=" * 80)
    print("[*] HMPA Microprotein Representation Supervisor Daemon")
    print(f"Total Target Catalog : {TOTAL_CATALOG:,} microproteins")
    print(f"Cycle Chunk Size     : {args.chunk_size} sequences/cycle")
    print(f"Concurrency Workers  : {args.workers} threads")
    print(f"Checkpoint Batch Size: {args.batch_size} sequences/file")
    print(f"Configured API Keys  : {len(configured_keys)} key(s) detected {[f'...{k[-6:]}' for k in configured_keys]}")
    print(f"Output Directory     : {HMPA_DIR}")
    print("=" * 80)

    consecutive_failures = 0
    cycle_num = 1
    session_start_time = time.time()
    initial_completed = get_current_completed_count()
    print(f"[*] Initial completed microproteins: {initial_completed:,} ({(initial_completed/TOTAL_CATALOG)*100:.2f}%)\n")

    while True:
        completed = get_current_completed_count()
        remaining = TOTAL_CATALOG - completed

        if remaining <= 0:
            print("\n" + "=" * 80)
            print("[SUCCESS] ALL 617,463 MICROPROTEINS HAVE BEEN FULLY EXTRACTED!")
            print("=" * 80)
            break

        # Calculate session statistics
        elapsed_sec = time.time() - session_start_time
        newly_done = max(completed - initial_completed, 0)
        rate = newly_done / max(elapsed_sec, 0.001)
        eta_str = format_duration(remaining / max(rate, 0.01)) if rate > 0.01 else "calculating..."

        now_str = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        print(f"[{now_str}] Cycle {cycle_num} | Done: {completed:,}/{TOTAL_CATALOG:,} ({(completed/TOTAL_CATALOG)*100:.2f}%) | Rate: {rate:.2f} seq/s | ETA: {eta_str}")

        chunk = min(args.chunk_size, remaining)
        cmd = [
            sys.executable,
            "-u",
            str(ROOT_DIR / "bindpred" / "extract_representations.py"),
            "--dataset", "hmpa",
            "--limit", str(chunk),
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
                time.sleep(1)  # Brief pause between healthy cycles
            elif res.returncode == 3:
                # All configured API keys exhausted daily credit limit
                consecutive_failures = 0
                now_utc = datetime.utcnow()
                today_reset = now_utc.replace(hour=0, minute=2, second=0, microsecond=0)
                reset_utc = today_reset if today_reset > now_utc else (today_reset + timedelta(days=1))
                wait_sec = max(int((reset_utc - now_utc).total_seconds()), 300)
                print("\n" + "!" * 80)
                print(f"[PAUSE] All Biohub API keys have exhausted daily credit limits.")
                print(f"[PAUSE] Quota refreshes daily at 00:00 UTC.")
                print(f"[PAUSE] Sleeping for {format_duration(wait_sec)} until {reset_utc.strftime('%Y-%m-%d %H:%M:%S UTC')}...")
                print("!" * 80 + "\n")
                time.sleep(wait_sec)
                cycle_num += 1
            else:
                consecutive_failures += 1
                backoff = min(15 * (2 ** (consecutive_failures - 1)), 300)
                print(f"[!] Cycle {cycle_num} exited with code {res.returncode}. Retry {consecutive_failures}/{args.max_retries} in {backoff}s...")
                time.sleep(backoff)

        except KeyboardInterrupt:
            print("\n[!] Supervisor interrupted by user. Graceful shutdown.")
            break
        except Exception as e:
            consecutive_failures += 1
            backoff = min(15 * (2 ** (consecutive_failures - 1)), 300)
            print(f"[!] Unexpected error in cycle {cycle_num}: {e}. Backing off for {backoff}s...")
            time.sleep(backoff)

        if consecutive_failures >= args.max_retries:
            print(f"\n[CRITICAL] Reached {args.max_retries} consecutive failures. Pausing supervisor.")
            break


if __name__ == "__main__":
    main()
