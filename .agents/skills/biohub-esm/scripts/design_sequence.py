#!/usr/bin/env python3
"""
Generative protein design & motif scaffolding script using Biohub ESM3.
Queries POST /api/v1/generate with iterative unmasking across tracks.
"""

import argparse
import json
import os
import sys
import urllib.error
import urllib.request
from pathlib import Path


def load_env_file(filepath: Path) -> dict:
    env_vars = {}
    if not filepath.is_file():
        return env_vars
    try:
        with open(filepath, "r", encoding="utf-8") as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith("#"):
                    continue
                if "=" in line:
                    k, v = line.split("=", 1)
                    env_vars[k.strip()] = v.strip().strip("'\"")
    except Exception:
        pass
    return env_vars


def get_api_key() -> str:
    key = os.environ.get("ESM_API_KEY")
    if key and key.strip():
        return key.strip()
    for p in [Path.cwd() / ".env", Path.home() / ".env"]:
        vars_map = load_env_file(p)
        if "ESM_API_KEY" in vars_map and vars_map["ESM_API_KEY"]:
            return vars_map["ESM_API_KEY"]
    print("[ERROR] ESM_API_KEY not found in environment, .env, or ~/.env", file=sys.stderr)
    sys.exit(1)


def main():
    parser = argparse.ArgumentParser(description="Generative sequence design and motif scaffolding with ESM3")
    parser.add_argument("prompt_sequence", help="Seed sequence pattern. Use '_' for unmasked residues to be designed (e.g. 'M___G__L___K')")
    parser.add_argument("-o", "--output", default="generated_design.fasta", help="Output FASTA filename")
    parser.add_argument("--model", default="esm3-open-2024-03", help="ESM3 model variant")
    parser.add_argument("--steps", type=int, default=20, help="Number of unmasking generation steps (default: 20)")
    parser.add_argument("--temperature", type=float, default=0.5, help="Sampling temperature (default: 0.5)")
    parser.add_argument("--top-p", type=float, default=1.0, help="Top-p nucleus sampling (default: 1.0)")
    parser.add_argument("--schedule", default="cosine", choices=["cosine", "linear"], help="Unmasking schedule")
    parser.add_argument("--strategy", default="random", choices=["random", "entropy"], help="Unmasking strategy")

    args = parser.parse_args()
    api_key = get_api_key()
    prompt = "".join(args.prompt_sequence.split()).upper()

    num_masked = prompt.count("_")
    print(f"Generating sequence with ESM3 ({len(prompt)} aa total, {num_masked} masked positions)...")
    print(f"Prompt: {prompt}")

    body = {
        "model": args.model,
        "track": "sequence",
        "inputs": {
            "sequence": prompt
        },
        "num_steps": min(args.steps, max(1, num_masked)),
        "temperature": args.temperature,
        "temperature_annealing": True,
        "top_p": args.top_p,
        "schedule": args.schedule,
        "strategy": args.strategy,
        "condition_on_coordinates_only": False,
        "potential_sequence_of_concern": False
    }

    req_data = json.dumps(body).encode("utf-8")
    headers = {
        "Authorization": f"Bearer {api_key}",
        "Content-Type": "application/json",
        "User-Agent": "Antigravity-BiohubESM-Skill/1.0"
    }

    req = urllib.request.Request("https://biohub.ai/api/v1/generate", data=req_data, headers=headers, method="POST")

    try:
        with urllib.request.urlopen(req) as resp:
            res = json.loads(resp.read().decode("utf-8"))
    except urllib.error.HTTPError as e:
        err = e.read().decode("utf-8", errors="replace")
        print(f"[ERROR] HTTP {e.code}: {err}", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"[ERROR] Generation failed: {e}", file=sys.stderr)
        sys.exit(1)

    tracks = res.get("tracks") or res.get("inputs") or res
    gen_seq = tracks.get("sequence") if isinstance(tracks, dict) else res.get("sequence")

    if not gen_seq:
        print("[ERROR] Response did not contain generated sequence.", file=sys.stderr)
        print(json.dumps(res, indent=2)[:500])
        sys.exit(1)

    print("\n--- Generation Completed ---")
    print(f"Designed Sequence: {gen_seq}")
    print(f"Length:            {len(gen_seq)} aa")

    out_path = Path(args.output).resolve()
    with open(out_path, "w", encoding="utf-8") as f:
        f.write(f">designed_esm3_temp_{args.temperature}\n{gen_seq}\n")
    print(f"Saved design to:   {out_path}")


if __name__ == "__main__":
    main()
