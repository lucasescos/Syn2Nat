#!/usr/bin/env python3
"""
Sequence embedding and representation helper script using Biohub ESMC models.
Queries https://biohub.ai/api/v1/logits to obtain sequence embeddings or mean hidden states.
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
                    key, val = line.split("=", 1)
                    key = key.strip()
                    val = val.strip().strip("'\"")
                    env_vars[key] = val
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
    print("Please set ESM_API_KEY before running embedding extraction.", file=sys.stderr)
    sys.exit(1)


def read_sequence(input_str: str) -> str:
    p = Path(input_str)
    if p.is_file():
        seq_lines = []
        with open(p, "r", encoding="utf-8") as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith(">"):
                    continue
                seq_lines.append(line)
        return "".join(seq_lines).upper()
    return "".join(input_str.split()).upper()


def main():
    parser = argparse.ArgumentParser(description="Extract protein representations using Biohub ESMC")
    parser.add_argument("sequence", help="Protein sequence string or FASTA file path")
    parser.add_argument("-o", "--output", default="embedding.json", help="Output filename (.json or .npy, default: embedding.json)")
    parser.add_argument("--model", default="esmc-600m-2024-12", choices=["esmc-300m-2024-12", "esmc-600m-2024-12", "esmc-6b-2024-12"], help="ESMC model variant")
    parser.add_argument("--layer", type=int, default=-2, help="Hidden layer index to extract (-2 is recommended default)")
    parser.add_argument("--per-residue", action="store_true", help="Extract per-residue embeddings instead of pooled sequence mean")
    parser.add_argument("--endpoint", default="https://biohub.ai/api/v1/logits", help="Biohub logits endpoint")

    args = parser.parse_args()
    api_key = get_api_key()
    sequence = read_sequence(args.sequence)

    if not sequence:
        print("[ERROR] Sequence is empty.", file=sys.stderr)
        sys.exit(1)

    print(f"Querying {args.model} for sequence of length {len(sequence)} aa (layer: {args.layer})...")

    # Construct request
    body = {
        "model": args.model,
        "protein": {
            "sequence": sequence
        },
        "config": {
            "sequence": True,
            "return_embeddings": args.per_residue,
            "return_mean_embedding": not args.per_residue,
            "return_hidden_states": False,
            "return_mean_hidden_states": True,
            "ith_hidden_layer": args.layer,
            "sae_config": None
        },
        "potential_sequence_of_concern": False
    }

    req_data = json.dumps(body).encode("utf-8")
    headers = {
        "Authorization": f"Bearer {api_key}",
        "Content-Type": "application/json",
        "User-Agent": "Antigravity-BiohubESM-Skill/1.0"
    }

    req = urllib.request.Request(args.endpoint, data=req_data, headers=headers, method="POST")

    try:
        with urllib.request.urlopen(req) as resp:
            result = json.loads(resp.read().decode("utf-8"))
    except urllib.error.HTTPError as e:
        err_msg = e.read().decode("utf-8", errors="replace")
        print(f"[ERROR] HTTP {e.code}: {err_msg}", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"[ERROR] Failed to query Biohub API: {e}", file=sys.stderr)
        sys.exit(1)

    # Extract embedding vector
    embeddings_data = None
    embedding_type = ""

    if args.per_residue:
        raw_emb = result.get("embeddings") or result.get("hidden_states")
        if raw_emb:
            # Slicing BOS/EOS if length is seq_len + 2
            if len(raw_emb) == len(sequence) + 2:
                embeddings_data = raw_emb[1:-1]
            else:
                embeddings_data = raw_emb
            embedding_type = f"per-residue ({len(embeddings_data)} x {len(embeddings_data[0]) if embeddings_data else 0})"
    else:
        raw_mean = result.get("mean_embedding") or result.get("mean_hidden_state")
        if raw_mean:
            embeddings_data = raw_mean
            embedding_type = f"sequence-mean (dim: {len(embeddings_data)})"

    if embeddings_data is None:
        print("[ERROR] Could not find embedding data in API response.", file=sys.stderr)
        print(json.dumps(result, indent=2)[:500])
        sys.exit(1)

    out_path = Path(args.output).resolve()
    out_path.parent.mkdir(parents=True, exist_ok=True)

    if out_path.suffix == ".npy":
        try:
            import numpy as np
            arr = np.array(embeddings_data, dtype=np.float32)
            np.save(out_path, arr)
            print(f"Saved NumPy array ({arr.shape}) to {out_path}")
        except ImportError:
            print("[WARN] numpy is not installed, saving as JSON instead.", file=sys.stderr)
            out_path = out_path.with_suffix(".json")
            with open(out_path, "w", encoding="utf-8") as f:
                json.dump({"model": args.model, "type": embedding_type, "embedding": embeddings_data}, f)
            print(f"Saved embedding data to {out_path}")
    else:
        with open(out_path, "w", encoding="utf-8") as f:
            json.dump({
                "model": args.model,
                "layer": args.layer,
                "type": embedding_type,
                "sequence_length": len(sequence),
                "embedding": embeddings_data
            }, f)
        print(f"Saved embedding JSON to {out_path}")

    print(f"Extraction successful: {embedding_type}")


if __name__ == "__main__":
    main()
