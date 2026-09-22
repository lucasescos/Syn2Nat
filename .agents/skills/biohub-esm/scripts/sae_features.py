#!/usr/bin/env python3
"""
Extract interpretable Sparse Autoencoder (SAE) features using Biohub ESMC.
Queries https://biohub.ai/api/v1/logits with sae_config to identify active biological motifs.
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
    for var in ["ESM_API_KEY", "ESM_API_KEY_ALT"]:
        key = os.environ.get(var)
        if key and key.strip():
            return key.strip()

    for p in [Path.cwd() / ".env", Path.home() / ".env"]:
        vars_map = load_env_file(p)
        for var in ["ESM_API_KEY", "ESM_API_KEY_ALT"]:
            if var in vars_map and vars_map[var]:
                return vars_map[var]

    print("[ERROR] ESM_API_KEY / ESM_API_KEY_ALT not found in environment, .env, or ~/.env", file=sys.stderr)
    print("Please set ESM_API_KEY or ESM_API_KEY_ALT before running SAE feature extraction.", file=sys.stderr)
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


def fetch_feature_dossier(feature_idx: int) -> dict | None:
    url = f"https://biohub.ai/esm/protein/api/v1alpha1/features/{feature_idx}"
    req = urllib.request.Request(url, headers={"User-Agent": "Antigravity-BiohubESM-Skill/1.0"})
    try:
        with urllib.request.urlopen(req, timeout=5) as resp:
            return json.loads(resp.read().decode("utf-8"))
    except Exception:
        return None


def main():
    parser = argparse.ArgumentParser(description="Extract and interpret ESMC Sparse Autoencoder (SAE) features")
    parser.add_argument("sequence", help="Protein sequence string or FASTA file path")
    parser.add_argument("-o", "--output", default="sae_features.json", help="Output JSON filename (default: sae_features.json)")
    parser.add_argument("--model", default="esmc-6b-2024-12", choices=["esmc-6b-2024-12", "esmc-600m-2024-12", "esmc-300m-2024-12"], help="ESMC model variant")
    parser.add_argument("--sae-id", default=None, help="Explicit SAE identifier (defaults to standard layer for selected model)")
    parser.add_argument("--topk", type=int, default=10, help="Number of top firing features to display")
    parser.add_argument("--annotate", action="store_true", help="Fetch natural language labels for top features from Atlas API")

    args = parser.parse_args()
    api_key = get_api_key()
    sequence = read_sequence(args.sequence)

    # Resolve default SAE ID if not specified
    sae_id = args.sae_id
    if not sae_id:
        if args.model == "esmc-6b-2024-12":
            sae_id = "esmc-6b-2024-12-sae-layer60-k64-codebook16384"
        elif args.model == "esmc-600m-2024-12":
            sae_id = "esmc-600m-2024-12-sae-layer27-k64-codebook16384"
        elif args.model == "esmc-300m-2024-12":
            sae_id = "esmc-300m-2024-12-sae-layer23-k64-codebook65536"

    print(f"Extracting SAE features using {args.model} -> {sae_id} (length: {len(sequence)} aa)...")

    body = {
        "model": args.model,
        "protein": {
            "sequence": sequence
        },
        "config": {
            "sequence": True,
            "return_embeddings": False,
            "return_mean_embedding": False,
            "return_hidden_states": False,
            "return_mean_hidden_states": False,
            "sae_config": {
                "models": [sae_id]
            }
        },
        "potential_sequence_of_concern": False
    }

    req_data = json.dumps(body).encode("utf-8")
    headers = {
        "Authorization": f"Bearer {api_key}",
        "Content-Type": "application/json",
        "User-Agent": "Antigravity-BiohubESM-Skill/1.0"
    }

    req = urllib.request.Request("https://biohub.ai/api/v1/logits", data=req_data, headers=headers, method="POST")

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

    sae_outputs = result.get("sae_outputs", {})
    features_payload = sae_outputs.get(sae_id)

    if not features_payload:
        print("[ERROR] SAE output not found in response.", file=sys.stderr)
        print(json.dumps(result, indent=2)[:500])
        sys.exit(1)

    # Aggregate active features across residues
    feature_max_activations = {}
    
    # Handle list of residue outputs or sparse dict
    if isinstance(features_payload, dict):
        indices = features_payload.get("indices", [])
        values = features_payload.get("values", [])
        for idx, val in zip(indices, values):
            feature_max_activations[idx] = max(feature_max_activations.get(idx, 0.0), float(val))
    elif isinstance(features_payload, list):
        # Could be per-residue list of sparse activations
        for res_entry in features_payload:
            if isinstance(res_entry, dict):
                for k, v in res_entry.items():
                    try:
                        ik = int(k)
                        feature_max_activations[ik] = max(feature_max_activations.get(ik, 0.0), float(v))
                    except ValueError:
                        pass

    sorted_features = sorted(feature_max_activations.items(), key=lambda x: x[1], reverse=True)
    top_features = sorted_features[:args.topk]

    print(f"\nTop {len(top_features)} Active SAE Features:")
    print("-" * 65)
    print(f"{'Feature Index':<15} {'Max Activation':<18} {'Biological Annotation'}")
    print("-" * 65)

    enriched_top = []
    for f_idx, act in top_features:
        label = "N/A"
        if args.annotate:
            dossier = fetch_feature_dossier(f_idx)
            if dossier:
                label = dossier.get("label") or dossier.get("summary") or "Uncharacterized"
        print(f"{f_idx:<15} {act:<18.4f} {label}")
        enriched_top.append({"feature_index": f_idx, "max_activation": act, "label": label})

    # Save output
    out_path = Path(args.output).resolve()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w", encoding="utf-8") as f:
        json.dump({
            "model": args.model,
            "sae_id": sae_id,
            "sequence_length": len(sequence),
            "top_features": enriched_top,
            "raw_sae_output": features_payload
        }, f, indent=2)

    print(f"\nFull SAE results saved to: {out_path}")


if __name__ == "__main__":
    main()
