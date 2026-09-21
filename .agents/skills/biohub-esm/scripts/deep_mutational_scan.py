#!/usr/bin/env python3
"""
In-silico deep mutational scanning (DMS) and load-bearing residue analysis using Biohub ESMC.
Queries POST /api/v1/logits to obtain sequence logits, computes per-position Shannon entropy,
and identifies mutationally constrained / load-bearing positions.
"""

import argparse
import json
import math
import os
import sys
import urllib.error
import urllib.request
from pathlib import Path

STANDARD_AAS = "ACDEFGHIKLMNPQRSTVWY"


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


def softmax(logits: list[float]) -> list[float]:
    max_l = max(logits)
    exp_l = [math.exp(x - max_l) for x in logits]
    sum_e = sum(exp_l)
    return [x / sum_e for x in exp_l]


def main():
    parser = argparse.ArgumentParser(description="In-silico mutational scan and entropy analysis using ESMC")
    parser.add_argument("sequence", help="Protein sequence string or FASTA file")
    parser.add_argument("-o", "--output", default="mutational_scan.tsv", help="Output TSV filename")
    parser.add_argument("--model", default="esmc-600m-2024-12", choices=["esmc-300m-2024-12", "esmc-600m-2024-12", "esmc-6b-2024-12"], help="ESMC model variant")
    parser.add_argument("--top-constrained", type=int, default=15, help="Number of most constrained residues to print")

    args = parser.parse_args()
    api_key = get_api_key()

    seq = args.sequence
    p = Path(seq)
    if p.is_file():
        lines = [line.strip() for line in open(p) if line.strip() and not line.startswith(">")]
        seq = "".join(lines).upper()
    else:
        seq = "".join(seq.split()).upper()

    print(f"Running mutational scan on {len(seq)} aa using {args.model}...")

    body = {
        "model": args.model,
        "protein": {"sequence": seq},
        "config": {
            "sequence": True,
            "return_embeddings": False,
            "return_mean_embedding": False,
            "return_hidden_states": False,
            "return_mean_hidden_states": False
        },
        "potential_sequence_of_concern": False
    }

    req = urllib.request.Request(
        "https://biohub.ai/api/v1/logits",
        data=json.dumps(body).encode("utf-8"),
        headers={
            "Authorization": f"Bearer {api_key}",
            "Content-Type": "application/json",
            "User-Agent": "Antigravity-BiohubESM-Skill/1.0"
        },
        method="POST"
    )

    try:
        with urllib.request.urlopen(req) as resp:
            res = json.loads(resp.read().decode("utf-8"))
    except urllib.error.HTTPError as e:
        err = e.read().decode("utf-8", errors="replace")
        print(f"[ERROR] HTTP {e.code}: {err}", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"[ERROR] Query failed: {e}", file=sys.stderr)
        sys.exit(1)

    logits = res.get("logits")
    if not logits:
        print("[ERROR] Response did not contain logits.", file=sys.stderr)
        sys.exit(1)

    # If logits includes BOS/EOS tokens, trim to match sequence length
    if len(logits) == len(seq) + 2:
        logits = logits[1:-1]

    results = []
    for pos, (wt_aa, pos_logits) in enumerate(zip(seq, logits)):
        # Normalize probabilities over the 20 standard amino acids (first 20 or matching vocabulary)
        # Using top standard amino acid logits
        aa_probs = softmax(pos_logits[:20])
        entropy = -sum(p * math.log2(p) for p in aa_probs if p > 1e-12)
        
        # Sort predicted preferences
        ranked_aas = sorted(zip(STANDARD_AAS, aa_probs), key=lambda x: x[1], reverse=True)
        top_alt = [f"{aa}({prob:.2f})" for aa, prob in ranked_aas[:3]]

        results.append({
            "position": pos + 1,
            "wt": wt_aa,
            "entropy": entropy,
            "top_predictions": ", ".join(top_alt)
        })

    # Output TSV
    out_path = Path(args.output).resolve()
    with open(out_path, "w", encoding="utf-8") as f:
        f.write("Position\tWT\tEntropy\tTop_Substitutions\n")
        for r in results:
            f.write(f"{r['position']}\t{r['wt']}\t{r['entropy']:.3f}\t{r['top_predictions']}\n")

    print(f"\nSaved mutational scan data to: {out_path}")

    # Display most mutationally constrained (lowest entropy / load-bearing positions)
    constrained = sorted(results, key=lambda x: x["entropy"])[:args.top_constrained]
    print(f"\nTop {len(constrained)} Load-Bearing / Constrained Positions (Lowest Shannon Entropy):")
    print("-" * 65)
    print(f"{'Position':<10} {'WT':<5} {'Entropy (bits)':<16} {'Model Preferences'}")
    print("-" * 65)
    for c in constrained:
        print(f"{c['position']:<10} {c['wt']:<5} {c['entropy']:<16.3f} {c['top_predictions']}")


if __name__ == "__main__":
    main()
