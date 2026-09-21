#!/usr/bin/env python3
"""
Tokenization & Track conversion script using Biohub ESM3/ESMC.
Queries POST /api/v1/encode and POST /api/v1/decode to convert between biological tracks
(sequence, coordinates, secondary structure, SASA, function) and discrete model tokens.
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
    parser = argparse.ArgumentParser(description="Encode biological tracks to tokens or decode tokens back with ESM")
    subparsers = parser.add_subparsers(dest="command", required=True)

    # Subcommand: encode
    enc_p = subparsers.add_parser("encode", help="Encode biological tracks to discrete tokens")
    enc_p.add_argument("sequence", help="Protein sequence or path to tracks JSON file")
    enc_p.add_argument("-o", "--output", default="tokens.json", help="Output JSON filename")
    enc_p.add_argument("--model", default="esm3-open-2024-03", help="Model variant")

    # Subcommand: decode
    dec_p = subparsers.add_parser("decode", help="Decode discrete tokens back to biological tracks")
    dec_p.add_argument("tokens_file", help="Path to tokens JSON file")
    dec_p.add_argument("-o", "--output", default="decoded_tracks.json", help="Output JSON filename")
    dec_p.add_argument("--model", default="esm3-open-2024-03", help="Model variant")

    args = parser.parse_args()
    api_key = get_api_key()

    headers = {
        "Authorization": f"Bearer {api_key}",
        "Content-Type": "application/json",
        "User-Agent": "Antigravity-BiohubESM-Skill/1.0"
    }

    if args.command == "encode":
        seq_input = args.sequence
        p = Path(seq_input)
        if p.is_file():
            tracks_input = json.loads(p.read_text(encoding="utf-8"))
        else:
            tracks_input = {"sequence": "".join(seq_input.split()).upper()}

        body = {
            "model": args.model,
            "inputs": tracks_input,
            "potential_sequence_of_concern": False
        }
        print(f"Encoding tracks to tokens using {args.model}...")
        endpoint = "https://biohub.ai/api/v1/encode"

    elif args.command == "decode":
        p = Path(args.tokens_file)
        if not p.is_file():
            print(f"[ERROR] Tokens file not found: {p}", file=sys.stderr)
            sys.exit(1)
        tokens_input = json.loads(p.read_text(encoding="utf-8"))
        if "tokens" in tokens_input:
            tokens_input = tokens_input["tokens"]

        body = {
            "model": args.model,
            "inputs": tokens_input,
            "potential_sequence_of_concern": False
        }
        print(f"Decoding tokens to tracks using {args.model}...")
        endpoint = "https://biohub.ai/api/v1/decode"

    req = urllib.request.Request(endpoint, data=json.dumps(body).encode("utf-8"), headers=headers, method="POST")

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

    out_path = Path(args.output).resolve()
    with open(out_path, "w", encoding="utf-8") as f:
        json.dump(res, f, indent=2)

    print(f"Saved conversion to: {out_path}")
    if args.command == "encode" and "sequence" in res.get("tokens", res):
        toks = res.get("tokens", res).get("sequence", [])
        print(f"Token count: {len(toks)} (preview: {toks[:10]}...)")
    elif args.command == "decode":
        tracks = res.get("tracks", res)
        if "sequence" in tracks:
            print(f"Decoded Sequence: {tracks['sequence']}")


if __name__ == "__main__":
    main()
