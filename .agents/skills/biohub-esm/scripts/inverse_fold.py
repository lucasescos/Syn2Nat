#!/usr/bin/env python3
"""
Inverse folding script using Biohub ESM3.
Reads backbone 3D coordinates from a PDB or mmCIF (.cif) file and queries POST /api/v1/inverse_fold
to design matching amino acid sequences.
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


def parse_pdb_coords(pdb_path: Path) -> tuple[list, str]:
    """Extracts backbone N, CA, C coordinates per residue from PDB."""
    residues = {}
    res_names = {}
    with open(pdb_path, "r", encoding="utf-8", errors="replace") as f:
        for line in f:
            if line.startswith(("ATOM", "HETATM")):
                atom_name = line[12:16].strip()
                if atom_name in ("N", "CA", "C"):
                    chain = line[21].strip()
                    res_num = line[22:26].strip()
                    key = (chain, res_num)
                    if key not in residues:
                        residues[key] = {}
                        res_names[key] = line[17:20].strip()
                    x = float(line[30:38].strip())
                    y = float(line[38:46].strip())
                    z = float(line[46:54].strip())
                    residues[key][atom_name] = [x, y, z]

    # Standard 37-atom or 3-atom backbone
    coords_list = []
    orig_seq = []
    aa_3to1 = {
        "ALA": "A", "CYS": "C", "ASP": "D", "GLU": "E", "PHE": "F",
        "GLY": "G", "HIS": "H", "ILE": "I", "LYS": "K", "LEU": "L",
        "MET": "M", "ASN": "N", "PRO": "P", "GLN": "Q", "ARG": "R",
        "SER": "S", "THR": "T", "VAL": "V", "TRP": "W", "TYR": "Y"
    }
    for key in sorted(residues.keys(), key=lambda k: (k[0], int(k[1]) if k[1].lstrip("-").isdigit() else k[1])):
        r = residues[key]
        if "CA" in r:
            ca = r["CA"]
            n = r.get("N", ca)
            c = r.get("C", ca)
            coords_list.append([n, ca, c])
            orig_seq.append(aa_3to1.get(res_names.get(key, ""), "X"))

    return coords_list, "".join(orig_seq)


def main():
    parser = argparse.ArgumentParser(description="Inverse fold backbone coordinates to sequence using ESM3")
    parser.add_argument("structure_file", help="Path to PDB or mmCIF structure file")
    parser.add_argument("-o", "--output", default="designed_sequence.fasta", help="Output FASTA filename")
    parser.add_argument("--temperature", type=float, default=0.1, help="Sampling temperature (default: 0.1, lower is more conservative)")
    parser.add_argument("--model", default="esm3-open-2024-03", help="ESM3 model variant")
    parser.add_argument("--num-samples", type=int, default=1, help="Number of candidate sequences to design")

    args = parser.parse_args()
    api_key = get_api_key()
    struct_path = Path(args.structure_file)

    if not struct_path.is_file():
        print(f"[ERROR] Structure file not found: {struct_path}", file=sys.stderr)
        sys.exit(1)

    print(f"Parsing backbone coordinates from {struct_path.name}...")
    coords, orig_seq = parse_pdb_coords(struct_path)
    if not coords:
        print("[ERROR] Could not extract valid backbone coordinates from structure file.", file=sys.stderr)
        sys.exit(1)

    print(f"Extracted backbone for {len(coords)} residues.")
    if orig_seq and "X" not in orig_seq:
        print(f"Original sequence (from structure): {orig_seq[:30]}... (len {len(orig_seq)})")

    headers = {
        "Authorization": f"Bearer {api_key}",
        "Content-Type": "application/json",
        "User-Agent": "Antigravity-BiohubESM-Skill/1.0"
    }

    designed_sequences = []
    for sample_idx in range(args.num_samples):
        body = {
            "model": args.model,
            "coordinates": coords,
            "inverse_folding_config": {
                "temperature": args.temperature,
                "invalid_ids": []
            },
            "potential_sequence_of_concern": False
        }

        req = urllib.request.Request("https://biohub.ai/api/v1/inverse_fold", data=json.dumps(body).encode("utf-8"), headers=headers, method="POST")
        try:
            with urllib.request.urlopen(req) as resp:
                res = json.loads(resp.read().decode("utf-8"))
                seq = res.get("sequence") or (res.get("tracks", {}).get("sequence") if isinstance(res.get("tracks"), dict) else None)
                if seq:
                    designed_sequences.append(seq)
                    print(f"Sample #{sample_idx + 1}: {seq}")
                else:
                    print(f"[WARN] No sequence field in response for sample {sample_idx + 1}: {res}")
        except urllib.error.HTTPError as e:
            err = e.read().decode("utf-8", errors="replace")
            print(f"[ERROR] HTTP {e.code}: {err}", file=sys.stderr)
            sys.exit(1)
        except Exception as e:
            print(f"[ERROR] Query failed: {e}", file=sys.stderr)
            sys.exit(1)

    if designed_sequences:
        out_path = Path(args.output).resolve()
        with open(out_path, "w", encoding="utf-8") as f:
            for i, s in enumerate(designed_sequences):
                f.write(f">designed_seq_{i+1}_temp_{args.temperature}\n{s}\n")
        print(f"\nSaved {len(designed_sequences)} designed sequence(s) to: {out_path}")


if __name__ == "__main__":
    main()
