#!/usr/bin/env python3
"""
Structure prediction helper script for Biohub ESMFold2.
Submits sequence(s) and complexes to https://biohub.ai/api/v1/fold or /fold_all_atom,
saving predicted mmCIF (.cif) structures and metrics (pLDDT, pTM, ipTM, PAE, distogram).
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
                    env_vars[key.strip()] = val.strip().strip("'\"")
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
    print("Please set ESM_API_KEY before running structure prediction.", file=sys.stderr)
    sys.exit(1)


def read_fasta_or_str(input_str: str) -> list[tuple[str, str]]:
    """Returns list of (id, sequence)."""
    p = Path(input_str)
    if p.is_file():
        entries = []
        cur_id = "A"
        cur_seq = []
        with open(p, "r", encoding="utf-8") as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                if line.startswith(">"):
                    if cur_seq:
                        entries.append((cur_id, "".join(cur_seq)))
                        cur_seq = []
                    cur_id = line[1:].split()[0]
                else:
                    cur_seq.append(line)
            if cur_seq:
                entries.append((cur_id, "".join(cur_seq)))
        return entries

    chains = input_str.split("|")
    chain_letters = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
    entries = []
    for i, seq in enumerate(chains):
        cid = chain_letters[i % len(chain_letters)]
        clean_seq = "".join(seq.split()).upper()
        entries.append((cid, clean_seq))
    return entries


def main():
    parser = argparse.ArgumentParser(description="Fold protein sequences and all-atom complexes with Biohub ESMFold2")
    parser.add_argument("sequence", help="Protein sequence, pipe-separated chains (e.g. 'SEQ1|SEQ2'), or path to FASTA file")
    parser.add_argument("-o", "--output", default="predicted_structure.cif", help="Output .cif or .pdb filename (default: predicted_structure.cif)")
    parser.add_argument("--model", default="esmfold2-fast-2026-05", choices=["esmfold2-fast-2026-05", "esmfold2-2026-05"], help="Model version")
    parser.add_argument("--num-loops", type=int, default=20, help="Number of folding loops (default: 20)")
    parser.add_argument("--num-steps", type=int, default=100, help="Diffusion sampling steps (default: 100)")
    parser.add_argument("--include-pae", action="store_true", help="Include predicted aligned error matrix")
    parser.add_argument("--include-distogram", action="store_true", help="Include distogram prediction")
    parser.add_argument("--include-embeddings", action="store_true", help="Include pair embeddings in response")
    parser.add_argument("--ligand", action="append", help="Small molecule ligand CCD code (e.g. 'SAH', 'ATP') or SMILES string")
    parser.add_argument("--dna", action="append", help="DNA sequence to co-fold into complex")
    parser.add_argument("--rna", action="append", help="RNA sequence to co-fold into complex")
    parser.add_argument("--pocket-chain", help="Binder chain ID for pocket conditioning")
    parser.add_argument("--pocket-contacts", help="JSON file or string describing target contact constraints")

    args = parser.parse_args()
    api_key = get_api_key()

    protein_chains = read_fasta_or_str(args.sequence)
    if not protein_chains:
        print("[ERROR] No valid sequence provided.", file=sys.stderr)
        sys.exit(1)

    is_all_atom = bool(args.ligand or args.dna or args.rna or args.pocket_chain)
    endpoint = "https://biohub.ai/api/v1/fold_all_atom" if is_all_atom else "https://biohub.ai/api/v1/fold"

    print(f"Folding with model '{args.model}' via {endpoint.split('/')[-1]}...")
    for cid, s in protein_chains:
        print(f"  Protein Chain {cid}: {len(s)} aa")

    config_payload = {
        "num_loops": args.num_loops,
        "num_sampling_steps": args.num_steps,
        "lm_dropout": 0.3,
        "include_pae": args.include_pae,
        "include_distogram": args.include_distogram,
        "include_pair_chains_iptm": len(protein_chains) > 1,
        "include_embeddings": args.include_embeddings
    }

    if is_all_atom:
        all_sequences = []
        for cid, s in protein_chains:
            all_sequences.append({"id": cid, "sequence": s, "type": "protein"})

        if args.dna:
            for i, dna_seq in enumerate(args.dna):
                did = f"DNA_{i+1}"
                all_sequences.append({"id": did, "sequence": dna_seq.upper(), "type": "dna"})
                print(f"  DNA Entity {did}: {len(dna_seq)} bp")

        if args.rna:
            for i, rna_seq in enumerate(args.rna):
                rid = f"RNA_{i+1}"
                all_sequences.append({"id": rid, "sequence": rna_seq.upper(), "type": "rna"})
                print(f"  RNA Entity {rid}: {len(rna_seq)} nt")

        if args.ligand:
            for i, lig in enumerate(args.ligand):
                lid = f"L_{i+1}"
                if len(lig) <= 5 and lig.isupper():
                    all_sequences.append({"id": lid, "ccd": [lig], "type": "ligand"})
                    print(f"  Ligand {lid}: CCD '{lig}'")
                else:
                    all_sequences.append({"id": lid, "smiles": lig, "type": "ligand"})
                    print(f"  Ligand {lid}: SMILES '{lig}'")

        all_atom_input = {"sequences": all_sequences}
        if args.pocket_chain and args.pocket_contacts:
            contacts = json.loads(args.pocket_contacts) if args.pocket_contacts.startswith("[") else json.loads(Path(args.pocket_contacts).read_text())
            all_atom_input["pocket"] = {
                "binder_chain_id": args.pocket_chain,
                "contacts": contacts
            }

        body = {
            "model": args.model,
            "all_atom_input": all_atom_input,
            "num_loops": args.num_loops,
            "num_sampling_steps": args.num_steps,
            "lm_dropout": 0.3,
            "include_pae": args.include_pae,
            "include_distogram": args.include_distogram,
            "include_pair_chains_iptm": True,
            "include_embeddings": args.include_embeddings,
            "potential_sequence_of_concern": False
        }
    else:
        sequences_payload = [{"id": cid, "sequence": s} for cid, s in protein_chains]
        body = {
            "model": args.model,
            "sequences": sequences_payload,
            "config": config_payload,
            "potential_sequence_of_concern": False
        }

    req_data = json.dumps(body).encode("utf-8")
    headers = {
        "Authorization": f"Bearer {api_key}",
        "Content-Type": "application/json",
        "User-Agent": "Antigravity-BiohubESM-Skill/1.0"
    }

    req = urllib.request.Request(endpoint, data=req_data, headers=headers, method="POST")

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

    complex_data = result.get("complex") or result.get("coordinates")
    if not complex_data:
        print("[ERROR] Response did not contain structure coordinate data.", file=sys.stderr)
        print(json.dumps(result, indent=2)[:500])
        sys.exit(1)

    out_path = Path(args.output).resolve()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w", encoding="utf-8") as f:
        f.write(complex_data if isinstance(complex_data, str) else json.dumps(complex_data))

    plddt_list = result.get("plddt", [])
    valid_plddt = [x for x in plddt_list if x is not None]
    mean_plddt = (sum(valid_plddt) / len(valid_plddt)) if valid_plddt else 0.0
    ptm = result.get("ptm")
    iptm = result.get("interface_ptm") or result.get("iptm")

    if mean_plddt >= 90:
        band = "Very high confidence (pLDDT >= 90)"
    elif mean_plddt >= 70:
        band = "Confident (70 <= pLDDT < 90)"
    elif mean_plddt >= 50:
        band = "Low confidence / flexible (50 <= pLDDT < 70)"
    else:
        band = "Very low confidence / disordered (pLDDT < 50)"

    print("\n--- ESMFold2 Prediction Results ---")
    print(f"Structure saved: {out_path}")
    print(f"Mean pLDDT:      {mean_plddt:.2f} ({band})")
    if ptm is not None:
        print(f"pTM score:       {float(ptm):.3f}")
    if iptm is not None:
        print(f"ipTM score:      {float(iptm):.3f}")

    if args.include_pae and "pae" in result:
        pae_path = out_path.with_suffix(".pae.json")
        with open(pae_path, "w", encoding="utf-8") as f:
            json.dump({"pae": result["pae"]}, f)
        print(f"PAE matrix saved:{pae_path}")

    if "pair_chains_iptm" in result and result["pair_chains_iptm"]:
        print("\nPairwise Chain ipTM Matrix:")
        for row in result["pair_chains_iptm"]:
            print(" ", [f"{v:.3f}" if v is not None else "---" for v in row])


if __name__ == "__main__":
    main()
