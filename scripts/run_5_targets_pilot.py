#!/usr/bin/env python3
"""
Run 5-Target Pilot Test for "Approach 1: Reverse Binder Lookup":
1. Canonical Target Pockets from AlphaFold 3D structures.
2. Conditioned Binder Design via Biohub ESM3 (/api/v1/generate).
3. Biophysical Representation Extraction via Biohub ESMC 6B (Layer 79, 2560-d).
4. Screening 7,264 Natural Microprotein Embeddings for Binder Mimetics.
5. Export Results to TSV & Parquet.
"""

import json
import math
import os
import sys
import time
import urllib.request
from pathlib import Path
import numpy as np
import pandas as pd
import torch

BASE_DIR = Path(__file__).resolve().parent.parent
if str(BASE_DIR) not in sys.path:
    sys.path.insert(0, str(BASE_DIR))

from esm.sdk.api import ESMProtein
from scripts.test_esm3_api import get_key
from bindpred.esmc_client import ESMCClient

TARGETS_CONFIG = [
    {
        "gene": "SNX11",
        "acc": "Q9Y5W9",
        "name": "Sorting nexin-11",
        "pdb": "structures/pilot_targets/SNX11_Q9Y5W9.pdb",
        "pocket_start": 80,
        "pocket_end": 105,
        "pocket_desc": "PX domain peptide/phosphoinositide binding groove"
    },
    {
        "gene": "PSMA1",
        "acc": "P25786",
        "name": "Proteasome subunit alpha type-1",
        "pdb": "structures/pilot_targets/PSMA1_P25786.pdb",
        "pocket_start": 60,
        "pocket_end": 85,
        "pocket_desc": "20S proteasome alpha-ring subunit interface/gate"
    },
    {
        "gene": "RPRD1A",
        "acc": "Q96P16",
        "name": "Regulation of nuclear pre-mRNA domain-containing protein 1A",
        "pdb": "structures/pilot_targets/RPRD1A_Q96P16.pdb",
        "pocket_start": 20,
        "pocket_end": 45,
        "pocket_desc": "CID (CTD-interacting domain) binding cleft"
    },
    {
        "gene": "CIB2",
        "acc": "O75838",
        "name": "Calcium and integrin-binding family member 2",
        "pdb": "structures/pilot_targets/CIB2_O75838.pdb",
        "pocket_start": 50,
        "pocket_end": 75,
        "pocket_desc": "EF-hand target engagement hydrophobic groove"
    },
    {
        "gene": "THAP2",
        "acc": "Q9H0W7",
        "name": "THAP domain-containing protein 2",
        "pdb": "structures/pilot_targets/THAP2_Q9H0W7.pdb",
        "pocket_start": 30,
        "pocket_end": 55,
        "pocket_desc": "THAP zinc-coordinating / protein-docking motif"
    }
]

def sanitize_coordinates(coords_list):
    """Replace all NaNs and missing atoms with None for standard JSON serialization."""
    clean_coords = []
    for res in coords_list:
        clean_res = []
        for atom in res:
            if atom is None:
                clean_res.append([None, None, None])
            else:
                clean_res.append([
                    None if (v is None or (isinstance(v, float) and math.isnan(v))) else float(v)
                    for v in atom
                ])
        clean_coords.append(clean_res)
    return clean_coords

def generate_binder_for_pocket(api_key: str, pocket_seq: str, pocket_coords: list, binder_len: int = 30) -> str:
    """Generate a binder sequence using ESM3 conditioned on pocket coordinates."""
    empty_res = [[None, None, None]] * 37
    prompt_seq = pocket_seq + ("_" * binder_len)
    prompt_coords = pocket_coords + ([empty_res] * binder_len)

    body = {
        'model': 'esm3-open-2024-03',
        'track': 'sequence',
        'inputs': {
            'sequence': prompt_seq,
            'coordinates': prompt_coords
        },
        'num_steps': 15,
        'temperature': 0.5,
        'temperature_annealing': True
    }
    req = urllib.request.Request(
        'https://biohub.ai/api/v1/generate',
        data=json.dumps(body).encode('utf-8'),
        headers={'Authorization': f'Bearer {api_key}', 'Content-Type': 'application/json'},
        method='POST'
    )
    with urllib.request.urlopen(req) as resp:
        res = json.loads(resp.read().decode('utf-8'))
        full_seq = res['outputs']['sequence']
        binder_seq = full_seq[len(pocket_seq):]
        return binder_seq

def main():
    print("=" * 80)
    print("5-TARGET PILOT: REVERSE BINDER LOOKUP VIA BIOHUB ESM3 & ESMC 6B")
    print("=" * 80)

    api_key = get_key()
    esmc_client = ESMCClient(token=api_key, default_model="esmc-6b-2024-12")

    # Load 7,264 pre-computed microprotein embeddings
    mp_parquet = BASE_DIR / "bindpred" / "representations" / "microproteins" / "microproteins_sequence_embeddings.parquet"
    print(f"Loading pre-computed microprotein embedding library from {mp_parquet}...")
    df_mp = pd.read_parquet(mp_parquet)
    print(f"Loaded {len(df_mp):,} microproteins (embedding matrix: {len(df_mp)} x 2560).")

    mp_matrix = np.stack(df_mp["vector"].values) # (7264, 2560)
    mp_norms = np.linalg.norm(mp_matrix, axis=1)

    results = []

    for idx, target in enumerate(TARGETS_CONFIG, 1):
        gene = target["gene"]
        acc = target["acc"]
        desc = target["pocket_desc"]
        print(f"\n[{idx}/5] Processing Target: {gene} ({acc}) - {target['name']}")
        print(f"    Pocket: Residues {target['pocket_start']}–{target['pocket_end']} ({desc})")

        # Load AlphaFold PDB structure
        pdb_path = BASE_DIR / target["pdb"]
        prot = ESMProtein.from_pdb(str(pdb_path))
        full_coords = sanitize_coordinates(prot.coordinates.tolist())
        full_seq = prot.sequence

        pocket_seq = full_seq[target["pocket_start"]:target["pocket_end"]]
        pocket_coords = full_coords[target["pocket_start"]:target["pocket_end"]]

        # Generate Binder
        print(f"    Generating 30-aa binder sequence with ESM3 conditioned on pocket coordinates...")
        t0 = time.time()
        try:
            binder_seq = generate_binder_for_pocket(api_key, pocket_seq, pocket_coords, binder_len=30)
            t_gen = time.time() - t0
            print(f"    [OK] Generated Binder ({len(binder_seq)} aa) in {t_gen:.2f}s: {binder_seq}")
        except Exception as e:
            print(f"    [ERROR] Binder generation failed: {e}")
            continue

        # Extract ESMC 6B Embedding
        print(f"    Extracting ESMC 6B Layer 79 embedding for binder...")
        t0 = time.time()
        rep_res = esmc_client.get_representation(binder_seq, layer_idx=-2, per_residue=False)
        binder_vec = rep_res["vector"]
        t_embed = time.time() - t0
        print(f"    [OK] Extracted vector (norm: {np.linalg.norm(binder_vec):.2f}) in {t_embed:.2f}s")

        # Cosine Similarity Search across 7,264 microproteins
        norm_b = np.linalg.norm(binder_vec)
        cos_sims = np.dot(mp_matrix, binder_vec) / (mp_norms * norm_b)

        # Retrieve Top 5 hits
        top5_indices = np.argsort(-cos_sims)[:5]
        print(f"    >>> Top 5 Microprotein Binders for {gene} Pocket:")
        for r_idx, hit_i in enumerate(top5_indices, 1):
            hit_row = df_mp.iloc[hit_i]
            sim = cos_sims[hit_i]
            print(f"        #{r_idx}: {hit_row['id']} (Gene: {hit_row['gene']}, {hit_row['length']} aa) | Cosine Sim: {sim:.4f}")

            results.append({
                "target_gene": gene,
                "target_uniprot": acc,
                "target_name": target["name"],
                "pocket_start": target["pocket_start"],
                "pocket_end": target["pocket_end"],
                "pocket_description": desc,
                "pocket_sequence": pocket_seq,
                "designed_binder_sequence": binder_seq,
                "hit_rank": r_idx,
                "microprotein_id": hit_row["id"],
                "microprotein_gene": hit_row["gene"],
                "microprotein_length": hit_row["length"],
                "cosine_similarity": round(float(sim), 4)
            })

    # Save consolidated results
    res_df = pd.DataFrame(results)
    tsv_out = BASE_DIR / "reports" / "pilot_5_targets_binder_search.tsv"
    parquet_out = BASE_DIR / "reports" / "pilot_5_targets_binder_search.parquet"
    res_df.to_csv(tsv_out, sep="\t", index=False)
    res_df.to_parquet(parquet_out, index=False)

    print("\n" + "=" * 80)
    print("PILOT SCREENING COMPLETE")
    print(f"Saved results to:\n  - {tsv_out}\n  - {parquet_out}")
    print("=" * 80)

if __name__ == "__main__":
    main()
