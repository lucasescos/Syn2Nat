#!/usr/bin/env python3
"""
Multi-Scale Size-Aware Binder Design & Microprotein Screening.

Tests 3 distinct binder size scales on canonical target pockets:
1. Short Linear Motif / Loop (L = 20 aa)
2. Intermediate Secondary Structure Element (L = 40 aa)
3. Structured Microprotein Mini-Domain (L = 70 aa)

Features:
- Full heavy-atom 3D coordinate conditioning via Biohub ESM3.
- Layer 79 biophysical representation extraction via Biohub ESMC 6B (2560-d).
- Strict microprotein catalog filtering (Length <= 100 aa).
- Size-stratified and global similarity retrieval across 7,264 natural smORFs.
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

BASE_DIR = Path(__file__).resolve().parent.parent
if str(BASE_DIR) not in sys.path:
    sys.path.insert(0, str(BASE_DIR))

from esm.sdk.api import ESMProtein
from scripts.test_esm3_api import get_key
from bindpred.esmc_client import ESMCClient

# Test targets
TARGETS = [
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
        "gene": "PSMA1",
        "acc": "P25786",
        "name": "Proteasome subunit alpha type-1",
        "pdb": "structures/pilot_targets/PSMA1_P25786.pdb",
        "pocket_start": 60,
        "pocket_end": 85,
        "pocket_desc": "20S proteasome alpha-ring subunit interface/gate"
    },
    {
        "gene": "SNX11",
        "acc": "Q9Y5W9",
        "name": "Sorting nexin-11",
        "pdb": "structures/pilot_targets/SNX11_Q9Y5W9.pdb",
        "pocket_start": 80,
        "pocket_end": 105,
        "pocket_desc": "PX domain peptide/phosphoinositide binding groove"
    }
]

BINDER_SCALES = [
    {"scale": "Short Motif (20 aa)", "len": 20, "search_min": 10, "search_max": 40},
    {"scale": "Intermediate Helix (40 aa)", "len": 40, "search_min": 25, "search_max": 60},
    {"scale": "Mini-Domain (70 aa)", "len": 70, "search_min": 50, "search_max": 100}
]

def sanitize_coordinates(coords_list):
    """Replace all NaNs and missing atoms with None for valid JSON serialization."""
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

def generate_binder(api_key: str, pocket_seq: str, pocket_coords: list, binder_len: int) -> str:
    """Condition on 3D pocket heavy atoms and unmask binder sequence."""
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
        'num_steps': min(20, max(10, binder_len // 2)),
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
        return full_seq[len(pocket_seq):]

def main():
    print("=" * 80)
    print("MULTI-SCALE BINDER DESIGN & SIZE-GATED MICROPROTEIN SCREENING")
    print("=" * 80)

    api_key = get_key()
    esmc = ESMCClient(token=api_key, default_model="esmc-6b-2024-12")

    # Load microprotein embeddings
    mp_parquet = BASE_DIR / "bindpred" / "representations" / "microproteins" / "microproteins_sequence_embeddings.parquet"
    print(f"Loading microprotein embedding database from {mp_parquet}...")
    df_mp = pd.read_parquet(mp_parquet)

    # Enforce strict microprotein size cap (<= 100 aa)
    df_mp_clean = df_mp[df_mp["length"] <= 100].copy().reset_index(drop=True)
    print(f"Loaded {len(df_mp):,} total entries -> Filtered to {len(df_mp_clean):,} bona fide smORFs (length <= 100 aa).")

    mp_matrix = np.stack(df_mp_clean["vector"].values) # (N, 2560)
    mp_norms = np.linalg.norm(mp_matrix, axis=1)

    all_records = []

    for t_idx, target in enumerate(TARGETS, 1):
        gene = target["gene"]
        acc = target["acc"]
        pdb_path = BASE_DIR / target["pdb"]
        print(f"\n================================================================================")
        print(f"TARGET [{t_idx}/{len(TARGETS)}]: {gene} ({acc}) - {target['name']}")
        print(f"Pocket: Residues {target['pocket_start']}–{target['pocket_end']} ({target['pocket_desc']})")
        print(f"================================================================================")

        prot = ESMProtein.from_pdb(str(pdb_path))
        full_coords = sanitize_coordinates(prot.coordinates.tolist())
        full_seq = prot.sequence

        pocket_seq = full_seq[target["pocket_start"]:target["pocket_end"]]
        pocket_coords = full_coords[target["pocket_start"]:target["pocket_end"]]
        print(f"Pocket Sequence ({len(pocket_seq)} aa): {pocket_seq}")

        for s_idx, scale_cfg in enumerate(BINDER_SCALES, 1):
            scale_name = scale_cfg["scale"]
            b_len = scale_cfg["len"]
            min_len = scale_cfg["search_min"]
            max_len = scale_cfg["search_max"]

            print(f"\n  [Scale {s_idx}/3] {scale_name} (Length: {b_len} aa)")

            # 1. Generate Binder with ESM3
            t0 = time.time()
            binder_seq = generate_binder(api_key, pocket_seq, pocket_coords, binder_len=b_len)
            t_gen = time.time() - t0
            print(f"    - Generated Binder ({len(binder_seq)} aa) in {t_gen:.2f}s:")
            print(f"      Sequence: {binder_seq}")

            # 2. Extract ESMC 6B Embedding
            t0 = time.time()
            rep = esmc.get_representation(binder_seq, layer_idx=-2, per_residue=False)
            b_vec = rep["vector"]
            t_emb = time.time() - t0
            print(f"    - Extracted ESMC 6B vector in {t_emb:.2f}s (norm: {np.linalg.norm(b_vec):.1f})")

            # 3. Compute Cosine Similarity against all clean microproteins (<= 100 aa)
            norm_b = np.linalg.norm(b_vec)
            cos_sims = np.dot(mp_matrix, b_vec) / (mp_norms * norm_b)
            df_mp_clean["current_sim"] = cos_sims

            # A. Size-Matched Top Hits (within targeted bracket)
            df_bracket = df_mp_clean[(df_mp_clean["length"] >= min_len) & (df_mp_clean["length"] <= max_len)]
            top_bracket = df_bracket.sort_values("current_sim", ascending=False).head(3)

            # B. Global smORF Top Hits (<= 100 aa)
            top_global = df_mp_clean.sort_values("current_sim", ascending=False).head(3)

            print(f"    >>> Top Size-Matched Hits ({min_len}–{max_len} aa):")
            for rank, (_, row) in enumerate(top_bracket.iterrows(), 1):
                print(f"        #{rank}: {row['id']:<18} Gene: {row['gene']:<15} Length: {row['length']:<3} aa | Sim: {row['current_sim']:.4f}")

            print(f"    >>> Top Global smORF Hits (<= 100 aa):")
            for rank, (_, row) in enumerate(top_global.iterrows(), 1):
                print(f"        #{rank}: {row['id']:<18} Gene: {row['gene']:<15} Length: {row['length']:<3} aa | Sim: {row['current_sim']:.4f}")
                
                all_records.append({
                    "target_gene": gene,
                    "target_uniprot": acc,
                    "scale_name": scale_name,
                    "binder_length": b_len,
                    "designed_binder": binder_seq,
                    "hit_rank": rank,
                    "microprotein_id": row["id"],
                    "microprotein_gene": row["gene"],
                    "microprotein_length": row["length"],
                    "cosine_similarity": round(float(row["current_sim"]), 4),
                    "is_size_matched": bool(min_len <= row["length"] <= max_len)
                })

    # Export results
    out_df = pd.DataFrame(all_records)
    tsv_out = BASE_DIR / "reports" / "multiscale_binder_pilot_results.tsv"
    parquet_out = BASE_DIR / "reports" / "multiscale_binder_pilot_results.parquet"
    out_df.to_csv(tsv_out, sep="\t", index=False)
    out_df.to_parquet(parquet_out, index=False)

    print("\n" + "=" * 80)
    print("MULTI-SCALE SCREENING COMPLETE")
    print(f"Exported to:\n  - {tsv_out}\n  - {parquet_out}")
    print("=" * 80)

if __name__ == "__main__":
    main()
