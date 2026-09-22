#!/usr/bin/env python3
"""
Test Approach 1 end-to-end on 1 target (CIB2).
1. Load PDB structure with ESMProtein.from_pdb.
2. Select pocket residues.
3. Call ESM3 generate with coordinates constraint.
4. Extract ESMC 6B embedding of designed binder.
5. Query 7,264 microprotein library and find top matching microproteins.
"""

import json
import os
import sys
import urllib.request
import numpy as np
import pandas as pd
import torch
from pathlib import Path

# Add project root to sys.path
BASE_DIR = Path(__file__).resolve().parent.parent
if str(BASE_DIR) not in sys.path:
    sys.path.insert(0, str(BASE_DIR))

from esm.sdk.api import ESMProtein
from scripts.test_esm3_api import get_key
from bindpred.esmc_client import ESMCClient

def run_cib2_test():
    api_key = get_key()
    pdb_path = "structures/pilot_targets/CIB2_O75838.pdb"
    print(f"Loading structure from {pdb_path}...")
    prot = ESMProtein.from_pdb(pdb_path)
    
    full_seq = prot.sequence
    full_coords = prot.coordinates.tolist()
    print(f"Target length: {len(full_seq)} aa")
    
    import math
    clean_coords = []
    for res in full_coords:
        clean_res = []
        for atom in res:
            if atom is None:
                clean_res.append([None, None, None])
            else:
                clean_res.append([None if (v is None or (isinstance(v, float) and math.isnan(v))) else float(v) for v in atom])
        clean_coords.append(clean_res)
    
    # Select a 25-residue pocket (e.g. residues 50 to 75, which forms a core EF-hand target binding cleft in CIB2)
    start_idx, end_idx = 50, 75
    pocket_seq = full_seq[start_idx:end_idx]
    pocket_coords = clean_coords[start_idx:end_idx]
    print(f"Pocket sequence ({len(pocket_seq)} aa): {pocket_seq}")
    
    # Binder length: 30 aa
    binder_len = 30
    prompt_seq = pocket_seq + ("_" * binder_len)
    
    # Unassigned coordinates for binder residues: 37 atoms of [None, None, None]
    empty_res = [[None, None, None]] * 37
    prompt_coords = pocket_coords + ([empty_res] * binder_len)
    
    print(f"Sending prompt to ESM3 generate (track=sequence, num_steps=15)...")
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
    try:
        with urllib.request.urlopen(req) as resp:
            res = json.loads(resp.read().decode('utf-8'))
            full_out_seq = res['outputs']['sequence']
            designed_binder = full_out_seq[len(pocket_seq):]
            print(f"[SUCCESS] Designed binder sequence ({len(designed_binder)} aa): {designed_binder}")
    except urllib.error.HTTPError as e:
        print("[ERROR] HTTP Error:", e.code, e.read().decode('utf-8', errors='replace'))
        return
        
    # Now extract ESMC 6B embedding of designed binder
    print("Extracting ESMC 6B layer 79 embedding for designed binder...")
    client = ESMCClient(token=api_key, default_model="esmc-6b-2024-12")
    rep_res = client.get_representation(designed_binder, layer_idx=-2, per_residue=False)
    binder_vec = rep_res["vector"]
    print(f"Binder vector shape: {binder_vec.shape}, norm: {np.linalg.norm(binder_vec):.2f}")
    
    # Load 7,264 microprotein sequence embeddings
    print("Loading 7,264 microprotein embeddings library...")
    df_mp = pd.read_parquet("bindpred/representations/microproteins/microproteins_sequence_embeddings.parquet")
    print(f"Loaded {len(df_mp)} microproteins.")
    
    # Compute cosine similarities
    mp_vectors = np.stack(df_mp["vector"].values) # Shape (7264, 2560)
    norm_mp = np.linalg.norm(mp_vectors, axis=1, keepdims=True)
    norm_binder = np.linalg.norm(binder_vec)
    
    cos_sims = np.dot(mp_vectors, binder_vec) / (norm_mp.squeeze() * norm_binder)
    df_mp["cosine_similarity"] = cos_sims
    
    top_hits = df_mp.sort_values("cosine_similarity", ascending=False).head(5)
    print("\n--- Top 5 Matching Microproteins for CIB2 Designed Binder ---")
    for rank, (_, row) in enumerate(top_hits.iterrows(), 1):
        print(f"#{rank}: {row['id']} (Gene: {row['gene']}, Length: {row['length']} aa) | Cosine Sim: {row['cosine_similarity']:.4f}")

if __name__ == "__main__":
    run_cib2_test()
