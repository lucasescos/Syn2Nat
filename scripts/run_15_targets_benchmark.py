#!/usr/bin/env python3
"""
Large-Scale Reverse Binder Lookup Benchmark across 15 Top Canonical Targets:
1. AlphaFold 3D Structures parsed into 37-atom heavy atom tensors.
2. High-confidence structural pockets (pLDDT > 90–98).
3. 35-aa Binder Generation via Biohub ESM3 (/api/v1/generate).
4. Biophysical Embedding via Biohub ESMC 6B (Layer 79, 2560-d).
5. Screening against 7,264 natural microproteins (filtered to Length <= 100 aa).
6. Comprehensive Tabular & Parquet Export.
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

TARGET_SPECS = [
    {"gene": "RPRD1A", "acc": "Q96P16", "name": "Pre-mRNA domain-containing protein 1A", "start": 25, "end": 50},
    {"gene": "MBTD1",  "acc": "Q05BQ5", "name": "MBT domain-containing protein 1", "start": 533, "end": 558},
    {"gene": "SNX11",  "acc": "Q9Y5W9", "name": "Sorting nexin-11", "start": 114, "end": 139},
    {"gene": "ABLIM1", "acc": "O14639", "name": "Actin-binding LIM protein 1", "start": 747, "end": 772},
    {"gene": "ATF7",   "acc": "P17544", "name": "Cyclic AMP-dependent TF ATF-7", "start": 356, "end": 381},
    {"gene": "PSMA1",  "acc": "P25786", "name": "Proteasome subunit alpha type-1", "start": 83, "end": 108},
    {"gene": "RNF217", "acc": "Q8TC41", "name": "E3 ubiquitin-protein ligase RNF217", "start": 502, "end": 527},
    {"gene": "ABI1",   "acc": "Q8IZP0", "name": "Abl interactor 1", "start": 87, "end": 112},
    {"gene": "PSMA5",  "acc": "P28066", "name": "Proteasome subunit alpha type-5", "start": 176, "end": 201},
    {"gene": "THAP2",  "acc": "Q9H0W7", "name": "THAP domain-containing protein 2", "start": 131, "end": 156},
    {"gene": "CIB2",   "acc": "O75838", "name": "Calcium and integrin-binding 2", "start": 143, "end": 168},
    {"gene": "RAB31",  "acc": "Q13636", "name": "Ras-related protein Rab-31", "start": 113, "end": 138},
    {"gene": "CDK15",  "acc": "Q96Q40", "name": "Cyclin-dependent kinase 15", "start": 362, "end": 387},
    {"gene": "BAP1",   "acc": "Q92560", "name": "Ubiquitin carboxyl-terminal hydrolase BAP1", "start": 106, "end": 131},
    {"gene": "DUSP16", "acc": "Q9BY84", "name": "Dual specificity phosphatase 16", "start": 237, "end": 262},
]

BINDER_LENGTH = 35

def sanitize_coordinates(coords_list):
    """Replace NaNs and missing coordinates with None for valid JSON serialization."""
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

def generate_binder(api_key: str, pocket_seq: str, pocket_coords: list, binder_len: int = 35) -> str:
    """Generate binder sequence conditioned on pocket 3D heavy atoms."""
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
        'num_steps': 18,
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
    print("15-TARGET REVERSE BINDER LOOKUP BENCHMARK")
    print("=" * 80)

    api_key = get_key()
    esmc = ESMCClient(token=api_key, default_model="esmc-6b-2024-12")

    # Load 7,264 microprotein embedding library
    mp_parquet = BASE_DIR / "bindpred" / "representations" / "microproteins" / "microproteins_sequence_embeddings.parquet"
    print(f"Loading microprotein embedding database from {mp_parquet}...")
    df_mp = pd.read_parquet(mp_parquet)

    # Filter to strict microprotein length (<= 100 aa)
    df_mp_smorfs = df_mp[df_mp["length"] <= 100].copy().reset_index(drop=True)
    print(f"Filtered catalog to {len(df_mp_smorfs):,} microproteins (Length <= 100 aa).")

    mp_matrix = np.stack(df_mp_smorfs["vector"].values) # (N, 2560)
    mp_norms = np.linalg.norm(mp_matrix, axis=1)

    all_hits = []

    for idx, spec in enumerate(TARGET_SPECS, 1):
        gene = spec["gene"]
        acc = spec["acc"]
        name = spec["name"]
        pdb_file = BASE_DIR / "structures" / "pilot_15_targets" / f"{gene}_{acc}.pdb"

        print(f"\n[{idx:2d}/15] Target: {gene} ({acc}) - {name}")
        print(f"       Pocket window: {spec['start']+1}–{spec['end']} (25 residues)")

        try:
            prot = ESMProtein.from_pdb(str(pdb_file))
            full_coords = sanitize_coordinates(prot.coordinates.tolist())
            full_seq = prot.sequence

            pocket_seq = full_seq[spec["start"]:spec["end"]]
            pocket_coords = full_coords[spec["start"]:spec["end"]]

            # 1. Generate Binder Sequence with ESM3
            t0 = time.time()
            binder_seq = generate_binder(api_key, pocket_seq, pocket_coords, binder_len=BINDER_LENGTH)
            t_gen = time.time() - t0
            print(f"       [ESM3] Binder ({len(binder_seq)} aa) generated in {t_gen:.2f}s: {binder_seq}")

            # 2. Extract ESMC 6B Layer 79 Embedding
            t0 = time.time()
            rep = esmc.get_representation(binder_seq, layer_idx=-2, per_residue=False)
            b_vec = rep["vector"]
            t_emb = time.time() - t0
            print(f"       [ESMC 6B] Embedding extracted in {t_emb:.2f}s (norm: {np.linalg.norm(b_vec):.1f})")

            # 3. Screen against Microprotein Library
            norm_b = np.linalg.norm(b_vec)
            sims = np.dot(mp_matrix, b_vec) / (mp_norms * norm_b)

            top5_idx = np.argsort(-sims)[:5]
            print(f"       >>> Top Natural Microprotein Mimetics:")
            for rank, hit_i in enumerate(top5_idx, 1):
                row = df_mp_smorfs.iloc[hit_i]
                sim_val = sims[hit_i]
                print(f"           #{rank}: {row['id']:<18} Gene: {row['gene']:<15} Len: {row['length']:<3} aa | Cosine Sim: {sim_val:.4f}")

                all_hits.append({
                    "target_rank": idx,
                    "target_gene": gene,
                    "target_uniprot": acc,
                    "target_name": name,
                    "pocket_start": spec["start"] + 1,
                    "pocket_end": spec["end"],
                    "pocket_sequence": pocket_seq,
                    "designed_binder": binder_seq,
                    "hit_rank": rank,
                    "microprotein_id": row["id"],
                    "microprotein_gene": row["gene"],
                    "microprotein_length": row["length"],
                    "cosine_similarity": round(float(sim_val), 4)
                })

            time.sleep(0.2) # Polite pacing

        except Exception as e:
            print(f"       [ERROR] Processing failed for {gene}: {e}")

    # Export consolidated benchmark dataset
    res_df = pd.DataFrame(all_hits)
    tsv_path = BASE_DIR / "reports" / "reverse_binder_lookup_15_targets.tsv"
    parquet_path = BASE_DIR / "reports" / "reverse_binder_lookup_15_targets.parquet"

    res_df.to_csv(tsv_path, sep="\t", index=False)
    res_df.to_parquet(parquet_path, index=False)

    print("\n" + "=" * 80)
    print("15-TARGET BENCHMARK COMPLETED")
    print(f"Total Matches Recorded: {len(res_df):,}")
    print(f"Mean Cosine Similarity: {res_df['cosine_similarity'].mean():.4f} (Max: {res_df['cosine_similarity'].max():.4f})")
    print(f"Saved artifacts to:")
    print(f"  - {tsv_path}")
    print(f"  - {parquet_path}")
    print("=" * 80)

if __name__ == "__main__":
    main()
