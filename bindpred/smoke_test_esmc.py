#!/usr/bin/env python3
"""
Smoke test comparing ESMC-600M and ESMC-6B representations on a microprotein and a human target.
Saves results to bindpred/smoke_test_output.json.
"""

import json
import os
import sys
from pathlib import Path
import numpy as np
import pandas as pd

# Add repo root to path
ROOT_DIR = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT_DIR))

from bindpred.esmc_client import ESMCClient


def cosine_similarity(u: np.ndarray, v: np.ndarray) -> float:
    return float(np.dot(u, v) / (np.linalg.norm(u) * np.linalg.norm(v) + 1e-9))


def run_smoke_test():
    print("=" * 80)
    print("ESMC Latent Biophysical Representation Extraction: Smoke Test")
    print("Comparing ESMC-600M vs ESMC-6B Frontier World Model")
    print("=" * 80)

    client = ESMCClient()
    results = {}

    # 1. Load microprotein target
    orfs_path = ROOT_DIR / "data" / "parsed" / "orfs.parquet"
    if not orfs_path.exists():
        raise FileNotFoundError(f"Missing {orfs_path}")

    df_orfs = pd.read_parquet(orfs_path)
    # Target microprotein: c11norep11 (IGF2-AS, 131 aa)
    micro_sample = df_orfs[df_orfs["orf_id"] == "c11norep11"]
    if micro_sample.empty:
        micro_sample = df_orfs.iloc[[0]]

    micro_row = micro_sample.iloc[0]
    micro_id = str(micro_row["orf_id"])
    micro_gene = str(micro_row.get("gene", "N/A"))
    micro_seq = str(micro_row["sequence"])

    print(f"\n[Microprotein Target] ID: {micro_id} | Gene: {micro_gene} | Length: {len(micro_seq)} aa")
    print(f"Sequence: {micro_seq[:50]}... (total {len(micro_seq)} aa)")

    # 2. Load human canonical target
    human_path = ROOT_DIR / "data" / "targets" / "human_canonical_proteome.parquet"
    if not human_path.exists():
        raise FileNotFoundError(f"Missing {human_path}")

    df_human = pd.read_parquet(human_path)
    # Target human protein: P01344 (IGF2, 180 aa, biological cognate of IGF2-AS!)
    human_sample = df_human[df_human["uniprot_id"] == "P01344"]
    if human_sample.empty:
        human_sample = df_human.iloc[[0]]

    human_row = human_sample.iloc[0]
    human_id = str(human_row["uniprot_id"])
    human_gene = str(human_row.get("gene_symbol", "N/A"))
    human_name = str(human_row.get("protein_name", "N/A"))
    human_seq = str(human_row["sequence"])

    print(f"\n[Human Canonical Target] ID: {human_id} | Gene: {human_gene} | Name: {human_name} | Length: {len(human_seq)} aa")
    print(f"Sequence: {human_seq[:50]}... (total {len(human_seq)} aa)")

    # Models to benchmark
    models = [
        ("esmc-600m-2024-12", "ESMC 600M (Standard High-Accuracy, 37 layers, 1152-dim)"),
        ("esmc-6b-2024-12", "ESMC 6B (Frontier World Model, 81 layers, 2560-dim)"),
    ]

    for model_name, desc in models:
        print(f"\n" + "-" * 80)
        print(f"Evaluating Model: {model_name}")
        print(f"Description:      {desc}")
        print("-" * 80)

        # Microprotein extraction
        print(f"-> Extracting representations for microprotein {micro_id}...")
        micro_res = client.get_representation(micro_seq, model=model_name, layer_idx=-2)
        print(
            f"   Done in {micro_res['latency_sec']}s | "
            f"Layers: {micro_res['num_layers']} | "
            f"Penultimate Layer Dim: {micro_res['hidden_dim']} | "
            f"L2 Norm: {micro_res['l2_norm']:.4f}"
        )

        # Human protein extraction
        print(f"-> Extracting representations for human protein {human_id} ({human_gene})...")
        human_res = client.get_representation(human_seq, model=model_name, layer_idx=-2)
        print(
            f"   Done in {human_res['latency_sec']}s | "
            f"Layers: {human_res['num_layers']} | "
            f"Penultimate Layer Dim: {human_res['hidden_dim']} | "
            f"L2 Norm: {human_res['l2_norm']:.4f}"
        )

        # Latent cosine similarity
        sim = cosine_similarity(micro_res["vector"], human_res["vector"])
        print(f"   Latent Cosine Similarity ({micro_id} <-> {human_id}): {sim:.4f}")

        results[model_name] = {
            "description": desc,
            "total_layers": micro_res["num_layers"],
            "hidden_dimension": micro_res["hidden_dim"],
            "selected_layer": micro_res["selected_layer"],
            "microprotein": {
                "id": micro_id,
                "gene": micro_gene,
                "length_aa": len(micro_seq),
                "l2_norm": micro_res["l2_norm"],
                "mean_activation": micro_res["mean_val"],
                "std_activation": micro_res["std_val"],
                "latency_sec": micro_res["latency_sec"],
                "vector_sample": [round(float(x), 5) for x in micro_res["vector"][:8]],
            },
            "human_protein": {
                "id": human_id,
                "gene": human_gene,
                "protein_name": human_name,
                "length_aa": len(human_seq),
                "l2_norm": human_res["l2_norm"],
                "mean_activation": human_res["mean_val"],
                "std_activation": human_res["std_val"],
                "latency_sec": human_res["latency_sec"],
                "vector_sample": [round(float(x), 5) for x in human_res["vector"][:8]],
            },
            "latent_cosine_similarity": round(sim, 5),
        }

    out_file = ROOT_DIR / "bindpred" / "smoke_test_output.json"
    with open(out_file, "w", encoding="utf-8") as f:
        json.dump(results, f, indent=2)

    print("\n" + "=" * 80)
    print(f"[SUCCESS] Smoke test completed cleanly. Results saved to:")
    print(f"          {out_file}")
    print("=" * 80)


if __name__ == "__main__":
    run_smoke_test()
