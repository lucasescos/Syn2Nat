import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
ids_path = ROOT / "bindpred" / "representations" / "human" / "completed_ids.json"
ids = json.load(open(ids_path))

ids_json_str = json.dumps(ids)

script = f'''# Full-scale A100 Extraction for 19,052 Human Proteins
import gc
import json
import os
import time
from pathlib import Path
import numpy as np
import pandas as pd
import torch
from esm.models.esmc.compatibility import _BatchedESMProteinTensor
from esm.sdk.api import ESMProtein

# 1. Load completed IDs (1,600 targets already extracted locally)
COMPLETED_IDS = set({ids_json_str})
print(f"Loaded {{len(COMPLETED_IDS):,}} completed IDs to skip.")

# 2. Filter remaining targets
df_remaining = df_all[~df_all["id"].isin(COMPLETED_IDS)].reset_index(drop=True)
print(f"Total human targets to process: {{len(df_remaining):,}}")

# Setup output directory
BATCHES_DIR = Path("human_batches")
BATCHES_DIR.mkdir(parents=True, exist_ok=True)

BATCH_SIZE = 200
START_BATCH_IDX = 8  # batches 0..7 already completed locally
total_remaining = len(df_remaining)
total_batches = (total_remaining + BATCH_SIZE - 1) // BATCH_SIZE

print(f"Processing {{total_remaining:,}} proteins in {{total_batches}} batches (size {{BATCH_SIZE}})...")
overall_t0 = time.time()
total_processed = 0

for b_i in range(total_batches):
    batch_idx = START_BATCH_IDX + b_i
    p_file = BATCHES_DIR / f"batch_{{batch_idx:04d}}_sequence.parquet"
    n_file = BATCHES_DIR / f"batch_{{batch_idx:04d}}_residues.npz"

    # Checkpoint: skip if already on disk
    if p_file.exists() and n_file.exists():
        print(f"[SKIP] Batch {{batch_idx:04d}} already exists.")
        continue

    start_idx = b_i * BATCH_SIZE
    end_idx = min(start_idx + BATCH_SIZE, total_remaining)
    batch_df = df_remaining.iloc[start_idx:end_idx]

    b_t0 = time.time()
    batch_results = []

    for _, row in batch_df.iterrows():
        uid = row["id"]
        gene = row["gene"]
        full_seq = row["sequence"]
        orig_len = len(full_seq)
        seq = full_seq[:2048]

        try:
            pt = model.encode(ESMProtein(sequence=seq))
            bpt = _BatchedESMProteinTensor.from_protein_tensor(pt)

            with torch.no_grad():
                res = model.model(input_ids=bpt.sequence, output_hidden_states=True)

            hs = res.hidden_states[79][0, 1:-1]
            mean_vec = hs.mean(dim=0).float().cpu().numpy()
            residues_mat = hs.float().cpu().numpy().astype(np.float16)

            batch_results.append({{
                "id": uid,
                "gene": gene,
                "length": len(seq),
                "orig_length": orig_len,
                "vector": mean_vec,
                "residues": residues_mat
            }})
        except Exception as e:
            print(f"[ERROR] Failed {{uid}}: {{e}}")

    if not batch_results:
        continue

    # Save batch
    npz_dict = {{r["id"]: r["residues"] for r in batch_results}}
    np.savez_compressed(n_file, **npz_dict)

    parquet_df = pd.DataFrame([{{
        "id": r["id"],
        "gene": r["gene"],
        "length": r["length"],
        "orig_length": r["orig_length"],
        "vector": r["vector"].tolist()
    }} for r in batch_results])
    parquet_df.to_parquet(p_file, index=False)

    total_processed += len(batch_results)
    b_time = time.time() - b_t0
    elapsed = time.time() - overall_t0
    rate = total_processed / (elapsed + 1e-6)
    rem_time = (total_remaining - total_processed) / (rate + 1e-6)

    print(f"[SAVED] Batch {{batch_idx:04d}} ({{len(batch_results)}} in {{b_time:.1f}}s) | Total: {{total_processed}}/{{total_remaining}} ({{total_processed/total_remaining*100:.1f}}%) | Rate: {{rate:.1f}} seq/s | ETA: {{rem_time/60:.1f}} min")

    del batch_results, npz_dict, parquet_df
    torch.cuda.empty_cache()
    gc.collect()

print(f"\\n[ALL COMPLETE] Finished {{total_processed:,}} proteins in {{(time.time()-overall_t0)/60:.1f}} minutes!")
'''

out_file = ROOT / "bindpred" / "representations" / "human" / "colab_run_code.py"
open(out_file, "w", encoding="utf-8").write(script)
print(f"Written {out_file} ({len(script):,} chars)")
