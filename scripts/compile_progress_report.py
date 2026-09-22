import json
from pathlib import Path
import pandas as pd
import numpy as np

base_dir = Path.cwd()
state_file = base_dir / "structures" / "cofolding_state.jsonl"
top50_path = base_dir / "reports" / "top50_microprotein_binding_partners_top10.parquet"
exp460_path = base_dir / "reports" / "expanded_ppi_candidates_460.parquet"
cont_path = base_dir / "reports" / "continuous_cofolding_candidates.parquet"
master_parquet = base_dir / "reports" / "master_microprotein_ppi_structural_atlas.parquet"
master_tsv = base_dir / "reports" / "master_microprotein_ppi_structural_atlas.tsv"

# Load state
records = []
with open(state_file, "r", encoding="utf-8") as f:
    for line in f:
        line = line.strip()
        if line:
            try:
                records.append(json.loads(line))
            except Exception:
                pass

df_state = pd.DataFrame(records).drop_duplicates(subset=["pair_id"], keep="last")

# Sync master atlas
cand_dfs = []
for p in [top50_path, exp460_path, cont_path]:
    if p.is_file():
        cand_dfs.append(pd.read_parquet(p))

if cand_dfs:
    all_cands = pd.concat(cand_dfs, ignore_index=True).drop_duplicates(subset=["pair_id"], keep="last")
    merged = pd.merge(all_cands, df_state, on="pair_id", how="inner", suffixes=("", "_state"))
    successful = merged[merged["status"] == "success"].copy()
    successful = successful.sort_values("iptm", ascending=False).reset_index(drop=True)
    successful.to_parquet(master_parquet, index=False)
    successful.to_csv(master_tsv, sep="\t", index=False)
    print(f"Master Atlas successfully updated: {len(successful)} complexes.")
else:
    successful = df_state[df_state["status"] == "success"].copy()

print(f"\n=======================================================")
print(f"  COFOLDING PROGRESS REPORT & ATLAS METRICS")
print(f"=======================================================")
print(f"Total Unique Pairs Evaluated in Ledger: {len(df_state)}")
print("Status Breakdown:")
for k, v in df_state["status"].value_counts().items():
    print(f"  - {k}: {v}")

print(f"\nTotal Successfully Folded Complexes: {len(successful)}")
print(f"Unique Microproteins in Successful Atlas: {successful['orf_id'].nunique() if 'orf_id' in successful.columns else 'N/A'}")
print(f"Unique Microprotein Genes in Atlas: {successful['microprotein_gene'].nunique() if 'microprotein_gene' in successful.columns else 'N/A'}")

iptm_s = successful["iptm"].dropna()
ptm_s = successful["ptm"].dropna()
up_plddt = successful["microprotein_plddt"].dropna()
p_plddt = successful["partner_plddt"].dropna()

print(f"\nInterface Predicted TM-score (ipTM):")
print(f"  Mean:   {iptm_s.mean():.4f}")
print(f"  Median: {iptm_s.median():.4f}")
print(f"  Max:    {iptm_s.max():.4f}")

print(f"\nComplex Predicted TM-score (pTM):")
print(f"  Mean:   {ptm_s.mean():.4f}")
print(f"  Median: {ptm_s.median():.4f}")
print(f"  Max:    {ptm_s.max():.4f}")

print(f"\nConfidence Stratification (ipTM thresholds):")
print(f"  - High Confidence (ipTM >= 0.75):       {(iptm_s >= 0.75).sum()}")
print(f"  - Plausible Interface (0.60 <= ipTM < 0.75): {((iptm_s >= 0.60) & (iptm_s < 0.75)).sum()}")
print(f"  - Weak / Transient (0.45 <= ipTM < 0.60):    {((iptm_s >= 0.45) & (iptm_s < 0.60)).sum()}")
print(f"  - Non-Interacting (ipTM < 0.45):             {(iptm_s < 0.45).sum()}")
print(f"  - Total ipTM >= 0.60 (High Quality):         {(iptm_s >= 0.60).sum()}")

print(f"\nTop 25 Complexes with Highest ipTM in the Entire Atlas:")
cols_show = ["pair_id", "microprotein_gene", "partner_gene_symbol", "iptm", "ptm", "microprotein_plddt", "partner_plddt", "binding_classification"]
available_cols = [c for c in cols_show if c in successful.columns]
print(successful[available_cols].head(25).to_string(index=False))

# Identify new high-confidence hits from continuous / AlphaGenome runs
print("\nNew High-Confidence Hits (ipTM >= 0.55) from Recent Continuous Runs:")
if "cohort" in successful.columns:
    new_runs = successful[successful["cohort"].str.contains("Continuous|AlphaGenome", na=False)]
    print(f"Total complexes folded in recent Continuous/AlphaGenome runs: {len(new_runs)}")
    top_new = new_runs[new_runs["iptm"] >= 0.50]
    print(f"Recent runs with ipTM >= 0.50: {len(top_new)}")
    if len(top_new) > 0:
        print(top_new[available_cols].to_string(index=False))
