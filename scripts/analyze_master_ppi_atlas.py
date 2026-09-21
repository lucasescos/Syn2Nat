import json
from pathlib import Path
import pandas as pd
import numpy as np

BASE_DIR = Path.cwd()
PDB_DIR = BASE_DIR / "structures" / "pdbs"
PAE_DIR = BASE_DIR / "structures" / "pae"
STATE_FILE = BASE_DIR / "structures" / "cofolding_state.jsonl"

TOP50_FILE = BASE_DIR / "reports" / "top50_microprotein_binding_partners_top10.parquet"
EXP_CAND_FILE = BASE_DIR / "reports" / "expanded_ppi_candidates_460.parquet"

MASTER_PARQUET = BASE_DIR / "reports" / "master_microprotein_ppi_structural_atlas.parquet"
MASTER_TSV = BASE_DIR / "reports" / "master_microprotein_ppi_structural_atlas.tsv"

print("Loading candidate definitions...")
top50_df = pd.read_parquet(TOP50_FILE)
exp_df = pd.read_parquet(EXP_CAND_FILE)

# Load all state records
records = []
with open(STATE_FILE, "r", encoding="utf-8") as f:
    for line in f:
        line = line.strip()
        if line:
            try:
                records.append(json.loads(line))
            except Exception:
                pass

state_df = pd.DataFrame(records).drop_duplicates(subset=["pair_id"], keep="last")
print(f"Total state records: {len(state_df)}")
print("Status breakdown in state:")
print(state_df["status"].value_counts())

# Build the master successfully folded dataset
# 1. Top 50 originally folded (345 pairs)
# 2. Track A backfilled pairs (155 pairs) -> replacing the 155 skipped in Top 50 to form a complete 500-complex Top 50 dataset!
# 3. Track B expanded Tier 1 peptidein pairs (305 pairs)
# Total successful: 345 + 155 + 305 = 805 pairs!

# Combine candidate definitions:
# For Top 50, let's create top50_complete (345 kept + 155 backfilled)
track_a_pairs = exp_df[exp_df["cohort"].isin(top50_df["cohort"].unique())].copy()
track_b_pairs = exp_df[exp_df["cohort"] == "5. Expanded Tier 1 Peptidein"].copy()
print(f"Track A candidates: {len(track_a_pairs)}, Track B candidates: {len(track_b_pairs)}")

# Map of pair_id to candidate metadata
all_cands = pd.concat([top50_df, exp_df], ignore_index=True).drop_duplicates(subset=["pair_id"], keep="last")

# Merge with state
master_df = pd.merge(all_cands, state_df, on="pair_id", how="inner", suffixes=("", "_state"))
# Filter for successful folds
successful_master = master_df[master_df["status"] == "success"].copy()
print(f"\nTotal successfully folded master complexes: {len(successful_master)}")
assert len(successful_master) == 805, f"Expected 805, got {len(successful_master)}"

# Verify PDB and PAE files exist on disk
pdb_exists = successful_master["pdb_path"].apply(lambda p: Path(p).is_file() and Path(p).stat().st_size > 0).all()
print(f"All 805 PDB files exist with valid size on disk: {pdb_exists}")

# Sort by ipTM descending
successful_master = successful_master.sort_values("iptm", ascending=False).reset_index(drop=True)

# Save Master Atlas
successful_master.to_parquet(MASTER_PARQUET, index=False)
successful_master.to_csv(MASTER_TSV, sep="\t", index=False)
print(f"Master atlas saved to {MASTER_PARQUET} and {MASTER_TSV}.")

# -------------------------------------------------------------
# Statistical Insights
# -------------------------------------------------------------
print("\n=== Master Structural Atlas Statistics (805 Complexes) ===")
print("Overall Metrics Summary:")
for col in ["iptm", "ptm", "mean_plddt", "microprotein_plddt", "partner_plddt", "esmfold_total_length", "elapsed_seconds"]:
    s = successful_master[col].dropna()
    print(f"  {col:22s}: Mean={s.mean():.3f}, Median={s.median():.3f}, Min={s.min():.3f}, Max={s.max():.3f}")

print("\nBinding Classification Breakdown:")
print(successful_master["binding_classification"].value_counts())

print("\nipTM Thresholds:")
print(f"  ipTM >= 0.70: {(successful_master['iptm'] >= 0.70).sum()}")
print(f"  ipTM >= 0.60: {(successful_master['iptm'] >= 0.60).sum()}")
print(f"  ipTM >= 0.50: {(successful_master['iptm'] >= 0.50).sum()}")
print(f"  ipTM >= 0.40: {(successful_master['iptm'] >= 0.40).sum()}")
print(f"  ipTM >= 0.30: {(successful_master['iptm'] >= 0.30).sum()}")

print("\nTop 30 Highest ipTM Complexes across the 805 Atlas:")
cols = ["microprotein_gene", "microprotein_biotype", "microprotein_length", "partner_gene_symbol", "partner_length", "esmfold_total_length", "iptm", "ptm", "microprotein_plddt", "partner_plddt", "cohort"]
print(successful_master[cols].head(30).to_string(index=True))
