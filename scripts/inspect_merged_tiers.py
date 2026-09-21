import pandas as pd

tiers = pd.read_parquet("data/parsed/orf_tiers.parquet")
orfs = pd.read_parquet("data/parsed/orfs.parquet")
top50 = pd.read_parquet("reports/top50_curated_candidates.parquet")

print("Tiers columns:", tiers.columns.tolist()[:10])
print("ORFs columns:", orfs.columns.tolist()[:10])

merged = pd.merge(tiers, orfs, on="orf_id", suffixes=("_tier", "_orf"))
t1_peps = merged[(merged["is_peptidein"] == True) & (merged["final_tier"].isin(["1A", "1B"]))]
rem_t1 = t1_peps[~t1_peps["orf_id"].isin(set(top50["orf_id"]))].copy()
print(f"Remaining Tier 1 peptideins count: {len(rem_t1)}")
len_col = "length_aa" if "length_aa" in rem_t1.columns else "length_orf"
print(f"Length col ({len_col}) stats:")
print(rem_t1[len_col].describe())
print("\nSample rows:")
print(rem_t1[["orf_id", "gene_tier", "orf_biotype_tier", "final_tier", len_col, "sequence"]].head(5).to_string(index=False))
