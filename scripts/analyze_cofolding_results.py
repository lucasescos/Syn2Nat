import pandas as pd
import numpy as np

# Load cofolding summary
df = pd.read_parquet("reports/esmfold2_cofolding_summary.parquet")
print(f"Total rows in summary: {len(df)}")
print(f"Columns: {df.columns.tolist()}")

# Status breakdown
print("\n--- Status Breakdown ---")
print(df["status"].value_counts(dropna=False))

# Filter successful cofolds
success_df = df[df["status"] == "success"].copy()
print(f"\nTotal successfully folded complexes: {len(success_df)}")

# Summary statistics of metrics
print("\n--- Metrics Summary (Successful Folds) ---")
for col in ["iptm", "ptm", "mean_plddt", "microprotein_plddt", "partner_plddt", "elapsed_seconds"]:
    if col in success_df.columns:
        s = success_df[col].dropna()
        print(f"  {col:20s}: Mean={s.mean():.3f}, Median={s.median():.3f}, Min={s.min():.3f}, Max={s.max():.3f}")

# Binding classification distribution
print("\n--- Binding Classification Breakdown ---")
print(success_df["binding_classification"].value_counts())

# Top hits by ipTM
print("\n--- Top 25 Highest ipTM Complexes ---")
cols_to_show = [
    "microprotein_rank", "microprotein_gene", "microprotein_biotype", "microprotein_length",
    "partner_rank", "partner_gene_symbol", "partner_length", "esmfold_total_length",
    "iptm", "ptm", "microprotein_plddt", "partner_plddt", "interaction_category", "biological_rationale"
]
top_hits = success_df.sort_values("iptm", ascending=False).head(25)
for idx, (_, r) in enumerate(top_hits.iterrows(), 1):
    print(f"\nHit #{idx:02d}: {r['microprotein_gene']} (Rank {r['microprotein_rank']}, {r['microprotein_biotype']}, {r['microprotein_length']} aa) + {r['partner_gene_symbol']} ({r['partner_length']} aa)")
    print(f"  ipTM: {r['iptm']:.4f} | pTM: {r['ptm']:.4f} | uProt pLDDT: {r['microprotein_plddt']:.2f} | Partner pLDDT: {r['partner_plddt']:.2f} | Total: {r['esmfold_total_length']} aa")
    print(f"  Category: {r['interaction_category']} | Rationale: {r['biological_rationale'][:100]}...")

# Top hits by induced folding (microprotein pLDDT >= 0.70)
print("\n--- Top Complexes with Induced Structural Folding of Microprotein (uProt pLDDT >= 0.70) ---")
structured_uprots = success_df[success_df["microprotein_plddt"] >= 0.70].sort_values(["iptm", "microprotein_plddt"], ascending=False)
print(f"Total complexes where microprotein adopts structured fold (pLDDT >= 70): {len(structured_uprots)}")
for idx, (_, r) in enumerate(structured_uprots.head(10).iterrows(), 1):
    print(f"  #{idx:02d}: {r['microprotein_gene']} + {r['partner_gene_symbol']} | ipTM: {r['iptm']:.3f}, uProt pLDDT: {r['microprotein_plddt']:.2f}, Partner pLDDT: {r['partner_plddt']:.2f}")

# Group by microprotein: which microproteins have the highest max ipTM?
print("\n--- Top Microproteins by Best Interface ipTM ---")
best_per_micro = success_df.groupby("microprotein_gene").apply(lambda g: g.loc[g["iptm"].idxmax()]).reset_index(drop=True)
best_per_micro = best_per_micro.sort_values("iptm", ascending=False)
for idx, (_, r) in enumerate(best_per_micro.head(15).iterrows(), 1):
    print(f"  #{idx:02d}: {r['microprotein_gene']:15s} (Rank {r['microprotein_rank']:2d}, {r['microprotein_biotype']:10s}) Best Partner: {r['partner_gene_symbol']:10s} | ipTM: {r['iptm']:.4f} | pTM: {r['ptm']:.4f}")
