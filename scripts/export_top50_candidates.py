#!/usr/bin/env python3
"""
Generate Curated Top 50 Human Microprotein Candidates Portfolio.
Balanced across:
  - Cohort 1: Autonomous Tier 1 Mass-Spec Peptideins (15)
  - Cohort 2: Autonomous lncRNA-ORFs (15)
  - Cohort 3: Autonomous Regulatory 5' UTR uORFs (12)
  - Cohort 4: Dual-Coding Tier 1 Peptideins (8)
"""

from pathlib import Path
import pandas as pd

BASE_DIR = Path(__file__).resolve().parent.parent
DATA_DIR = BASE_DIR / "data"
REPORTS_DIR = BASE_DIR / "reports"

df_avi = pd.read_parquet(DATA_DIR / "screening" / "microprotein_avi_scores.parquet")
df_tiers = pd.read_parquet(DATA_DIR / "parsed" / "orf_tiers.parquet")
df_dis = pd.read_parquet(DATA_DIR / "screening" / "microprotein_single_variant_disruptions.parquet")

tier_cols = [c for c in ["peptide_support", "hla_detected", "nonhla_detected", "n_hla_peptides", "n_nonhla_peptides"] if c in df_tiers.columns]
df = df_avi.merge(df_tiers[["orf_id"] + tier_cols], on="orf_id", how="left")
dis_cols = [c for c in ["primary_mechanism", "is_high_impact_hotspot", "convergent_modalities_count"] if c in df_dis.columns]
df = df.merge(df_dis[["orf_id"] + dis_cols], on="orf_id", how="left")

# Cohort 1: Autonomous Tier 1 Peptideins (15)
c1 = df[
    df["is_peptidein"] &
    df["final_tier"].str.startswith("1") &
    df["orf_biotype"].isin(["lncRNA-ORF", "uORF", "dORF"])
].sort_values("median_avi_phred", ascending=False).head(15).copy()
c1["cohort"] = "1. Autonomous Tier 1 Peptidein"

# Cohort 2: Autonomous lncRNA-ORFs (15) - exclude c1
c2 = df[
    (df["orf_biotype"] == "lncRNA-ORF") &
    (~df["orf_id"].isin(c1["orf_id"])) &
    (df["median_avi_phred"] >= 15.0)
].sort_values("max_avi_phred", ascending=False).head(15).copy()
c2["cohort"] = "2. Autonomous lncRNA-ORF"

# Cohort 3: Autonomous 5' UTR uORFs (12) - exclude c1
c3 = df[
    (df["orf_biotype"] == "uORF") &
    (~df["orf_id"].isin(c1["orf_id"])) &
    (df["median_avi_phred"] >= 18.0)
].sort_values("max_avi_phred", ascending=False).head(12).copy()
c3["cohort"] = "3. Autonomous Regulatory uORF"

# Cohort 4: CDS-Overlapping Tier 1 Peptideins (8)
c4 = df[
    df["is_peptidein"] &
    df["final_tier"].str.startswith("1") &
    (~df["orf_biotype"].isin(["lncRNA-ORF", "uORF", "dORF"]))
].sort_values("median_avi_phred", ascending=False).head(8).copy()
c4["cohort"] = "4. Dual-Coding Tier 1 Peptidein"

top50 = pd.concat([c1, c2, c3, c4], ignore_index=True)
top50["rank"] = range(1, len(top50) + 1)

cols_export = [
    "rank", "cohort", "orf_id", "gene", "orf_biotype", "length_aa", "final_tier",
    "is_peptidein", "peptide_support", "median_avi_phred", "max_avi_phred",
    "pct_phred_ge_20", "peak_variant", "primary_mechanism",
    "convergent_modalities_count", "atlas_url"
]
top50_out = top50[cols_export]

top50_out.to_csv(REPORTS_DIR / "top50_curated_candidates.tsv", sep="\t", index=False)
top50_out.to_parquet(REPORTS_DIR / "top50_curated_candidates.parquet", index=False)
print(f"Successfully exported {len(top50_out)} curated candidates to {REPORTS_DIR / 'top50_curated_candidates.tsv'}")
