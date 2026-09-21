#!/usr/bin/env python3
"""
Multi-Omics Comparative Analysis of the Human Microprotein Catalog (7,264 smORFs).
Integrates:
  - orfs.parquet & orf_tiers.parquet (Biotypes, mass-spec tiers)
  - microprotein_intervals_all_7264.parquet (Genomic intervals, host genes, coordinates)
  - microprotein_avi_scores.parquet (AlphaGenome Variant Impact constraint scores)
  - microprotein_single_variant_disruptions.parquet (Differential functional disruptions across 10 scorers)
  - microprotein_single_variant_significant_tracks.parquet (Tissue/biosample track impacts across 54 GTEx tissues)

Generates:
  - Statistical hypothesis tests (Kruskal-Wallis, Mann-Whitney U, Chi-Square, Spearman)
  - 4 publication-grade figures (300 DPI)
  - Top 20 candidate outlier table (Median Phred >= 20 & High-Impact Hotspot)
  - Summary metrics JSON & TSV outputs
"""

import json
import os
import shutil
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np
import pandas as pd
from scipy import stats
import seaborn as sns

# Set style
plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.sans-serif"] = ["DejaVu Sans", "Arial", "Helvetica"]
plt.rcParams["axes.edgecolor"] = "#2c3e50"
plt.rcParams["axes.linewidth"] = 1.0
plt.rcParams["xtick.color"] = "#2c3e50"
plt.rcParams["ytick.color"] = "#2c3e50"

BASE_DIR = Path(__file__).resolve().parent.parent
DATA_DIR = BASE_DIR / "data"
REPORTS_DIR = BASE_DIR / "reports"
FIGURES_DIR = REPORTS_DIR / "figures"
ARTIFACTS_DIR = Path(r"C:\Users\Lucas.Escosteguy\.gemini\antigravity\brain\b4954ee9-d996-4b97-9d18-0ab5359d6b32")
ARTIFACTS_FIG_DIR = ARTIFACTS_DIR / "figures"

REPORTS_DIR.mkdir(parents=True, exist_ok=True)
FIGURES_DIR.mkdir(parents=True, exist_ok=True)
ARTIFACTS_FIG_DIR.mkdir(parents=True, exist_ok=True)


def load_and_merge_data():
    print("[1/5] Loading and merging multi-omics datasets...")
    df_orfs = pd.read_parquet(DATA_DIR / "parsed" / "orfs.parquet")
    df_tiers = pd.read_parquet(DATA_DIR / "parsed" / "orf_tiers.parquet")
    df_intervals = pd.read_parquet(DATA_DIR / "curated" / "microprotein_intervals_all_7264.parquet")
    df_avi = pd.read_parquet(DATA_DIR / "screening" / "microprotein_avi_scores.parquet")
    df_disrupt = pd.read_parquet(DATA_DIR / "screening" / "microprotein_single_variant_disruptions.parquet")

    print(f"  - orfs: {df_orfs.shape}")
    print(f"  - tiers: {df_tiers.shape}")
    print(f"  - intervals: {df_intervals.shape}")
    print(f"  - avi scores: {df_avi.shape}")
    print(f"  - disruptions: {df_disrupt.shape}")

    # Merge core datasets (7,264 rows)
    df_core = df_avi.copy()

    # Add extra metadata from df_tiers & df_orfs
    tier_meta_cols = [c for c in ["peptide_support", "annotation_decision_final", "hla_detected", "nonhla_detected", "n_hla_peptides", "n_nonhla_peptides", "hpp_category"] if c in df_tiers.columns]
    df_core = df_core.merge(df_tiers[["orf_id"] + tier_meta_cols], on="orf_id", how="left")

    orf_meta_cols = [c for c in ["ORBLv", "ORBLq"] if c in df_orfs.columns]
    if orf_meta_cols:
        df_core = df_core.merge(df_orfs[["orf_id"] + orf_meta_cols], on="orf_id", how="left")

    # Add disruption metadata (6,964 rows)
    disrupt_cols = [c for c in df_disrupt.columns if c not in df_core.columns or c == "orf_id"]
    df_full = df_core.merge(df_disrupt[disrupt_cols], on="orf_id", how="left")

    print(f"  -> Master merged dataset: {df_full.shape[0]} microproteins, {df_full.shape[1]} columns")
    return df_full


def perform_statistical_analyses(df):
    print("\n[2/5] Performing statistical hypothesis testing...")
    results = {}

    # 1. Biotype Stratification
    biotypes = [b for b in df["orf_biotype"].value_counts().index if len(df[df["orf_biotype"] == b]) >= 20]
    groups_median = [df[df["orf_biotype"] == b]["median_avi_phred"].dropna() for b in biotypes]
    groups_max = [df[df["orf_biotype"] == b]["max_avi_phred"].dropna() for b in biotypes]
    groups_pct20 = [df[df["orf_biotype"] == b]["pct_phred_ge_20"].dropna() for b in biotypes]

    # Kruskal-Wallis tests
    kw_median = stats.kruskal(*groups_median)
    kw_max = stats.kruskal(*groups_max)
    kw_pct20 = stats.kruskal(*groups_pct20)

    results["biotype_kruskal_wallis"] = {
        "median_avi_phred": {"stat": float(kw_median.statistic), "p_value": float(kw_median.pvalue)},
        "max_avi_phred": {"stat": float(kw_max.statistic), "p_value": float(kw_max.pvalue)},
        "pct_phred_ge_20": {"stat": float(kw_pct20.statistic), "p_value": float(kw_pct20.pvalue)},
        "biotypes_tested": biotypes,
    }

    # Pairwise Mann-Whitney U: lncRNA-ORF vs uORF
    lncrna_med = df[df["orf_biotype"] == "lncRNA-ORF"]["median_avi_phred"].dropna()
    uorf_med = df[df["orf_biotype"] == "uORF"]["median_avi_phred"].dropna()
    mwu_lncrna_vs_uorf = stats.mannwhitneyu(lncrna_med, uorf_med, alternative="two-sided")

    # Pairwise: intORF vs lncRNA-ORF
    intorf_med = df[df["orf_biotype"] == "intORF"]["median_avi_phred"].dropna()
    mwu_intorf_vs_lncrna = stats.mannwhitneyu(intorf_med, lncrna_med, alternative="two-sided")

    # Pairwise: intORF vs uORF
    mwu_intorf_vs_uorf = stats.mannwhitneyu(intorf_med, uorf_med, alternative="two-sided")

    results["biotype_pairwise_mwu"] = {
        "lncRNA_vs_uORF": {
            "u_stat": float(mwu_lncrna_vs_uorf.statistic),
            "p_value": float(mwu_lncrna_vs_uorf.pvalue),
            "median_lncrna": float(lncrna_med.median()),
            "median_uorf": float(uorf_med.median()),
            "mean_lncrna": float(lncrna_med.mean()),
            "mean_uorf": float(uorf_med.mean()),
        },
        "intORF_vs_lncRNA": {
            "u_stat": float(mwu_intorf_vs_lncrna.statistic),
            "p_value": float(mwu_intorf_vs_lncrna.pvalue),
            "median_intorf": float(intorf_med.median()),
            "median_lncrna": float(lncrna_med.median()),
        },
        "intORF_vs_uORF": {
            "u_stat": float(mwu_intorf_vs_uorf.statistic),
            "p_value": float(mwu_intorf_vs_uorf.pvalue),
            "median_intorf": float(intorf_med.median()),
            "median_uorf": float(uorf_med.median()),
        }
    }

    # Biotype breakdown summary table
    biotype_summary = []
    for b in df["orf_biotype"].unique():
        sub = df[df["orf_biotype"] == b]
        biotype_summary.append({
            "biotype": b,
            "count": len(sub),
            "pct_catalog": len(sub) / len(df) * 100,
            "median_avi_phred_median": sub["median_avi_phred"].median(),
            "median_avi_phred_mean": sub["median_avi_phred"].mean(),
            "max_avi_phred_median": sub["max_avi_phred"].median(),
            "pct_ge_15_mean": sub["pct_phred_ge_15"].mean(),
            "pct_ge_20_mean": sub["pct_phred_ge_20"].mean(),
            "count_ge_20": (sub["median_avi_phred"] >= 20).sum(),
            "pct_with_median_ge_20": (sub["median_avi_phred"] >= 20).mean() * 100,
        })
    df_biotype_summary = pd.DataFrame(biotype_summary).sort_values("count", ascending=False)
    results["biotype_summary"] = df_biotype_summary.to_dict(orient="records")

    # Primary disruption mechanism by biotype contingency test (Chi-Square)
    df_disrupt_valid = df[df["primary_mechanism"].notna() & (df["primary_mechanism"] != "")]
    contingency = pd.crosstab(df_disrupt_valid["orf_biotype"], df_disrupt_valid["primary_mechanism"])
    chi2_res = stats.chi2_contingency(contingency)
    results["mechanism_by_biotype_chi2"] = {
        "chi2_stat": float(chi2_res.statistic),
        "p_value": float(chi2_res.pvalue),
        "dof": int(chi2_res.dof),
    }

    # 2. Experimental Tier Validation: Tier 1 Peptideins (84) vs Remainder
    t1_pep = df[(df["is_peptidein"] == True) & (df["final_tier"].str.startswith("1", na=False))]
    non_t1_pep = df[~df["orf_id"].isin(t1_pep["orf_id"])]

    mwu_t1_median = stats.mannwhitneyu(t1_pep["median_avi_phred"].dropna(), non_t1_pep["median_avi_phred"].dropna(), alternative="two-sided")
    mwu_t1_max = stats.mannwhitneyu(t1_pep["max_avi_phred"].dropna(), non_t1_pep["max_avi_phred"].dropna(), alternative="two-sided")
    mwu_t1_pct20 = stats.mannwhitneyu(t1_pep["pct_phred_ge_20"].dropna(), non_t1_pep["pct_phred_ge_20"].dropna(), alternative="two-sided")

    # High-impact hotspot proportion test
    t1_hotspots = t1_pep["is_high_impact_hotspot"].sum()
    t1_total = t1_pep["is_high_impact_hotspot"].notna().sum()
    nont1_hotspots = non_t1_pep["is_high_impact_hotspot"].sum()
    nont1_total = non_t1_pep["is_high_impact_hotspot"].notna().sum()

    table_hotspot = [[t1_hotspots, t1_total - t1_hotspots], [nont1_hotspots, nont1_total - nont1_hotspots]]
    fisher_hotspot = stats.fisher_exact(table_hotspot)

    results["tier1_peptidein_validation"] = {
        "n_tier1_peptideins": len(t1_pep),
        "n_remainder": len(non_t1_pep),
        "tier1_median_phred_median": float(t1_pep["median_avi_phred"].median()),
        "remainder_median_phred_median": float(non_t1_pep["median_avi_phred"].median()),
        "mwu_median_phred": {"u_stat": float(mwu_t1_median.statistic), "p_value": float(mwu_t1_median.pvalue)},
        "tier1_max_phred_median": float(t1_pep["max_avi_phred"].median()),
        "remainder_max_phred_median": float(non_t1_pep["max_avi_phred"].median()),
        "mwu_max_phred": {"u_stat": float(mwu_t1_max.statistic), "p_value": float(mwu_t1_max.pvalue)},
        "tier1_pct_ge_20_mean": float(t1_pep["pct_phred_ge_20"].mean()),
        "remainder_pct_ge_20_mean": float(non_t1_pep["pct_phred_ge_20"].mean()),
        "mwu_pct_ge_20": {"u_stat": float(mwu_t1_pct20.statistic), "p_value": float(mwu_t1_pct20.pvalue)},
        "tier1_hotspot_rate": float(t1_hotspots / t1_total) if t1_total > 0 else 0,
        "remainder_hotspot_rate": float(nont1_hotspots / nont1_total) if nont1_total > 0 else 0,
        "fisher_hotspot": {"odds_ratio": float(fisher_hotspot.statistic), "p_value": float(fisher_hotspot.pvalue)},
    }

    # 3. Length vs Constraint Correlation
    spearman_len_med = stats.spearmanr(df["length_aa"], df["median_avi_phred"])
    pearson_len_med = stats.pearsonr(df["length_aa"], df["median_avi_phred"])
    spearman_len_max = stats.spearmanr(df["length_aa"], df["max_avi_phred"])
    spearman_len_nvar = stats.spearmanr(df["length_aa"], df["n_variants_scored"])

    results["length_correlation"] = {
        "spearman_length_vs_median_phred": {"rho": float(spearman_len_med.statistic), "p_value": float(spearman_len_med.pvalue)},
        "pearson_length_vs_median_phred": {"r": float(pearson_len_med.statistic), "p_value": float(pearson_len_med.pvalue)},
        "spearman_length_vs_max_phred": {"rho": float(spearman_len_max.statistic), "p_value": float(spearman_len_max.pvalue)},
        "spearman_length_vs_n_variants": {"rho": float(spearman_len_nvar.statistic), "p_value": float(spearman_len_nvar.pvalue)},
        "median_length_aa": float(df["length_aa"].median()),
        "iqr_length_aa": [float(df["length_aa"].quantile(0.25)), float(df["length_aa"].quantile(0.75))],
    }

    print("  -> Hypothesis testing completed successfully.")
    return results, df_biotype_summary


def analyze_tissue_tracks():
    print("\n[3/5] Querying and aggregating tissue & track impacts (9.9M records)...")
    tracks_file = DATA_DIR / "screening" / "microprotein_single_variant_significant_tracks.parquet"
    
    # Read relevant columns
    cols = ["orf_id", "output_type", "gtex_tissue", "raw_score", "quantile_score", "transcription_factor", "histone_mark"]
    df_tracks = pd.read_parquet(tracks_file, columns=cols)
    
    # Filter GTEx tracks
    df_gtex = df_tracks[df_tracks["gtex_tissue"] != ""].copy()
    print(f"  - Total GTEx significant tracks: {len(df_gtex):,} across {df_gtex['gtex_tissue'].nunique()} tissues")

    # Aggregate by tissue and modality
    gtex_tissue_modality = pd.crosstab(df_gtex["gtex_tissue"], df_gtex["output_type"])
    gtex_tissue_modality["Total"] = gtex_tissue_modality.sum(axis=1)
    gtex_top15 = gtex_tissue_modality.sort_values("Total", ascending=False).head(15)

    # Unique ORFs per tissue
    orfs_per_tissue = df_gtex.groupby("gtex_tissue")["orf_id"].nunique().sort_values(ascending=False)

    # ATAC and DNASE cell lines
    df_cell = df_tracks[df_tracks["gtex_tissue"] == ""].copy()
    tf_counts = df_cell["transcription_factor"].value_counts().head(10)
    histone_counts = df_cell["histone_mark"].value_counts().head(10)

    tissue_results = {
        "top_15_gtex_summary": gtex_top15.reset_index().to_dict(orient="records"),
        "top_tissues_unique_orfs": orfs_per_tissue.head(15).to_dict(),
        "top_transcription_factors": tf_counts.to_dict(),
        "top_histone_marks": histone_counts.to_dict(),
    }
    
    return gtex_top15, tissue_results


def generate_figures(df, gtex_top15, stats_results):
    print("\n[4/5] Generating publication-quality figures (300 DPI)...")
    sns.set_theme(style="whitegrid", font="sans-serif")

    # -------------------------------------------------------------
    # FIGURE 1: Biotype Selection Distribution (Violin / Box Plot)
    # -------------------------------------------------------------
    fig, axes = plt.subplots(1, 2, figsize=(16, 7), sharey=False)
    
    biotype_order = (
        df.groupby("orf_biotype")["median_avi_phred"]
        .median()
        .sort_values(ascending=False)
        .index.tolist()
    )

    # 1A: Median AVI Phred (Negative selection / coding constraint)
    sns.violinplot(
        data=df,
        x="orf_biotype",
        y="median_avi_phred",
        order=biotype_order,
        ax=axes[0],
        palette="viridis",
        inner="box",
        cut=0,
        linewidth=1.2
    )
    axes[0].axhline(20, color="#d90429", linestyle="--", linewidth=1.5, label="Top 1% Constraint (Phred ≥ 20)")
    axes[0].axhline(15, color="#f77f00", linestyle=":", linewidth=1.5, label="Top 3% Constraint (Phred ≥ 15)")
    axes[0].set_title("A. Background Purifying Selection (Median AVI Phred)", fontsize=14, fontweight="bold", pad=12)
    axes[0].set_xlabel("smORF Biotype", fontsize=12, fontweight="bold")
    axes[0].set_ylabel("Median AVI Phred Score", fontsize=12, fontweight="bold")
    axes[0].set_xticklabels(axes[0].get_xticklabels(), rotation=30, ha="right", fontsize=10)
    axes[0].legend(loc="upper right", frameon=True, framealpha=0.9)
    axes[0].set_ylim(-1, 35)

    for idx, b in enumerate(biotype_order):
        sub = df[df["orf_biotype"] == b]
        pct20 = (sub["median_avi_phred"] >= 20).mean() * 100
        axes[0].text(
            idx,
            sub["median_avi_phred"].quantile(0.95) + 1.2,
            f"n={len(sub)}\n{pct20:.1f}%",
            ha="center",
            va="bottom",
            fontsize=8.5,
            fontweight="bold",
            color="#1d3557"
        )

    # 1B: Peak Hotspot Phred
    sns.boxplot(
        data=df,
        x="orf_biotype",
        y="max_avi_phred",
        order=biotype_order,
        ax=axes[1],
        palette="viridis",
        fliersize=2.5,
        linewidth=1.2,
        boxprops=dict(alpha=0.85)
    )
    axes[1].axhline(30, color="#7209b7", linestyle="--", linewidth=1.5, label="Top 0.1% Hotspot (Phred ≥ 30)")
    axes[1].axhline(20, color="#d90429", linestyle=":", linewidth=1.5, label="Top 1% Hotspot (Phred ≥ 20)")
    axes[1].set_title("B. Peak Deleterious Hotspot (Max AVI Phred)", fontsize=14, fontweight="bold", pad=12)
    axes[1].set_xlabel("smORF Biotype", fontsize=12, fontweight="bold")
    axes[1].set_ylabel("Peak Hotspot AVI Phred Score", fontsize=12, fontweight="bold")
    axes[1].set_xticklabels(axes[1].get_xticklabels(), rotation=30, ha="right", fontsize=10)
    axes[1].legend(loc="lower right", frameon=True, framealpha=0.9)
    axes[1].set_ylim(0, 75)

    plt.tight_layout()
    fig1_path = FIGURES_DIR / "fig1_biotype_selection_distribution.png"
    plt.savefig(fig1_path, dpi=300)
    plt.close()
    print(f"  [+] Saved Figure 1: {fig1_path}")

    # -------------------------------------------------------------
    # FIGURE 2: Disruption Mechanism Breakdown (Stacked Bar Chart)
    # -------------------------------------------------------------
    fig, ax = plt.subplots(figsize=(14, 7))
    df_mech = df[df["primary_mechanism"].notna() & (df["primary_mechanism"] != "")].copy()
    
    t1_pep = df_mech[(df_mech["is_peptidein"] == True) & (df_mech["final_tier"].str.startswith("1", na=False))].copy()
    t1_pep["cohort"] = f"Tier 1 Peptideins (n={len(t1_pep)})"
    
    df_mech["cohort"] = df_mech["orf_biotype"] + " (n=" + df_mech["orf_biotype"].map(df_mech["orf_biotype"].value_counts()).astype(str) + ")"
    
    combined_mech = pd.concat([df_mech[["cohort", "primary_mechanism"]], t1_pep[["cohort", "primary_mechanism"]]], ignore_index=True)
    
    ct = pd.crosstab(combined_mech["cohort"], combined_mech["primary_mechanism"], normalize="index") * 100
    
    mech_cols_ordered = [
        "Severe Splicing Disruption",
        "Splicing Disruption",
        "Chromatin Opening",
        "Chromatin Closing / Repression",
        "TF Binding Motif Disruption",
        "Transcriptional Activation",
        "Transcriptional Silencing",
        "Subtle / Regulatory Shift",
        "Benign / Tolerated",
        "Epigenetic Mark Alteration"
    ]
    mech_cols_present = [c for c in mech_cols_ordered if c in ct.columns]
    ct = ct[mech_cols_present]
    
    row_order = [r for r in ct.index if "Tier 1" in r] + [r for r in ct.index if "Tier 1" not in r]
    ct = ct.loc[row_order]

    colors = [
        "#9d0208", # Severe Splicing Disruption
        "#e85d04", # Splicing Disruption
        "#2a9d8f", # Chromatin Opening
        "#264653", # Chromatin Closing
        "#e76f51", # TF Binding
        "#0077b6", # Transcriptional Activation
        "#03045e", # Transcriptional Silencing
        "#a8dadc", # Subtle / Regulatory Shift
        "#adb5bd", # Benign / Tolerated
        "#6a0572"  # Epigenetic Mark Alteration
    ]
    color_map = {m: c for m, c in zip(mech_cols_ordered, colors)}

    ct.plot(
        kind="barh",
        stacked=True,
        ax=ax,
        color=[color_map[c] for c in ct.columns],
        edgecolor="#222222",
        linewidth=0.7,
        width=0.72
    )

    ax.set_title("Primary Functional Disruption Mechanisms of Peak Hotspots Across smORF Biotypes", fontsize=14, fontweight="bold", pad=12)
    ax.set_xlabel("Percentage of Evaluated Candidates (%)", fontsize=12, fontweight="bold")
    ax.set_ylabel("Cohort / Biotype", fontsize=12, fontweight="bold")
    ax.set_xlim(0, 100)
    ax.xaxis.set_major_formatter(ticker.PercentFormatter(100))
    ax.legend(title="Primary Disruption Modality", bbox_to_anchor=(1.02, 1), loc="upper left", frameon=True)
    ax.grid(axis="x", linestyle="--", alpha=0.5)

    plt.tight_layout()
    fig2_path = FIGURES_DIR / "fig2_disruption_mechanism_breakdown.png"
    plt.savefig(fig2_path, dpi=300)
    plt.close()
    print(f"  [+] Saved Figure 2: {fig2_path}")

    # -------------------------------------------------------------
    # FIGURE 3: Constraint vs ORF Length Landscape
    # -------------------------------------------------------------
    fig, ax = plt.subplots(figsize=(12, 7))
    
    hb = ax.hexbin(
        df["length_aa"],
        df["median_avi_phred"],
        gridsize=45,
        cmap="Blues",
        mincnt=1,
        bins="log",
        edgecolors="none",
        alpha=0.75
    )
    cb = fig.colorbar(hb, ax=ax, orientation="vertical", pad=0.02)
    cb.set_label("Log10(Candidate Density)", fontsize=11, fontweight="bold")

    t1_subset = df[(df["is_peptidein"] == True) & (df["final_tier"].str.startswith("1", na=False))]
    ax.scatter(
        t1_subset["length_aa"],
        t1_subset["median_avi_phred"],
        color="#d90429",
        edgecolor="black",
        s=48,
        alpha=0.9,
        label=f"Tier 1 Peptideins (n={len(t1_subset)})",
        zorder=5
    )

    ax.axhline(20, color="#d90429", linestyle="--", linewidth=1.4, alpha=0.8, label="Phred ≥ 20 (Top 1%)")
    ax.axhline(15, color="#f77f00", linestyle=":", linewidth=1.4, alpha=0.8, label="Phred ≥ 15 (Top 3.2%)")

    top_t1 = t1_subset.sort_values("median_avi_phred", ascending=False).head(3)
    for _, row in top_t1.iterrows():
        ax.annotate(
            f"{row['gene']}\n({row['orf_id']})",
            xy=(row["length_aa"], row["median_avi_phred"]),
            xytext=(row["length_aa"] + 5, row["median_avi_phred"] + 1.2),
            arrowprops=dict(facecolor="#2b2d42", arrowstyle="->", lw=1.0),
            fontsize=9,
            fontweight="bold",
            color="#2b2d42",
            zorder=10
        )

    rho = stats_results["length_correlation"]["spearman_length_vs_median_phred"]["rho"]
    p_val = stats_results["length_correlation"]["spearman_length_vs_median_phred"]["p_value"]
    
    ax.text(
        0.03, 0.94,
        f"Spearman's $\\rho$ = {rho:.3f} (p = {p_val:.2e})\nPearson's r = {stats_results['length_correlation']['pearson_length_vs_median_phred']['r']:.3f}",
        transform=ax.transAxes,
        fontsize=11,
        fontweight="bold",
        verticalalignment="top",
        bbox=dict(boxstyle="round,pad=0.5", facecolor="#ffffff", edgecolor="#2c3e50", alpha=0.9)
    )

    ax.set_title("Human smORF Length (aa) vs. AlphaGenome Purifying Selection", fontsize=14, fontweight="bold", pad=12)
    ax.set_xlabel("smORF Peptide Length (Amino Acids)", fontsize=12, fontweight="bold")
    ax.set_ylabel("Median AVI Phred Score (Background Constraint)", fontsize=12, fontweight="bold")
    ax.set_xlim(0, 220)
    ax.set_ylim(0, 32)
    ax.legend(loc="lower right", frameon=True, framealpha=0.9)

    plt.tight_layout()
    fig3_path = FIGURES_DIR / "fig3_constraint_vs_length.png"
    plt.savefig(fig3_path, dpi=300)
    plt.close()
    print(f"  [+] Saved Figure 3: {fig3_path}")

    # -------------------------------------------------------------
    # FIGURE 4: GTEx Tissue Impact Matrix
    # -------------------------------------------------------------
    fig, ax = plt.subplots(figsize=(14, 8))
    
    gtex_matrix = gtex_top15[["SPLICE_JUNCTIONS", "SPLICE_SITE_USAGE", "RNA_SEQ"]].copy()
    clean_labels = [idx.replace("Cells_", "").replace("_", " ") for idx in gtex_matrix.index]
    gtex_matrix.index = clean_labels

    sns.heatmap(
        gtex_matrix,
        annot=True,
        fmt=",d",
        cmap="YlOrRd",
        linewidths=1.0,
        linecolor="#ffffff",
        cbar_kws={"label": "Significant Hotspot Disruption Track Count"},
        ax=ax
    )

    ax.set_title("Distribution of High-Impact Hotspot Disruptions Across Top 15 GTEx Tissues", fontsize=14, fontweight="bold", pad=12)
    ax.set_xlabel("AlphaGenome Modality Track", fontsize=12, fontweight="bold")
    ax.set_ylabel("GTEx Tissue Site / Cell Line", fontsize=12, fontweight="bold")
    ax.set_xticklabels(["Splice Junctions", "Splice Site Usage", "RNA Expression (Log2FC)"], fontsize=11, fontweight="bold")
    ax.set_yticklabels(ax.get_yticklabels(), fontsize=10, rotation=0)

    plt.tight_layout()
    fig4_path = FIGURES_DIR / "fig4_gtex_tissue_impact_matrix.png"
    plt.savefig(fig4_path, dpi=300)
    plt.close()
    print(f"  [+] Saved Figure 4: {fig4_path}")

    # Copy all figures to artifacts directory
    for f in [fig1_path, fig2_path, fig3_path, fig4_path]:
        shutil.copy(f, ARTIFACTS_FIG_DIR / f.name)
    print(f"  [+] Copied figures to artifacts directory: {ARTIFACTS_FIG_DIR}")


def build_top20_outliers_table(df):
    print("\n[5/5] Extracting ranked top 20 candidate outliers...")
    subset = df[(df["median_avi_phred"] >= 20.0) & (df["is_high_impact_hotspot"] == True)].copy()
    subset_sorted = subset.sort_values("max_avi_phred", ascending=False).reset_index(drop=True)
    top20 = subset_sorted.head(20).copy()

    top20["rank"] = range(1, len(top20) + 1)
    cols = [
        "rank", "orf_id", "gene", "orf_biotype", "length_aa", "final_tier", "is_peptidein",
        "median_avi_phred", "max_avi_phred", "peak_variant", "primary_mechanism",
        "convergent_modalities_count", "atlas_url"
    ]
    top20_out = top20[cols]

    top20_out.to_csv(REPORTS_DIR / "top20_outliers.tsv", sep="\t", index=False)
    top20_out.to_parquet(REPORTS_DIR / "top20_outliers.parquet", index=False)
    print(f"  [+] Saved Top 20 Outliers to {REPORTS_DIR / 'top20_outliers.tsv'}")

    return top20_out


def main():
    print("=" * 80)
    print("STARTING MULTI-OMICS EXPLORATORY ANALYSIS OF HUMAN MICROPROTEINS")
    print("=" * 80)

    df = load_and_merge_data()
    stats_results, df_biotype_summary = perform_statistical_analyses(df)
    gtex_top15, tissue_results = analyze_tissue_tracks()
    generate_figures(df, gtex_top15, stats_results)
    top20_df = build_top20_outliers_table(df)

    all_summary = {
        "statistical_tests": stats_results,
        "tissue_and_tracks": tissue_results,
        "total_catalog_size": len(df),
        "total_variants_scored": int(df["n_variants_scored"].sum()),
    }
    with open(REPORTS_DIR / "multiomics_statistical_summary.json", "w") as f:
        json.dump(all_summary, f, indent=2)

    print("\n" + "=" * 80)
    print("MULTI-OMICS PIPELINE COMPLETED SUCCESSFULLY!")
    print("=" * 80)


if __name__ == "__main__":
    main()
