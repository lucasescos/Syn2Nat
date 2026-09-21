#!/usr/bin/env python3
"""
De-confounded Multi-Omics Analysis of Human Microproteins.
Separates Autonomous smORFs (lncRNA-ORFs, uORFs, dORFs) from Canonical CDS-Overlapping smORFs (intORFs, uoORFs, doORFs, mixed)
to eliminate the 'Canonical CDS Hitchhiking' confounder and identify genuine microprotein-driven functional selection.
"""

import json
import os
import shutil
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np
import pandas as pd
from scipy import stats
import seaborn as sns

BASE_DIR = Path(__file__).resolve().parent.parent
DATA_DIR = BASE_DIR / "data"
REPORTS_DIR = BASE_DIR / "reports"
FIGURES_DIR = REPORTS_DIR / "figures"
ARTIFACTS_DIR = Path(r"C:\Users\Lucas.Escosteguy\.gemini\antigravity\brain\b4954ee9-d996-4b97-9d18-0ab5359d6b32")
ARTIFACTS_FIG_DIR = ARTIFACTS_DIR / "figures"

REPORTS_DIR.mkdir(parents=True, exist_ok=True)
FIGURES_DIR.mkdir(parents=True, exist_ok=True)
ARTIFACTS_FIG_DIR.mkdir(parents=True, exist_ok=True)


def load_dataset():
    print("[1/5] Loading and classifying datasets by CDS autonomy...")
    df_avi = pd.read_parquet(DATA_DIR / "screening" / "microprotein_avi_scores.parquet")
    df_tiers = pd.read_parquet(DATA_DIR / "parsed" / "orf_tiers.parquet")
    df_disrupt = pd.read_parquet(DATA_DIR / "screening" / "microprotein_single_variant_disruptions.parquet")

    # Core merge
    df = df_avi.copy()
    tier_cols = [c for c in ["peptide_support", "final_tier_raw", "hla_detected", "nonhla_detected"] if c in df_tiers.columns]
    df = df.merge(df_tiers[["orf_id"] + tier_cols], on="orf_id", how="left")

    disrupt_cols = [c for c in df_disrupt.columns if c not in df.columns or c == "orf_id"]
    df = df.merge(df_disrupt[disrupt_cols], on="orf_id", how="left")

    # Define Genomic Autonomy / CDS Overlap Status
    # Autonomous: lncRNA-ORF (no CDS on transcript), uORF (5' UTR non-coding), dORF (3' UTR non-coding)
    # CDS-Overlapping: intORF (internal CDS), uoORF (spans start codon into CDS), doORF (spans stop codon from CDS), mixed
    autonomous_biotypes = ["lncRNA-ORF", "uORF", "dORF"]
    df["autonomy_group"] = df["orf_biotype"].apply(
        lambda b: "Autonomous (Non-CDS)" if b in autonomous_biotypes else "CDS-Overlapping"
    )
    df["is_autonomous"] = df["orf_biotype"].isin(autonomous_biotypes)

    print(f"  - Total Microproteins: {len(df):,}")
    print(f"  - Autonomous smORFs: {(df['is_autonomous']).sum():,} ({df['is_autonomous'].mean()*100:.1f}%)")
    print(f"  - CDS-Overlapping smORFs: {(~df['is_autonomous']).sum():,} ({(~df['is_autonomous']).mean()*100:.1f}%)")
    return df


def compute_deconfounded_statistics(df):
    print("\n[2/5] Computing de-confounded statistical comparisons...")
    results = {}

    auto = df[df["is_autonomous"]].copy()
    overlap = df[~df["is_autonomous"]].copy()

    # 1. Autonomous vs CDS-Overlapping constraint tests
    mwu_median = stats.mannwhitneyu(overlap["median_avi_phred"].dropna(), auto["median_avi_phred"].dropna(), alternative="two-sided")
    mwu_max = stats.mannwhitneyu(overlap["max_avi_phred"].dropna(), auto["max_avi_phred"].dropna(), alternative="two-sided")
    mwu_pct20 = stats.mannwhitneyu(overlap["pct_phred_ge_20"].dropna(), auto["pct_phred_ge_20"].dropna(), alternative="two-sided")

    # Extreme constraint count (Phred >= 20)
    auto_ge20 = (auto["median_avi_phred"] >= 20).sum()
    overlap_ge20 = (overlap["median_avi_phred"] >= 20).sum()
    table_ge20 = [[overlap_ge20, len(overlap) - overlap_ge20], [auto_ge20, len(auto) - auto_ge20]]
    fisher_ge20 = stats.fisher_exact(table_ge20)

    # Moderate constraint count (Phred >= 15)
    auto_ge15 = (auto["median_avi_phred"] >= 15).sum()
    overlap_ge15 = (overlap["median_avi_phred"] >= 15).sum()

    results["hitchhiking_comparison"] = {
        "autonomous_count": len(auto),
        "overlapping_count": len(overlap),
        "autonomous_median_phred_median": float(auto["median_avi_phred"].median()),
        "autonomous_median_phred_mean": float(auto["median_avi_phred"].mean()),
        "overlapping_median_phred_median": float(overlap["median_avi_phred"].median()),
        "overlapping_median_phred_mean": float(overlap["median_avi_phred"].mean()),
        "mwu_median_phred": {"u_stat": float(mwu_median.statistic), "p_value": float(mwu_median.pvalue)},
        "mwu_max_phred": {"u_stat": float(mwu_max.statistic), "p_value": float(mwu_max.pvalue)},
        "mwu_pct_ge_20": {"u_stat": float(mwu_pct20.statistic), "p_value": float(mwu_pct20.pvalue)},
        "autonomous_pct_with_median_ge_20": float(auto_ge20 / len(auto) * 100),
        "overlapping_pct_with_median_ge_20": float(overlap_ge20 / len(overlap) * 100),
        "autonomous_count_ge_20": int(auto_ge20),
        "overlapping_count_ge_20": int(overlap_ge20),
        "autonomous_count_ge_15": int(auto_ge15),
        "overlapping_count_ge_15": int(overlap_ge15),
        "fisher_extreme_constraint": {"odds_ratio": float(fisher_ge20.statistic), "p_value": float(fisher_ge20.pvalue)},
    }

    # 2. lncRNA-ORF Pure Autonomous Benchmark
    lnc = df[df["orf_biotype"] == "lncRNA-ORF"].copy()
    lnc_ge20 = (lnc["median_avi_phred"] >= 20).sum()
    lnc_ge15 = (lnc["median_avi_phred"] >= 15).sum()
    lnc_ge10 = (lnc["median_avi_phred"] >= 10).sum()

    results["lncrna_orf_benchmark"] = {
        "total_lncrna_orfs": len(lnc),
        "median_phred_median": float(lnc["median_avi_phred"].median()),
        "median_phred_mean": float(lnc["median_avi_phred"].mean()),
        "max_phred_median": float(lnc["max_avi_phred"].median()),
        "count_ge_20": int(lnc_ge20),
        "pct_ge_20": float(lnc_ge20 / len(lnc) * 100),
        "count_ge_15": int(lnc_ge15),
        "pct_ge_15": float(lnc_ge15 / len(lnc) * 100),
        "count_ge_10": int(lnc_ge10),
        "pct_ge_10": float(lnc_ge10 / len(lnc) * 100),
    }

    # 3. Mechanism shift Chi-Square test
    df_mech = df[df["primary_mechanism"].notna() & (df["primary_mechanism"] != "")]
    ct_mech = pd.crosstab(df_mech["is_autonomous"], df_mech["primary_mechanism"])
    chi2_mech = stats.chi2_contingency(ct_mech)

    results["mechanism_shift_chi2"] = {
        "chi2_stat": float(chi2_mech.statistic),
        "p_value": float(chi2_mech.pvalue),
        "dof": int(chi2_mech.dof),
    }

    # 4. Tier 1 Peptideins stratified by Autonomy
    t1 = df[df["is_peptidein"] & df["final_tier"].str.startswith("1")].copy()
    t1_auto = t1[t1["is_autonomous"]]
    t1_over = t1[~t1["is_autonomous"]]

    mwu_t1_auto_vs_rest = stats.mannwhitneyu(
        t1_auto["median_avi_phred"].dropna(),
        auto[~auto["orf_id"].isin(t1_auto["orf_id"])]["median_avi_phred"].dropna(),
        alternative="two-sided"
    )

    results["tier1_peptideins_deconfounded"] = {
        "total_tier1": len(t1),
        "autonomous_tier1_count": len(t1_auto),
        "overlapping_tier1_count": len(t1_over),
        "autonomous_tier1_median_phred": float(t1_auto["median_avi_phred"].median()),
        "autonomous_tier1_mean_phred": float(t1_auto["median_avi_phred"].mean()),
        "autonomous_tier1_ge_15_count": int((t1_auto["median_avi_phred"] >= 15).sum()),
        "mwu_t1_auto_vs_autonomous_remainder": {
            "u_stat": float(mwu_t1_auto_vs_rest.statistic),
            "p_value": float(mwu_t1_auto_vs_rest.pvalue)
        }
    }

    print("  -> Deconfounded statistical tests complete.")
    return results


def generate_deconfounded_figures(df, stats_results):
    print("\n[3/5] Generating de-confounded publication figures (300 DPI)...")
    sns.set_theme(style="whitegrid", font="sans-serif")

    # -------------------------------------------------------------------------
    # FIGURE 5: Autonomous vs. CDS-Overlapping Constraint Shift (Violin & CDF)
    # -------------------------------------------------------------------------
    fig, axes = plt.subplots(1, 2, figsize=(16, 6.5))

    # 5A: Grouped violin plot by autonomy group and detailed biotype
    biotype_order = ["lncRNA-ORF", "uORF", "dORF", "intORF", "uoORF", "mixed", "doORF"]
    palette = {
        "lncRNA-ORF": "#457b9d", "uORF": "#1d3557", "dORF": "#a8dadc",
        "intORF": "#e63946", "uoORF": "#f77f00", "mixed": "#d62828", "doORF": "#fcbf49"
    }

    sns.violinplot(
        data=df,
        x="orf_biotype",
        y="median_avi_phred",
        order=biotype_order,
        palette=palette,
        inner="quartile",
        cut=0,
        linewidth=1.2,
        ax=axes[0]
    )
    axes[0].axvline(2.5, color="#2c3e50", linestyle="-", linewidth=2.0, alpha=0.7)
    axes[0].text(1.0, 31, "AUTONOMOUS smORFs\n(Unconfounded, n=5,387)", ha="center", fontsize=11, fontweight="bold", color="#1d3557", bbox=dict(boxstyle="round,pad=0.3", facecolor="#f1faee", edgecolor="#1d3557", alpha=0.9))
    axes[0].text(4.5, 31, "CDS-OVERLAPPING smORFs\n(Host Hitchhiking, n=1,877)", ha="center", fontsize=11, fontweight="bold", color="#e63946", bbox=dict(boxstyle="round,pad=0.3", facecolor="#ffebee", edgecolor="#e63946", alpha=0.9))

    axes[0].axhline(20, color="#d90429", linestyle="--", linewidth=1.3, alpha=0.8, label="Phred ≥ 20 (Top 1%)")
    axes[0].axhline(15, color="#f77f00", linestyle=":", linewidth=1.3, alpha=0.8, label="Phred ≥ 15 (Top 3.2%)")

    axes[0].set_title("A. Purifying Selection by Genomic CDS Autonomy", fontsize=13, fontweight="bold", pad=10)
    axes[0].set_xlabel("smORF Biotype", fontsize=11, fontweight="bold")
    axes[0].set_ylabel("Median AVI Phred Score", fontsize=11, fontweight="bold")
    axes[0].set_ylim(-1, 35)
    axes[0].legend(loc="upper right", frameon=True)

    # 5B: Cumulative Distribution Function (CDF)
    auto_phreds = np.sort(df[df["is_autonomous"]]["median_avi_phred"].dropna())
    overlap_phreds = np.sort(df[~df["is_autonomous"]]["median_avi_phred"].dropna())
    lnc_phreds = np.sort(df[df["orf_biotype"] == "lncRNA-ORF"]["median_avi_phred"].dropna())

    axes[1].plot(auto_phreds, np.linspace(0, 1, len(auto_phreds)), label="All Autonomous smORFs (n=5,387)", color="#1d3557", linewidth=2.5)
    axes[1].plot(lnc_phreds, np.linspace(0, 1, len(lnc_phreds)), label="lncRNA-ORFs (Pure Autonomy, n=2,012)", color="#457b9d", linewidth=2.0, linestyle="--")
    axes[1].plot(overlap_phreds, np.linspace(0, 1, len(overlap_phreds)), label="CDS-Overlapping (n=1,877)", color="#e63946", linewidth=2.5)

    axes[1].axvline(20, color="#d90429", linestyle="--", linewidth=1.3, alpha=0.8)
    axes[1].axvline(15, color="#f77f00", linestyle=":", linewidth=1.3, alpha=0.8)

    odds_r = stats_results["hitchhiking_comparison"]["fisher_extreme_constraint"]["odds_ratio"]
    axes[1].text(0.48, 0.25, f"Odds Ratio for Phred ≥ 20:\nOR = {odds_r:.1f}x in Overlapping\n(p < 1e-300, Fisher's Exact)", transform=axes[1].transAxes, fontsize=10, fontweight="bold", bbox=dict(boxstyle="round,pad=0.4", facecolor="#ffffff", edgecolor="#2c3e50", alpha=0.9))

    axes[1].set_title("B. Empirical Cumulative Distribution Functions (CDF)", fontsize=13, fontweight="bold", pad=10)
    axes[1].set_xlabel("Median AVI Phred Score", fontsize=11, fontweight="bold")
    axes[1].set_ylabel("Cumulative Proportion", fontsize=11, fontweight="bold")
    axes[1].legend(loc="lower right", frameon=True)
    axes[1].set_xlim(0, 30)
    axes[1].set_ylim(0, 1.02)

    plt.tight_layout()
    fig5_path = FIGURES_DIR / "fig5_autonomous_vs_overlapping_selection.png"
    plt.savefig(fig5_path, dpi=300)
    plt.close()
    print(f"  [+] Saved Figure 5: {fig5_path}")

    # -------------------------------------------------------------------------
    # FIGURE 6: Autonomous vs CDS-Overlapping Disruption Modality Divergence
    # -------------------------------------------------------------------------
    fig, ax = plt.subplots(figsize=(13, 6.5))
    df_mech = df[df["primary_mechanism"].notna() & (df["primary_mechanism"] != "")].copy()

    df_mech["display_group"] = df_mech["orf_biotype"].map({
        "lncRNA-ORF": "lncRNA-ORFs (Pure Autonomous, n=1,927)",
        "uORF": "uORFs (5' UTR Autonomous, n=2,801)",
        "dORF": "dORFs (3' UTR Autonomous, n=439)",
        "intORF": "intORFs (Canonical CDS, n=687)",
        "uoORF": "uoORFs (CDS-Overlapping, n=611)",
        "mixed": "mixed (CDS-Overlapping, n=426)"
    })
    df_mech_valid = df_mech[df_mech["display_group"].notna()]

    ct = pd.crosstab(df_mech_valid["display_group"], df_mech_valid["primary_mechanism"], normalize="index") * 100

    order_groups = [
        "lncRNA-ORFs (Pure Autonomous, n=1,927)",
        "uORFs (5' UTR Autonomous, n=2,801)",
        "dORFs (3' UTR Autonomous, n=439)",
        "intORFs (Canonical CDS, n=687)",
        "uoORFs (CDS-Overlapping, n=611)",
        "mixed (CDS-Overlapping, n=426)"
    ]
    ct = ct.loc[order_groups]

    colors = {
        "Severe Splicing Disruption": "#9d0208",
        "Splicing Disruption": "#e85d04",
        "Chromatin Opening": "#2a9d8f",
        "Chromatin Closing / Repression": "#264653",
        "TF Binding Motif Disruption": "#e76f51",
        "Transcriptional Activation": "#0077b6",
        "Transcriptional Silencing": "#03045e",
        "Subtle / Regulatory Shift": "#a8dadc",
        "Benign / Tolerated": "#adb5bd",
        "Epigenetic Mark Alteration": "#6a0572"
    }
    cols_order = [c for c in colors.keys() if c in ct.columns]
    ct = ct[cols_order]

    ct.plot(
        kind="barh",
        stacked=True,
        ax=ax,
        color=[colors[c] for c in ct.columns],
        edgecolor="#222222",
        linewidth=0.7,
        width=0.70
    )

    ax.set_title("Shift in Functional Disruption Mechanisms: Autonomous vs. CDS-Overlapping smORFs", fontsize=13, fontweight="bold", pad=12)
    ax.set_xlabel("Proportion of Evaluated Hotspots (%)", fontsize=11, fontweight="bold")
    ax.set_ylabel("Cohort & Genomic Context", fontsize=11, fontweight="bold")
    ax.set_xlim(0, 100)
    ax.xaxis.set_major_formatter(ticker.PercentFormatter(100))
    ax.legend(title="Primary Disruption Modality", bbox_to_anchor=(1.02, 1), loc="upper left", frameon=True)
    ax.grid(axis="x", linestyle="--", alpha=0.5)

    plt.tight_layout()
    fig6_path = FIGURES_DIR / "fig6_autonomous_disruption_mechanisms.png"
    plt.savefig(fig6_path, dpi=300)
    plt.close()
    print(f"  [+] Saved Figure 6: {fig6_path}")

    # -------------------------------------------------------------------------
    # FIGURE 7: The Unconfounded Autonomous Standout Landscape
    # -------------------------------------------------------------------------
    fig, ax = plt.subplots(figsize=(13, 7.5))

    auto_df = df[df["is_autonomous"]].copy()

    # Hexbin of Autonomous background
    hb = ax.hexbin(
        auto_df["length_aa"],
        auto_df["median_avi_phred"],
        gridsize=40,
        cmap="YlGnBu",
        mincnt=1,
        bins="log",
        edgecolors="none",
        alpha=0.65
    )
    cb = fig.colorbar(hb, ax=ax, orientation="vertical", pad=0.02)
    cb.set_label("Log10(Autonomous smORF Density)", fontsize=11, fontweight="bold")

    # Overlay Autonomous Tier 1 Peptideins
    t1_auto = auto_df[auto_df["is_peptidein"] & auto_df["final_tier"].str.startswith("1")]
    ax.scatter(
        t1_auto["length_aa"],
        t1_auto["median_avi_phred"],
        color="#d90429",
        edgecolor="black",
        s=60,
        alpha=0.9,
        label=f"Autonomous Tier 1 Peptideins (n={len(t1_auto)})",
        zorder=6
    )

    # Highlight top unconfounded lncRNA-ORFs with Phred >= 20
    top_lnc = auto_df[(auto_df["orf_biotype"] == "lncRNA-ORF") & (auto_df["median_avi_phred"] >= 20.0)]
    ax.scatter(
        top_lnc["length_aa"],
        top_lnc["median_avi_phred"],
        facecolors="none",
        edgecolors="#7209b7",
        linewidths=1.8,
        s=120,
        label=f"Top Autonomous lncRNA-ORFs (Phred ≥ 20, n={len(top_lnc)})",
        zorder=7
    )

    # Annotate key autonomous loci
    standouts_to_annotate = [
        ("c3riboseqorf106", "ZBTB11-AS1 (Tier 1B Pep)"),
        ("c16norep119", "ZFHX3-AS1 (lncRNA)"),
        ("c2riboseqorf55", "WBP1 (Tier 1B uORF)"),
        ("c3riboseqorf183", "BCL6 (uORF)"),
        ("c2norep207", "TTN-AS1 (Tier 1B Pep)"),
        ("c12riboseqorf146", "PXN-AS1 (lncRNA)"),
        ("c7norep47", "HOXA-AS3 (lncRNA)")
    ]

    for orf_id, label in standouts_to_annotate:
        match = auto_df[auto_df["orf_id"] == orf_id]
        if not match.empty:
            r = match.iloc[0]
            ax.annotate(
                label,
                xy=(r["length_aa"], r["median_avi_phred"]),
                xytext=(r["length_aa"] + 4, r["median_avi_phred"] + 1.1),
                arrowprops=dict(facecolor="#2b2d42", arrowstyle="->", lw=1.0),
                fontsize=9,
                fontweight="bold",
                color="#1d3557",
                zorder=10
            )

    ax.axhline(20, color="#d90429", linestyle="--", linewidth=1.3, alpha=0.8, label="Phred ≥ 20 (Top 1%)")
    ax.axhline(15, color="#f77f00", linestyle=":", linewidth=1.3, alpha=0.8, label="Phred ≥ 15 (Top 3.2%)")

    ax.set_title("Unconfounded Autonomous smORFs (lncRNA-ORFs & 5' UTR uORFs) Under Purifying Selection", fontsize=13, fontweight="bold", pad=12)
    ax.set_xlabel("smORF Peptide Length (Amino Acids)", fontsize=11, fontweight="bold")
    ax.set_ylabel("Median AVI Phred Score (Background Constraint)", fontsize=11, fontweight="bold")
    ax.set_xlim(0, 180)
    ax.set_ylim(0, 30)
    ax.legend(loc="upper right", frameon=True, framealpha=0.9)

    plt.tight_layout()
    fig7_path = FIGURES_DIR / "fig7_top_autonomous_candidates_landscape.png"
    plt.savefig(fig7_path, dpi=300)
    plt.close()
    print(f"  [+] Saved Figure 7: {fig7_path}")

    # Copy to artifacts
    for f in [fig5_path, fig6_path, fig7_path]:
        shutil.copy(f, ARTIFACTS_FIG_DIR / f.name)
    print(f"  [+] Copied Figures 5, 6, 7 to artifacts dir: {ARTIFACTS_FIG_DIR}")


def build_deconfounded_leaderboards(df):
    print("\n[4/5] Building de-confounded standout leaderboards...")

    # LEADERBOARD 1: Top Autonomous lncRNA-ORFs
    # Pure non-coding transcript, no canonical CDS anywhere
    lnc_hits = df[
        (df["orf_biotype"] == "lncRNA-ORF") &
        (df["median_avi_phred"] >= 15.0) &
        (df["is_high_impact_hotspot"] == True)
    ].sort_values("max_avi_phred", ascending=False).reset_index(drop=True)
    lnc_hits["rank"] = range(1, len(lnc_hits) + 1)
    
    top20_lnc = lnc_hits.head(20)[[
        "rank", "orf_id", "gene", "length_aa", "final_tier", "is_peptidein",
        "median_avi_phred", "max_avi_phred", "peak_variant", "primary_mechanism",
        "convergent_modalities_count", "atlas_url"
    ]]
    top20_lnc.to_csv(REPORTS_DIR / "top20_autonomous_lncrna_orfs.tsv", sep="\t", index=False)
    print(f"  [+] Saved Top 20 Autonomous lncRNA-ORFs: {REPORTS_DIR / 'top20_autonomous_lncrna_orfs.tsv'}")

    # LEADERBOARD 2: Top Autonomous uORFs (Strictly 5' UTR, non-CDS overlapping)
    uorf_hits = df[
        (df["orf_biotype"] == "uORF") &
        (df["median_avi_phred"] >= 18.0) &
        (df["is_high_impact_hotspot"] == True)
    ].sort_values("max_avi_phred", ascending=False).reset_index(drop=True)
    uorf_hits["rank"] = range(1, len(uorf_hits) + 1)

    top20_uorf = uorf_hits.head(20)[[
        "rank", "orf_id", "gene", "length_aa", "final_tier", "is_peptidein",
        "median_avi_phred", "max_avi_phred", "peak_variant", "primary_mechanism",
        "convergent_modalities_count", "atlas_url"
    ]]
    top20_uorf.to_csv(REPORTS_DIR / "top20_autonomous_uorfs.tsv", sep="\t", index=False)
    print(f"  [+] Saved Top 20 Autonomous uORFs: {REPORTS_DIR / 'top20_autonomous_uorfs.tsv'}")

    # LEADERBOARD 3: Top Autonomous Mass-Spec Tier 1 Peptideins
    t1_auto_hits = df[
        (df["is_autonomous"] == True) &
        (df["is_peptidein"] == True) &
        (df["final_tier"].str.startswith("1", na=False))
    ].sort_values("median_avi_phred", ascending=False).reset_index(drop=True)
    t1_auto_hits["rank"] = range(1, len(t1_auto_hits) + 1)

    top_t1_auto = t1_auto_hits.head(20)[[
        "rank", "orf_id", "gene", "orf_biotype", "length_aa", "final_tier",
        "median_avi_phred", "max_avi_phred", "peak_variant", "primary_mechanism",
        "convergent_modalities_count", "atlas_url"
    ]]
    top_t1_auto.to_csv(REPORTS_DIR / "top_autonomous_tier1_peptideins.tsv", sep="\t", index=False)
    print(f"  [+] Saved Top Autonomous Tier 1 Peptideins: {REPORTS_DIR / 'top_autonomous_tier1_peptideins.tsv'}")

    return top20_lnc, top20_uorf, top_t1_auto


def main():
    print("=" * 80)
    print("DE-CONFOUNDED MULTI-OMICS ANALYSIS (AUTONOMOUS vs CDS-OVERLAPPING smORFs)")
    print("=" * 80)

    df = load_dataset()
    stats_results = compute_deconfounded_statistics(df)
    generate_deconfounded_figures(df, stats_results)
    top20_lnc, top20_uorf, top_t1_auto = build_deconfounded_leaderboards(df)

    with open(REPORTS_DIR / "deconfounded_statistical_summary.json", "w") as f:
        json.dump(stats_results, f, indent=2)

    print("\n" + "=" * 80)
    print("DE-CONFOUNDED PIPELINE COMPLETED SUCCESSFULLY!")
    print("=" * 80)


if __name__ == "__main__":
    main()
