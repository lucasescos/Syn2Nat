# /// script
# dependencies = [
#   "alphagenome",
#   "python-dotenv",
#   "pandas",
#   "pyarrow",
#   "anndata",
#   "numpy",
#   "tqdm",
# ]
# ///
#!/usr/bin/env python3
"""
Step 4: AlphaGenome Single Variant Functional Disruption Analysis.

Executes DeepMind's AlphaGenome DNA model (`alphagenome.models.dna_client`) across
the non-canonical microprotein catalog. For each microprotein, analyzes its peak
deleterious hotspot mutation (identified from AVI screening) to assess functional
disruption on gene expression, splicing, chromatin accessibility, and TF binding.

Features:
- Multi-threaded batching with ThreadPoolExecutor
- Automatic exponential backoff on gRPC quota / rate limits
- Two-tier output: Master Disruption Summary + Filtered Significant Tracks
- Resumable: scans completed batches and only queries uncompleted ORFs
- Automatic master consolidation to Parquet & TSV
"""

import argparse
import os
import sys
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import dotenv
import numpy as np
import pandas as pd
from tqdm import tqdm

from alphagenome.data import genome
from alphagenome.models import dna_client, variant_scorers

BASE_DIR = Path(__file__).resolve().parent.parent
DATA_DIR = BASE_DIR / "data"
SCREENING_DIR = DATA_DIR / "screening"
BATCHES_DIR = SCREENING_DIR / "single_variant_batches"
CONSOLIDATED_SUMMARY_PARQUET = SCREENING_DIR / "microprotein_single_variant_disruptions.parquet"
CONSOLIDATED_SUMMARY_TSV = SCREENING_DIR / "microprotein_single_variant_disruptions.tsv"
CONSOLIDATED_TRACKS_PARQUET = SCREENING_DIR / "microprotein_single_variant_significant_tracks.parquet"

SEQ_LENGTH = 2**20  # 1,048,576 bp (recommended length)


def get_recommended_scorers():
    """Build the 10 recommended differential scorers."""
    return [
        variant_scorers.RECOMMENDED_VARIANT_SCORERS[m]
        for m in variant_scorers.RECOMMENDED_VARIANT_SCORERS
        if "ACTIVE" not in m and "CAGE" not in m and "PROCAP" not in m
    ]


def score_single_variant(client, scorers, row):
    """
    Score a single microprotein hotspot variant across all 10 differential scorers.
    Returns (summary_dict, significant_tracks_df).
    """
    orf_id = row["orf_id"]
    variant_str = row["peak_variant"]
    host_gene = str(row.get("gene", ""))

    if not variant_str or pd.isna(variant_str):
        return None, None

    # Parse chr:pos:ref>alt (1-based)
    chrom, pos_str, ref_alt = variant_str.split(":")
    ref, alt = ref_alt.split(">")
    pos = int(pos_str)

    variant = genome.Variant(chrom, pos, ref, alt)
    interval = variant.reference_interval.resize(SEQ_LENGTH)

    # gRPC query with retry on quota / rate exhaustion
    max_retries = 5
    scores_list = None
    for attempt in range(max_retries):
        try:
            scores_list = client.score_variant(
                interval=interval,
                variant=variant,
                variant_scorers=scorers,
            )
            break
        except Exception as e:
            err_str = str(e)
            if "RESOURCE_EXHAUSTED" in err_str or "Quota exceeded" in err_str:
                wait_time = 30 * (attempt + 1)
                print(f"\n[Quota limit] Pausing {wait_time}s before retry for {orf_id} ({attempt+1}/{max_retries})...")
                time.sleep(wait_time)
            else:
                raise e

    if scores_list is None:
        raise RuntimeError(f"Failed to score {orf_id} ({variant_str}) after {max_retries} retries.")

    # Tidy into tabular DataFrame
    dfs = []
    for score_adata in scores_list:
        try:
            df = variant_scorers.tidy_scores([score_adata], match_gene_strand=True)
            if df is not None and not df.empty:
                dfs.append(df)
        except Exception as e:
            pass

    if not dfs:
        return None, None

    tidy_df = pd.concat(dfs, ignore_index=True)

    # -------------------------------------------------------------
    # 1. Modality Extraction & Summary Metrics
    # -------------------------------------------------------------
    summary = {
        "orf_id": orf_id,
        "gene": host_gene,
        "transcript": row.get("transcript", ""),
        "orf_biotype": row.get("orf_biotype", ""),
        "final_tier": row.get("final_tier", ""),
        "is_peptidein": row.get("is_peptidein", False),
        "strand": row.get("strand", ""),
        "chrom": chrom,
        "length_aa": row.get("length_aa", 0),
        "peak_variant": variant_str,
        "max_avi_phred": row.get("max_avi_phred", np.nan),
        "median_avi_phred": row.get("median_avi_phred", np.nan),
    }

    # Helper to extract max stats for a subset
    def extract_stats(sub_df, prefix):
        if sub_df.empty:
            return {
                f"{prefix}_max_abs_raw": 0.0,
                f"{prefix}_raw_at_max": 0.0,
                f"{prefix}_max_abs_quantile": 0.0,
                f"{prefix}_top_feature": "",
                f"{prefix}_top_biosample": "",
            }
        # Rank by absolute quantile first, then absolute raw
        abs_raw = sub_df["raw_score"].abs()
        abs_q = sub_df["quantile_score"].abs()
        idx_max = abs_q.idxmax() if abs_q.max() > 0.9 else abs_raw.idxmax()
        top_row = sub_df.loc[idx_max]

        feature = str(top_row.get("gene_name", top_row.get("transcription_factor", top_row.get("histone_mark", ""))))
        biosample = str(top_row.get("biosample_name", top_row.get("gtex_tissue", "")))

        return {
            f"{prefix}_max_abs_raw": round(float(abs_raw.max()), 4),
            f"{prefix}_raw_at_max": round(float(top_row["raw_score"]), 4),
            f"{prefix}_max_abs_quantile": round(float(abs_q.max()), 5),
            f"{prefix}_top_feature": feature if feature != "nan" else "",
            f"{prefix}_top_biosample": biosample if biosample != "nan" else "",
        }

    # RNA_SEQ (Gene expression)
    rna_df = tidy_df[tidy_df["output_type"] == "RNA_SEQ"]
    summary.update(extract_stats(rna_df, "rna"))
    # Host-gene specific RNA effect
    if host_gene and not rna_df.empty:
        host_rna = rna_df[rna_df["gene_name"] == host_gene]
        if not host_rna.empty:
            h_idx = host_rna["quantile_score"].abs().idxmax()
            summary["rna_host_raw"] = round(float(host_rna.loc[h_idx, "raw_score"]), 4)
            summary["rna_host_quantile"] = round(float(host_rna.loc[h_idx, "quantile_score"]), 5)
        else:
            summary["rna_host_raw"] = np.nan
            summary["rna_host_quantile"] = np.nan
    else:
        summary["rna_host_raw"] = np.nan
        summary["rna_host_quantile"] = np.nan

    # SPLICING (Splice Junctions, Splice Sites, Splice Site Usage)
    splice_types = ["SPLICE_JUNCTIONS", "SPLICE_SITES", "SPLICE_SITE_USAGE"]
    splicing_df = tidy_df[tidy_df["output_type"].isin(splice_types)]
    summary.update(extract_stats(splicing_df, "splicing"))

    # Splicing breakdowns
    junc_df = tidy_df[tidy_df["output_type"] == "SPLICE_JUNCTIONS"]
    summary["splice_junc_max_raw"] = round(float(junc_df["raw_score"].abs().max()), 4) if not junc_df.empty else 0.0
    summary["splice_junc_max_q"] = round(float(junc_df["quantile_score"].abs().max()), 5) if not junc_df.empty else 0.0

    site_df = tidy_df[tidy_df["output_type"] == "SPLICE_SITES"]
    summary["splice_site_max_raw"] = round(float(site_df["raw_score"].abs().max()), 4) if not site_df.empty else 0.0

    usage_df = tidy_df[tidy_df["output_type"] == "SPLICE_SITE_USAGE"]
    summary["splice_usage_max_raw"] = round(float(usage_df["raw_score"].abs().max()), 4) if not usage_df.empty else 0.0

    # Host-gene specific splicing effect
    if host_gene and not splicing_df.empty:
        host_splice = splicing_df[splicing_df["gene_name"] == host_gene]
        if not host_splice.empty:
            hs_idx = host_splice["quantile_score"].abs().idxmax()
            summary["splicing_host_raw"] = round(float(host_splice.loc[hs_idx, "raw_score"]), 4)
            summary["splicing_host_quantile"] = round(float(host_splice.loc[hs_idx, "quantile_score"]), 5)
        else:
            summary["splicing_host_raw"] = np.nan
            summary["splicing_host_quantile"] = np.nan
    else:
        summary["splicing_host_raw"] = np.nan
        summary["splicing_host_quantile"] = np.nan

    # CHROMATIN ACCESSIBILITY (DNASE & ATAC)
    dnase_df = tidy_df[tidy_df["output_type"] == "DNASE"]
    summary.update(extract_stats(dnase_df, "dnase"))

    atac_df = tidy_df[tidy_df["output_type"] == "ATAC"]
    summary.update(extract_stats(atac_df, "atac"))

    chromatin_raw = max(summary["dnase_max_abs_raw"], summary["atac_max_abs_raw"])
    chromatin_q = max(summary["dnase_max_abs_quantile"], summary["atac_max_abs_quantile"])
    summary["chromatin_max_abs_raw"] = chromatin_raw
    summary["chromatin_max_abs_quantile"] = chromatin_q

    # TRANSCRIPTION FACTORS (CHIP_TF)
    tf_df = tidy_df[tidy_df["output_type"] == "CHIP_TF"]
    summary.update(extract_stats(tf_df, "chip_tf"))

    # HISTONE MARKS (CHIP_HISTONE)
    histone_df = tidy_df[tidy_df["output_type"] == "CHIP_HISTONE"]
    summary.update(extract_stats(histone_df, "chip_histone"))

    # CONTACT MAPS
    contact_df = tidy_df[tidy_df["output_type"] == "CONTACT_MAPS"]
    summary["contact_map_max_abs_raw"] = round(float(contact_df["raw_score"].abs().max()), 4) if not contact_df.empty else 0.0

    # -------------------------------------------------------------
    # 2. Primary Mechanism Classification (per interpretation guide)
    # -------------------------------------------------------------
    # Rule: Raw score < 0.1 is noise/benign even if quantile is high
    # Strong: raw >= 0.5 & quantile >= 0.99
    # Splicing priority: Exon skipping / cryptic splice disruption
    mech = "Benign / Low Molecular Effect"
    convergences = 0

    splicing_sig = summary["splicing_max_abs_raw"] >= 0.5 and summary["splicing_max_abs_quantile"] >= 0.99
    rna_sig = summary["rna_max_abs_raw"] >= 0.5 and summary["rna_max_abs_quantile"] >= 0.99
    chromatin_sig = summary["chromatin_max_abs_raw"] >= 0.5 and summary["chromatin_max_abs_quantile"] >= 0.99
    tf_sig = summary["chip_tf_max_abs_raw"] >= 0.5 and summary["chip_tf_max_abs_quantile"] >= 0.99
    histone_sig = summary["chip_histone_max_abs_raw"] >= 0.5 and summary["chip_histone_max_abs_quantile"] >= 0.99

    convergences = sum([splicing_sig, rna_sig, chromatin_sig, tf_sig, histone_sig])

    if splicing_sig:
        mech = "Severe Splicing Disruption" if summary["splicing_max_abs_raw"] >= 1.0 else "Splicing Disruption"
    elif rna_sig:
        mech = "Transcriptional Silencing" if summary["rna_raw_at_max"] < 0 else "Transcriptional Activation"
    elif chromatin_sig:
        mech = "Chromatin Closing / Repression" if summary["dnase_raw_at_max"] < 0 else "Chromatin Opening"
    elif tf_sig:
        mech = "TF Binding Motif Disruption"
    elif histone_sig:
        mech = "Epigenetic Mark Alteration"
    elif any(summary[f"{p}_max_abs_raw"] >= 0.1 and summary[f"{p}_max_abs_quantile"] >= 0.99 for p in ["splicing", "rna", "dnase", "atac", "chip_tf"]):
        mech = "Subtle / Regulatory Shift"
    else:
        mech = "Benign / Tolerated"

    summary["primary_mechanism"] = mech
    summary["convergent_modalities_count"] = convergences
    summary["is_high_impact_hotspot"] = (convergences >= 1 and mech != "Subtle / Regulatory Shift" and mech != "Benign / Tolerated")

    # -------------------------------------------------------------
    # 3. Filter Significant Tracks (|q| >= 0.99 & |raw| >= 0.1)
    # -------------------------------------------------------------
    sig_mask = (tidy_df["quantile_score"].abs() >= 0.99) & (tidy_df["raw_score"].abs() >= 0.10)
    sig_tracks = tidy_df[sig_mask].copy()

    if not sig_tracks.empty:
        sig_tracks["orf_id"] = orf_id
        if "variant_id" in sig_tracks.columns:
            sig_tracks["variant_id"] = sig_tracks["variant_id"].apply(
                lambda v: f"{v.chromosome}:{v.position}:{v.reference_bases}>{v.alternate_bases}"
                if hasattr(v, "chromosome")
                else str(v)
            )
        keep_cols = [
            "orf_id",
            "variant_id",
            "output_type",
            "gene_name",
            "biosample_name",
            "raw_score",
            "quantile_score",
            "transcription_factor",
            "histone_mark",
            "gtex_tissue",
        ]
        available_cols = [c for c in keep_cols if c in sig_tracks.columns]
        sig_tracks = sig_tracks[available_cols].copy()
        for col in available_cols:
            if col not in ["raw_score", "quantile_score"]:
                sig_tracks[col] = sig_tracks[col].fillna("").astype(str)

    return summary, sig_tracks


def consolidate_results():
    """Consolidate all single variant batch files into master outputs."""
    summary_files = sorted(BATCHES_DIR.glob("batch_*_summary.parquet"))
    tracks_files = sorted(BATCHES_DIR.glob("batch_*_tracks.parquet"))

    if not summary_files:
        print("No summary batch files found to consolidate.")
        return

    print(f"\nConsolidating {len(summary_files)} summary batch files...")
    dfs_summary = [pd.read_parquet(f) for f in summary_files]
    df_master = pd.concat(dfs_summary, ignore_index=True).drop_duplicates(subset=["orf_id"])
    df_master = df_master.sort_values("max_avi_phred", ascending=False)

    df_master.to_parquet(CONSOLIDATED_SUMMARY_PARQUET, index=False)
    df_master.to_csv(CONSOLIDATED_SUMMARY_TSV, sep="\t", index=False)
    print(f" - Master Summary Parquet: {CONSOLIDATED_SUMMARY_PARQUET} ({len(df_master):,} ORFs)")
    print(f" - Master Summary TSV:     {CONSOLIDATED_SUMMARY_TSV}")

    if tracks_files:
        print(f"Consolidating {len(tracks_files)} significant tracks batch files...")
        dfs_tracks = [pd.read_parquet(f) for f in tracks_files]
        df_tracks_master = pd.concat(dfs_tracks, ignore_index=True)
        df_tracks_master.to_parquet(CONSOLIDATED_TRACKS_PARQUET, index=False)
        print(f" - Significant Tracks Parquet: {CONSOLIDATED_TRACKS_PARQUET} ({len(df_tracks_master):,} tracks)")

    # Print summary insights
    print("\n" + "=" * 70)
    print("SINGLE VARIANT FUNCTIONAL DISRUPTION SUMMARY")
    print("=" * 70)
    print(f"Total Microproteins Analyzed: {len(df_master):,}")
    print(f"High-Impact Hotspots (|raw| >= 0.5, |q| >= 0.99): {df_master['is_high_impact_hotspot'].sum():,} ({df_master['is_high_impact_hotspot'].mean()*100:.1f}%)")

    print("\nPrimary Molecular Mechanisms:")
    print(df_master["primary_mechanism"].value_counts().to_string())

    print("\nTop 5 Disruptive Splicing Hotspots:")
    splicing_top = df_master.sort_values("splicing_max_abs_raw", ascending=False).head(5)
    print(splicing_top[["orf_id", "gene", "peak_variant", "splicing_max_abs_raw", "splicing_top_feature", "primary_mechanism"]].to_string(index=False))

    print("\nTop 5 Disruptive Transcriptional Silencing Hotspots:")
    rna_top = df_master[df_master["rna_raw_at_max"] < 0].sort_values("rna_max_abs_raw", ascending=False).head(5)
    if not rna_top.empty:
        print(rna_top[["orf_id", "gene", "peak_variant", "rna_raw_at_max", "rna_top_feature", "rna_top_biosample"]].to_string(index=False))


def main():
    parser = argparse.ArgumentParser(description="AlphaGenome Single Variant Functional Disruption Analysis.")
    parser.add_argument("--batch-size", type=int, default=50, help="Number of ORFs per checkpoint batch.")
    parser.add_argument("--limit", type=int, default=None, help="Optional max ORFs to query (for quick runs).")
    parser.add_argument("--workers", type=int, default=6, help="Parallel worker threads for gRPC queries.")
    parser.add_argument("--start-idx", type=int, default=0, help="Starting index in the catalog.")
    parser.add_argument("--input-file", type=str, default=str(SCREENING_DIR / "microprotein_avi_scores.parquet"), help="Path to input parquet containing peak_variant.")
    args = parser.parse_args()

    print("=" * 70)
    print("AlphaGenome Single Variant Functional Disruption Pipeline")
    print("=" * 70)

    # 1. Load credentials
    dotenv.load_dotenv(BASE_DIR / ".env")
    dotenv.load_dotenv(Path.home() / ".env")
    api_key = os.environ.get("ALPHAGENOME_API_KEY")
    if not api_key:
        print("ERROR: ALPHAGENOME_API_KEY is missing from .env!")
        sys.exit(1)

    print("Connecting to AlphaGenome DNA Model API...")
    client = dna_client.create(
        api_key=api_key,
        address="dns:///gdmscience.googleapis.com:443",
    )
    scorers = get_recommended_scorers()
    print(f"Client authenticated. Using {len(scorers)} differential scorers:")
    print(" - " + "\n - ".join([s.name for s in scorers]))

    BATCHES_DIR.mkdir(parents=True, exist_ok=True)

    input_path = Path(args.input_file)
    print(f"\nLoading input dataset: {input_path}...")
    df_all = pd.read_parquet(input_path)
    print(f"Total targets in input file: {len(df_all):,}")

    if "peak_variant" not in df_all.columns:
        print("ERROR: 'peak_variant' column is missing from input file!")
        sys.exit(1)

    # Filter out rows with null/empty peak_variant
    df_all = df_all[df_all["peak_variant"].notna() & (df_all["peak_variant"] != "")].copy()
    print(f"Targets with valid peak_variant: {len(df_all):,}")

    # Check already scored ORFs across existing batch checkpoints
    scored_ids = set()
    existing_batches = sorted(BATCHES_DIR.glob("batch_*_summary.parquet"))
    if existing_batches:
        print(f"Found {len(existing_batches)} existing batch files in {BATCHES_DIR}. Checking completed ORFs...")
        for b in existing_batches:
            try:
                b_df = pd.read_parquet(b)
                scored_ids.update(b_df["orf_id"].tolist())
            except Exception as e:
                print(f"Warning reading {b.name}: {e}")
        print(f"Already scored: {len(scored_ids):,} ORFs.")

    # Filter to remaining uncompleted ORFs
    df_todo = df_all[~df_all["orf_id"].isin(scored_ids)].copy()
    if args.start_idx > 0:
        df_todo = df_todo.iloc[args.start_idx:]
    if args.limit:
        df_todo = df_todo.head(args.limit)

    print(f"ORFs to query in this execution run: {len(df_todo):,}")
    if len(df_todo) == 0:
        print("All requested ORFs are already scored! Proceeding to consolidation...")
        consolidate_results()
        return

    # Process in batches
    batch_num = len(existing_batches)
    rows = [row for _, row in df_todo.iterrows()]
    total_processed = 0
    t0 = time.time()

    def process_row(r):
        try:
            return score_single_variant(client, scorers, r)
        except Exception as err:
            return {"orf_id": r["orf_id"], "error": str(err)}, None

    for i in range(0, len(rows), args.batch_size):
        chunk = rows[i : i + args.batch_size]
        batch_num += 1
        batch_t0 = time.time()

        print(f"\n--- Querying Batch {batch_num} ({len(chunk)} ORFs, workers={args.workers}) ---")
        chunk_summaries = []
        chunk_tracks = []

        with ThreadPoolExecutor(max_workers=args.workers) as executor:
            future_to_orf = {executor.submit(process_row, r): r["orf_id"] for r in chunk}
            for fut in as_completed(future_to_orf):
                summary, tracks = fut.result()
                if summary and "error" not in summary:
                    chunk_summaries.append(summary)
                    if tracks is not None and not tracks.empty:
                        chunk_tracks.append(tracks)
                elif summary and "error" in summary:
                    print(f"  [!] Error on {summary['orf_id']}: {summary['error']}")

        batch_elapsed = time.time() - batch_t0
        total_processed += len(chunk_summaries)

        # Save batch checkpoints
        if chunk_summaries:
            batch_df = pd.DataFrame(chunk_summaries)
            summary_path = BATCHES_DIR / f"batch_{batch_num:05d}_summary.parquet"
            batch_df.to_parquet(summary_path, index=False)

            if chunk_tracks:
                tracks_df = pd.concat(chunk_tracks, ignore_index=True)
                tracks_path = BATCHES_DIR / f"batch_{batch_num:05d}_tracks.parquet"
                tracks_df.to_parquet(tracks_path, index=False)
            else:
                tracks_df = pd.DataFrame()

            print(f"Saved checkpoint {summary_path.name}: {len(batch_df)} ORFs, {len(tracks_df)} significant tracks in {batch_elapsed:.1f}s ({batch_elapsed/max(1, len(chunk)):.2f}s/ORF)")

            # Preview top mechanism in batch
            high_impact = batch_df[batch_df["is_high_impact_hotspot"]]
            print(f"  Batch High-Impact Hotspots: {len(high_impact)}/{len(batch_df)}")
            if not high_impact.empty:
                top_hit = high_impact.sort_values("splicing_max_abs_raw", ascending=False).iloc[0]
                print(f"  Top Disruption: {top_hit['orf_id']} ({top_hit['gene']}) -> {top_hit['primary_mechanism']} (Splicing Raw: {top_hit['splicing_max_abs_raw']})")

    total_elapsed = time.time() - t0
    print("\n" + "=" * 70)
    print(f"Finished querying run: {total_processed} ORFs in {total_elapsed:.1f} s ({total_elapsed/max(1, total_processed):.2f} s/ORF)")
    print("=" * 70)

    consolidate_results()


if __name__ == "__main__":
    main()
