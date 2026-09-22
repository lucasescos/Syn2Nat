#!/usr/bin/env python3
"""
Generate Prioritized List of 1,000 Human Canonical Proteins using AlphaGenome Data.

Integrates:
1. Reference Human Canonical Proteome (UniProt UP000005640, 20,652 entries).
2. DeepMind AlphaGenome Variant Impact (AVI) purifying selection scores across 2.9M dense SNVs.
3. DeepMind AlphaGenome 1-Mb Multimodal DNA Model track disruption predictions:
   - Splicing disruption (cryptic donors/acceptors, junction usage, exon skipping).
   - Gene expression (RNA-seq).
   - Chromatin accessibility (DNase I & ATAC-seq).
   - Transcription factor binding ablation (ChIP-seq).
   - Convergent modality count and high-impact regulatory hotspots.

Outputs:
- reports/canonical_proteome_alphagenome_top1000.tsv
- reports/canonical_proteome_alphagenome_top1000.parquet
- data/targets/canonical_proteome_alphagenome_top1000.tsv
"""

import sys
from pathlib import Path
import numpy as np
import pandas as pd

BASE_DIR = Path(__file__).resolve().parent.parent
DATA_DIR = BASE_DIR / "data"
REPORTS_DIR = BASE_DIR / "reports"
TARGETS_DIR = DATA_DIR / "targets"
SCREENING_DIR = DATA_DIR / "screening"

def main():
    print("=" * 80)
    print("Prioritizing 1,000 Canonical Proteins Using DeepMind AlphaGenome Data")
    print("=" * 80)

    # 1. Load Canonical Proteome
    proteome_path = TARGETS_DIR / "human_canonical_proteome.parquet"
    print(f"Loading canonical proteome from {proteome_path}...")
    df_canon = pd.read_parquet(proteome_path)
    print(f"Loaded {len(df_canon):,} canonical proteins.")

    # 2. Load AlphaGenome Screening Data
    avi_path = SCREENING_DIR / "microprotein_avi_scores.parquet"
    dis_path = SCREENING_DIR / "microprotein_single_variant_disruptions.parquet"
    print(f"Loading AlphaGenome AVI scores from {avi_path}...")
    df_avi = pd.read_parquet(avi_path)
    print(f"Loading AlphaGenome multimodal disruptions from {dis_path}...")
    df_dis = pd.read_parquet(dis_path)

    # Merge AVI and Disruption features
    dis_cols = [
        'orf_id', 'convergent_modalities_count', 'is_high_impact_hotspot', 
        'primary_mechanism', 'splicing_max_abs_raw', 'rna_max_abs_raw', 
        'chromatin_max_abs_raw', 'chip_tf_max_abs_raw'
    ]
    df_ag = df_avi.merge(
        df_dis[[c for c in dis_cols if c in df_dis.columns]],
        on='orf_id',
        how='left'
    )
    print(f"Merged AlphaGenome screening features: {len(df_ag):,} loci.")

    # 3. Merge with Canonical Proteome by Gene Symbol
    merged = df_canon.merge(df_ag, left_on='gene_symbol', right_on='gene', how='inner')
    print(f"Mapped {len(merged):,} AlphaGenome loci across {merged['uniprot_id'].nunique():,} unique canonical proteins.")

    # 4. Aggregate AlphaGenome Metrics to Canonical Protein Level
    records = []
    for (uniprot_id, entry_name, gene_symbol, protein_name, length, seq), group in merged.groupby(
        ['uniprot_id', 'entry_name', 'gene_symbol', 'protein_name', 'length', 'sequence']
    ):
        med_phred = float(group['median_avi_phred'].max())
        mean_med_phred = float(group['median_avi_phred'].mean())
        max_phred = float(group['max_avi_phred'].max())
        p90_phred = float(group['p90_avi_phred'].max())
        pct20 = float(group['pct_phred_ge_20'].max())
        pct15 = float(group['pct_phred_ge_15'].max())
        total_variants = int(group['n_variants_scored'].sum())

        conv_max = group['convergent_modalities_count'].max()
        conv_val = int(conv_max) if pd.notna(conv_max) else 0
        high_impact = bool(group['is_high_impact_hotspot'].fillna(False).any())

        splice_max = float(group['splicing_max_abs_raw'].max()) if pd.notna(group['splicing_max_abs_raw'].max()) else 0.0
        rna_max = float(group['rna_max_abs_raw'].max()) if pd.notna(group['rna_max_abs_raw'].max()) else 0.0
        chrom_max = float(group['chromatin_max_abs_raw'].max()) if pd.notna(group['chromatin_max_abs_raw'].max()) else 0.0
        tf_max = float(group['chip_tf_max_abs_raw'].max()) if pd.notna(group['chip_tf_max_abs_raw'].max()) else 0.0

        # Best row for variant & atlas link
        best_row = group.loc[group['median_avi_phred'].idxmax()]
        biotypes = set(group['orf_biotype'].dropna())
        has_cds_overlap = any(b in ['intORF', 'uoORF', 'doORF', 'mixed'] for b in biotypes)
        has_utr = any(b in ['uORF', 'dORF'] for b in biotypes)

        records.append({
            'uniprot_id': uniprot_id,
            'entry_name': entry_name,
            'gene_symbol': gene_symbol,
            'protein_name': protein_name,
            'length': length,
            'sequence': seq,
            'smorf_loci_count': len(group),
            'biotypes_present': '; '.join(sorted(biotypes)),
            'has_cds_overlap': has_cds_overlap,
            'has_utr_smorf': has_utr,
            'max_median_avi_phred': round(med_phred, 3),
            'mean_median_avi_phred': round(mean_med_phred, 3),
            'peak_hotspot_phred': round(max_phred, 3),
            'p90_avi_phred': round(p90_phred, 3),
            'pct_variants_phred_ge_20': round(pct20, 2),
            'pct_variants_phred_ge_15': round(pct15, 2),
            'total_variants_scored': total_variants,
            'convergent_modalities_count': conv_val,
            'is_high_impact_hotspot': high_impact,
            'splicing_disruption_raw': round(splice_max, 4),
            'rna_disruption_raw': round(rna_max, 4),
            'chromatin_disruption_raw': round(chrom_max, 4),
            'chip_tf_disruption_raw': round(tf_max, 4),
            'primary_mechanism': str(best_row.get('primary_mechanism', 'Unknown')),
            'peak_variant': str(best_row.get('peak_variant', '')),
            'atlas_url': str(best_row.get('atlas_url', ''))
        })

    agg_df = pd.DataFrame(records)
    print(f"Consolidated into {len(agg_df):,} unique canonical proteins with AlphaGenome annotations.")

    # 5. Calculate AlphaGenome Priority Score (AGPS)
    # A. Purifying Selection Pillar: median Phred + 0.3 * clipped peak Phred
    agg_df['avi_selection_score'] = agg_df['max_median_avi_phred'] + 0.3 * np.clip(agg_df['peak_hotspot_phred'], 0, 50)

    # B. Multimodal Disruption Pillar: convergent modalities (x3.0) + high impact hotspot flag (5.0) + splicing impact (x2.5)
    agg_df['multimodal_score'] = (
        agg_df['convergent_modalities_count'] * 3.0 + 
        agg_df['is_high_impact_hotspot'].astype(float) * 5.0 +
        np.clip(agg_df['splicing_disruption_raw'], 0, 2.0) * 2.5
    )

    # C. Locus Architecture Pillar: 1.25x weighting for verified coding sequence constraint
    agg_df['locus_weight'] = np.where(agg_df['has_cds_overlap'], 1.25, 1.0)

    # D. Composite AlphaGenome Priority Score
    agg_df['alphagenome_priority_score'] = np.round(
        (agg_df['avi_selection_score'] * 1.5 + agg_df['multimodal_score']) * agg_df['locus_weight'],
        3
    )

    # 6. Rank and Filter to Top 1,000
    agg_df = agg_df.sort_values('alphagenome_priority_score', ascending=False).reset_index(drop=True)
    top1000 = agg_df.head(1000).copy()
    top1000['priority_rank'] = range(1, 1001)

    # Rearrange columns logically
    ordered_cols = [
        'priority_rank',
        'uniprot_id',
        'entry_name',
        'gene_symbol',
        'protein_name',
        'length',
        'alphagenome_priority_score',
        'max_median_avi_phred',
        'mean_median_avi_phred',
        'peak_hotspot_phred',
        'p90_avi_phred',
        'pct_variants_phred_ge_20',
        'pct_variants_phred_ge_15',
        'total_variants_scored',
        'convergent_modalities_count',
        'is_high_impact_hotspot',
        'primary_mechanism',
        'splicing_disruption_raw',
        'rna_disruption_raw',
        'chromatin_disruption_raw',
        'chip_tf_disruption_raw',
        'has_cds_overlap',
        'has_utr_smorf',
        'smorf_loci_count',
        'biotypes_present',
        'peak_variant',
        'atlas_url',
        'sequence'
    ]
    top1000 = top1000[ordered_cols]

    # 7. Export Datasets
    tsv_report_path = REPORTS_DIR / "canonical_proteome_alphagenome_top1000.tsv"
    parquet_report_path = REPORTS_DIR / "canonical_proteome_alphagenome_top1000.parquet"
    tsv_targets_path = TARGETS_DIR / "canonical_proteome_alphagenome_top1000.tsv"

    top1000.to_csv(tsv_report_path, sep="\t", index=False)
    top1000.to_parquet(parquet_report_path, index=False)
    top1000.to_csv(tsv_targets_path, sep="\t", index=False)

    print("\n" + "=" * 80)
    print("SUCCESSFULLY GENERATED TOP 1,000 CANONICAL PROTEOME DATASET")
    print("=" * 80)
    print(f"TSV Export (Reports):  {tsv_report_path} ({tsv_report_path.stat().st_size:,} bytes)")
    print(f"Parquet Export:        {parquet_report_path} ({parquet_report_path.stat().st_size:,} bytes)")
    print(f"TSV Export (Targets):  {tsv_targets_path} ({tsv_targets_path.stat().st_size:,} bytes)")
    
    # Statistical Summary
    print("\n--- Summary Statistics for Top 1,000 Canonical Proteins ---")
    print(f"Mean AlphaGenome Priority Score: {top1000['alphagenome_priority_score'].mean():.2f} (Range: {top1000['alphagenome_priority_score'].min():.2f} - {top1000['alphagenome_priority_score'].max():.2f})")
    print(f"Mean Median AVI Phred:           {top1000['max_median_avi_phred'].mean():.2f} (Max: {top1000['max_median_avi_phred'].max():.2f})")
    print(f"Mean Peak Hotspot Phred:         {top1000['peak_hotspot_phred'].mean():.2f} (Max: {top1000['peak_hotspot_phred'].max():.2f})")
    print(f"Proteins with Phred >= 20 (Top 1% constraint): {(top1000['max_median_avi_phred'] >= 20).sum()} ({((top1000['max_median_avi_phred'] >= 20).sum() / 1000)*100:.1f}%)")
    print(f"Proteins with High-Impact Hotspot:             {top1000['is_high_impact_hotspot'].sum()} ({(top1000['is_high_impact_hotspot'].sum() / 1000)*100:.1f}%)")
    print(f"Mean Convergent Modalities Count:              {top1000['convergent_modalities_count'].mean():.2f}")
    print(f"Mean Protein Length:                           {top1000['length'].mean():.1f} aa (Median: {top1000['length'].median():.0f} aa, Range: {top1000['length'].min()} - {top1000['length'].max()} aa)")
    print("\nTop Primary Disruption Mechanisms:")
    print(top1000['primary_mechanism'].value_counts())

if __name__ == "__main__":
    main()
