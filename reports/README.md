# Microprotein Empirical Reports & Milestone Index

This directory serves as the centralized repository for all computed scientific reports, statistical validations, and empirical structural milestones generated across the project.

---

## Chronological & Methodological Pipeline Index

The reports follow a progressive experimental workflow that moves from raw non-canonical genomic coordinates to 3D atomic-resolution macromolecular complexes:

```mermaid
flowchart TD
    R1["Report 1: Multi-Omics Selection & Hitchhiking<br/>(2.9M SNVs, Proved 49.1x Hitchhiking Confounder)"] --> R2["Report 2: Curated Top 50 Portfolio<br/>(50 Candidates across 4 Balanced Cohorts)"]
    R2 --> R3["Report 3: Binding Partner Prioritization<br/>(500 Physiological Target Complexes)"]
    R3 --> R4["Report 4: Initial ESMFold2 Cofolding (Phase 1)<br/>(345 Folded / 155 Exceeded 768 aa Limit)"]
    R4 --> R5["Report 5: Master Structural Atlas (Authoritative)<br/>(805 Complexes, 100% Top 50 Coverage + 61 Tier 1 Interactomes)"]

    style R1 fill:#e1f5fe,stroke:#0288d1,stroke-width:1px
    style R2 fill:#e1f5fe,stroke:#0288d1,stroke-width:1px
    style R3 fill:#e1f5fe,stroke:#0288d1,stroke-width:1px
    style R4 fill:#fff3e0,stroke:#f57c00,stroke-width:1px
    style R5 fill:#e8f5e9,stroke:#388e3c,stroke-width:2px
```

---

## Detailed Report Summaries

| Report File | Title & Scope | Status | Key Deliverables & Artifacts |
| :--- | :--- | :--- | :--- |
| [**`microprotein_multiomics_analysis_report.md`**](microprotein_multiomics_analysis_report.md) | **Comprehensive Multi-Omics Characterization & Genomic Selection**<br>• Scope: 7,264 human smORFs<br>• 2,900,430 SNVs scored via AlphaGenome AVI | `COMPLETED` *(Exploratory)* | • Statistical proof of 49.1x CDS-overlapping hitchhiking confounder ($p < 10^{-300}$)<br>• Definition of the Autonomous smORF benchmark ($n = 5,387$)<br>• [`data/screening/microprotein_avi_scores.parquet`](../data/screening/microprotein_avi_scores.parquet) |
| [**`top50_curated_microprotein_portfolio.md`**](top50_curated_microprotein_portfolio.md) | **Curated Top 50 Human Microprotein Portfolio**<br>• Scope: 50 prioritized candidates across 4 balanced cohorts | `COMPLETED` *(Baseline)* | • Stratified portfolio (Tier 1 autonomous, lncRNA-ORFs, 5' UTR uORFs, dual-coding)<br>• Standardized CDS coordinates and biotype assignments<br>• [`top50_curated_candidates.parquet`](top50_curated_candidates.parquet) |
| [**`top50_microprotein_binding_partners_summary.md`**](top50_microprotein_binding_partners_summary.md) | **Top 50 Microprotein Binding Partners Summary**<br>• Scope: 500 candidate complexes (Top 50 $\times$ 10 partners) | `COMPLETED` *(Target Mapping)* | • 4-tier biological evidence curation (HIPPIE, BioGRID, IntAct, STRING)<br>• Complete sequence pairs prepared for cofolding<br>• [`top50_microprotein_binding_partners_top10.parquet`](top50_microprotein_binding_partners_top10.parquet) |
| [**`esmfold2_cofolding_report.md`**](esmfold2_cofolding_report.md) | **Initial Biohub ESMFold2 Cofolding (Phase 1)**<br>• Scope: 500 candidate complexes screened on Biohub API | `INTERMEDIATE` *(Superseded)* | • 345 complexes successfully folded ($\le 768\text{ aa}$)<br>• Identified 155 partner size bottlenecks exceeding interactive limit<br>• *Note: See Report 5 for the complete remediated atlas* |
| [**`expanded_ppi_atlas_report.md`**](expanded_ppi_atlas_report.md) | **Master Microprotein–Protein Structural Atlas (805 Complexes)**<br>• Scope: 805 full-atom cofolded complexes | `AUTHORITATIVE` *(Benchmark)* | • **100% Top 50 Portfolio coverage** (500/500 complexes)<br>• Full Tier 1 peptidein interactome (305 complexes)<br>• Discovery of nanomolar-like interfaces (`SNX13`–`SERPINE2` $\text{ipTM} = 0.8525$)<br>• [`master_microprotein_ppi_structural_atlas.parquet`](master_microprotein_ppi_structural_atlas.parquet)<br>• 805 PDBs in [`structures/pdbs/`](../structures/pdbs/) |

---

## Quick Reference: Authoritative Datasets & Outputs

* **Master Structural Table**: [`reports/master_microprotein_ppi_structural_atlas.parquet`](master_microprotein_ppi_structural_atlas.parquet) (805 rows, 35 metrics including ipTM, pTM, pLDDT, contacts, PAE).
* **Tab-Separated Table**: [`reports/master_microprotein_ppi_structural_atlas.tsv`](master_microprotein_ppi_structural_atlas.tsv).
* **3D Atomic Structures**: [`structures/pdbs/`](../structures/pdbs/) (805 full-atom `.pdb` files).
* **Inter-Chain PAE Matrices**: [`structures/pae/`](../structures/pae/) (805 `.json` matrices).
