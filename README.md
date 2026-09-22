# Syn2Nat: Binder-Guided Latent-Space Retrieval of Natural Microproteins for Protein–Protein Binding Prediction

> **Authoritative Master Reference**: See [**`PROJECT_MASTER_DOCUMENT.md`**](PROJECT_MASTER_DOCUMENT.md) for the single-file master reference detailing genomic selection datasets, hitchhiking proofs, the Top 50 portfolio, and the **Syn2Nat** discovery framework.
> **Methodology Specification**: See [**`docs/latent_space_binder_matching_methodology.md`**](docs/latent_space_binder_matching_methodology.md) for full architectural blueprints, complexity comparisons, and workflows.

---

## 🎯 The Syn2Nat Paradigm

> *"Rather than testing every candidate microprotein against every human target ($\mathcal{O}(N \times M)$ combinatorial explosion), Syn2Nat inverts the discovery direction: we generate unconstrained de novo binders for prioritized targets, project them into ESMC latent space, and mine natural microproteins whose representations align with these target-specific fingerprints."*

### Core Workflow:
1. **Target Prioritization**: Prioritize high-value human targets from the **AlphaGenome Top 1,000 list** based on dense purifying selection (AVI Phred scores), tissue-specific epigenetic regulation, and disease relevance.
2. **De Novo Binder Generation (Unconstrained Biologically Honest Scaffolds)**: Run the gradient-guided optimization protocol ([`notebooks/syn2nat_colab.ipynb`](notebooks/syn2nat_colab.ipynb)) directly on Google Colab or local GPU clusters. Binder geometries are **never** artificially forced to match microprotein lengths—their scaffolds emerge naturally from target epitope biophysics.
3. **ESMC Latent Space Embedding**: Extract representations using the same **ESMC-6B** / **ESMC-600M** pipeline utilized for the natural microprotein atlas.
4. **Multi-Modal Latent & Sequence Search ("Computationally Free" Space)**: Screen thousands of microproteins in seconds using vector similarity (cosine distance, $k$-NN) combined with primary sequence motif alignment and physicochemical profiling.
5. **Targeted Structural Validation**: Validate only the top-filtered microprotein candidates using high-resolution **ESMFold2** co-folding ($\text{ipTM} \ge 0.60$, inter-chain $\text{PAE} < 10\text{ \AA}$).

---

## 🧭 Repository Navigation

```
microproteinproject/
├── PROJECT_MASTER_DOCUMENT.md             # Authoritative scientific monolith (all phases, proofs & metrics)
├── README.md                              # This repository portal
│
├── reports/                               # Empirical results, validation reports & structural atlas
│   ├── README.md                          # Progressive report index (Reports 1 to 5)
│   ├── master_microprotein_ppi_structural_atlas.parquet # Master database of 805 cofolded complexes
│   └── ...                                # Individual markdown milestone reports
│
├── docs/                                  # Technical specifications, methods & API manuals
│   ├── README.md                          # Documentation index and domain taxonomy
│   ├── mappie.md                          # MAPPIE 128-d latent space & functional discovery guide
│   ├── representation_pipeline_roadmap.md # Large-scale representation engine roadmap (HMPA 617k)
│   ├── ppi_evaluation_hierarchy.md        # 4-pillar biophysical evaluation framework
│   ├── esmfold2_latent_interface_screening.md # ESMC 6B latent interface screening architecture
│   ├── hippie.md, biogrid.md, intact.md, string.md, pepbind.md # Canonical interactome guides
│   ├── biohub.md                          # Biohub ESM / ESMFold2 platform documentation
│   └── alphagenome/                       # DeepMind AlphaGenome multimodal DNA model guide
│
├── bindpred/                              # Active ESMC 6B representation extraction engine
│   ├── extract_representations.py         # Batch extraction of layer 79 sequence & residue embeddings
│   ├── esmc_client.py                     # Biohub ESMC API client
│   └── representations/                   # Checkpointed Parquet & NPZ representation archives
│
├── mappie/                                # Cloned MAPPIE core framework (Uni Mainz CBDM)
│   ├── core_algorithm/                    # Latent projection & hypergeometric enrichment scripts
│   └── data_processed/                    # Preprocessed runtime indices
│
├── data/                                  # Project datasets & reference interactome indices
│   ├── mappie/                            # 199,138 reference PPI latent index (latent_index.npz)
│   ├── curated/                           # Top 50 portfolio, binding partners, query intervals
│   ├── parsed/                            # Parsed smORF catalogs (orfs.parquet, orf_tiers.parquet)
│   ├── screening/                         # AlphaGenome AVI purifying selection scores
│   ├── targets/                           # Canonical human proteome (UP000005640)
│   └── pfam/                              # Pfam-A profile HMMs
│
├── structures/                            # 3D macromolecular structures & confidence matrices
│   ├── pdbs/                              # 805 full-atom cofolded PDB files
│   └── pae/                               # 805 inter-chain PAE JSON matrices
│
├── scripts/                               # Executable clients and pipeline scripts
├── plans/                                 # Multi-agent virtual lab protocols & Paperclip orchestration
└── literature/                            # Published studies, preprints & supplementary data
```

---

## 📊 Empirical Milestones & Reports Hub

Detailed findings are organized into five progressive reports in [**`reports/`**](reports/README.md):

| # | Milestone Report | Scientific Focus & Deliverable | Status |
| :-: | :--- | :--- | :--- |
| **1** | [**`microprotein_multiomics_analysis_report.md`**](reports/microprotein_multiomics_analysis_report.md) | **Multi-Omics Selection & Hitchhiking Resolution**<br>Screened 2.9M SNVs via AlphaGenome AVI. Proved that CDS-overlapping smORFs (`intORFs`) exhibit a 49.1x false constraint artifact from host gene hitchhiking ($p < 10^{-300}$). Established the autonomous benchmark ($n = 5,387$). | `COMPLETED` |
| **2** | [**`top50_curated_microprotein_portfolio.md`**](reports/top50_curated_microprotein_portfolio.md) | **Curated Top 50 Microprotein Portfolio**<br>Constructed a de-confounded portfolio across 4 balanced cohorts (Tier 1 Autonomous, lncRNA-ORFs, 5' UTR uORFs, MS-validated dual-coding). | `COMPLETED` |
| **3** | [**`top50_microprotein_binding_partners_summary.md`**](reports/top50_microprotein_binding_partners_summary.md) | **Prioritized Binding Partner Interactome**<br>Curated 500 candidate complexes (Top 50 $\times$ 10 partners) across a 4-tier biological evidence hierarchy (HIPPIE, BioGRID, IntAct, STRING). | `COMPLETED` |
| **4** | [**`esmfold2_cofolding_report.md`**](reports/esmfold2_cofolding_report.md) | **Initial ESMFold2 Batch Screening (Phase 1)**<br>Successfully folded 345 complexes $\le 768\text{ aa}$. Flagged 155 partner size bottlenecks exceeding interactive API limits. *(Superseded by Report 5)*. | `INTERMEDIATE` |
| **5** | [**`expanded_ppi_atlas_report.md`**](reports/expanded_ppi_atlas_report.md) | **Master Structural Atlas (805 Complexes)**<br>**100% complete coverage of Top 50** (500 complexes) + full Tier 1 peptidein interactome (305 complexes). Identified breakthrough complexes (`SNX13`–`SERPINE2` $\text{ipTM} = 0.8525$, `AP001372.2`–`PCNA` $\text{ipTM} = 0.8350$). | `AUTHORITATIVE` |

---

## 🛠️ Integrated Interactomes & Documentation

Complete technical manuals and Python API clients are indexed in [**`docs/`**](docs/README.md):

* **Latent Space & Deorphanization**:
  * [**`docs/mappie.md`**](docs/mappie.md): MAPPIE continuous latent space, autoencoder projection, and functional discovery.
  * [**`docs/representation_pipeline_roadmap.md`**](docs/representation_pipeline_roadmap.md): Strategic roadmap for HMPA 617k smORF representation banking.
  * [**`docs/ppi_evaluation_hierarchy.md`**](docs/ppi_evaluation_hierarchy.md): 4-pillar biophysical evaluation criteria.
  * [**`docs/esmfold2_latent_interface_screening.md`**](docs/esmfold2_latent_interface_screening.md): ESMC 6B trunk decoupling architecture.
* **Scored Interactomes & Python Clients**:
  * [**HIPPIE v2.3**](docs/hippie.md) ([`scripts/hippie_client.py`](scripts/hippie_client.py)): Scored physical interactome ($0.0 \le S \le 1.0$), GTEx tissue filtering, MAPPIE reference base.
  * [**BioGRID v4.4**](docs/biogrid.md) ([`scripts/biogrid_client.py`](scripts/biogrid_client.py)): Experimental evidence codes (AP-MS, Two-Hybrid).
  * [**IntAct / IMEx**](docs/intact.md) ([`scripts/intact_client.py`](scripts/intact_client.py)): Deep molecular curation, contact mutations, and stoichiometry.
  * [**STRING v12.0**](docs/string.md) ([`scripts/string_client.py`](scripts/string_client.py)): Multi-evidence functional network integration.
  * [**PepBind**](docs/pepbind.md) ([`scripts/pepbind_client.py`](scripts/pepbind_client.py)): Structurally resolved peptide–protein interfaces.
* **Platforms & Foundation Models**:
  * [**Biohub ESM Platform**](docs/biohub.md): `esmc-6b-2024-12`, `esmfold2-fast-2026-05`, and ESM Atlas.
  * [**DeepMind AlphaGenome**](docs/alphagenome/README.md): 14-part guide to the 1-Mb multimodal genomic sequence model.

---

## 📚 Foundational Literature

1. **Primary Study**: *"Expanding the human proteome with microproteins and peptideins"*, *Nature* 2026 (DOI: [10.1038/s41586-026-10459-x](https://doi.org/10.1038/s41586-026-10459-x)) &mdash; Stored under [`literature/s41586-026-10459-x.pdf`](literature/s41586-026-10459-x.pdf).
2. **MAPPIE Framework**: *"A map of human protein-protein interaction embeddings for functional discovery"*, *bioRxiv* 2026 (DOI: [10.64898/2026.08.07.743440](https://doi.org/10.64898/2026.08.07.743440)) &mdash; Stored under [`literature/2026.08.07.743440v1.full.pdf`](literature/2026.08.07.743440v1.full.pdf).
3. **HMPA Platform**: *"Comprehensive annotation and analysis of human microproteins by human microprotein atlas platform"*, *Communications Chemistry* 2026 (DOI: [10.1038/s42004-026-01588-y](https://doi.org/10.1038/s42004-026-01588-y)).
