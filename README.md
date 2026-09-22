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
Syn2Nat/
├── README.md                              # Repository portal & quickstart
├── PROJECT_MASTER_DOCUMENT.md             # Authoritative scientific & architectural monograph
│
├── notebooks/                             # Interactive GPU workflows
│   └── syn2nat_colab.ipynb                # End-to-end Google Colab binder generation & retrieval pipeline
│
├── bindpred/                              # ESMC representation extraction engine
│   ├── esmc_client.py                     # Biohub ESMC 6B/600M inference client
│   └── extract_representations.py         # Batch extraction of sequence vectors & residue matrices
│
├── scripts/                               # Executable clients and pipeline scripts
│   ├── run_5_targets_pilot.py             # Pilot binder retrieval execution
│   ├── run_15_targets_benchmark.py        # 15-target reverse lookup benchmark
│   ├── annotate_interface_sites.py        # Contact chemistry, active site proximity & ΔSASA extraction
│   └── ...                                # Interactome & foundation model clients
│
├── docs/                                  # Technical specifications & API manuals
│   ├── README.md                          # Documentation index and domain taxonomy
│   ├── latent_space_binder_matching_methodology.md # Syn2Nat methodology specification
│   ├── biohub.md                          # Biohub ESM / ESMFold2 platform documentation
│   ├── alphagenome/                       # DeepMind AlphaGenome multimodal DNA model guide
│   └── ...                                # Scored interactome reference guides
│
├── reports/                               # Active Syn2Nat benchmarks & reports
│   ├── README.md                          # Syn2Nat benchmark & reports portal
│   ├── interface_site_characterization_report.md # Master interface & active site mechanism report
│   ├── canonical_proteome_alphagenome_top1000.* # AlphaGenome Top 1,000 prioritized human targets
│   ├── pilot_5_targets_binder_search.*    # Pilot binder search benchmarks
│   ├── reverse_binder_lookup_15_targets.* # 15-target reverse binder matching benchmark
│   └── authentic_pockets_benchmark_results.* # Authentic pocket validation benchmark
│
└── data/                                  # Reference catalogs & curated indices
    ├── parsed/                            # Parsed smORF catalogs (orfs.parquet, orf_tiers.parquet)
    ├── screening/                         # AlphaGenome AVI purifying selection scores
    └── targets/                           # Canonical human proteome (UP000005640)
```

---

## 📊 Active Syn2Nat Benchmarks & Reports Hub

Active computational deliverables and structural benchmarks are maintained in [**`reports/`**](reports/README.md):

| Category | Benchmark / Report | Focus & Deliverable | Status |
| :--- | :--- | :--- | :--- |
| **Target Prioritization** | [**`canonical_proteome_alphagenome_top1000`**](reports/canonical_proteome_alphagenome_top1000.tsv) | **AlphaGenome Top 1,000 Target Selection**<br>Curated and prioritized high-value human targets from the AlphaGenome atlas based on dense purifying selection (AVI Phred scores), tissue epigenetics, and disease relevance. | `ACTIVE` |
| **Syn2Nat Retrieval** | [**`pilot_5_targets_binder_search`**](reports/pilot_5_targets_binder_search.tsv) | **Pilot De Novo Binder Retrieval**<br>End-to-end pilot validating latent-space and sequence-level matching of unconstrained synthetic binders against the microprotein catalog. | `BENCHMARK` |
| **Reverse Lookup** | [**`reverse_binder_lookup_15_targets`**](reports/reverse_binder_lookup_15_targets.tsv) | **15-Target Reverse Binder Matching**<br>Reverse lookup evaluation evaluating cross-target specificity and multi-modal ranking across 15 high-priority human targets. | `BENCHMARK` |
| **Multi-Scale Binders** | [**`multiscale_binder_pilot_results`**](reports/multiscale_binder_pilot_results.tsv) | **Multi-Scale Binder Pilot Evaluation**<br>Evaluation across varying unconstrained binder length regimes without artificial size forcing. | `BENCHMARK` |
| **Pocket Verification** | [**`authentic_pockets_benchmark_results`**](reports/authentic_pockets_benchmark_results.tsv) | **Authentic Pocket Structural Benchmark**<br>Rigorous structural pocket validation comparing binder-predicted interaction footprints with native binding clefts. | `BENCHMARK` |
| **Interface Characterization** | [**`interface_site_characterization_report.md`**](reports/interface_site_characterization_report.md) | **Master Interface Site & Mechanism Atlas (1,812 Complexes)**<br>Biophysical classification across 1,812 complexes: catalytic active site insertion ($\le 5.0\text{ \AA}$), pocket occlusion, buried surface area ($\Delta\text{SASA}$), and PDBe-KB experimental overlap. | `ACTIVE` |

> [!NOTE]
> **Exploratory Phase 1 Archival**: Early exploratory analyses (the initial 5 milestone reports on genomic selection, hitchhiking proofs, Top 50 curation, and preliminary 805-complex tables) have been uncommitted from git and safely archived locally under [`archive/exploratory_phase1/`](archive/) to maintain a lightweight, production-grade repository focused on **Syn2Nat**. Zero data was deleted.

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
