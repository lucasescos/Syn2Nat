# Project Technical Documentation & API Manuals

This directory contains technical specifications, computational frameworks, database integration manuals, and client guides for all tools and resources integrated into this research project.

---

## Documentation Architecture

```
docs/
├── README.md                                 # This navigation hub
│
├── [Domain 1: Latent Space & Methodological Frameworks]
│   ├── mappie.md                             # 128-d latent PPI space & functional discovery engine
│   ├── esmfold2_latent_interface_screening.md# ESMC 6B pairwise screening & deorphanization architecture
│   ├── latent_space_binder_matching_methodology.md # De novo binder-guided latent search for target deorphanization
│   ├── ppi_evaluation_hierarchy.md           # 4-pillar biophysical evaluation (ipTM, PAE, contacts)
│   └── representation_pipeline_roadmap.md    # HMPA 617k smORF curation & quota-managed representation engine
│
├── [Domain 2: Curated Physical & Functional Interactomes]
│   ├── hippie.md                             # HIPPIE scored human interactome (MAPPIE training foundation)
│   ├── biogrid.md                            # BioGRID HTP/LTP experimental evidence & REST API client
│   ├── intact.md                             # IntAct / IMEx molecular curation, stoichiometry & interfaces
│   ├── string.md                             # STRING v12.0 functional association network & client
│   └── pepbind.md                            # PepBind curated peptide–protein interface database
│
├── [Domain 3: Foundation Models & AI Platforms]
│   ├── biohub.md                             # Biohub ESM platform: ESMC 6B, ESMFold2 & ESM Atlas guide
│   └── alphagenome/                          # DeepMind AlphaGenome 1-Mb multimodal model (14-part guide)
│       └── README.md                         # AlphaGenome index & navigation
│
└── [Archive & Raw Scrapes]
    ├── archive/                              # Superseded exploratory proposals
    │   └── alphagenome_exploratory_idea.md   # Initial 3-step exploratory screening scratchpad
    └── hippie_reference/                     # Raw reference documentation from HIPPIE web portal
        ├── download.md                       # Download specifications
        └── information.md                    # Database scoring methodologies
```

---

## Domain 1: Latent Space & Methodological Frameworks

| Guide | Description | Key Modules / Implementation |
| :--- | :--- | :--- |
| [**`latent_space_binder_matching_methodology.md`**](latent_space_binder_matching_methodology.md) | **Syn2Nat: De Novo Binder-Guided Latent Search for Target Deorphanization**<br>Inverted high-throughput discovery architecture: prioritizes targets from AlphaGenome Top 1,000, generates unconstrained *de novo* binders via ESMFold2, projects into ESMC latent space, conducts multi-modal sequence/embedding searches, and validates top microproteins. | • Pipeline: `Syn2Nat` (Synthetic-to-Natural)<br>• Design: `binder_design.py`<br>• Representation: `ESMC-6B` / `ESMC-600M`<br>• Search: Latent + Sequence multi-modal |
| [**`mappie.md`**](mappie.md) | **MAPPIE Latent PPI Manifold & Functional Discovery**<br>Details the 128-dimensional continuous latent space learned by deep autoencoder compression of ESM-2 pairwise embeddings ($e_A \odot e_B$) across 199,138 reference human PPIs. Enables unsupervised functional discovery without network topology. | • Local code: [`mappie/`](../mappie/)<br>• Latent index: [`data/mappie/latent_index.npz`](../data/mappie/latent_index.npz)<br>• Web: [Uni-Mainz MAPPIE](https://cbdm-01.zdv.uni-mainz.de/~mcihan/mappie/) |
| [**`esmfold2_latent_interface_screening.md`**](esmfold2_latent_interface_screening.md) | **GENCODE Peptidein Deorphanization Architecture**<br>Outlines the 6-phase pipeline that decouples the frozen 80-layer ESMC 6B language model trunk from generative diffusion to achieve scalable cross-chain interface scoring before full 3D rollout. | • Concept: Semantic $\mathcal{O}(1)$ interface filter<br>• Trunk: `esmc-6b-2024-12`<br>• Fold: `esmfold2-fast-2026-05` |
| [**`ppi_evaluation_hierarchy.md`**](ppi_evaluation_hierarchy.md) | **The 4-Pillar Biophysical Interaction Hierarchy**<br>Establishes objective biophysical criteria for classifying non-canonical complexes: Interface TM-score ($\text{ipTM} \ge 0.60$), inter-chain PAE certainty ($< 10\text{ Å}$), induced folding footprint, and atomic contact chemistry. | • Benchmarks: Nanomolar-like complexes<br>• Validated candidates: `SNX13`–`SERPINE2`, `AP001372.2`–`PCNA` |
| [**`representation_pipeline_roadmap.md`**](representation_pipeline_roadmap.md) | **Large-Scale Microprotein Representation Engine**<br>Defines the active project roadmap: curating 617,462 smORFs from the Human Microprotein Atlas (HMPA) and extracting penultimate layer 79 ESMC 6B embeddings under automated, daily quota-managed batches. | • Data source: HMPA (Zenodo 19547725)<br>• Pipeline: [`bindpred/`](../bindpred/) |

---

## Domain 2: Curated Physical & Functional Interactomes

Each interactome database has a corresponding programmatic client located under [`scripts/`](../scripts/):

| Database | Documentation | Client Script | Primary Role & Unique Capability |
| :--- | :--- | :--- | :--- |
| **HIPPIE v2.3** | [**`hippie.md`**](hippie.md) | [`scripts/hippie_client.py`](../scripts/hippie_client.py) | **Experimental confidence scoring** ($0.0 \le S \le 1.0$) and GTEx tissue specificity. Serves as the ground-truth corpus behind MAPPIE's 199,138 indexed PPIs. |
| **BioGRID v4.4** | [**`biogrid.md`**](biogrid.md) | [`scripts/biogrid_client.py`](../scripts/biogrid_client.py) | **Experimental evidence breadth** across HTP and LTP systems (Affinity Capture-MS, Two-Hybrid, Co-fractionation). |
| **IntAct / IMEx** | [**`intact.md`**](intact.md) | [`scripts/intact_client.py`](../scripts/intact_client.py) | **Deep molecular curation**: PSI-MI compliance, binding interfaces, mutation effects, and complex stoichiometry. |
| **STRING v12.0** | [**`string.md`**](string.md) | [`scripts/string_client.py`](../scripts/string_client.py) | **Broad functional association networks** integrating physical interaction evidence with text-mining, co-expression, and genomic neighborhood. |
| **PepBind** | [**`pepbind.md`**](pepbind.md) | [`scripts/pepbind_client.py`](../scripts/pepbind_client.py) | **Peptide–protein interfaces**: Structurally resolved and experimentally tested complexes focused on short linear motifs (SLiMs). |

---

## Domain 3: Foundation Models & AI Platforms

| Platform / Model | Documentation | Primary Application |
| :--- | :--- | :--- |
| **Biohub ESM & ESMFold2** | [**`biohub.md`**](biohub.md) | Reference for Biohub's API endpoints (`fold_all_atom`, `representations`), SDK patterns, and ESM Atlas (>1.1B structures). |
| **DeepMind AlphaGenome** | [**`alphagenome/README.md`**](alphagenome/README.md) | Comprehensive 14-chapter developer guide for Google DeepMind's 1-Mb multimodal genomic sequence model. Used in Phase 1 to calculate dense purifying selection (AVI Phred scores) and prove the hitchhiking effect. |

---

## Archive & Reference Scrapes

* [**`archive/alphagenome_exploratory_idea.md`**](archive/alphagenome_exploratory_idea.md): Historical scratchpad of the initial 3-step genomic screening proposal (completed and superseded by [`reports/microprotein_multiomics_analysis_report.md`](../reports/microprotein_multiomics_analysis_report.md)).
* [**`hippie_reference/`**](hippie_reference/): Direct documentation and download specifications scraped from the HIPPIE web portal (`download.md`, `information.md`).
