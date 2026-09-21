# Human Microprotein Characterization & ESM Analysis

> **Master Project Reference**: See [**`PROJECT_MASTER_DOCUMENT.md`**](PROJECT_MASTER_DOCUMENT.md) for the comprehensive, single-file master document covering all datasets, computational screening pipelines, hitchhiking confounder proofs, the de-confounded Top 50 portfolio, 805 cofolded complexes, and biophysical evaluation metrics.

Research repository for exploring, cataloging, and functionally characterizing human microproteins and peptideins derived from small open reading frames (smORFs) using Evolutionary Scale Modeling (ESM), ESM Atlas clustering, and Pfam domain annotations.

---

## Reference Literature

* **Primary Study**: *"Expanding the human proteome with microproteins and peptideins"*
  * **Citation**: Nature 2026, DOI: [10.1038/s41586-026-10459-x](https://doi.org/10.1038/s41586-026-10459-x)
  * **PDF & Supplements**: Stored under [`literature/`](literature/)
* **MAPPIE Framework**: *"A map of human protein-protein interaction embeddings for functional discovery"*
  * **Citation**: bioRxiv 2026.08.07.743440, DOI: [10.64898/2026.08.07.743440](https://doi.org/10.64898/2026.08.07.743440)
  * **PDF**: Stored under [`literature/2026.08.07.743440v1.full.pdf`](literature/2026.08.07.743440v1.full.pdf)
  * **Web Server**: [cbdm-01.zdv.uni-mainz.de/~mcihan/mappie/](https://cbdm-01.zdv.uni-mainz.de/~mcihan/mappie/)
  * **Code & Guide**: Cloned under [`mappie/`](mappie/) | Documentation: [`docs/mappie.md`](docs/mappie.md)

---

## Directory Layout

```
microprotein/
├── data/                                  # Working datasets, parsed tables, and profile HMMs
│   ├── mappie/                            # MAPPIE datasets (UMAP table, 128-d latent space, latent_index.npz)
│   ├── parsed/                            # Parsed parquet datasets (orfs.parquet, orf_tiers.parquet)
│   ├── pfam/                              # Pfam-A profile HMMs (Pfam-A.hmm.gz, Pfam-A.hmm.h3f)
│   ├── queries/                           # Curated microprotein query sets (Tier 1 fastas & metadata)
│   ├── targets/                           # Canonical human proteome (UP000005640 parquet & fasta)
│   ├── screening/                         # AlphaGenome variant impact (AVI) scores & batch files
│   └── raw/                               # Raw supplementary tables (supp_table11_orbl.xlsx)
│
├── mappie/                                # Cloned upstream MAPPIE framework (github.com/mcihan0bioinf/mappie)
│   ├── core_algorithm/                    # Novel PPI projection & hypergeometric enrichment
│   ├── dark_interactome/                  # Dark hub candidate prediction scripts
│   ├── benchmark/                         # Benchmark & validation suites
│   ├── training/                          # Autoencoder & UMAP training pipelines
│   └── data_processed/                    # Preprocessed runtime indices (latent_index.npz)
│
├── docs/                                  # Technical and API documentation
│   ├── mappie.md                          # MAPPIE framework guide: continuous latent space & functional discovery
│   ├── esmfold2_latent_interface_screening.md # Framework for deorphanizing GENCODE peptideins
│   ├── biohub.md                          # Biohub ESM API, SDK cheat sheet & ESM Atlas search reference
│   └── alphagenome/                       # AlphaGenome regulatory variant effect prediction documentation
│
├── literature/                            # Research papers and supplementary data
│   ├── 2026.08.07.743440v1.full.pdf       # MAPPIE preprint (bioRxiv 2026)
│   ├── s41586-026-10459-x.pdf             # Nature 2026 primary research article
│   └── supplementary/                     # Supplementary tables & figures (MOESM1 through MOESM18)
│
├── plans/                                 # Research frameworks and multi-agent virtual lab protocols
│   └── paperclip_agent_framework.md       # Virtual Lab & Paperclip orchestration framework
│
├── scripts/                               # Pipeline execution scripts
│   ├── phase1_curate_datasets.py          # Phase 1: Query & canonical target curation
│   └── step1_extract_intervals.py ...     # AlphaGenome genomic variant screening steps
│
├── .agents/skills/biohub-esm/             # Antigravity agent skill & Python scripts for ESM operations
│   ├── SKILL.md                           # Skill definition and instructions
│   ├── references/api_reference.md        # Technical API specification
│   └── scripts/                           # Tool scripts (atlas_search.py, fold.py, embed.py, etc.)
│
├── .env                                   # Local environment variables (e.g. ESM_API_KEY)
└── .gitignore                             # Git configuration
```

---

## Documentation & Tooling

* **BioGRID Interactome Web Service**: See [`docs/biogrid.md`](docs/biogrid.md) for full BioGRID REST API endpoints, parameter references, TAB2/JSON schemas, and Python client integration ([`scripts/biogrid_client.py`](scripts/biogrid_client.py)).
* **HIPPIE Human Scored Interactome**: See [`docs/hippie.md`](docs/hippie.md) for experiment-based confidence scoring ($0.0 \le S \le 1.0$), GTEx tissue filtering, and programmatic client ([`scripts/hippie_client.py`](scripts/hippie_client.py)).
* **MAPPIE Functional Annotation Engine**: See [`docs/mappie.md`](docs/mappie.md) for latent space projection of microprotein PPIs, cosine nearest-neighbor search, and hypergeometric enrichment over reference human interactomes.
* **Deorphanization Framework**: See [`docs/esmfold2_latent_interface_screening.md`](docs/esmfold2_latent_interface_screening.md) for the 6-phase screening pipeline using ESMC 6B 2D pairwise latent representations, ESMFold2-Fast structural validation, and AlphaGenome Atlas genomic variant integration.
* **Biohub ESM & ESM Atlas Platform Guide**: See [`docs/biohub.md`](docs/biohub.md) for endpoints, `curl` recipes, Python SDK patterns, SAE feature analysis, and ESMFold2 configuration.
* **AlphaGenome Regulatory Model**: See [`docs/alphagenome/README.md`](docs/alphagenome/README.md) for DeepMind's multimodal genomic variant effect model.
* **Agent Automation**: Use the `.agents/skills/biohub-esm/` skill for automated ESM embedding, folding, and cluster search workflows.

