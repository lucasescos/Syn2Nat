# MAPPIE: Map of Protein–Protein Interaction Embeddings

> **Local Path**: [`mappie/`](../mappie/)  
> **Datasets Path**: [`data/mappie/`](../data/mappie/)  
> **Web Server**: [https://cbdm-01.zdv.uni-mainz.de/~mcihan/mappie/](https://cbdm-01.zdv.uni-mainz.de/~mcihan/mappie/)  
> **Preprint**: Mert Cihan, Ute Distler, Miguel A Andrade-Navarro, *"A map of human protein-protein interaction embeddings for functional discovery"*, bioRxiv 2026.08.07.743440 ([DOI: 10.64898/2026.08.07.743440](https://doi.org/10.64898/2026.08.07.743440))  
> **Preprint PDF**: [`literature/2026.08.07.743440v1.full.pdf`](../literature/2026.08.07.743440v1.full.pdf)  
> **Underlying PPI Database**: [HIPPIE Documentation (`docs/hippie.md`)](hippie.md)

---

## 1. Overview & Theoretical Framework

MAPPIE is a deep learning framework designed to infer the function of protein–protein interactions (PPIs) independent of network connectivity. 

### The Network Sparsity Problem
Traditional guilt-by-association methods infer the function of a protein by examining its direct neighbors in a physical PPI network. However:
1. **The "Dark Interactome"**: Thousands of uncharacterized proteins (including novel smORF-derived microproteins and peptideins) have zero or very few documented interaction partners.
2. **Topology Bias**: Well-studied hubs (e.g., TP53, UBC) dominate graph-based algorithms, introducing extreme degree-related biases.
3. **Set-Level Limitation**: Standard functional enrichment looks at sets of isolated proteins rather than the specific functional context of an individual interaction interface.

### MAPPIE's Architecture
MAPPIE circumvents network topology by learning a continuous latent space directly from sequence representations:
1. **Sequence Embeddings**: Proteins are encoded using the **ESM-2** protein language model (`esm2_t33_650M_UR50D`, 1280-dimensional mean-pooled representations).
2. **Pairwise Merging**: Given proteins $A$ and $B$, their embeddings $e_A$ and $e_B$ are merged element-wise via multiplication ($e_A \odot e_B$).
3. **Autoencoder Compression**: A symmetrical deep autoencoder projects the merged vector into a compact **128-dimensional latent space** ($z \in \mathbb{R}^{128}$).
4. **Metric Space Proximity**: In this latent space, PPIs that share functional mechanisms cluster closely together, even if the proteins never interact with each other in nature.
5. **Hypergeometric Functional Enrichment**: Functional terms (GO, Reactome pathways, Pfam domains, InterPro families, Rhea reactions, and 3did domain–domain interactions) are evaluated for statistical enrichment across the $k$-nearest neighbor interactions in the latent manifold.

```
Protein A (Seq) ──> [ ESM-2 650M ] ──> e_A (1280-d) ┐
                                                     ├─> Element-wise Multiply ──> [ Autoencoder ] ──> Latent Vector z (128-d) ──> UMAP (2D)
Protein B (Seq) ──> [ ESM-2 650M ] ──> e_B (1280-d) ┘                                                                         │
                                                                                                                               └──> k-NN Search
                                                                                                                                     └──> Functional Enrichment
```

---

## 2. Directory & Dataset Organization

All code and data files are organized locally within this project:

```
microproteinproject/
├── mappie/                                # Cloned upstream GitHub repository
│   ├── annotation_db/                     # SQLite database creation scripts
│   ├── benchmark/                         # CORUM complex recovery, degree stratification, STRING comparisons
│   ├── core_algorithm/                    # Core inference & enrichment scripts
│   │   ├── project_ppi.py                 # Novel PPI projection into latent & UMAP space
│   │   ├── run_enrichment.py              # Hypergeometric term enrichment engine
│   │   └── enrichment_plots.py            # Visualizations (dotplots, heatmaps, bar charts)
│   ├── dark_interactome/                  # Screening pipelines for dark hub candidate discovery
│   ├── embedding/                         # ESM-2 embedding generation & merging utilities
│   ├── model_selection/                   # DDI density matrix evaluation & latent space optimization
│   ├── training/                          # Autoencoder & UMAP training pipelines
│   └── data_processed/                    # Preprocessed binary runtime indices
│       └── latent_index.npz               # Binary vector index (keys & vecs) for script execution
│
├── data/
│   └── mappie/                            # Bulk reference datasets downloaded from web server
│       ├── mappie_umap_table.csv          # 14 MB: 2D UMAP coordinates & HIPPIE scores for 199,138 PPIs
│       ├── mappie_latent_space.csv.gz     # 124 MB: Full 128-d latent space representations (gzipped)
│       └── latent_index.npz               # Compiled NumPy archive (keys, vecs) for zero-latency queries
│
└── docs/
    └── mappie.md                          # This reference manual
```

### Dataset Schemas

#### 1. UMAP Coordinates Table (`data/mappie/mappie_umap_table.csv`)
* **Size**: 14 MB (199,139 rows, header + 199,138 PPIs)
* **Columns**:
  * `ppikey`: Unique interaction identifier formatted as `<PROTEIN_A>_<PROTEIN_B>` (e.g. `MDM2_HUMAN_P53_HUMAN`).
  * `proteinA`: UniProt mnemonic for partner A.
  * `proteinB`: UniProt mnemonic for partner B.
  * `hippie_score`: Experimental confidence score from HIPPIE v2.3 ($0.65 \le s \le 1.0$).
  * `UMAP1`: First 2D projection coordinate.
  * `UMAP2`: Second 2D projection coordinate.

#### 2. Latent Space Table (`data/mappie/mappie_latent_space.csv.gz`)
* **Size**: 124 MB gzipped (199,139 rows)
* **Columns**:
  * `ppikey`: Interaction identifier.
  * `dim_0` to `dim_127`: 128 float values representing the compressed interaction latent vector.

#### 3. Binary Vector Index (`data/mappie/latent_index.npz`)
* **Format**: Compressed NumPy archive (`.npz`)
* **Arrays**:
  * `keys`: 1D array of strings (`shape=(199138,)`), interaction keys.
  * `vecs`: 2D `float32` array (`shape=(199138, 128)`), L2-normalized latent vectors ready for cosine nearest-neighbor search.

---

## 3. Practical Usage & Code Recipes

### Recipe 1: Quick Data Loading with Pandas
```python
import pandas as pd

# Load UMAP coordinates
df_umap = pd.read_csv("data/mappie/mappie_umap_table.csv")
print(f"Loaded {len(df_umap):,} interactions with UMAP coordinates.")

# Load 128-D latent embeddings
df_latent = pd.read_csv("data/mappie/mappie_latent_space.csv.gz", compression="gzip")
print(f"Latent matrix shape: {df_latent.shape}")
```

### Recipe 2: Zero-Latency Nearest-Neighbor Queries
```python
import numpy as np
from sklearn.neighbors import NearestNeighbors

# Load precompiled index
data = np.load("data/mappie/latent_index.npz", allow_pickle=True)
vecs, keys = data["vecs"], data["keys"]
key_to_idx = {k: i for i, k in enumerate(keys)}

# Fit cosine neighbor search engine
nn = NearestNeighbors(n_neighbors=10, metric="cosine", algorithm="brute")
nn.fit(vecs)

# Query neighbors for a known interaction (e.g., P53-MDM2)
query_key = "MDM2_HUMAN_P53_HUMAN"
if query_key not in key_to_idx:
    query_key = "P53_HUMAN_MDM2_HUMAN"

q_idx = key_to_idx[query_key]
distances, indices = nn.kneighbors(vecs[q_idx : q_idx + 1])

print(f"Top 10 Latent Neighbors for {query_key}:")
for rank, (idx, dist) in enumerate(zip(indices[0], distances[0]), start=1):
    print(f" {rank:2d}. {keys[idx]:<35} (Cosine distance: {dist:.4f})")
```

### Recipe 3: Projecting Novel Microprotein Interactions
To project a novel candidate pair (e.g., a microprotein and a canonical human target) into MAPPIE's latent space:

```python
import sys
from pathlib import Path
sys.path.insert(0, str(Path("mappie").resolve()))

from core_algorithm.project_ppi import project_single_ppi

# Option A: From UniProt IDs (for canonical partners)
# Option B: From raw amino acid sequences (for novel microproteins)
microprotein_seq = "MGNKLSIKDLFRR..." # Novel smORF translation
target_uniprot = "P04637"           # Human p53

# Computes ESM-2 embedding, multiplies, feeds through Autoencoder and UMAP
result = project_single_ppi(
    prot_a=microprotein_seq,
    prot_b=target_uniprot,
    is_seq_a=True,
    is_seq_b=False
)

print("Latent Vector (128-d):", result["latent_vector"])
print("UMAP Coordinates:", result["umap_x"], result["umap_y"])
```

---

## 4. Synergy with Human Microprotein Screening

In this workspace, we deorphanize smORF-derived microproteins using:
1. **ESMC 6B pairwise latent representations** (`bindpred/`) to detect interaction interfaces.
2. **ESMFold2** multi-chain structural co-folding to score physical interface confidence ($pLDDT$, $pTM$, $ipTM$).
3. **AlphaGenome Variant Impact (AVI)** scores to prioritize disease-associated non-coding regulatory smORFs.

### Where MAPPIE Fits:
Once a high-confidence microprotein–target pair is identified by ESMFold2 ($ipTM \ge 0.75$), MAPPIE provides **interaction-level functional annotation**:
* By projecting the microprotein–target complex into MAPPIE's latent space, we find the reference human PPIs that occupy the same structural/functional niche.
* Hypergeometric enrichment over its $k$-nearest neighbors reveals the likely pathway (e.g. *Spliceosome assembly*, *Ribosome biogenesis*, or *Ubiquitin ligase complex*) without requiring previous biochemical studies of the microprotein.
