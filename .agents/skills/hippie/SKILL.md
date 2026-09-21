---
name: hippie
description: >-
  Query and retrieve human protein-protein interactions (PPIs) with experimental confidence scores
  from the HIPPIE database (Human Integrated Protein-Protein Interaction rEference v2.4, CBDM Group / JGU Mainz).
  Use to validate physical interactions, check interaction scores between protein pairs, discover binding partners,
  construct Layer 0/1 subnetworks, and filter by confidence thresholds or GTEx tissue specificity.
---

# HIPPIE: Human Integrated Protein-Protein Interaction rEference

This skill equips agents with direct access to the **HIPPIE (v2.4)** database and REST API to query confidence-scored human protein-protein interactions.

---

## 1. When to Use HIPPIE

* **PPI Validation**: Determine whether two human proteins (e.g. a microprotein candidate target and a cellular receptor/kinase) have experimentally verified physical interactions.
* **Confidence Scoring**: Retrieve evidence-weighted confidence scores ($0.00$ to $1.00$) calculated from experimental techniques (X-ray, NMR, SPR, Y2H, etc.), independent publication counts, and orthology.
* **Neighborhood Discovery**: Find all high-confidence interactors for one or more proteins of interest.
* **Subnetwork Extraction**: Construct internal subgraphs (Layer 0) among a candidate set of genes to discover interactome modules or functional complexes.
* **Benchmarking & Validation**: Serve as experimental ground truth for computational interaction predictions (e.g., ESMFold2 interface screening, BindPred, MAPPIE embeddings).

---

## 2. Confidence Threshold Guidelines

HIPPIE confidence scores range from $0.00$ to $1.00$:

| Score Range | Tier | Quartile | Recommended Usage |
| :--- | :--- | :--- | :--- |
| $\ge 0.72$ (or $0.73$) | **High Confidence** | Q3 (top 25%) | Stringent validation, structural modeling, training positive benchmarks. |
| $\ge 0.63$ | **Medium Confidence** | Q2 (median) | Exploratory screening, candidate target identification. |
| $\ge 0.65$ | **MAPPIE Benchmark** | – | Cutoff utilized by MAPPIE autoencoder training. |
| $< 0.63$ | **Low / Exploratory** | Q1 | Single low-throughput screen or high-throughput Y2H with single evidence. |

---

## 3. Pre-built CLI & Python Tools

The client script is available at [`scripts/hippie_client.py`](../../../scripts/hippie_client.py).

### CLI Recipes

```powershell
# 1. Query whether two proteins interact and get their confidence score
python scripts/hippie_client.py --proteins TP53,MDM2 --pair

# 2. Retrieve high-confidence interactors for a protein (score >= 0.72)
python scripts/hippie_client.py --proteins TP53 --conf 0.72

# 3. Construct an internal subnetwork (Layer 0) among candidate genes
python scripts/hippie_client.py --proteins TP53,MDM2,CDKN1A,ATM --layer 0 --conf 0.63

# 4. Export query results to TSV or JSON
python scripts/hippie_client.py --proteins BRCA1,BARD1 --format mitab --output data/targets/brca_network.tsv
```

### Python API Recipes

```python
from scripts.hippie_client import HippieClient, query_pair, get_interactors

client = HippieClient()

# Query a pair
result = query_pair("TP53", "MDM2")
if result:
    print(f"Interaction Score: {result['score']}")

# Query interactors with DataFrame output
df = client.query("CDKN1A", layer=1, conf_thres=0.63)
print(df[["gene_1", "gene_2", "score"]].head())

# Query Layer 0 internal connections
genes = ["TP53", "MDM2", "CDKN1A", "BAX", "PML"]
subnetwork = client.query_subnetwork(genes, conf_thres=0.63)
```

---

## 4. Reference Resources

* Comprehensive Documentation: [`docs/hippie.md`](../../../docs/hippie.md)
* Experimental Scoring Weights: [`data/hippie/experimental_scores.tsv`](../../../data/hippie/experimental_scores.tsv)
* Offline Java Tools Manuals: [`data/hippie/README_NC.txt`](../../../data/hippie/README_NC.txt) and [`data/hippie/README_RS.txt`](../../../data/hippie/README_RS.txt)
* Upstream Service: `https://cbdm-01.zdv.uni-mainz.de/~mschaefer/hippie/`
