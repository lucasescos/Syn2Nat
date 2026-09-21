# ESM Atlas Public API Reference

The ESM Atlas provides open, keyless access to a public database indexing over **1.1 billion clustered protein structures** and 16,384 Sparse Autoencoder (SAE) functional features.

---

## 1. Base URL & Access Limits

* **Base URL**: `https://biohub.ai/esm/protein/api/v1alpha1`
* **Authentication**: None required (open public alpha).
* **Constraints**:
  * Sequence length $\le 2048$ aa.
  * Maximum 100 results per similarity search query.
  * Maximum 500 protein hashes per batch submission.
  * On-the-fly folding (`fold_on_miss=true`) available for sequences $<700$ aa.

---

## 2. Endpoints Overview

### 2.1 Similarity & Homology Search: `GET /similarity-search`
Searches the 1.1B database for structural and sequence homologs.
* **Query Parameters**:
  * `sequence` *(string, required)*: Protein sequence ($\le 2048$ aa).
  * `topk_results` *(int, default: 10, max: 100)*: Number of nearest matches to return.
  * `include_cluster_info` *(bool, default: false)*: Returns cluster size and representative hash.
* **Returns**:
  * `similar_proteins`: List of matches (accession, length, similarity score, pLDDT, name, cluster size).
  * `top_features_across_results`: Most frequent SAE features shared across top hits.

### 2.2 Protein Hash Lookup: `GET /proteins/{hash}`
Retrieves pre-computed structure, sequence, metadata, and active SAE features using sequence MD5 hash.
* **Query Parameters**:
  * `fold_on_miss` *(bool, default: false)*: Triggers on-the-fly folding if the structure was not previously cached ($<700$ aa).

### 2.3 UniProt Identifier Resolver: `GET /uniprot/{uniprot_id}`
Maps UniProt accession ID (e.g. `P00520`) to:
* Sequence & length
* Gene name & organism taxonomy
* Functional description
* Corresponding Atlas sequence MD5 hash

### 2.4 SAE Feature Dossier: `GET /features/{index}`
Fetches the biological dossier for any of the 16,384 SAE codebook features ($0 \le \text{index} \le 16383$):
* `label`: Human-readable functional title.
* `category`: Broad classification (e.g. enzymatic, structural, binding).
* `summary`: Curated biological description.
* `top_swissprot_activations`: Known UniProt proteins with the highest feature activation.
* `decoder_nearest_neighbors`: Semantically adjacent feature indices.

### 2.5 Structure Thumbnails: `GET /proteins/{hash}/thumbnail/{type}`
Returns a PNG thumbnail of the protein structure.
* `type`:
  * `"plddt"`: Colored by local structural confidence (blue = high, orange = low).
  * `"pct-characterized"`: Colored by sequence novelty / characterized residue ratio.

### 2.6 Batch Queries: `POST /proteins/batch`
Submits up to 500 hashes for asynchronous batch retrieval.
* **Status**: Immediate HTTP 200 (zip) or HTTP 202 (`job_id`).
* Polled via `GET /proteins/batch/jobs/{id}`.

---

## 3. CLI Script Usage (`scripts/atlas_search.py`)

```bash
# 1. Homology & similarity search
python .agents/skills/biohub-esm/scripts/atlas_search.py similarity "MSHHWGYGKH..." --topk 5

# 2. UniProt lookup to Atlas hash and annotation
python .agents/skills/biohub-esm/scripts/atlas_search.py uniprot P00520

# 3. Download structure thumbnail
python .agents/skills/biohub-esm/scripts/atlas_search.py thumbnail <protein_hash> --type plddt -o thumb.png

# 4. Inspect SAE functional feature dossier
python .agents/skills/biohub-esm/scripts/atlas_search.py feature 100

# 5. Batch lookup for multiple hashes
python .agents/skills/biohub-esm/scripts/atlas_search.py batch hashes.txt -o batch_results.zip
```
