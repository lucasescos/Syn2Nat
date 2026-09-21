# Comprehensive Biohub ESM API & SDK Reference

Complete technical specification for interacting with Biohub's Evolutionary Scale Modeling (ESM) platform (`https://biohub.ai`).

> [!TIP]
> **Modular Reference Guides**: For dedicated, in-depth documentation by workflow, see:
> * [**Python SDK & `ESMProtein`**](python_sdk_reference.md): Official `pip install esm`, `ESM3ForgeInferenceClient`, `GenerationConfig`, and PDB methods.
> * [**ESM3 Generative Design**](esm3_generation_and_design.md): Multi-track generation (`sequence`, `secondary_structure`, `structure`, `function`), sampling schedules.
> * [**Structure Prediction & Complexes**](structure_prediction_and_complexes.md): ESMFold2, multimers, all-atom complexes with DNA/RNA/ligands.
> * [**Representations & Interpretability**](representations_and_interpretability.md): ESMC embeddings, logits, SAE codebook features, DMS entropy.
> * [**Inverse Folding & Tokenizers**](inverse_folding_and_tokenizers.md): Coordinates to sequence, track tokenization, raw tensor generation.
> * [**ESM Atlas Public API**](esm_atlas_public_api.md): Keyless public database (>1.1B structures, homology search, UniProt lookups).

---

## 1. Authentication & Base URLs

### Platform Endpoints (Requires API Key)
* **Base URL**: `https://biohub.ai/api/v1`
* **Auth Header**: `Authorization: Bearer <ESM_API_KEY>`
* **Key Source**: [Biohub Developer Console](https://biohub.ai/developer-console/api-keys)
* **Environment Variable**: `ESM_API_KEY` (in `~/.env` or `.env`)

### ESM Atlas Endpoints (Public, Alpha, Keyless)
* **Base URL**: `https://biohub.ai/esm/protein/api/v1alpha1`
* **Auth**: None required (open alpha API)
* **Caps**: Sequences $\le 2048$ aa, up to 100 results per similarity query, up to 500 hashes per batch, on-the-fly folding for $<700$ aa.

---

## 2. Supported Models

| Model ID | Family | Purpose | Layers | Dimension |
| :--- | :--- | :--- | :--- | :--- |
| `esmc-300m-2024-12` | ESMC | Fast sequence embedding & scoring | 30 | 960 |
| `esmc-600m-2024-12` | ESMC | High-accuracy sequence representation | 36 | 1152 |
| `esmc-6b-2024-12` | ESMC | Frontier PLM, world model, SAE base | 80 | 2560 |
| `esmfold2-fast-2026-05`| ESMFold2 | Default structure prediction (fast) | — | — |
| `esmfold2-2026-05` | ESMFold2 | Full structure prediction (supports MSAs) | — | — |
| `esmfold2-2026-05-cutoff-2025` | ESMFold2 | Server-side batch client folding | — | — |
| `esm3-open-2024-03` | ESM3 | Generative design, unmasking, inverse fold | — | — |

---

## 3. Platform Endpoints (`POST https://biohub.ai/api/v1/…`)

### 3.1 `POST /fold` — Protein Sequence & Backbone Structure
Folds single or multi-chain protein sequences into 3D coordinates.
* **Payload**:
  * `model`: `"esmfold2-fast-2026-05"` or `"esmfold2-2026-05"`
  * `sequences`: `[{"id": "A", "sequence": "..."}]` or `{"sequence": "CHAIN_A|CHAIN_B"}`
  * `config`:
    * `num_loops`: `20` (0–20)
    * `num_sampling_steps`: `100` (1–100)
    * `lm_dropout`: `0.3`
    * `lm_mask_pct`: `0.0`
    * `include_pae`: `bool`
    * `include_distogram`: `bool`
    * `include_pair_chains_iptm`: `bool`
    * `include_embeddings`: `bool`
  * `potential_sequence_of_concern`: `false`

### 3.2 `POST /fold_all_atom` — Multi-Entity Macromolecular Complexes
Predicts structures of proteins co-folded with DNA, RNA, Small-Molecule Ligands, and Chemical Modifications.
* **Payload**:
  * `all_atom_input`:
    * `sequences`: Array of entities:
      * `ProteinInput`: `{"id": "A", "sequence": "...", "type": "protein"}`
      * `DNAInput`: `{"id": "B", "sequence": "GATCGATC", "type": "dna"}`
      * `RNAInput`: `{"id": "C", "sequence": "GAUCGAUC", "type": "rna"}`
      * `LigandInput`: `{"id": "L", "ccd": ["SAH"], "type": "ligand"}` or `{"id": "L", "smiles": "...", "type": "ligand"}`
    * `pocket`: `PocketConditioning` (`binder_chain_id`, `contacts`) for epitope targeting.
    * `covalent_bonds`: `CovalentBond` (`chain_id1`, `res_idx1`, `atom_idx1`, `chain_id2`, `res_idx2`, `atom_idx2`).
    * `distogram_conditioning`: Distogram distance constraints.
  * Parameters match `/fold`.

### 3.3 `POST /logits` — Representations, Hidden States & SAEs
Extracts sequence logits, residue embeddings, layer-wise representations, and Sparse Autoencoder (SAE) features.
* **Payload**:
  * `model`: `"esmc-600m-2024-12"`, `"esmc-6b-2024-12"`, etc.
  * `protein`: `{"sequence": "..."}`
  * `config`:
    * `sequence`: `true`
    * `return_embeddings`: `bool` (per-residue $L \times D$)
    * `return_mean_embedding`: `bool` (pooled $1 \times D$)
    * `return_hidden_states`: `bool`
    * `return_mean_hidden_states`: `bool`
    * `ith_hidden_layer`: Integer layer index (`-2` recommended default; `-1` for all layers on 300M/600M)
    * `sae_config`: `{"models": ["esmc-6b-2024-12-sae-layer60-k64-codebook16384"]}`

### 3.4 `POST /generate` — Generative Multi-Track Design
Generative protein design via iterative unmasking using ESM3.
* **Tracks**: `sequence`, `structure`, `secondary_structure`, `sasa`, `function`.
* **Payload**:
  * `model`: `"esm3-open-2024-03"`
  * `track`: Track name to generate (e.g. `"sequence"`)
  * `inputs`: Seed `Tracks` object (e.g. `{"sequence": "M___G__L___K"}`)
  * `num_steps`: Iterative unmasking steps (default `20`)
  * `temperature`: `0.5`
  * `temperature_annealing`: `true`
  * `top_p`: `1.0`
  * `schedule`: `"cosine"` or `"linear"`
  * `strategy`: `"random"` or `"entropy"`
  * `condition_on_coordinates_only`: `bool`

### 3.5 `POST /inverse_fold` — Backbone to Sequence
Recovers or designs amino acid sequences matching 3D backbone coordinates.
* **Payload**:
  * `model`: `"esm3-open-2024-03"`
  * `coordinates`: Array of backbone 3D coordinates ($N \times 3$ or $N \times 37 \times 3$)
  * `inverse_folding_config`:
    * `temperature`: Sampling temperature (default `0.1`)
    * `invalid_ids`: Token IDs to disallow

### 3.6 `POST /forward_and_sample` — Single Pass Sampling
Executes a single model forward pass with per-track multinomial sampling.

### 3.7 `POST /encode` & `POST /decode` — Tokenizer Transformations
Converts between biological continuous/string tracks (`Tracks`) and discrete token indices (`Tokens`).

### 3.8 `POST /generate_tensor` — Raw Token Generation
Operates directly on discrete token tensors for custom unmasking protocols.

---

## 4. ESM Atlas Endpoints (`https://biohub.ai/esm/protein/api/v1alpha1/…`)

Public endpoints indexing >1.1 billion clustered structures:

1. `GET /similarity-search`:
   * Params: `sequence` (req, $\le 2048$ aa), `topk_results` (10, max 100), `include_cluster_info` (bool).
   * Returns: `similar_proteins` (accession, length, similarity score, pLDDT, name, cluster size) and `top_features_across_results`.
2. `GET /proteins/{hash}`:
   * Lookup by sequence MD5 hash. `fold_on_miss=true` folds on-the-fly for sequences $<700$ aa.
3. `POST /proteins/batch`:
   * Batch lookup for up to 500 hashes: `include_structure`, `include_cluster_info`, `include_sequence`.
   * Returns immediate zip file (200) or job status (202).
4. `GET /proteins/batch/jobs/{id}` & `DELETE /proteins/batch/jobs/{id}`:
   * Polls or cancels queued batch jobs.
5. `GET /proteins/{hash}/thumbnail/{type}`:
   * Returns structure PNG thumbnail. `type`: `plddt` or `pct-characterized`.
6. `GET /clusters/{protein_hash}`:
   * Cluster size, representative hash, and characterization ratio.
7. `GET /features`:
   * Bulk unpaginated list of all 16,384 SAE features.
8. `GET /features/{index}`:
   * Feature dossier: label, category, summary, threshold, top SwissProt activations, decoder nearest neighbors.
9. `GET /uniprot/{uniprot_id}`:
   * Maps UniProt Accession to sequence, gene, organism, functional description, and Atlas protein hash.

---

## 5. Python SDK Classes (`esm.sdk`)

* `esm.sdk.client(model, url, token)`: ESM3 client.
* `esm.sdk.esmc_client(model, url, token)`: ESMC client.
* `esm.sdk.esmfold2_client(model, url, token)`: ESMFold2 client.
* `esm.sdk.parallel_executor()`: Context manager for parallel batch job processing.
* `esm.sdk.ForgeBatchClient`: High-throughput asynchronous batch submission client.
* `esm.utils.structure.input_builder`: `ProteinInput`, `DNAInput`, `RNAInput`, `LigandInput`, `Modification`, `CovalentBond`, `PocketConditioning`.
