# Biohub ESM platform — API + SDK cheat sheet

Compiled 2026-09-01 against the live `https://biohub.ai/openapi.json` and `esm` SDK v3.4.0.
Full formatted guide: https://claude.ai/code/artifact/ba16e8e3-8969-4dd5-8c27-1885ada082c1

SOURCING — main API endpoints/fields/defaults come from the machine-readable OpenAPI 3.1
schema; SDK signatures parsed from esm v3.4.0 source; tutorial hyperparameters read from the
cookbook notebooks. The **Atlas section is weaker**: no OpenAPI spec exists
(`/esm/protein/api/v1alpha1/openapi.json` → 404), so those params are transcribed from prose
docs, its response field lists are labelled "key selections" upstream (non-exhaustive), and it
is alpha with an empty changelog. **Rate limits are undocumented for both APIs.**

## Setup

```bash
pip install esm            # Python >=3.12, torch 2.11.x
export ESM_API_KEY=...     # from https://biohub.ai/developer-console/api-keys
```

```python
from esm.sdk import client, esmc_client, esmfold2_client, parallel_executor, ForgeBatchClient
esm3   = client(model="esm3-open-2024-03", url="https://biohub.ai")
esmc   = esmc_client(model="esmc-600m-2024-12", url="https://biohub.ai")
folder = esmfold2_client(model="esmfold2-2026-05", url="https://biohub.ai")
```
Defaults: `esm3-sm-open-v1`, `esmc-300m-2024-12`, `esmfold2-fast-2026-05`.
Classes still named "Forge" (API moved from forge.evolutionaryscale.ai → biohub.ai).

## Models

| ID | Family | Layers | d_model |
|---|---|---|---|
| esmc-300m-2024-12 | ESMC | 30 | 960 |
| esmc-600m-2024-12 | ESMC | 36 | 1152 |
| esmc-6b-2024-12 | ESMC | 80 | 2560 |
| esmfold2-fast-2026-05 | ESMFold2 | — | — |
| esmfold2-2026-05 | ESMFold2 | — | — |
| esmfold2-2026-05-cutoff-2025 | ESMFold2 | — | — |
| esm3-open-2024-03 | ESM3 | — | — |

NOTE: `esmfold2-2026-05-cutoff-2025` is NOT accepted by `/fold` or `/fold_all_atom`. It is
hardcoded by `FoldMaxAccuracyHandler._prepare_request` (esm/sdk/forge.py) as the model for the
batch `fold_max_accuracy` job, so you reach it via `ForgeBatchClient.fold_max_accuracy` and never
pass it yourself. Local equivalents: HF `ESMFold2-Experimental-*-Cutoff2025` (used by binder_design.py).

`max_ith_hidden_layer`: esmc-300m 30, esmc-600m 36, esmc-6b 80, esm3-small 48, esm3-medium 96.
`ith_hidden_layer=-1` (all layers) NOT supported on esmc-6b or any ESM3.

SAE ids (ESMC only, all Top-K k=64):
- esmc-300m-2024-12-sae-layer23-k64-codebook65536  ← requires normalize_features=False
- esmc-600m-2024-12-sae-layer27-k64-codebook{16384,65536}
- esmc-6b-2024-12-sae-layer60-k64-codebook{16384,65536}  ← layer60/16384 is the Atlas/paper SAE

## Endpoints — POST https://biohub.ai/api/v1/…

| Endpoint | Models | Purpose |
|---|---|---|
| /logits | ESM3, ESMC | logits, embeddings, hidden states, SAE features |
| /encode, /decode | ESM3, ESMC | Tracks ⇄ Tokens tokenizer round-trip |
| /generate | ESM3 | iterative unmasking of one track |
| /generate_tensor | ESM3 | same, raw token in/out |
| /forward_and_sample | ESM3 | one forward pass, per-track sampling |
| /fold | ESMFold2, ESM3 | sequence → backbone (`\|` separates chains) |
| /fold_all_atom | ESMFold2 | protein + DNA + RNA + ligand complexes |
| /inverse_fold | ESM3 | backbone → sequence |

All requests take `potential_sequence_of_concern: bool = false`.
All responses carry `model`, `created`, `warning_messages`.

### /generate defaults
`num_steps=20` (1–100, ≤ seq len), `temperature=0.5`, `temperature_annealing=true`,
`top_p=1.0`, `schedule=cosine`, `strategy=random`, `condition_on_coordinates_only=true`,
`invalid_ids=[]`, `only_compute_backbone_rmsd=false`.
track ∈ {sequence, structure, secondary_structure, sasa, function}.

### /fold + /fold_all_atom defaults
`num_loops=20` (0–20), `num_sampling_steps=100` (1–100), `lm_dropout=0.3`,
`lm_mask_pct=0.0` (unset ⇒ 0.1 for -fast, 0.0 for full), `msa_max_depth=1024`,
`msa_column_mask_rate=0.1`, `include_{distogram,pae,pair_chains_iptm,embeddings}=false`.
Response: coordinates (N×37×3), complex, plddt, ptm, interface_ptm, pae (L×L),
pair_chains_iptm (C×C), output_embedding_sequence, output_embedding_pair_pooled,
residue_index, entity_id.

### /inverse_fold
`inverse_folding_config = {temperature: 0.1, invalid_ids: []}`; optional `sequence` conditioning.
For diversity change the seed, not the temperature.

## LogitsConfig

sequence, structure*, secondary_structure*, sasa*, function*, residue_annotations,
return_embeddings, return_mean_embedding, return_hidden_states, return_mean_hidden_states,
ith_hidden_layer=-1, sae_config.
(* not served on the platform — local weights only.)

## SDK objects

- `ESMProtein`: sequence / secondary_structure / sasa / function_annotations / coordinates
  in; plddt, ptm, pae, crmsd, globularity, interface_ptm, pair_chains_iptm,
  output_embedding_*, residue_index, entity_id out.
- `ESMProteinTensor`: token form of the same tracks.
- `ESMProteinError(error_code, error_msg)` is **returned, not raised** — always isinstance-check.
- `esm.utils.structure.input_builder`: ProteinInput, RNAInput, DNAInput, LigandInput,
  Modification(position, ccd), CovalentBond(chain_id1,res_idx1,atom_idx1,…),
  PocketConditioning(binder_chain_id, contacts), DistogramConditioning,
  StructurePredictionInput(sequences, pocket, distogram_conditioning, covalent_bonds).
  `id=["A","B"]` on one sequence = homomer.

## Core patterns

```python
# ESMC embedding, all layers at once
CFG = LogitsConfig(sequence=True, return_mean_hidden_states=True)
out = esmc.logits(esmc.encode(ESMProtein(sequence=s)), CFG)
out.mean_hidden_state            # [n_layers+1, D]

# Leave-one-out masking → entropy / LLR / pseudo-perplexity
seqs = [s[:i] + "_" + s[i+1:] for i in range(len(s))]
with parallel_executor() as ex:
    outs = ex.execute_batch(user_func=get_logits, client=esmc, sequence=seqs)
# read position i+1 (BOS offset)

# SAE features
sae = "esmc-6b-2024-12-sae-layer60-k64-codebook16384"
out = esmc6b.logits(t, LogitsConfig(sae_config=SAEConfig(models=[sae])), return_bytes=False)
feats = out.sae_outputs[sae].to_dense().numpy()[1:-1]     # [L, 16384]

# ESMFold2 all-atom
spi = StructurePredictionInput(sequences=[ProteinInput(id="A", sequence=SEQ)])
res = folder.fold_all_atom(spi, config=FoldingConfig(num_loops=10, include_pae=True))
res.plddt, res.ptm, res.iptm, res.complex.to_mmcif()

# ESM3 motif scaffolding
prompt = ESMProtein(sequence="".join(seq_prompt), coordinates=structure_prompt)
out = esm3.generate(prompt, GenerationConfig(track="sequence", num_steps=n, temperature=0.5))
```

Layer choice heuristics from the tutorials: `[-2]` = safe default and best for property
regression; ~3/4 depth best for functional/EC/GO similarity; mid-network (~1/3–1/2)
best for structural-class separation.

## ESM Atlas API — SEPARATE API, no key

Base: `https://biohub.ai/esm/protein/api/v1alpha1/` — fully public, no auth, **alpha**
(schemas may change without notice). ~1.1B clustered proteins (70% id), ~6.5B full.
Sources: UniProt, JGI IMG, MGnify, SPIRE. Same 16,384 SAE features as
`esmc-6b-2024-12-sae-layer60-k64-codebook16384`.

Find Atlas clusters similar to a query protein by **ESMC SAE feature similarity**
(not BLAST / sequence alignment). Each result = one cluster (its representative).

### Similarity Search Endpoint

```
GET /esm/protein/api/v1alpha1/similarity-search
```

| Parameter | Type / Default | Meaning |
|---|---|---|
| `sequence` | *string (required)* | Query protein sequence (≤2048 aa) |
| `topk_results` | *int (10, max 100)* | Number of similar clusters to return |
| `include_cluster_info` | *bool (false)* | Set `true` to return `cluster_size` and `protein_name` (annotation) |
| `topk_features` | *int (20, max 100)* | Top shared SAE features to return |
| `min_similarity` | *float (optional)* | Minimum similarity score threshold |
| `cluster_pct_characterized_max` | *float (optional)* | Upper bound on % characterized. Set low (e.g. `0.1`) for **uncharacterized / divergent (dark proteome)** neighbors |

> [!WARNING]
> **Payload Size Warning**: The API response inlines a full **PDB structure per hit** in the `pdb` field. Strip or ignore this field client-side if not needed, as it quickly dominates memory and payload size.

### Quick Terminal Inspection (curl + Python)

```bash
curl -s "https://biohub.ai/esm/protein/api/v1alpha1/similarity-search\
?sequence=<SEQUENCE>\
&topk_results=10&include_cluster_info=true&topk_features=10" \
| python3 -c '
import json, sys
d = json.load(sys.stdin)
print("cluster_id  members  simil.  annotation")
for p in d["similar_proteins"]:
    cid = p["protein_hash"][:8]                 # cluster id is first 8 hex chars
    sim = round(p["similarity_score"], 3)
    print(cid, p["cluster_size"], sim, p["protein_name"])
'
```

### Response Schema

* **Top level**: `query_sequence`, `protein_hash`, `similar_proteins[]`, `top_features_across_results`, `restricted_count`.
  * `protein_hash`: Present at top-level *only* if the exact query sequence already exists in the Atlas.
  * `restricted_count`: Hits withheld by biosecurity filters (non-zero means incomplete results).
* **Each `similar_proteins[]` entry**:
  ```
  protein_hash        # cluster id = first 8 hex chars
  protein_accession
  sequence_length
  similarity_score
  cluster_size        # member count (when include_cluster_info=true)
  protein_name        # annotation text (when include_cluster_info=true)
  ptm, mean_plddt, residues_plddt[], pdb
  ```
* **Each `top_features_across_results[]` entry**:
  `feature_index`, `occurrence_count`, `min_activation`, `max_activation`, `mean_activation`.

### Companion Endpoints (Drill-Down)

| Endpoint | Returns / Purpose |
|---|---|
| `GET /esm/protein/api/v1alpha1/proteins/{protein_hash}?topk_features=10&normalize_features=true` | One protein's top SAE features (`fold_on_miss=true` folds on-the-fly if <700 aa) |
| `GET /esm/protein/api/v1alpha1/clusters/{protein_hash}?topk_features=10` | Cluster members, representative hash, and cluster's top SAE features |
| `GET /esm/protein/api/v1alpha1/features/{idx}` | Full SAE feature dossier (label, summary, category, threshold, nearest neighbors) |
| `GET /esm/protein/api/v1alpha1/features` | All 16,384 SAE features in one bulk call |
| `POST /esm/protein/api/v1alpha1/proteins/batch` | Batch lookup (≤500 hashes): structure, cluster info, features → 200 zip or 202 job |
| `GET /esm/protein/api/v1alpha1/proteins/batch/jobs/{id}` | Poll batch job (200 done · 202 running · 410 expired) |
| `DELETE /esm/protein/api/v1alpha1/proteins/batch/jobs/{id}` | Cancel batch job (204) |
| `GET /esm/protein/api/v1alpha1/proteins/{hash}/thumbnail/{type}` | Structure PNG thumbnail (`plddt` or `pct-characterized`) |
| `GET /esm/protein/api/v1alpha1/uniprot/{uniprot_id}` | Map UniProt Accession → Atlas protein hash |

**Full Discovery Flow**: Query sequence → Top clusters (`/similarity-search`) → Per-cluster members (`/clusters/{hash}`) → SAE feature dossier (`/features/{idx}`).

### Python API Query Pattern

```python
import hashlib, httpx
h = hashlib.md5(b"FVNQHLCGSHLVEALYLVCGERGFFYTPKT").hexdigest()   # protein hash = MD5 of sequence
p = httpx.get(f"https://biohub.ai/esm/protein/api/v1alpha1/proteins/{h}",
              params={"topk_features": 5}).json()
c = httpx.get(f"https://biohub.ai/esm/protein/api/v1alpha1/clusters/{p['cluster_rep_protein_hash']}").json()

r = httpx.get("https://biohub.ai/esm/protein/api/v1alpha1/similarity-search",
              params={"sequence": SEQ, "topk_results": 5, "include_cluster_info": True}).json()
```

### Atlas Characteristics & Rules of Thumb
- **On-the-fly folding**: Miss on `/proteins/{hash}` with `fold_on_miss=true` folds with ESMFold2 on the fly (for <700 aa) and computes SAE features in real time. Zero-setup fold+annotate.
- **Feature space clustering**: Clustering runs in SAE feature space (Linclust-inspired, MinHash+LSH, exact Jaccard verify); members share ≥60% feature overlap (Jaccard ≥ 0.6). Similarity = cosine over feature activations.
- **Characterization status**:
  - *Characterized*: Own Pfam domain.
  - *Partially characterized*: Cluster-mate has a Pfam domain.
  - *Uncharacterized*: No known domain in cluster.
  - Low `cluster_pct_characterized_max` filters toward dark proteome.
- **Feature activation**: <1% of features fire per protein. Normalized score = activation ÷ max over 208M UniRef90, weighted by feature rarity/selectivity — raw activations are not comparable.
- **Published pLDDT bands**: 90–100 very high · 70–90 confident · 50–70 low/possible disorder · <50 likely intrinsically disordered. pTM > ~0.5 = confident.

### Bulk Offline Data (AWS S3 Open Data)
For offline large-scale cluster processing without API rate limits, access the public open dataset bucket:
- **S3 Bucket**: `s3://esm-protein-atlas/v1/` (`--no-sign-request`)
- **Key paths**:
  - `clusters/data/representative_proteins.parquet`
  - `clusters/indexes/secondary/`
  - `sae/data_shards/`
- **Registry info**: https://registry.opendata.aws/biohub-esm-atlas/

## Interaction-propensity signals (relevant to this project)

| Signal | How |
|---|---|
| ipTM scalar | `fold_all_atom` → `result.iptm` |
| per-chain-pair ipTM | `include_pair_chains_iptm=True` → C×C |
| residue-level interface confidence | `include_pae=True` AND `include_embeddings=True`, off-diagonal block by `entity_id` (entity_id/residue_index are only returned when embeddings are requested) |
| pair embeddings for a learned head | `include_embeddings=True` → `output_embedding_pair_pooled` |
| epitope-directed scoring | `PocketConditioning(binder_chain_id, contacts)` |
| cheap sequence-only screen | ESMC `mean_hidden_state` at a swept layer |
| load-bearing residues | leave-one-out entropy + LLR |
| functional neighbours / dark-proteome check | Atlas `/similarity-search`, keyless |

Binder design protocol (paper, nanomolar hits): `cookbook/tutorials/binder_design.py`,
Modal H100, inverts ESMFold2. `#` = design, letter = fix. AdamW lr 0.1, 150 steps,
losses intra_contact 0.5 / inter_contact 0.5 / glob 0.2 / epitope 0.5.
Selection = 0.5·ipTM + 0.5·distogram-ipTM-proxy. ≈$2/job, ≈$150 for a 256-job sweep.

## Gotchas

- BOS/EOS: slice `[1:-1]` or index `+1` for per-residue outputs.
- `GenerationConfig` SDK default temperature is 1.0; HTTP default is 0.5. Set it explicitly.
- `plddt` has a **null slot per chain separator** — drop nulls before `mean()`.
- `esmfold2-fast` silently ignores MSAs (no alignment reader). Use `esmfold2-2026-05`.
- `all_atom_input` present ⇒ top-level `sequence`/`msa` ignored.
- SASA "inf" encoding changed 2025-03-18 from `-1` to `1000`.
- Function annotations = `(interpro_tag, start, end)`, 1-indexed inclusive, InterPro 95.0.
- SAE naming: schema uses `…-sae-layer27-k64-codebook16384`; one cookbook snippet still
  shows the stale underscore form. Trust the schema.
- `SAEConfig(model=...)` singular is DEPRECATED (warns); passing both model and models raises
  ValueError. 300M SAEs raise ValueError unless `normalize_features=False` (default is True).
- `batch_executor()` deprecated → `parallel_executor()`.
- Rate limits & private model entitlements live behind auth in the Developer Console.
- Atlas is a DIFFERENT API: different base path, no key, no SDK, alpha stability, its own
  undocumented rate limits. Caps: 2048 aa search, 100 results, 500 hashes/batch, <700 aa fold.
- Platform applies keyword/sequence guardrails for controlled pathogens and toxins;
  local HF weights (MIT) do not.

## Local execution

```python
from esm.models.esmc import EsmcForMaskedLM, EsmcTokenizer, EsmcSaeModel, EsmcForSequenceClassification
from esm.models.esmfold2 import EsmFold2Model, ESMFold2InputBuilder
EsmcForMaskedLM.from_pretrained("biohub/ESMC-6B", device="cuda")
EsmFold2Model.from_pretrained("biohub/ESMFold2-Fast", device="cuda").set_kernel_backend("fused")
```
PEFT fine-tune: `LoraConfig(r=8, lora_alpha=16, target_modules=["out_proj"],
target_parameters=["layernorm_qkv.weight","ffn.fc1_weight","ffn.fc2_weight"],
modules_to_save=["classifier"])` — fused params aren't nn.Linear, hence `target_parameters`.

Apple Silicon: MLX port via `mlx_lm.models.esmfold2.ESMFold2Model`; 32 GB min, ~25 GB weights.

## Sources

- https://biohub.ai/openapi.json (ground truth for endpoints/fields)
- https://github.com/Biohub/esm (SDK 3.4.0 + cookbook/tutorials)
- https://biohub.ai/api-reference, /models/esmc, /models/esmfold2, /models/esm3
- https://biohub.ai/esm/protein/atlas/api-docs/ (api_reference, concepts, examples, faqs)
- https://registry.opendata.aws/biohub-esm-atlas/ (Atlas S3 registry)
- Candido et al. 2026, bioRxiv 10.64898/2026.06.03.729735
- Hayes et al. 2025, Science 10.1126/science.ads0018
