# Representations, Hidden States & Sparse Autoencoders (SAEs)

ESM Cambrian (ESMC) models are frontier protein language world models trained on billions of biological sequences. They provide dense vector representations, per-residue embeddings, and interpretable biological feature dictionaries via Sparse Autoencoders.

---

## 1. Supported ESMC Models

| Model ID | Layers | Embedding Dim | Context Length | Primary Use Case |
| :--- | :--- | :--- | :--- | :--- |
| `esmc-300m-2024-12` | 30 | 960 | 2048 | Ultra-fast sequence embedding and scoring |
| `esmc-600m-2024-12` | 36 | 1152 | 2048 | Standard high-accuracy sequence representation |
| `esmc-6b-2024-12` | 80 | 2560 | 2048 | Frontier world model & base for SAE codebooks |

---

## 2. Platform Endpoint: `POST /api/v1/logits`

Extracts logits, hidden states, pooled representations, or SAE activations.

### Request Payload Schema
```json
{
  "model": "esmc-600m-2024-12",
  "protein": {
    "sequence": "MSHHWGYGKH..."
  },
  "config": {
    "sequence": true,
    "return_embeddings": true,
    "return_mean_embedding": true,
    "return_hidden_states": false,
    "return_mean_hidden_states": false,
    "ith_hidden_layer": -2,
    "sae_config": {
      "models": ["esmc-6b-2024-12-sae-layer60-k64-codebook16384"]
    }
  }
}
```

### Parameter Guidelines
* `ith_hidden_layer` *(int, default: -2)*: Penultimate layer (-2) is optimal for downstream classification, clustering, and property prediction. Use `-1` for the final layer.
* `return_mean_embedding` *(bool)*: Returns pooled sequence-level vector ($1 \times D$) with BOS/EOS tokens removed.
* `return_embeddings` *(bool)*: Returns full per-residue matrix ($L \times D$).

---

## 3. Sparse Autoencoders (SAEs) for Interpretability

SAEs decompose dense transformer activations into ~16,384 sparse, interpretable biological features (e.g. catalytic dyads, metal binding sites, transmembrane helices, post-translational modification motifs).

### SAE Model Identifier
* `"esmc-6b-2024-12-sae-layer60-k64-codebook16384"`: Extracts Top-$k$ (64 active features) across 16k codebook units at layer 60.

### SAE Response Schema
```json
{
  "sae_activations": {
    "indices": [100, 452, 1205],
    "values": [3.42, 1.88, 0.95]
  }
}
```

Each feature index corresponds to an annotated biological dossier accessible via the ESM Atlas (`GET /features/{index}`).

---

## 4. In-Silico Deep Mutational Scanning (DMS)

By computing position-by-position Shannon entropy from output logits:
$$H_i = -\sum_{a \in \mathcal{A}} p_{i,a} \log_2 p_{i,a}$$
Low entropy positions indicate mutationally constrained, load-bearing residues essential for folding or catalytic activity.

---

## 5. CLI Script Usage

```bash
# 1. Extract pooled sequence embedding (penultimate layer)
python .agents/skills/biohub-esm/scripts/embed.py "MSHHWGYGKH..." -o results/embed.json --model esmc-600m-2024-12 --layer -2

# 2. Extract per-residue embeddings
python .agents/skills/biohub-esm/scripts/embed.py sequence.fasta -o results/residue_embed.json --per-residue

# 3. Extract active SAE features and fetch Atlas biological annotations
python .agents/skills/biohub-esm/scripts/sae_features.py "MSHHWGYGKH..." --topk 10 --annotate

# 4. Perform in-silico DMS scan to identify constrained sites
python .agents/skills/biohub-esm/scripts/deep_mutational_scan.py sequence.fasta -o results/dms.tsv --top-constrained 15
```
