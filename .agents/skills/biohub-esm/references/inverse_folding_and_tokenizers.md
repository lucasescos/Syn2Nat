# Inverse Folding, Tokenizers & Track Transformations

This guide details methods for converting 3D coordinates into sequence candidates (inverse folding), transforming biological tracks into discrete tokens, and operating directly on token tensors.

---

## 1. Inverse Folding: `POST /api/v1/inverse_fold`

Recovers amino acid sequences compatible with specified 3D backbone coordinates.

### Request Payload Schema
```json
{
  "model": "esm3-open-2024-03",
  "coordinates": [
    [[-7.499, 2.045, 10.377], [-6.562, 2.812, 9.562], [-5.800, 1.901, 8.606]],
    [[-5.594, 0.853, 8.760], [-4.375, 0.233, 8.250], [-4.120, 1.250, 7.150]]
  ],
  "inverse_folding_config": {
    "temperature": 0.1,
    "invalid_ids": []
  }
}
```

### Parameters
* `coordinates` *(list of floats, required)*: Backbone 3D coordinates array ($N \times 3$ or $N \times 37 \times 3$).
* `temperature` *(float, default: 0.1)*: Sampling temperature. Lower values ($0.1$) produce conservative, consensus-like sequences; higher values ($0.5$) introduce diverse candidate designs.
* `invalid_ids` *(list of ints, optional)*: Specific token IDs to mask out/disallow from emission.

### CLI Script Usage
```bash
# Inverse fold from PDB or mmCIF coordinate file
python .agents/skills/biohub-esm/scripts/inverse_fold.py structure.pdb -o results/designed_seqs.fasta --temperature 0.1 --num-samples 3
```

---

## 2. Track Tokenization: `POST /api/v1/encode` & `POST /api/v1/decode`

ESM3 represents all biological modalities as discrete tokens. These endpoints convert between continuous / human-readable representations (`Tracks`) and discrete integer indices (`Tokens`).

### `POST /api/v1/encode`
Encodes continuous/string tracks (sequence, secondary structure, coordinates) into discrete token arrays:
```json
{
  "model": "esm3-open-2024-03",
  "tracks": {
    "sequence": "MSHHWGYGKH..."
  }
}
```

### `POST /api/v1/decode`
Decodes discrete token integer sequences back into biological tracks:
```json
{
  "model": "esm3-open-2024-03",
  "tokens": {
    "sequence_tokens": [1, 23, 45, 12, 2]
  }
}
```

### CLI Script Usage
```bash
# Encode sequence/tracks into discrete tokens
python .agents/skills/biohub-esm/scripts/tokenize_tracks.py encode "MSHHWGYGKH..." -o tokens.json

# Decode discrete tokens back into tracks
python .agents/skills/biohub-esm/scripts/tokenize_tracks.py decode tokens.json -o decoded.json
```

---

## 3. Low-Level Generation & Sampling Endpoints

* **`POST /api/v1/generate_tensor`**: Performs iterative unmasking directly on discrete token tensors for custom sampling schedules, masked diffusion, or guided generation.
* **`POST /api/v1/forward_and_sample`**: Executes a single model forward pass with per-track multinomial sampling without unmasking loops.
