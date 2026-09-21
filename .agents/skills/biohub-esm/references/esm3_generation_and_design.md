# ESM3 Generative Protein Design & Motif Scaffolding

ESM3 is EvolutionaryScale's multimodal generative biology model. It reasons simultaneously across biological modalities—sequence, 3D structure, secondary structure, SASA, and function tokens—using discrete representation tracks.

---

## 1. Supported ESM3 Models

| Model Name | Parameters | Capabilities | Best For |
| :--- | :--- | :--- | :--- |
| `esm3-open-2024-03` | 1.4B | Sequence infilling, structure prediction, secondary structure | Rapid prototyping, local/open weights compatible |
| `esm3-medium-2024-08` | 7B | Multi-track generation and refinement | Balanced quality and inference latency |
| `esm3-large-2024-03` | 98B | Frontier reasoning, complex scaffolding, de novo generation | High-fidelity protein design & complex functional motifs |
| `esm3-sm-open-v1` | 1.4B | Compact open checkpoint | Quick generation tests |

---

## 2. Generation Modalities (Tracks)

ESM3 models operate on parallel discrete token tracks:
* **`sequence`**: Amino acid residues (20 standard amino acids + masking tokens `_`).
* **`structure`**: 3D atomic backbone and heavy atom coordinates.
* **`secondary_structure`**: 8-state or 3-state secondary structure strings (`H` = helix, `E` = sheet, `C` = coil/loop, `T` = turn, `S` = bend, `G` = 3_10 helix).
* **`sasa`**: Solvent Accessible Surface Area tokens (burial / exposure profiling).
* **`function`**: InterPro / Gene Ontology function keywords and annotation tokens.

---

## 3. Platform Endpoint: `POST https://biohub.ai/api/v1/generate`

Executes iterative unmasking and multi-track sampling.

### Headers
```http
POST /api/v1/generate HTTP/1.1
Host: biohub.ai
Authorization: Bearer <ESM_API_KEY>
Content-Type: application/json
```

### Request Payload Schema
```json
{
  "model": "esm3-open-2024-03",
  "track": "sequence",
  "inputs": {
    "sequence": "M___G__L___K",
    "secondary_structure": null,
    "coordinates": null
  },
  "num_steps": 20,
  "temperature": 0.5,
  "temperature_annealing": true,
  "top_p": 1.0,
  "schedule": "cosine",
  "strategy": "random",
  "condition_on_coordinates_only": false,
  "potential_sequence_of_concern": false
}
```

### Parameter Reference
* `track` *(string, required)*: The target track to generate or infill (`"sequence"`, `"structure"`, `"secondary_structure"`, `"sasa"`, `"function"`).
* `inputs` *(object, required)*: Partial constraints/prompt. Masked sequence positions should be represented by `_`.
* `num_steps` *(integer, default: 20)*: Number of iterative unmasking steps. Typically set to $\min(\text{steps}, \text{masked\_residues})$.
* `temperature` *(float, default: 0.5)*: Sampling temperature.
  * For **structure**: Use `0.0` for deterministic, reproducible coordinates.
  * For **sequence**: Use `0.2`–`0.7` depending on desired sequence diversity.
* `temperature_annealing` *(bool, default: true)*: Linearly cools temperature as steps proceed.
* `top_p` *(float, default: 1.0)*: Nucleus sampling probability cutoff.
* `schedule` *(string, choices: `["cosine", "linear"]`)*: Token unmasking schedule over steps.
* `strategy` *(string, choices: `["random", "entropy"]`)*:
  * `"random"`: Unmasks randomly sampled masked positions.
  * `"entropy"`: Unmasks the lowest-entropy (highest confidence) positions first.

---

## 4. Response Payload Schema

```json
{
  "model": "esm3-open-2024-03",
  "created": "2026-09-16T13:13:09.877056",
  "potential_sequence_of_concern": false,
  "tokens_used": 10,
  "warning_messages": null,
  "outputs": {
    "sequence": "MLILGLLLLLLK",
    "secondary_structure": "CCCCCCCCCCCSGGGCCTTCCCCTT...",
    "sasa": null,
    "function": null,
    "coordinates": [[-7.499, 2.045, 10.377], ...],
    "plddt": [0.535, 0.549, 0.570, ...],
    "ptm": 0.012,
    "crmsd": null,
    "globularity": null,
    "interface": null,
    "interface_ptm": null,
    "pae": null
  }
}
```

### Output Fields
* `potential_sequence_of_concern` *(bool)*: Automated biosecurity screening flag against regulated pathogen sequences.
* `tokens_used` *(int)*: Number of tokens consumed in generation.
* `outputs.sequence` *(string)*: Completed designed amino acid sequence.
* `outputs.secondary_structure` *(string)*: Secondary structure state for each residue.
* `outputs.coordinates` *(list)*: 3D atomic coordinates.
* `outputs.plddt` *(list of floats)*: Per-residue confidence metric (0.0 to 1.0).
* `outputs.ptm` *(float)*: Predicted TM-score (0.0 to 1.0).

---

## 5. CLI Script Usage

Run the pre-configured CLI script:
```bash
# Design unmasked positions ('_' indicates residues to generate)
python .agents/skills/biohub-esm/scripts/design_sequence.py "M___G__L___K" -o results/design.fasta --steps 20 --temperature 0.5
```
