# ESM Python SDK Reference (`esm` package)

The official `esm` Python SDK enables direct programmatic interaction with the Biohub ESM platform (`https://biohub.ai`) using native typed objects.

---

## 1. Installation

```bash
pip install esm
```

Requirements: Python $\ge 3.10$, `torch`, `biotite`.

---

## 2. ESM3 Client Initialization

```python
from esm.sdk.forge import ESM3ForgeInferenceClient
from esm.sdk.api import ESMProtein, GenerationConfig

client = ESM3ForgeInferenceClient(
    model="esm3-open-2024-03",   # or "esm3-medium-2024-08", "esm3-large-2024-03"
    url="https://biohub.ai",
    token="<ESM_API_KEY>"
)
```

---

## 3. The `ESMProtein` Object

The `ESMProtein` class is the universal data structure representing proteins across all multimodal tracks.

### Attributes
* `sequence` *(str | None)*: Amino acid sequence string. Use `_` for masked positions.
* `coordinates` *(torch.Tensor | None)*: Atomic coordinates of shape `[Length, 37, 3]` covering all heavy atoms per residue in standard PDB order.
* `secondary_structure` *(str | None)*: Secondary structure string.
* `sasa` *(torch.Tensor | None)*: Solvent accessible surface area values.
* `function_annotations` *(list | None)*: Functional keywords / ontology annotations.
* `plddt` *(torch.Tensor | None)*: Per-residue pLDDT confidence scores.
* `ptm` *(torch.Tensor | None)*: Global predicted TM-score.
* `potential_sequence_of_concern` *(bool)*: Automated biosecurity check.

### Built-in Methods
* `to_pdb(filepath: str)`: Exports the 3D structure directly to a standard `.pdb` file.
* `to_pdb_string() -> str`: Returns the protein structure as a PDB-formatted string.
* `from_pdb(filepath: str) -> ESMProtein`: Loads an `ESMProtein` from an existing `.pdb` file.
* `to_protein_chain()`: Converts to an individual chain representation.
* `to_protein_complex()`: Converts to a multi-entity complex representation.
* `copy() -> ESMProtein`: Deep copy of the protein object.

---

## 4. `GenerationConfig` Reference

```python
config = GenerationConfig(
    track="sequence",             # Target track: "sequence", "structure", "secondary_structure"
    num_steps=8,                  # Number of unmasking steps
    temperature=0.7,              # Sampling temperature (0.0 for structure, 0.2-0.7 for sequence)
    schedule="cosine"             # "cosine" or "linear"
)
```

---

## 5. End-to-End Examples

### Sequence Design / Infilling
```python
from esm.sdk.forge import ESM3ForgeInferenceClient
from esm.sdk.api import ESMProtein, GenerationConfig

client = ESM3ForgeInferenceClient(model="esm3-open-2024-03", url="https://biohub.ai", token=api_key)

# Prompt with masked residues
prompt = ESMProtein(sequence="M___G__L___K")
config = GenerationConfig(track="sequence", num_steps=8, temperature=0.7)

designed_protein = client.generate(prompt, config)
print("Designed sequence:", designed_protein.sequence)
```

### Structure Folding & PDB Export
```python
# Predict 3D coordinates from sequence
protein = ESMProtein(sequence="MSHHWGYGKH")
config = GenerationConfig(track="structure", num_steps=1, temperature=0.0)

folded = client.generate(protein, config)
print("Coordinates tensor shape:", folded.coordinates.shape)  # torch.Size([10, 37, 3])
print("Mean pLDDT:", folded.plddt.mean().item())

# Save directly to PDB
folded.to_pdb("predicted_structure.pdb")
```

---

## 6. Other SDK Client Constructors

* `esm.sdk.esmc_client(model, url, token)`: Dedicated client for ESMC representations & embeddings.
* `esm.sdk.esmfold2_client(model, url, token)`: Dedicated client for ESMFold2 fast folding.
* `esm.sdk.ForgeBatchClient`: High-throughput asynchronous batch submission client.
