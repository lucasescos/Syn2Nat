# Structure Prediction & All-Atom Complex Modeling

Biohub ESM provides state-of-the-art structure prediction for single proteins, multimers, and all-atom complexes through ESMFold2 (`POST /api/v1/fold` and `POST /api/v1/fold_all_atom`).

---

## 1. Supported Folding Models

| Model ID | Purpose | Supported Inputs |
| :--- | :--- | :--- |
| `esmfold2-fast-2026-05` | Default, high-speed single-chain & multimer folding | Protein sequences |
| `esmfold2-2026-05` | Full structure prediction with high accuracy | Protein, DNA, RNA, Small-molecule ligands |
| `esmfold2-2026-05-cutoff-2025` | Temporal benchmark / batch folding | Single & multi-chain proteins |

---

## 2. Standard Structure Prediction: `POST /api/v1/fold`

Predicts 3D atomic coordinates for single chains or protein-protein complexes.

### Request Payload Schema
```json
{
  "model": "esmfold2-fast-2026-05",
  "sequences": [
    {"id": "A", "sequence": "MSHHWGYGKH..."},
    {"id": "B", "sequence": "GATAGCGCTATC..."}
  ],
  "config": {
    "num_loops": 20,
    "num_sampling_steps": 100,
    "lm_dropout": 0.3,
    "lm_mask_pct": 0.0,
    "include_pae": true,
    "include_distogram": false,
    "include_pair_chains_iptm": true,
    "include_embeddings": false
  },
  "potential_sequence_of_concern": false
}
```

*Multi-chain syntax shorthand:* Sequences can also be provided as a pipe-delimited string `CHAIN_A_SEQ|CHAIN_B_SEQ`.

---

## 3. All-Atom Complex Prediction: `POST /api/v1/fold_all_atom`

Predicts structures of proteins co-folded with Nucleic Acids, Small-Molecule Ligands, and Chemical Modifications.

### Entity Types
1. **Protein**: `{"id": "A", "sequence": "...", "type": "protein"}`
2. **DNA**: `{"id": "B", "sequence": "GATCGATC", "type": "dna"}`
3. **RNA**: `{"id": "C", "sequence": "GAUCGAUC", "type": "rna"}`
4. **Ligand (CCD Code)**: `{"id": "L", "ccd": ["SAH"], "type": "ligand"}`
5. **Ligand (SMILES)**: `{"id": "L", "smiles": "CC(=O)NC1...", "type": "ligand"}`

### Epitope / Pocket Conditioning
Direct binder generation or complex conformation by specifying pocket residues:
```json
"pocket": {
  "binder_chain_id": "B",
  "contacts": [
    {"target_res_idx": 45, "binder_res_idx": 12, "distance_threshold": 4.5}
  ]
}
```

### Covalent Bonds
Specify inter-chain or residue-ligand covalent links:
```json
"covalent_bonds": [
  {
    "chain_id1": "A", "res_idx1": 25, "atom_idx1": "SG",
    "chain_id2": "L", "res_idx2": 1, "atom_idx2": "C1"
  }
]
```

---

## 4. Metrics & Confidence Outputs

* **`plddt`**: Per-residue confidence (0 to 100 or 0.0 to 1.0). $>90$ indicates very high confidence; $<50$ indicates disordered regions.
* **`ptm`**: Predicted TM-score (0.0 to 1.0) assessing overall fold topology.
* **`iptm`**: Interface predicted TM-score assessing inter-chain complex orientation.
* **`pae`**: Predicted Aligned Error matrix ($L \times L$) measuring relative position uncertainty between residue pairs.

---

## 5. CLI Script Usage (`scripts/fold.py`)

```bash
# Single protein sequence
python .agents/skills/biohub-esm/scripts/fold.py "MSHHWGYGKH..." -o results/ca2.cif

# Multi-chain protein complex (pipe-separated)
python .agents/skills/biohub-esm/scripts/fold.py "CHAIN_A_SEQ|CHAIN_B_SEQ" -o results/complex.cif --include-pae

# All-Atom Complex: Protein + Ligand + DNA
python .agents/skills/biohub-esm/scripts/fold.py "MSHHWGYGKH..." --ligand SAH --dna "GATAGCGCTATC" -o results/all_atom.cif
```
