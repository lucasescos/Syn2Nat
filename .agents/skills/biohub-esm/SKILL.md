---
name: biohub-esm
description: >-
  Predict protein and complex 3D structures with ESMFold2, extract sequence embeddings
  and representations with ESMC (300M, 600M, 6B), generate and design sequences with ESM3,
  perform inverse folding and mutational scanning, analyze functional motifs with Sparse Autoencoders (SAEs),
  and query the public ESM Atlas (>1.1B structures). Requires ESM_API_KEY for platform endpoints.
---

# Biohub ESM: Complete Suite for Structure Prediction, Design, and World Models

This skill equips the agent with a complete toolset to interact with Biohub's Evolutionary Scale Modeling (ESM) platform (`https://biohub.ai`), covering all authenticated platform endpoints (`/fold`, `/fold_all_atom`, `/logits`, `/generate`, `/inverse_fold`, `/encode`, `/decode`) and public Atlas endpoints.

---

## Prerequisites & Credentials Verification

Platform endpoints require an `ESM_API_KEY`.

### Verify Credentials Before Running Workflows
Before executing any script that queries authenticated platform endpoints, **ALWAYS** verify that `ESM_API_KEY` is present:

```bash
python .agents/skills/biohub-esm/scripts/verify_credentials.py
```

* If verification **succeeds** (exit code 0), proceed with the requested workflow.
* If verification **fails** (exit code 1), prompt the user to store their key into `~/.env` using the masked terminal command:
  * **Windows (PowerShell)**:
    ```powershell
    $key = Read-Host -AsSecureString 'Enter ESM_API_KEY'
    $BSTR = [System.Runtime.InteropServices.Marshal]::SecureStringToBSTR($key)
    $plain = [System.Runtime.InteropServices.Marshal]::PtrToStringAuto($BSTR)
    Add-Content -Path ~/.env -Value "ESM_API_KEY=$plain"
    ```
  * **Linux / macOS**:
    ```bash
    printf "Enter ESM_API_KEY (hidden): " && read -s k && echo && echo "ESM_API_KEY=$k" >> ~/.env
    ```
  * Key registration: [Biohub Developer Console](https://biohub.ai/developer-console/api-keys).

---

## Tool Catalog & CLI Scripts

All scripts are located in `.agents/skills/biohub-esm/scripts/` and can be run with `python` or `uv run`. You can also use the official `esm` Python SDK (`from esm.sdk.forge import ESM3ForgeInferenceClient`) directly.

### Quick Routing: Where to Look by Task

| Task / Biological Question | Modular Reference Guide | Primary Endpoint / SDK Class | Pre-built CLI Script |
| :--- | :--- | :--- | :--- |
| **Generative Protein Design & Infilling** | [esm3_generation_and_design.md](references/esm3_generation_and_design.md) | `POST /api/v1/generate` <br> `ESM3ForgeInferenceClient` | `scripts/design_sequence.py` |
| **Python SDK (`esm` package) & PDB Export** | [python_sdk_reference.md](references/python_sdk_reference.md) | `esm.sdk.forge` <br> `esm.sdk.api.ESMProtein` | `test_esm3_sdk.py` |
| **Structure Folding & All-Atom Complexes** | [structure_prediction_and_complexes.md](references/structure_prediction_and_complexes.md) | `POST /api/v1/fold` <br> `POST /api/v1/fold_all_atom` | `scripts/fold.py` |
| **Embeddings, Logits & Interpretability (SAEs)** | [representations_and_interpretability.md](references/representations_and_interpretability.md) | `POST /api/v1/logits` <br> ESMC models | `scripts/embed.py` <br> `scripts/sae_features.py` |
| **Deep Mutational Scanning (DMS Entropy)** | [representations_and_interpretability.md](references/representations_and_interpretability.md) | `POST /api/v1/logits` | `scripts/deep_mutational_scan.py` |
| **Inverse Folding & Tokenizer Tracks** | [inverse_folding_and_tokenizers.md](references/inverse_folding_and_tokenizers.md) | `POST /api/v1/inverse_fold` <br> `POST /api/v1/encode`, `/decode` | `scripts/inverse_fold.py` <br> `scripts/tokenize_tracks.py` |
| **Public ESM Atlas (>1.1B Structures)** | [esm_atlas_public_api.md](references/esm_atlas_public_api.md) | Keyless public Atlas API | `scripts/atlas_search.py` |

---

### 1. Structure Prediction & Complex Modeling (`fold.py`)
Predicts atomic coordinates, pLDDT confidence, pTM, and ipTM for single chains, multimers, and all-atom complexes (`POST /api/v1/fold` and `POST /api/v1/fold_all_atom`):

```bash
# Fold single protein sequence
python .agents/skills/biohub-esm/scripts/fold.py "MSHHWGYGKH..." -o results/ca2.cif

# Fold multi-chain protein complex (pipe-separated)
python .agents/skills/biohub-esm/scripts/fold.py "CHAIN_A_SEQ|CHAIN_B_SEQ" -o results/complex.cif --include-pae

# All-Atom Complex: Protein + Ligand (CCD code or SMILES) + DNA
python .agents/skills/biohub-esm/scripts/fold.py "MSHHWGYGKH..." --ligand SAH --dna "GATAGCGCTATC" -o results/all_atom.cif

# Epitope/Pocket-conditioned folding
python .agents/skills/biohub-esm/scripts/fold.py "TARGET_SEQ|BINDER_SEQ" --pocket-chain B --pocket-contacts pocket.json
```

### 2. Generative Protein Design & Motif Scaffolding (`design_sequence.py`)
Iteratively unmasks sequences and scaffolds functional motifs using ESM3 (`POST /api/v1/generate`):

```bash
# Design unmasked positions ('_' indicates residues to generate)
python .agents/skills/biohub-esm/scripts/design_sequence.py "M___G__L___K" -o results/design.fasta --steps 20 --temperature 0.5
```

### 3. Inverse Folding (`inverse_fold.py`)
Recovers candidate amino acid sequences matching 3D backbone coordinates (`POST /api/v1/inverse_fold`):

```bash
# Inverse fold from PDB or mmCIF coordinate file
python .agents/skills/biohub-esm/scripts/inverse_fold.py structure.pdb -o results/designed_seqs.fasta --temperature 0.1 --num-samples 3
```

### 4. Deep Mutational Scanning & Entropy (`deep_mutational_scan.py`)
Computes position-by-position Shannon entropy and identifies mutationally constrained / load-bearing positions via ESMC logits:

```bash
# In-silico DMS scan to identify critical residues
python .agents/skills/biohub-esm/scripts/deep_mutational_scan.py sequence.fasta -o results/dms.tsv --top-constrained 15
```

### 5. Sequence Embeddings & Representations (`embed.py`)
Extracts pooled sequence representations or residue-level embeddings across transformer layers using ESMC (`POST /api/v1/logits`):

```bash
# Pooled mean embedding from penultimate layer (-2)
python .agents/skills/biohub-esm/scripts/embed.py "MSHHWGYGKH..." -o results/ca2_embed.json --model esmc-600m-2024-12 --layer -2

# Per-residue embeddings (dim L x D, BOS/EOS trimmed)
python .agents/skills/biohub-esm/scripts/embed.py sequence.fasta -o results/residue_embed.json --per-residue
```

### 6. Sparse Autoencoders (SAEs) for Interpretability (`sae_features.py`)
Extracts interpretable biological functional features (~16k codebook units) decomposed by ESMC SAEs:

```bash
# Extract active SAE features and fetch natural language annotations from Atlas
python .agents/skills/biohub-esm/scripts/sae_features.py "MSHHWGYGKH..." --topk 10 --annotate
```

### 7. Tokenizer Transformations (`tokenize_tracks.py`)
Converts between biological tracks and discrete token representations (`POST /api/v1/encode` and `POST /api/v1/decode`):

```bash
# Encode sequence/tracks into discrete tokens
python .agents/skills/biohub-esm/scripts/tokenize_tracks.py encode "MSHHWGYGKH..." -o tokens.json

# Decode discrete tokens back into tracks
python .agents/skills/biohub-esm/scripts/tokenize_tracks.py decode tokens.json -o decoded.json
```

### 8. ESM Atlas Public Suite (`atlas_search.py`)
Queries the public keyless Atlas endpoint (`https://biohub.ai/esm/protein/api/v1alpha1/`) indexing >1.1B structures:

```bash
# Homology & similarity search (returns top matches and shared SAE features)
python .agents/skills/biohub-esm/scripts/atlas_search.py similarity "MSHHWGYGKH..." --topk 5

# Look up UniProt accession -> protein hash, sequence, gene, function
python .agents/skills/biohub-esm/scripts/atlas_search.py uniprot P00520

# Download structure thumbnail image (plddt or pct-characterized)
python .agents/skills/biohub-esm/scripts/atlas_search.py thumbnail <protein_hash> --type plddt -o thumb.png

# Query SAE Feature dossier (0 to 16,383)
python .agents/skills/biohub-esm/scripts/atlas_search.py feature 100

# Submit batch lookup for multiple protein hashes (up to 500)
python .agents/skills/biohub-esm/scripts/atlas_search.py batch hashes.txt -o batch_results.zip
```

---

## Detailed API & SDK Documentation
For in-depth schemas, parameters, and full SDK methods, consult the dedicated reference files:

* [**references/README.md**](references/README.md): Master index, authentication, and endpoint summary.
* [**references/python_sdk_reference.md**](references/python_sdk_reference.md): Official `esm` package, `ESM3ForgeInferenceClient`, `ESMProtein`, and PDB generation.
* [**references/esm3_generation_and_design.md**](references/esm3_generation_and_design.md): ESM3 multimodal generation, iterative unmasking, and track prompts.
* [**references/structure_prediction_and_complexes.md**](references/structure_prediction_and_complexes.md): ESMFold2, multimer complexes, and all-atom co-folding with DNA/RNA/ligands.
* [**references/representations_and_interpretability.md**](references/representations_and_interpretability.md): ESMC embeddings, hidden states, Sparse Autoencoders (SAEs), and in-silico DMS.
* [**references/inverse_folding_and_tokenizers.md**](references/inverse_folding_and_tokenizers.md): 3D coordinates to sequence, track tokenization, and discrete tensor generation.
* [**references/esm_atlas_public_api.md**](references/esm_atlas_public_api.md): Public keyless ESM Atlas API (>1.1B structures, homology search, UniProt resolver).
* [**references/api_reference.md**](references/api_reference.md): High-level monolithic technical reference overview.
