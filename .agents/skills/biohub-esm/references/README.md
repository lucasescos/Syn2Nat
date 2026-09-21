# Biohub ESM Reference Documentation Index

This directory contains modular technical references for interacting with the Biohub Evolutionary Scale Modeling (ESM) platform (`https://biohub.ai`) and the official `esm` Python SDK.

---

## Quick Routing Matrix: Where to Look

When a specific biological or modeling task is requested, refer directly to the corresponding documentation and script:

| Task / Intent | Dedicated Reference File | Primary API Endpoint / SDK Class | CLI Script |
| :--- | :--- | :--- | :--- |
| **Generative Protein Design & Motif Scaffolding** | [`esm3_generation_and_design.md`](esm3_generation_and_design.md) | `POST /api/v1/generate` <br> `ESM3ForgeInferenceClient` | `scripts/design_sequence.py` |
| **Python SDK (`esm` package) & `ESMProtein`** | [`python_sdk_reference.md`](python_sdk_reference.md) | `esm.sdk.forge` <br> `esm.sdk.api.ESMProtein` | `test_esm3_sdk.py` |
| **Structure Prediction & All-Atom Complexes** | [`structure_prediction_and_complexes.md`](structure_prediction_and_complexes.md) | `POST /api/v1/fold` <br> `POST /api/v1/fold_all_atom` | `scripts/fold.py` |
| **Embeddings, Logits & Interpretability (SAEs)** | [`representations_and_interpretability.md`](representations_and_interpretability.md) | `POST /api/v1/logits` <br> ESMC models | `scripts/embed.py` <br> `scripts/sae_features.py` <br> `scripts/deep_mutational_scan.py` |
| **Inverse Folding & Tokenizer Tracks** | [`inverse_folding_and_tokenizers.md`](inverse_folding_and_tokenizers.md) | `POST /api/v1/inverse_fold` <br> `POST /api/v1/encode`, `/decode` | `scripts/inverse_fold.py` <br> `scripts/tokenize_tracks.py` |
| **Public ESM Atlas (>1.1B Structures)** | [`esm_atlas_public_api.md`](esm_atlas_public_api.md) | `GET https://biohub.ai/esm/protein/api/v1alpha1/...` (Keyless) | `scripts/atlas_search.py` |

---

## Authentication & Base URLs

### Platform Endpoints (Authenticated)
* **Base URL**: `https://biohub.ai/api/v1`
* **Authorization**: `Bearer <ESM_API_KEY>`
* **Credential Source**: Environment variable `ESM_API_KEY` (in `~/.env` or `.env`).
* **Verification Script**: `python .agents/skills/biohub-esm/scripts/verify_credentials.py`

### ESM Atlas Public API (Keyless)
* **Base URL**: `https://biohub.ai/esm/protein/api/v1alpha1`
* **Authorization**: None required (open public alpha).
