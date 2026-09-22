#!/usr/bin/env python3
"""
ESMC Client for Sequence-Derived Biophysical Representation Extraction.
Built on the official Biohub ESM SDK (`esm.sdk.esmc_client`).

Supports:
- esmc-6b-2024-12 (Flagship 6B parameter frontier model, 81 layers, 2560-dim)
- esmc-600m-2024-12 (Standard 600M parameter model, 37 layers, 1152-dim)
- esmc-300m-2024-12 (Lightweight 300M parameter model, 31 layers, 960-dim)
"""

import os
import sys
import time
from pathlib import Path
from typing import Any, Dict, List, Optional
import numpy as np
import torch
from esm.sdk import esmc_client
from esm.sdk.api import ESMProtein, LogitsConfig


def load_env_file(filepath: Path) -> dict:
    env_vars = {}
    if not filepath.is_file():
        return env_vars
    try:
        with open(filepath, "r", encoding="utf-8") as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith("#"):
                    continue
                if "=" in line:
                    key, val = line.split("=", 1)
                    env_vars[key.strip()] = val.strip().strip("'\"")
    except Exception:
        pass
    return env_vars


def get_api_keys() -> list[str]:
    """Return all configured Biohub API keys (ESM_API_KEY, ESM_API_KEY_ALT, ESM_API_KEY_3, etc.)."""
    keys = []
    # 1. Environment variables
    for k, v in os.environ.items():
        if k.startswith("ESM_API_KEY") and v.strip() and v.strip() not in keys:
            keys.append(v.strip())

    # 2. .env files (workspace root and user home)
    for p in [Path.cwd() / ".env", Path.home() / ".env"]:
        vars_map = load_env_file(p)
        for k, v in vars_map.items():
            if k.startswith("ESM_API_KEY") and v.strip() and v.strip() not in keys:
                keys.append(v.strip())
    return keys


def get_api_key(prefer_alt: bool = False) -> str:
    """Return an active API key, optionally preferring alternative keys."""
    keys = get_api_keys()
    if not keys:
        raise ValueError("No Biohub ESM API keys (ESM_API_KEY*) found in environment, .env, or ~/.env")
    if prefer_alt and len(keys) > 1:
        return keys[1]
    return keys[0]



class ESMCClient:
    """Wrapper around Biohub's ESMC client for biophysical representation extraction."""

    def __init__(
        self,
        token: Optional[str] = None,
        default_model: str = "esmc-6b-2024-12",
        url: str = "https://biohub.ai",
    ):
        self.token = token or get_api_key()
        self.url = url
        self.default_model = default_model
        self._cached_models: Dict[str, Any] = {}

    def _get_model(self, model_name: str):
        if model_name not in self._cached_models:
            self._cached_models[model_name] = esmc_client(
                model=model_name, url=self.url, token=self.token
            )
        return self._cached_models[model_name]

    MAX_LAYERS = {
        "esmc-300m-2024-12": 30,
        "esmc-600m-2024-12": 36,
        "esmc-6b-2024-12": 80,
    }

    def get_representation(
        self,
        sequence: str,
        model: Optional[str] = None,
        layer_idx: int = -2,
        per_residue: bool = False,
    ) -> Dict[str, Any]:
        """
        Extract sequence representations.

        Args:
            sequence: Amino acid sequence string.
            model: Model name ('esmc-6b-2024-12', 'esmc-600m-2024-12', 'esmc-300m-2024-12').
            layer_idx: Layer index (-2 for penultimate layer, -1 for final layer).
            per_residue: If True, returns (L, D) array where every amino acid has a D-dim vector.
                         If False, returns (D,) array mean-pooled across all amino acids.

        Returns:
            Dict containing 'vector' or 'matrix' (numpy array), 'shape', 'dim', 'latency_sec', etc.
        """
        seq_clean = "".join(sequence.split()).upper()
        if not seq_clean:
            raise ValueError("Sequence cannot be empty.")

        model_name = model or self.default_model
        client_model = self._get_model(model_name)
        max_layer = self.MAX_LAYERS.get(model_name, 80)

        # Resolve negative layer indices
        if layer_idx < 0:
            resolved_layer = max_layer + 1 + layer_idx
        else:
            resolved_layer = layer_idx

        t0 = time.time()
        protein = ESMProtein(sequence=seq_clean)
        protein_tensor = client_model.encode(protein)

        if per_residue:
            # Request per-residue hidden states for the resolved layer
            config = LogitsConfig(
                sequence=True,
                return_hidden_states=True,
                ith_hidden_layer=resolved_layer,
            )
            output = client_model.logits(protein_tensor, config)
            latency = time.time() - t0

            # Shape: [1, 1, L + 2, D] -> slice off batch and <BOS>/<EOS>
            tensor_l_d = output.hidden_states.float().squeeze(0).squeeze(0)[1:-1]
            arr = tensor_l_d.detach().cpu().numpy()  # (L, D)

            return {
                "model": model_name,
                "sequence_length": len(seq_clean),
                "granularity": "per_residue",
                "layer_idx": resolved_layer,
                "shape": list(arr.shape),
                "dim": int(arr.shape[-1]),
                "array": arr,
                "latency_sec": round(latency, 3),
            }
        else:
            # Request mean-pooled across sequence length for all layers
            config = LogitsConfig(sequence=True, return_mean_hidden_states=True)
            output = client_model.logits(protein_tensor, config)
            latency = time.time() - t0

            mean_hidden = output.mean_hidden_state.float().squeeze(0)  # [num_layers, hidden_dim]
            num_layers, hidden_dim = mean_hidden.shape

            selected_tensor = mean_hidden[layer_idx]
            vector_np = selected_tensor.detach().cpu().numpy()

            return {
                "model": model_name,
                "sequence_length": len(seq_clean),
                "granularity": "sequence_pooled",
                "num_layers": int(num_layers),
                "hidden_dim": int(hidden_dim),
                "selected_layer": layer_idx,
                "resolved_layer": resolved_layer,
                "shape": list(vector_np.shape),
                "dim": int(hidden_dim),
                "vector": vector_np,
                "l2_norm": float(np.linalg.norm(vector_np)),
                "mean_val": float(np.mean(vector_np)),
                "std_val": float(np.std(vector_np)),
                "latency_sec": round(latency, 3),
            }

