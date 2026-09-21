#!/usr/bin/env python3
"""
Safe credential verification script for Biohub ESM.
Checks for the presence of ESM_API_KEY without exposing or printing the key value.
"""

import os
import sys
from pathlib import Path


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
                    key = key.strip()
                    val = val.strip().strip("'\"")
                    env_vars[key] = val
    except Exception:
        pass
    return env_vars


def get_api_key() -> str | None:
    # 1. Environment variable
    if "ESM_API_KEY" in os.environ and os.environ["ESM_API_KEY"].strip():
        return os.environ["ESM_API_KEY"].strip()

    # 2. Local workspace .env
    local_env = Path.cwd() / ".env"
    vars_local = load_env_file(local_env)
    if "ESM_API_KEY" in vars_local and vars_local["ESM_API_KEY"]:
        return vars_local["ESM_API_KEY"]

    # 3. User home directory ~/.env
    home_env = Path.home() / ".env"
    vars_home = load_env_file(home_env)
    if "ESM_API_KEY" in vars_home and vars_home["ESM_API_KEY"]:
        return vars_home["ESM_API_KEY"]

    return None


def main():
    key = get_api_key()
    if key:
        print("[OK] ESM_API_KEY is configured.")
        sys.exit(0)
    else:
        print("[ERROR] ESM_API_KEY is not set in environment, .env, or ~/.env", file=sys.stderr)
        print("\nTo set it safely without exposing typing in history, run:", file=sys.stderr)
        print("  Windows (PowerShell):", file=sys.stderr)
        print("    $key = Read-Host -AsSecureString 'Enter ESM_API_KEY'", file=sys.stderr)
        print("    $BSTR = [System.Runtime.InteropServices.Marshal]::SecureStringToBSTR($key)", file=sys.stderr)
        print("    $plain = [System.Runtime.InteropServices.Marshal]::PtrToStringAuto($BSTR)", file=sys.stderr)
        print("    Add-Content -Path ~/.env -Value \"ESM_API_KEY=$plain\"", file=sys.stderr)
        print("  Unix / macOS:", file=sys.stderr)
        print("    printf 'Enter ESM_API_KEY: ' && read -s k && echo && echo \"ESM_API_KEY=$k\" >> ~/.env", file=sys.stderr)
        print("\nGet your API key at: https://biohub.ai/developer-console/api-keys", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
