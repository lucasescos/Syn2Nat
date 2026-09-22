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


def get_credentials_status() -> dict:
    status = {"primary": False, "alt": False}
    # 1. Environment variable
    if "ESM_API_KEY" in os.environ and os.environ["ESM_API_KEY"].strip():
        status["primary"] = True
    if "ESM_API_KEY_ALT" in os.environ and os.environ["ESM_API_KEY_ALT"].strip():
        status["alt"] = True

    # 2. Local workspace .env and user home .env
    for p in [Path.cwd() / ".env", Path.home() / ".env"]:
        vars_map = load_env_file(p)
        if "ESM_API_KEY" in vars_map and vars_map["ESM_API_KEY"]:
            status["primary"] = True
        if "ESM_API_KEY_ALT" in vars_map and vars_map["ESM_API_KEY_ALT"]:
            status["alt"] = True

    return status


def get_api_key() -> str | None:
    # 1. Environment variable
    for var in ["ESM_API_KEY", "ESM_API_KEY_ALT"]:
        if var in os.environ and os.environ[var].strip():
            return os.environ[var].strip()

    # 2. Local workspace .env & user home ~/.env
    for p in [Path.cwd() / ".env", Path.home() / ".env"]:
        vars_map = load_env_file(p)
        for var in ["ESM_API_KEY", "ESM_API_KEY_ALT"]:
            if var in vars_map and vars_map[var]:
                return vars_map[var]

    return None


def main():
    status = get_credentials_status()
    if status["primary"] and status["alt"]:
        print("[OK] ESM_API_KEY is configured (Primary active + Alternative fallback key available).")
        sys.exit(0)
    elif status["primary"] or status["alt"]:
        active_type = "Primary" if status["primary"] else "Alternative"
        print(f"[OK] ESM_API_KEY is configured ({active_type} key available).")
        sys.exit(0)
    else:
        print("[ERROR] Neither ESM_API_KEY nor ESM_API_KEY_ALT is set in environment, .env, or ~/.env", file=sys.stderr)
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
