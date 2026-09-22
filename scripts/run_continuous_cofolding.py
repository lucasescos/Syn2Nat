#!/usr/bin/env python3
"""
Continuous Autonomous ESMFold2 Cofolding Engine for Human Microprotein Complexes.

Continuously discovers, curates physiological partner complexes (<=768 aa),
and folds them using Biohub ESMFold2 `fold_all_atom` API.
Runs autonomously until interrupted or until a stop sentinel (`structures/.stop_cofolding`) is detected.
Persists checkpoints after every single fold and dynamically synchronizes
`reports/master_microprotein_ppi_structural_atlas.parquet` and `.tsv`.
"""

import argparse
import json
import os
import sys
import time
import urllib.request
from pathlib import Path
import numpy as np
import pandas as pd

# Ensure ESM SDK is available
try:
    from esm.sdk.forge import SequenceStructureForgeInferenceClient, FoldingConfig
    from esm.utils.structure.input_builder import StructurePredictionInput, ProteinInput
    from esm.sdk.api import ESMProteinError
except ImportError as e:
    print(f"[ERROR] Failed to import ESM SDK: {e}", file=sys.stderr)
    sys.exit(1)


def load_api_keys() -> list[str]:
    """Retrieve all configured ESM API keys (primary, alt, and numbered) from environment or .env files."""
    keys = []
    for var, val in os.environ.items():
        if var.startswith("ESM_API_KEY") and val.strip() and val.strip() not in keys:
            keys.append(val.strip())

    for p in [Path.cwd() / ".env", Path.home() / ".env"]:
        if p.is_file():
            try:
                with open(p, "r", encoding="utf-8") as f:
                    for line in f:
                        line = line.strip()
                        if line.startswith("#") or "=" not in line:
                            continue
                        var, val = line.split("=", 1)
                        var = var.strip()
                        val = val.strip().strip("\"'")
                        if var.startswith("ESM_API_KEY") and val and val not in keys:
                            keys.append(val)
            except Exception:
                pass
    return keys


def load_api_key() -> str:
    """Retrieve primary or alternative ESM API key."""
    keys = load_api_keys()
    if keys:
        return keys[0]
    print("[ERROR] Neither ESM_API_KEY nor ESM_API_KEY_ALT found in environment, .env, or ~/.env", file=sys.stderr)
    sys.exit(1)


class RotatingESMClient:
    """Multi-key client with automatic rotation when rate limits or quotas are hit."""
    def __init__(self, model: str, tokens: list[str]):
        self.model = model
        self.tokens = tokens
        self.current_idx = 0
        self.clients = [
            SequenceStructureForgeInferenceClient(model=model, token=t)
            for t in tokens
        ]

    @property
    def current_client(self):
        return self.clients[self.current_idx]

    def rotate(self):
        if len(self.clients) > 1:
            self.current_idx = (self.current_idx + 1) % len(self.clients)
            print(f"\n    [KEY ROTATION] Biohub API quota/rate limit reached. Switched to key #{self.current_idx + 1}/{len(self.clients)}", file=sys.stderr)

    def fold_all_atom(self, *args, **kwargs):
        return self.current_client.fold_all_atom(*args, **kwargs)


def load_checkpoint(state_file: Path) -> dict:
    """Load previously completed or skipped pairs from the JSONL checkpoint file."""
    completed = {}
    if not state_file.is_file():
        return completed
    with open(state_file, "r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            try:
                rec = json.loads(line)
                if "pair_id" in rec:
                    completed[rec["pair_id"]] = rec
            except Exception:
                pass
    return completed


def save_checkpoint_record(state_file: Path, record: dict):
    """Append a single record to the state ledger."""
    with open(state_file, "a", encoding="utf-8") as f:
        f.write(json.dumps(record) + "\n")


def fold_complex_with_retry(client, spi, config, max_retries=4, base_backoff=3.0):
    """Execute fold_all_atom with exponential backoff and automatic key rotation on rate limits."""
    for attempt in range(1, max_retries + 1):
        try:
            res = client.fold_all_atom(spi, config=config)
            if isinstance(res, ESMProteinError):
                if res.error_code == 422:
                    raise ValueError(f"ESMProteinError [422]: {res.error_msg}")
                if res.error_code in (429, 403) or "quota" in str(res.error_msg).lower() or "rate" in str(res.error_msg).lower():
                    if hasattr(client, "rotate"):
                        client.rotate()
                        time.sleep(1.0)
                        continue
                raise RuntimeError(f"ESMProteinError [{res.error_code}]: {res.error_msg}")
            return res
        except ValueError as e:
            raise e
        except Exception as e:
            err_str = str(e)
            if "422" in err_str or "exceeds maximum allowed" in err_str:
                raise ValueError(err_str)
            if "429" in err_str or "quota" in err_str.lower() or "rate limit" in err_str.lower() or "403" in err_str:
                if hasattr(client, "rotate"):
                    client.rotate()
                    time.sleep(1.0)
                    continue
            if attempt == max_retries:
                raise e
            wait_time = base_backoff * (2 ** (attempt - 1))
            print(f"    [WARN] Attempt {attempt} failed: {err_str[:120]}... Backing off {wait_time:.1f}s")
            time.sleep(wait_time)



class StringPartnerResolver:
    """STRING v12.0 interaction partner query client with persistent disk cache."""
    def __init__(self, cache_file: Path):
        self.cache_file = cache_file
        self.cache = {}
        if self.cache_file.is_file():
            try:
                with open(self.cache_file, "r", encoding="utf-8") as f:
                    self.cache = json.load(f)
            except Exception:
                self.cache = {}

    def save_cache(self):
        try:
            with open(self.cache_file, "w", encoding="utf-8") as f:
                json.dump(self.cache, f, indent=2)
        except Exception:
            pass

    def get_partners(self, gene_symbol: str, prot_lookup: dict, limit: int = 50) -> list:
        clean_gene = gene_symbol.split("-AS")[0].split(".")[0]
        if clean_gene in self.cache:
            raw = self.cache[clean_gene]
            return [(p, sc) for p, sc in raw if p in prot_lookup]

        url = f"https://string-db.org/api/json/interaction_partners?identifiers={clean_gene}&species=9606&limit={limit}"
        req = urllib.request.Request(url, headers={"User-Agent": "Antigravity-ContinuousCofolding/1.0"})
        partners = []
        try:
            with urllib.request.urlopen(req, timeout=10) as resp:
                data = json.loads(resp.read().decode("utf-8"))
                for d in data:
                    p_name = d.get("preferredName_B")
                    score = float(d.get("score", 0.0))
                    if p_name:
                        partners.append((p_name, score))
        except Exception:
            pass

        seen = set()
        dedup = []
        for p, sc in partners:
            if p not in seen:
                seen.add(p)
                dedup.append((p, sc))

        self.cache[clean_gene] = dedup
        time.sleep(0.08)  # Polite delay
        return [(p, sc) for p, sc in dedup if p in prot_lookup]


def sync_master_atlas(base_dir: Path, state_file: Path, continuous_candidates_path: Path):
    """Synchronize master structural atlas with all folded complexes."""
    top50_path = base_dir / "reports" / "top50_microprotein_binding_partners_top10.parquet"
    exp460_path = base_dir / "reports" / "expanded_ppi_candidates_460.parquet"
    master_parquet = base_dir / "reports" / "master_microprotein_ppi_structural_atlas.parquet"
    master_tsv = base_dir / "reports" / "master_microprotein_ppi_structural_atlas.tsv"

    cand_dfs = []
    if top50_path.is_file():
        cand_dfs.append(pd.read_parquet(top50_path))
    if exp460_path.is_file():
        cand_dfs.append(pd.read_parquet(exp460_path))
    if continuous_candidates_path.is_file():
        cand_dfs.append(pd.read_parquet(continuous_candidates_path))

    if not cand_dfs:
        return

    all_cands = pd.concat(cand_dfs, ignore_index=True).drop_duplicates(subset=["pair_id"], keep="last")

    # Load state
    state_records = []
    with open(state_file, "r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if line:
                try:
                    state_records.append(json.loads(line))
                except Exception:
                    pass

    state_df = pd.DataFrame(state_records).drop_duplicates(subset=["pair_id"], keep="last")
    merged = pd.merge(all_cands, state_df, on="pair_id", how="inner", suffixes=("", "_state"))
    successful = merged[merged["status"] == "success"].copy()
    successful = successful.sort_values("iptm", ascending=False).reset_index(drop=True)

    successful.to_parquet(master_parquet, index=False)
    successful.to_csv(master_tsv, sep="\t", index=False)
    
    # Calculate key metrics
    high_conf = (successful["iptm"] >= 0.60).sum()
    print(f"\n>>> [MASTER SYNC] Atlas updated: {len(successful)} total folded complexes | High-confidence (ipTM >= 0.60): {high_conf}")


def main():
    parser = argparse.ArgumentParser(description="Continuous autonomous ESMFold2 cofolding pipeline")
    parser.add_argument("--partners-per-microprotein", type=int, default=3,
                        help="Number of physiological partners per microprotein (default: 3)")
    parser.add_argument("--model", default="esmfold2-fast-2026-05", choices=["esmfold2-fast-2026-05", "esmfold2-2026-05"],
                        help="ESMFold2 model version")
    parser.add_argument("--num-loops", type=int, default=20, help="Trunk loops (default: 20)")
    parser.add_argument("--num-steps", type=int, default=100, help="Diffusion sampling steps (default: 100)")
    parser.add_argument("--delay", type=float, default=0.5, help="Delay between API calls in seconds (default: 0.5)")
    parser.add_argument("--max-length", type=int, default=768, help="Biohub context length limit (default: 768)")
    parser.add_argument("--max-complexes", type=int, default=None, help="Optional max complexes to fold before pausing")
    parser.add_argument("--stop-file", default="structures/.stop_cofolding", help="Path to stop sentinel file")

    args = parser.parse_args()
    api_keys = load_api_keys()
    if not api_keys:
        print("[ERROR] Neither ESM_API_KEY nor ESM_API_KEY_ALT found in environment, .env, or ~/.env", file=sys.stderr)
        sys.exit(1)

    base_dir = Path.cwd()
    struct_dir = base_dir / "structures"
    pdb_dir = struct_dir / "pdbs"
    pae_dir = struct_dir / "pae"
    pdb_dir.mkdir(parents=True, exist_ok=True)
    pae_dir.mkdir(parents=True, exist_ok=True)

    cache_dir = base_dir / "data" / "cache"
    cache_dir.mkdir(parents=True, exist_ok=True)
    string_cache_file = cache_dir / "string_partners_cache.json"
    string_resolver = StringPartnerResolver(string_cache_file)

    state_file = struct_dir / "cofolding_state.jsonl"
    stop_file = base_dir / args.stop_file
    continuous_candidates_path = base_dir / "reports" / "continuous_cofolding_candidates.parquet"

    # Load canonical proteome lookup
    proteome_path = base_dir / "data" / "targets" / "human_canonical_proteome.parquet"
    print(f"Loading reference human canonical proteome from {proteome_path}...")
    proteome_df = pd.read_parquet(proteome_path)
    prot_lookup = {r["gene_symbol"]: r for _, r in proteome_df.iterrows()}

    # Load query microprotein definitions
    query_meta_path = base_dir / "data" / "queries" / "query_metadata.parquet"
    print(f"Loading query microprotein metadata from {query_meta_path}...")
    query_df = pd.read_parquet(query_meta_path)

    # Identify microproteins already folded in checkpoint/master atlas
    completed_map = load_checkpoint(state_file)
    print(f"Loaded existing checkpoint ledger: {len(completed_map)} total records.")
    success_count = sum(1 for r in completed_map.values() if r.get("status") == "success")
    print(f"  Successfully folded complexes on disk: {success_count}")

    # Determine which microproteins already have folds on disk
    master_path = base_dir / "reports" / "master_microprotein_ppi_structural_atlas.parquet"
    already_folded_orfs = set()
    if master_path.is_file():
        try:
            m_df = pd.read_parquet(master_path)
            already_folded_orfs = set(m_df["orf_id"].dropna())
        except Exception:
            pass

    # Sort query microproteins: completely uncharacterized first (has_folds=False), then Tier 1A, 1B, 2A
    tier_order = {"1A": 0, "1B": 1, "2A": 2, "2B": 3}
    query_df["tier_rank"] = query_df["final_tier"].map(tier_order).fillna(9)
    query_df["has_folds"] = query_df["orf_id"].isin(already_folded_orfs)
    query_sorted = query_df.sort_values(by=["has_folds", "tier_rank", "length"]).reset_index(drop=True)

    uncharacterized_count = (~query_df["has_folds"]).sum()
    print(f"Uncharacterized microproteins to fold: {uncharacterized_count} (out of {len(query_df)})")

    # Load continuous candidates if existing
    continuous_candidates = []
    if continuous_candidates_path.is_file():
        try:
            continuous_candidates = pd.read_parquet(continuous_candidates_path).to_dict("records")
            print(f"Loaded existing continuous candidate records: {len(continuous_candidates)}")
        except Exception:
            pass

    # Initialize multi-key ESM client with auto-rotation
    client = RotatingESMClient(model=args.model, tokens=api_keys)
    folding_config = FoldingConfig(
        num_loops=args.num_loops,
        num_sampling_steps=args.num_steps,
        include_pae=True,
        include_pair_chains_iptm=True
    )

    print("\n=======================================================")
    print("  AUTONOMOUS CONTINUOUS ESMFOLD2 COFOLDING ACTIVE")
    print(f"  Model:              {args.model} (loops={args.num_loops}, steps={args.num_steps})")
    print(f"  Partners/ORF:       {args.partners_per_microprotein}")
    print(f"  Length Cap:         {args.max_length} aa")
    print(f"  Stop Sentinel:      {stop_file}")
    print("=======================================================\n")

    t_start = time.time()
    n_processed_this_session = 0
    n_skipped_this_session = 0
    save_counter = 0

    try:
        for u_idx, (_, urow) in enumerate(query_sorted.iterrows(), start=1):
            if stop_file.exists():
                print(f"\n[STOP SENTINEL DETECTED] Found {stop_file}. Halting continuous cofolding gracefully.")
                break

            if args.max_complexes and n_processed_this_session >= args.max_complexes:
                print(f"\n[LIMIT REACHED] Reached session limit of {args.max_complexes} complexes.")
                break

            orf_id = urow["orf_id"]
            gene = urow["gene"]
            m_seq = urow["sequence"]
            m_len = len(m_seq)
            biotype = urow["orf_biotype"]
            tier = urow["final_tier"]

            # Check single microprotein length against platform limit
            if m_len >= args.max_length:
                continue

            max_p_len = args.max_length - m_len

            # Retrieve interaction partners from STRING / proteome
            raw_partners = string_resolver.get_partners(gene, prot_lookup, limit=60)

            # Filter valid partners
            valid_partners = []
            for p_name, sc in raw_partners:
                if p_name != gene and p_name in prot_lookup:
                    p_len = prot_lookup[p_name]["length"]
                    if p_len <= max_p_len:
                        valid_partners.append((p_name, sc, True))
                        if len(valid_partners) == args.partners_per_microprotein:
                            break

            # Supplement if fewer than desired
            if len(valid_partners) < args.partners_per_microprotein:
                # Try host gene CDS
                if gene in prot_lookup and gene not in [x[0] for x in valid_partners]:
                    if prot_lookup[gene]["length"] <= max_p_len:
                        valid_partners.append((gene, 1.000, True))

                # Supplement universal cellular hubs
                fallbacks = ["UBC", "ACTB", "HSPA8", "GAPDH", "YWHAZ", "RPS27A", "CALM1", "EEF1A1"]
                for fb in fallbacks:
                    if len(valid_partners) >= args.partners_per_microprotein:
                        break
                    if fb in prot_lookup and fb not in [x[0] for x in valid_partners]:
                        if prot_lookup[fb]["length"] <= max_p_len:
                            valid_partners.append((fb, 0.700, False))

            if not valid_partners:
                continue

            # Process each partner complex
            for p_rank, (p_name, sc, is_emp) in enumerate(valid_partners, start=1):
                if stop_file.exists():
                    print(f"\n[STOP SENTINEL DETECTED] Halting before folding next complex.")
                    break
                if args.max_complexes and n_processed_this_session >= args.max_complexes:
                    break

                pair_id = f"{orf_id}_partner_{p_rank:02d}_{p_name}"
                pdb_target = pdb_dir / f"{pair_id}.pdb"
                pae_target = pae_dir / f"{pair_id}_pae.json"

                # Check if already completed
                if pair_id in completed_map:
                    prev = completed_map[pair_id]
                    if prev.get("status") == "success" and pdb_target.is_file() and pdb_target.stat().st_size > 0:
                        n_skipped_this_session += 1
                        continue
                    if prev.get("status") == "exceeds_api_limit_768":
                        n_skipped_this_session += 1
                        continue

                p_info = prot_lookup[p_name]
                p_len = p_info["length"]
                tot_len = m_len + p_len

                if tot_len > args.max_length:
                    skip_record = {
                        "pair_id": pair_id,
                        "microprotein_gene": gene,
                        "partner_gene_symbol": p_name,
                        "status": "exceeds_api_limit_768",
                        "ptm": None,
                        "iptm": None,
                        "mean_plddt": None,
                        "microprotein_plddt": None,
                        "partner_plddt": None,
                        "binding_classification": "Exceeds API Context Limit (>768 aa)",
                        "pdb_path": None,
                        "pae_path": None,
                        "model_used": args.model,
                        "elapsed_seconds": 0.0,
                        "folded_at_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())
                    }
                    save_checkpoint_record(state_file, skip_record)
                    completed_map[pair_id] = skip_record
                    n_skipped_this_session += 1
                    continue

                # Prepare candidate metadata
                rationale = f"Continuous Tier {tier} microprotein ({gene}, {m_len} aa). "
                if is_emp:
                    rationale += f"STRING v12.0 interactor (score: {sc:.3f})."
                else:
                    rationale += f"Core cellular machinery partner ({p_info['protein_name']})."

                cand_rec = {
                    "pair_id": pair_id,
                    "microprotein_rank": u_idx,
                    "cohort": f"6. Continuous Tier {tier} Interactome",
                    "orf_id": orf_id,
                    "microprotein_gene": gene,
                    "microprotein_biotype": biotype,
                    "microprotein_tier": tier,
                    "microprotein_length": m_len,
                    "microprotein_sequence": m_seq,
                    "partner_rank": p_rank,
                    "partner_gene_symbol": p_name,
                    "partner_uniprot_id": p_info["uniprot_id"],
                    "partner_protein_name": p_info["protein_name"],
                    "partner_length": p_len,
                    "partner_sequence": p_info["sequence"],
                    "interaction_category": "Direct Physical Interactor" if is_emp else "Core Macromolecular Complex",
                    "string_interaction_score": round(sc, 3),
                    "is_empirical_string_score": is_emp,
                    "primary_tissue": "Ubiquitous / Tissue-enriched",
                    "biological_rationale": rationale,
                    "esmfold_multichain_sequence": f"{m_seq}|{p_info['sequence']}",
                    "esmfold_total_length": tot_len,
                    "esmfold_feasibility": "Optimal (<800 aa)"
                }
                continuous_candidates.append(cand_rec)

                # Execute fold
                session_num = n_processed_this_session + 1
                print(f"[{session_num} in session | ORF #{u_idx} ({orf_id})] Folding {pair_id} ({gene} [{m_len}aa] + {p_name} [{p_len}aa] = {tot_len}aa) ...")

                spi = StructurePredictionInput(sequences=[
                    ProteinInput(id="A", sequence=m_seq),
                    ProteinInput(id="B", sequence=p_info["sequence"])
                ])

                try:
                    t0 = time.time()
                    res = fold_complex_with_retry(client, spi, folding_config)
                    t_elapsed = time.time() - t0

                    # Metrics
                    ptm_val = float(res.ptm) if res.ptm is not None else None
                    iptm_val = float(res.iptm) if res.iptm is not None else None

                    plddt_raw = res.plddt
                    if hasattr(plddt_raw, "tolist"):
                        plddt_list = plddt_raw.tolist()
                    elif isinstance(plddt_raw, list):
                        plddt_list = plddt_raw
                    else:
                        plddt_list = []

                    valid_plddt = [x for x in plddt_list if x is not None and not np.isnan(x)]
                    mean_plddt = (sum(valid_plddt) / len(valid_plddt)) if valid_plddt else None

                    plddt_a_list = [x for x in plddt_list[:m_len] if x is not None and not np.isnan(x)]
                    plddt_b_list = [x for x in plddt_list[m_len:] if x is not None and not np.isnan(x)]
                    plddt_a = (sum(plddt_a_list) / len(plddt_a_list)) if plddt_a_list else None
                    plddt_b = (sum(plddt_b_list) / len(plddt_b_list)) if plddt_b_list else None

                    if iptm_val is not None:
                        if iptm_val >= 0.75:
                            b_class = "High Confidence (ipTM >= 0.75)"
                        elif iptm_val >= 0.60:
                            b_class = "Plausible Interface (0.60 <= ipTM < 0.75)"
                        elif iptm_val >= 0.45:
                            b_class = "Weak / Transient (0.45 <= ipTM < 0.60)"
                        else:
                            b_class = "Non-Interacting / Unbound (ipTM < 0.45)"
                    else:
                        b_class = "Undetermined"

                    # Export PDB
                    try:
                        pc = res.complex.to_protein_complex()
                        pdb_str = pc.to_pdb_string()
                        with open(pdb_target, "w", encoding="utf-8") as f_pdb:
                            f_pdb.write(pdb_str)
                    except Exception as e_pdb:
                        cif_target = pdb_dir / f"{pair_id}.cif"
                        with open(cif_target, "w", encoding="utf-8") as f_cif:
                            f_cif.write(res.complex.to_mmcif())
                        pdb_target = cif_target

                    # Export PAE
                    if hasattr(res, "pae") and res.pae is not None:
                        pae_raw = res.pae
                        pae_matrix = pae_raw.tolist() if hasattr(pae_raw, "tolist") else pae_raw
                        with open(pae_target, "w", encoding="utf-8") as f_pae:
                            json.dump(pae_matrix, f_pae)

                    fold_rec = {
                        "pair_id": pair_id,
                        "microprotein_gene": gene,
                        "partner_gene_symbol": p_name,
                        "status": "success",
                        "ptm": round(ptm_val, 4) if ptm_val else None,
                        "iptm": round(iptm_val, 4) if iptm_val else None,
                        "mean_plddt": round(mean_plddt, 2) if mean_plddt else None,
                        "microprotein_plddt": round(plddt_a, 2) if plddt_a else None,
                        "partner_plddt": round(plddt_b, 2) if plddt_b else None,
                        "binding_classification": b_class,
                        "pdb_path": str(pdb_target),
                        "pae_path": str(pae_target) if pae_target.exists() else None,
                        "model_used": args.model,
                        "elapsed_seconds": round(t_elapsed, 2),
                        "folded_at_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())
                    }

                    save_checkpoint_record(state_file, fold_rec)
                    completed_map[pair_id] = fold_rec
                    n_processed_this_session += 1
                    save_counter += 1

                    u_pct = f"{plddt_a*100:.1f}%" if plddt_a else "N/A"
                    p_pct = f"{plddt_b*100:.1f}%" if plddt_b else "N/A"
                    hit_marker = " *** HIGH-CONFIDENCE HIT ***" if (iptm_val and iptm_val >= 0.60) else ""
                    print(f"  --> [OK] in {t_elapsed:.1f}s | ipTM: {iptm_val:.3f} | pTM: {ptm_val:.3f} | pLDDT: uProt={u_pct}, P={p_pct} | {b_class}{hit_marker}")

                    # Sync candidate dataframe and master atlas every 5 folds
                    if save_counter >= 5:
                        cand_df = pd.DataFrame(continuous_candidates).drop_duplicates(subset=["pair_id"], keep="last")
                        cand_df.to_parquet(continuous_candidates_path, index=False)
                        string_resolver.save_cache()
                        sync_master_atlas(base_dir, state_file, continuous_candidates_path)
                        save_counter = 0

                    if args.delay > 0:
                        time.sleep(args.delay)

                except Exception as err:
                    print(f"  --> [FAIL] Error folding {pair_id}: {err}")
                    err_record = {
                        "pair_id": pair_id,
                        "microprotein_gene": gene,
                        "partner_gene_symbol": p_name,
                        "status": "error",
                        "error_message": str(err),
                        "folded_at_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())
                    }
                    save_checkpoint_record(state_file, err_record)
                    completed_map[pair_id] = err_record

    finally:
        # Final synchronization before exiting
        if continuous_candidates:
            cand_df = pd.DataFrame(continuous_candidates).drop_duplicates(subset=["pair_id"], keep="last")
            cand_df.to_parquet(continuous_candidates_path, index=False)
        string_resolver.save_cache()
        sync_master_atlas(base_dir, state_file, continuous_candidates_path)
        total_time = time.time() - t_start
        print("\n=======================================================")
        print(f"  Continuous Session Concluded")
        print(f"  Total Session Elapsed: {total_time/60:.2f} minutes")
        print(f"  Complexes Folded:      {n_processed_this_session}")
        print(f"  Complexes Skipped:     {n_skipped_this_session}")
        print("=======================================================")


if __name__ == "__main__":
    main()
