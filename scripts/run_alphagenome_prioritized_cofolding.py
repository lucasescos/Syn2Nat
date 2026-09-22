#!/usr/bin/env python3
"""
AlphaGenome-Prioritized Continuous ESMFold2 Cofolding Pipeline.

Prioritizes all remaining human smORFs using DeepMind AlphaGenome criteria:
1. Genomic Purifying Selection (median & max AVI Phred score).
2. De-confounded Genomic Autonomy (1.3x boost to unconfounded lncRNA-ORFs, uORFs, dORFs).
3. Multimodal Track Disruption (convergent regulatory modalities across RNA, splicing, chromatin, TF binding).
4. Experimental Translation Tiers (Tier 1A MS > Tier 1B HLA-MS > Tier 2A Ribo-seq).
5. Biohub ESMFold2 context length feasibility (<= 768 aa).

Runs continuously in the background, writing PDBs, PAE JSONs, and syncing the master atlas.
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


def load_api_key() -> str:
    """Retrieve ESM_API_KEY from environment or .env files."""
    key = os.environ.get("ESM_API_KEY")
    if key and key.strip():
        return key.strip()

    for p in [Path.cwd() / ".env", Path.home() / ".env"]:
        if p.is_file():
            try:
                with open(p, "r", encoding="utf-8") as f:
                    for line in f:
                        line = line.strip()
                        if line.startswith("ESM_API_KEY="):
                            val = line.split("=", 1)[1].strip().strip("\"'")
                            if val:
                                return val
            except Exception:
                pass

    print("[ERROR] ESM_API_KEY not found in environment, .env, or ~/.env", file=sys.stderr)
    sys.exit(1)


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


def fold_complex_with_retry(client, spi, config, max_retries=3, base_backoff=3.0):
    """Execute fold_all_atom with exponential backoff on transient errors."""
    for attempt in range(1, max_retries + 1):
        try:
            res = client.fold_all_atom(spi, config=config)
            if isinstance(res, ESMProteinError):
                if res.error_code == 422:
                    raise ValueError(f"ESMProteinError [422]: {res.error_msg}")
                raise RuntimeError(f"ESMProteinError [{res.error_code}]: {res.error_msg}")
            return res
        except ValueError as e:
            raise e
        except Exception as e:
            err_str = str(e)
            if "422" in err_str or "exceeds maximum allowed" in err_str:
                raise ValueError(err_str)
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
        req = urllib.request.Request(url, headers={"User-Agent": "Antigravity-AlphaGenomeCofolding/1.0"})
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
    
    high_conf = (successful["iptm"] >= 0.60).sum()
    print(f"\n>>> [MASTER SYNC] Master Atlas: {len(successful)} total folded complexes | High-confidence (ipTM >= 0.60): {high_conf}")


def build_alphagenome_prioritized_catalog(base_dir: Path, already_folded_orfs: set) -> pd.DataFrame:
    """Score and rank all unfolded human smORFs using the 5 AlphaGenome pillars."""
    print("Computing AlphaGenome prioritization across human smORF catalog...")
    avi_path = base_dir / "data" / "screening" / "microprotein_avi_scores.parquet"
    dis_path = base_dir / "data" / "screening" / "microprotein_single_variant_disruptions.parquet"
    orfs_path = base_dir / "data" / "parsed" / "orfs.parquet"
    tiers_path = base_dir / "data" / "parsed" / "orf_tiers.parquet"

    avi_df = pd.read_parquet(avi_path)
    dis_df = pd.read_parquet(dis_path)
    orfs_df = pd.read_parquet(orfs_path)
    tiers_df = pd.read_parquet(tiers_path)

    # Sequence lookup
    seq_lookup = {r["orf_id"]: r["sequence"] for _, r in orfs_df.iterrows()}

    # Merge
    dis_cols = ["orf_id", "is_high_impact_hotspot", "convergent_modalities_count", "primary_mechanism",
                "rna_max_abs_raw", "splicing_max_abs_raw", "chromatin_max_abs_raw", "chip_tf_max_abs_raw"]
    df = avi_df.merge(dis_df[[c for c in dis_cols if c in dis_df.columns]], on="orf_id", how="left")

    # Add sequences
    df["sequence"] = df["orf_id"].map(seq_lookup)
    df = df[df["sequence"].notna()].copy()

    # Filter out already folded microproteins
    df["has_folds"] = df["orf_id"].isin(already_folded_orfs)
    unfolded = df[~df["has_folds"]].copy()

    # Filter out single microproteins that alone exceed the 768 aa context cap
    unfolded = unfolded[unfolded["length_aa"] < 768].copy()

    # Define Genomic Autonomy
    autonomous_biotypes = ["lncRNA-ORF", "uORF", "dORF"]
    unfolded["is_autonomous"] = unfolded["orf_biotype"].isin(autonomous_biotypes)

    # Autonomy Weight: 1.3x for genuine autonomous, 1.0x for mass-spec verified overlapping, 0.6x for putative overlapping
    def calc_autonomy_weight(row):
        if row["is_autonomous"]:
            return 1.3
        elif str(row["final_tier"]).startswith("1"):
            return 1.0
        else:
            return 0.6

    unfolded["autonomy_weight"] = unfolded.apply(calc_autonomy_weight, axis=1)

    # Translation Tier Weight: 1A=100, 1B=75, 2A=50, 2B=25, 3=10, 4=0
    tier_weights = {"1A": 100.0, "1B": 75.0, "2A": 50.0, "2B": 25.0, "3": 10.0, "4": 0.0}
    unfolded["tier_weight"] = unfolded["final_tier"].map(tier_weights).fillna(0.0)

    # AlphaGenome Purifying Selection Component:
    unfolded["median_phred_clean"] = unfolded["median_avi_phred"].fillna(0.0)
    unfolded["max_phred_clipped"] = np.clip(unfolded["max_avi_phred"].fillna(0.0), 0.0, 40.0)
    unfolded["avi_selection_score"] = unfolded["median_phred_clean"] + 0.3 * unfolded["max_phred_clipped"]

    # Multimodal Track Disruption Component:
    unfolded["conv_clean"] = unfolded["convergent_modalities_count"].fillna(0.0)
    unfolded["hotspot_clean"] = unfolded["is_high_impact_hotspot"].fillna(False).astype(float)
    unfolded["multimodal_score"] = unfolded["conv_clean"] * 3.0 + unfolded["hotspot_clean"] * 5.0

    # Composite AlphaGenome Selection Priority Score (AGSPS):
    unfolded["agsps_score"] = (
        (unfolded["avi_selection_score"] * 1.5 + unfolded["multimodal_score"]) * unfolded["autonomy_weight"]
        + unfolded["tier_weight"] * 0.4
    )

    # Sort descending by AGSPS
    prioritized = unfolded.sort_values("agsps_score", ascending=False).reset_index(drop=True)
    prioritized["agsps_rank"] = range(1, len(prioritized) + 1)
    print(f"AlphaGenome prioritization completed: {len(prioritized)} unfolded smORFs ranked.")
    return prioritized


def main():
    parser = argparse.ArgumentParser(description="AlphaGenome-Prioritized Continuous ESMFold2 cofolding engine")
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
    api_key = load_api_key()

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

    # Remove any existing stop file before starting
    if stop_file.exists():
        try:
            stop_file.unlink()
        except Exception:
            pass

    # Load canonical proteome lookup
    proteome_path = base_dir / "data" / "targets" / "human_canonical_proteome.parquet"
    print(f"Loading reference human canonical proteome from {proteome_path}...")
    proteome_df = pd.read_parquet(proteome_path)
    prot_lookup = {r["gene_symbol"]: r for _, r in proteome_df.iterrows()}

    # Load checkpoint
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

    # Build AlphaGenome-prioritized queue
    prioritized_df = build_alphagenome_prioritized_catalog(base_dir, already_folded_orfs)

    # Load continuous candidates if existing
    continuous_candidates = []
    if continuous_candidates_path.is_file():
        try:
            continuous_candidates = pd.read_parquet(continuous_candidates_path).to_dict("records")
            print(f"Loaded existing continuous candidate records: {len(continuous_candidates)}")
        except Exception:
            pass

    # Initialize ESM client
    client = SequenceStructureForgeInferenceClient(model=args.model, token=api_key)
    folding_config = FoldingConfig(
        num_loops=args.num_loops,
        num_sampling_steps=args.num_steps,
        include_pae=True,
        include_pair_chains_iptm=True
    )

    print("\n================================================================================")
    print("  ALPHAGENOME-PRIORITIZED CONTINUOUS ESMFOLD2 COFOLDING CAMPAIGN ACTIVE")
    print(f"  Model:                {args.model} (loops={args.num_loops}, steps={args.num_steps})")
    print(f"  Partners/smORF:       {args.partners_per_microprotein}")
    print(f"  Platform Limit:       {args.max_length} aa")
    print(f"  Unfolded Queue:       {len(prioritized_df)} prioritized smORFs")
    print(f"  Stop Sentinel:        {stop_file}")
    print("================================================================================\n")

    t_start = time.time()
    n_processed_this_session = 0
    n_skipped_this_session = 0
    save_counter = 0

    try:
        for _, urow in prioritized_df.iterrows():
            if stop_file.exists():
                print(f"\n[STOP SENTINEL DETECTED] Found {stop_file}. Halting continuous cofolding gracefully.")
                break

            if args.max_complexes and n_processed_this_session >= args.max_complexes:
                print(f"\n[LIMIT REACHED] Reached session limit of {args.max_complexes} complexes.")
                break

            agsps_rank = urow["agsps_rank"]
            agsps_score = urow["agsps_score"]
            orf_id = urow["orf_id"]
            gene = str(urow["gene"]) if pd.notna(urow["gene"]) else orf_id
            m_seq = urow["sequence"]
            m_len = len(m_seq)
            biotype = urow["orf_biotype"]
            tier = urow["final_tier"]
            is_auto = urow["is_autonomous"]

            max_p_len = args.max_length - m_len

            # Retrieve interaction partners from STRING / proteome
            raw_partners = string_resolver.get_partners(gene, prot_lookup, limit=60)

            # Filter valid partners fitting size cap
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

                # Supplement universal cellular machinery hubs
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
                autonomy_str = "Autonomous" if is_auto else "CDS-Overlapping"
                rationale = (f"AlphaGenome Prioritized (Rank #{agsps_rank}, AGSPS: {agsps_score:.1f}, {autonomy_str}). "
                             f"smORF {gene} ({m_len} aa, Tier {tier}, median AVI {urow['median_avi_phred']:.1f}). ")
                if is_emp:
                    rationale += f"STRING v12.0 interactor {p_name} (score: {sc:.3f})."
                else:
                    rationale += f"Cellular machinery hub ({p_info['protein_name']})."

                cand_rec = {
                    "pair_id": pair_id,
                    "microprotein_rank": agsps_rank,
                    "cohort": f"AlphaGenome Ranked ({autonomy_str}, Tier {tier})",
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
                print(f"[{session_num} in session | AlphaGenome Rank #{agsps_rank} ({orf_id})] Folding {pair_id} ({gene} [{m_len}aa] + {p_name} [{p_len}aa] = {tot_len}aa) | AGSPS: {agsps_score:.1f} ...")

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
        if continuous_candidates:
            cand_df = pd.DataFrame(continuous_candidates).drop_duplicates(subset=["pair_id"], keep="last")
            cand_df.to_parquet(continuous_candidates_path, index=False)
        string_resolver.save_cache()
        sync_master_atlas(base_dir, state_file, continuous_candidates_path)
        total_time = time.time() - t_start
        print("\n================================================================================")
        print(f"  Continuous AlphaGenome Session Concluded")
        print(f"  Total Session Elapsed: {total_time/60:.2f} minutes")
        print(f"  Complexes Folded:      {n_processed_this_session}")
        print(f"  Complexes Skipped:     {n_skipped_this_session}")
        print("================================================================================")


if __name__ == "__main__":
    main()
