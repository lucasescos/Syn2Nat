#!/usr/bin/env python3
"""
Automated Protein-Protein Interface Site & Binding Epitope Annotation Pipeline

Evaluates WHERE microproteins join canonical binding partners:
1. Residue-level atomic contact extraction and chemical classification.
2. Buried Surface Area (Delta-SASA in A^2) and solvent exposure via Biopython ShrakeRupley.
3. PRODIGY contact chemistry profiling (IC distribution, link density, predicted dG).
4. UniProt catalytic/active site and domain proximity mapping (ACT_SITE, BINDING, DOMAIN).
5. EMBL-EBI PDBe-KB experimental macromolecular interface overlap scoring.
6. Automated molecular mechanism of action classification.
7. High-throughput checkpointing and multi-threading across the entire atlas.
"""

import argparse
import concurrent.futures
import json
import math
import os
import sys
import threading
import time
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple

import numpy as np
import pandas as pd
import requests

from Bio.PDB import PDBParser
from Bio.PDB.SASA import ShrakeRupley

# Project Directories
BASE_DIR = Path(__file__).resolve().parent.parent
STRUCTURES_DIR = BASE_DIR / "structures"
PDBS_DIR = STRUCTURES_DIR / "pdbs"
PAES_DIR = STRUCTURES_DIR / "pae"
CACHE_DIR = BASE_DIR / "data" / "cache"
UNIPROT_CACHE_DIR = CACHE_DIR / "uniprot"
PDBEKB_CACHE_DIR = CACHE_DIR / "pdbekb"
REPORTS_DIR = BASE_DIR / "reports"
INTERFACES_DIR = REPORTS_DIR / "interfaces"

# Ensure directories exist
UNIPROT_CACHE_DIR.mkdir(parents=True, exist_ok=True)
PDBEKB_CACHE_DIR.mkdir(parents=True, exist_ok=True)
INTERFACES_DIR.mkdir(parents=True, exist_ok=True)

# Residue Chemistry Definitions
ACIDIC = {"ASP", "GLU"}
BASIC = {"ARG", "LYS", "HIS"}
HYDROPHOBIC = {"ALA", "VAL", "LEU", "ILE", "MET", "PHE", "TRP", "PRO"}
POLAR = {"SER", "THR", "CYS", "ASN", "GLN", "TYR"}
AROMATIC = {"PHE", "TYR", "TRP", "HIS"}

# PRODIGY classification mapping (A: Apolar, C: Charged, P: Polar)
PRODIGY_CLASS_MAP = {
    "ALA": "A", "VAL": "A", "LEU": "A", "ILE": "A", "MET": "A", "PHE": "A", "PRO": "A", "TRP": "A", "CYS": "A", "GLY": "A",
    "GLU": "C", "ASP": "C", "ARG": "C", "LYS": "C", "HIS": "C",
    "SER": "P", "THR": "P", "ASN": "P", "GLN": "P", "TYR": "P"
}

# Thread lock for file writes
CHECKPOINT_LOCK = threading.Lock()
PRINT_LOCK = threading.Lock()


# ==============================================================================
# 1. Coordinate Parsing & Contact Geometry
# ==============================================================================

def parse_pdb_coordinates(pdb_path) -> Dict[str, Dict[int, dict]]:
    """Parse ATOM records grouped by chain (A: microprotein, B: partner) and residue number."""
    pdb_path = Path(pdb_path)
    chains: Dict[str, Dict[int, dict]] = {"A": {}, "B": {}}
    if not pdb_path.is_file():
        raise FileNotFoundError(f"PDB file not found: {pdb_path}")

    with open(pdb_path, "r", encoding="utf-8") as f:
        for line in f:
            if line.startswith(("ATOM", "HETATM")):
                atom_name = line[12:16].strip()
                res_name = line[17:20].strip()
                chain_id = line[21].strip()
                res_num = int(line[22:26].strip())
                x = float(line[30:38].strip())
                y = float(line[38:46].strip())
                z = float(line[46:54].strip())
                bfactor = float(line[60:66].strip())

                if chain_id in chains:
                    if res_num not in chains[chain_id]:
                        chains[chain_id][res_num] = {
                            "res_num": res_num,
                            "res_name": res_name,
                            "atoms": [],
                            "coords": [],
                            "bfactors": []
                        }
                    chains[chain_id][res_num]["atoms"].append(atom_name)
                    chains[chain_id][res_num]["coords"].append(np.array([x, y, z]))
                    chains[chain_id][res_num]["bfactors"].append(bfactor)
    return chains


def classify_contact_type(res_a_name: str, atom_a: str, res_b_name: str, atom_b: str, dist: float) -> str:
    """Classify physical nature of atomic contact."""
    elem_a = atom_a[0]
    elem_b = atom_b[0]

    # Salt bridge check (<= 4.0 A)
    if (res_a_name in ACIDIC and res_b_name in BASIC) or (res_a_name in BASIC and res_b_name in ACIDIC):
        acidic_atoms = {"OD1", "OD2", "OE1", "OE2", "OXT"}
        basic_atoms = {"NZ", "NH1", "NH2", "NE", "ND1", "NE2"}
        if ((atom_a in acidic_atoms and atom_b in basic_atoms) or (atom_a in basic_atoms and atom_b in acidic_atoms)) and dist <= 4.0:
            return "Salt Bridge"

    # Hydrogen bond check (<= 3.5 A)
    if ((elem_a == "N" and elem_b == "O") or (elem_a == "O" and elem_b == "N") or (elem_a == "O" and elem_b == "O")) and dist <= 3.5:
        return "Hydrogen Bond"

    # Hydrophobic packing (<= 4.5 A)
    if res_a_name in HYDROPHOBIC and res_b_name in HYDROPHOBIC:
        sidechain_a = atom_a not in {"N", "CA", "C", "O"}
        sidechain_b = atom_b not in {"N", "CA", "C", "O"}
        if sidechain_a and sidechain_b and dist <= 4.5:
            return "Hydrophobic Packing"

    # Aromatic contact (<= 4.8 A)
    if res_a_name in AROMATIC and res_b_name in AROMATIC:
        if atom_a not in {"N", "CA", "C", "O"} and atom_b not in {"N", "CA", "C", "O"} and dist <= 4.8:
            return "Aromatic Contact"

    # Van der Waals contact (<= 4.0 A)
    if dist <= 4.0:
        return "Van der Waals"
    elif dist <= 4.5:
        return "Peripheral Contact"
    return "Proximal"


def extract_contacts(chains: dict, pae_path: Optional[Path] = None, dist_cutoff: float = 4.5) -> Tuple[pd.DataFrame, dict]:
    """Calculate atomic contacts between Chain A and Chain B."""
    res_a = chains["A"]
    res_b = chains["B"]
    len_a = len(res_a)
    len_b = len(res_b)

    # Load PAE matrix if available
    pae_matrix = None
    if pae_path and pae_path.is_file():
        try:
            with open(pae_path, "r", encoding="utf-8") as f:
                pae_matrix = np.array(json.load(f))
        except Exception:
            pass

    contact_rows = []
    interacting_a = set()
    interacting_b = set()

    # PRODIGY IC bins
    ic_bins = {"CC": 0, "CP": 0, "AC": 0, "PP": 0, "AP": 0, "AA": 0}

    for r_a_num, data_a in res_a.items():
        coords_a = np.array(data_a["coords"])
        atoms_a = data_a["atoms"]
        name_a = data_a["res_name"]
        plddt_a = float(np.mean(data_a["bfactors"]))

        for r_b_num, data_b in res_b.items():
            coords_b = np.array(data_b["coords"])
            atoms_b = data_b["atoms"]
            name_b = data_b["res_name"]
            plddt_b = float(np.mean(data_b["bfactors"]))

            dists = np.linalg.norm(coords_a[:, np.newaxis, :] - coords_b[np.newaxis, :, :], axis=2)
            min_idx = np.unravel_index(np.argmin(dists), dists.shape)
            min_dist = float(dists[min_idx])

            if min_dist <= dist_cutoff:
                atom_a_best = atoms_a[min_idx[0]]
                atom_b_best = atoms_b[min_idx[1]]

                # PAE lookup
                pae_val = None
                if pae_matrix is not None:
                    idx_a = r_a_num - 1
                    idx_b = len_a + (r_b_num - 1)
                    if idx_a < pae_matrix.shape[0] and idx_b < pae_matrix.shape[1]:
                        pae_val = round(float((pae_matrix[idx_a, idx_b] + pae_matrix[idx_b, idx_a]) / 2.0), 2)

                itype = classify_contact_type(name_a, atom_a_best, name_b, atom_b_best, min_dist)
                interacting_a.add(r_a_num)
                interacting_b.add(r_b_num)

                # PRODIGY contact bin
                c_a = PRODIGY_CLASS_MAP.get(name_a, "A")
                c_b = PRODIGY_CLASS_MAP.get(name_b, "A")
                bin_key = "".join(sorted([c_a, c_b]))
                if bin_key in ic_bins:
                    ic_bins[bin_key] += 1

                contact_rows.append({
                    "uProt_res_num": r_a_num,
                    "uProt_res_name": name_a,
                    "uProt_atom": atom_a_best,
                    "uProt_plddt": round(plddt_a, 1),
                    "partner_res_num": r_b_num,
                    "partner_res_name": name_b,
                    "partner_atom": atom_b_best,
                    "partner_plddt": round(plddt_b, 1),
                    "min_distance_A": round(min_dist, 2),
                    "interaction_type": itype,
                    "pae_A": pae_val
                })

    contacts_df = pd.DataFrame(contact_rows).sort_values("min_distance_A").reset_index(drop=True) if contact_rows else pd.DataFrame()

    # Metrics
    uprot_contact_pct = (len(interacting_a) / len_a) * 100 if len_a > 0 else 0.0
    mean_pae = float(contacts_df["pae_A"].dropna().mean()) if not contacts_df.empty and "pae_A" in contacts_df and contacts_df["pae_A"].notnull().sum() > 0 else None
    
    # Link density
    max_possible_contacts = max(1, len(interacting_a) * len(interacting_b))
    link_density = len(contacts_df) / max_possible_contacts if not contacts_df.empty else 0.0

    # Estimated binding free energy dG
    est_dg = (0.368 * ic_bins["CC"] + 0.477 * ic_bins["AC"] - 0.906 * ic_bins["PP"] - 0.501 * ic_bins["AP"] - 10.0) if len(contacts_df) > 0 else 0.0

    stats = {
        "microprotein_len": len_a,
        "partner_len": len_b,
        "total_contacts": len(contacts_df),
        "microprotein_contact_residues": sorted(list(interacting_a)),
        "partner_contact_residues": sorted(list(interacting_b)),
        "microprotein_contact_count": len(interacting_a),
        "partner_contact_count": len(interacting_b),
        "microprotein_surface_involvement_pct": round(uprot_contact_pct, 1),
        "mean_interface_pae": round(mean_pae, 2) if mean_pae is not None else None,
        "hbond_count": int((contacts_df["interaction_type"] == "Hydrogen Bond").sum()) if not contacts_df.empty else 0,
        "salt_bridge_count": int((contacts_df["interaction_type"] == "Salt Bridge").sum()) if not contacts_df.empty else 0,
        "hydrophobic_count": int((contacts_df["interaction_type"] == "Hydrophobic Packing").sum()) if not contacts_df.empty else 0,
        "aromatic_count": int((contacts_df["interaction_type"] == "Aromatic Contact").sum()) if not contacts_df.empty else 0,
        "vdw_count": int((contacts_df["interaction_type"] == "Van der Waals").sum()) if not contacts_df.empty else 0,
        "prodigy_bins": ic_bins,
        "link_density": round(link_density, 3),
        "predicted_dg_kcal": round(est_dg, 2)
    }
    return contacts_df, stats


# ==============================================================================
# 2. Solvent Accessibility & Buried Surface Area (Delta-SASA)
# ==============================================================================

def compute_buried_sasa(pdb_path) -> dict:
    """Compute Lee-Richards / ShrakeRupley SASA for complex, chain A, chain B, and Buried SASA."""
    pdb_path = Path(pdb_path)
    parser = PDBParser(QUIET=True)
    struct = parser.get_structure("complex", str(pdb_path))
    model = struct[0]

    sr = ShrakeRupley()
    # 1. Complex SASA
    sr.compute(model, level="C")
    sasa_complex = sum(c.sasa for c in model)

    # 2. Isolated Chain A & Chain B SASA
    ch_a = model["A"]
    ch_b = model["B"]
    sr.compute(ch_a, level="C")
    sr.compute(ch_b, level="C")
    sasa_a = ch_a.sasa
    sasa_b = ch_b.sasa

    # Buried Surface Area (Delta-SASA)
    delta_sasa = max(0.0, (sasa_a + sasa_b) - sasa_complex)

    return {
        "sasa_microprotein_A2": round(sasa_a, 1),
        "sasa_partner_A2": round(sasa_b, 1),
        "sasa_complex_A2": round(sasa_complex, 1),
        "buried_sasa_A2": round(delta_sasa, 1)
    }


# ==============================================================================
# 3. UniProt REST API: Active Site & Domain Proximity
# ==============================================================================

def fetch_uniprot_features(uniprot_id: str) -> Optional[dict]:
    """Fetch structured UniProt features with disk caching."""
    if not uniprot_id or uniprot_id.upper() in ["NONE", "NAN", ""]:
        return None

    cache_file = UNIPROT_CACHE_DIR / f"{uniprot_id}.json"
    if cache_file.is_file():
        try:
            with open(cache_file, "r", encoding="utf-8") as f:
                return json.load(f)
        except Exception:
            pass

    url = f"https://rest.uniprot.org/uniprotkb/{uniprot_id}.json"
    try:
        resp = requests.get(url, timeout=15)
        if resp.status_code == 200:
            data = resp.json()
            with open(cache_file, "w", encoding="utf-8") as f:
                json.dump(data, f)
            return data
    except Exception as e:
        pass
    return None


def evaluate_uniprot_proximity(uniprot_data: Optional[dict], chains: dict) -> dict:
    """Evaluate 3D proximity to active sites, binding sites, domains, and PTMs."""
    if not uniprot_data or "features" not in uniprot_data:
        return {
            "has_uniprot_annotation": False,
            "catalytic_sites_count": 0,
            "active_site_min_dist_A": None,
            "is_active_site_cleft": False,
            "active_site_descriptions": [],
            "binding_sites_count": 0,
            "binding_site_min_dist_A": None,
            "is_binding_site_overlap": False,
            "binding_site_descriptions": [],
            "uniprot_domains_raw": []
        }

    feats = uniprot_data.get("features", [])
    act_sites = []
    binding_sites = []
    domains = []
    ptms = []

    for f in feats:
        ftype = f.get("type")
        desc = f.get("description", "")
        loc = f.get("location", {})
        start = loc.get("start", {}).get("value")
        end = loc.get("end", {}).get("value")

        if start is None:
            continue
        end = end if end is not None else start

        if ftype == "Active site":
            act_sites.append((start, end, desc))
        elif ftype == "Binding site":
            binding_sites.append((start, end, desc))
        elif ftype in ["Domain", "Region"]:
            domains.append((start, end, desc))
        elif ftype == "Modified residue":
            ptms.append((start, end, desc))

    res_a = chains["A"]
    res_b = chains["B"]
    
    all_coords_a = np.concatenate([data["coords"] for data in res_a.values()]) if res_a else np.empty((0, 3))

    def min_dist_to_feature_residues(site_list):
        min_d = float("inf")
        hit_descs = []
        for start, end, desc in site_list:
            for pos in range(start, end + 1):
                if pos in res_b:
                    b_coords = np.array(res_b[pos]["coords"])
                    dists = np.linalg.norm(all_coords_a[:, np.newaxis, :] - b_coords[np.newaxis, :, :], axis=2)
                    d_cur = float(np.min(dists))
                    if d_cur < min_d:
                        min_d = d_cur
                    if d_cur <= 5.0:
                        hit_descs.append(f"{desc} (res {pos})")
        return (round(min_d, 2) if min_d != float("inf") else None), hit_descs

    act_min_dist, act_hit_descs = min_dist_to_feature_residues(act_sites)
    bind_min_dist, bind_hit_descs = min_dist_to_feature_residues(binding_sites)

    return {
        "has_uniprot_annotation": True,
        "catalytic_sites_count": len(act_sites),
        "active_site_min_dist_A": act_min_dist,
        "is_active_site_cleft": bool(act_min_dist is not None and act_min_dist <= 5.0),
        "active_site_descriptions": act_hit_descs,
        "binding_sites_count": len(binding_sites),
        "binding_site_min_dist_A": bind_min_dist,
        "is_binding_site_overlap": bool(bind_min_dist is not None and bind_min_dist <= 5.0),
        "binding_site_descriptions": bind_hit_descs,
        "uniprot_domains_raw": domains
    }


# ==============================================================================
# 4. EMBL-EBI PDBe-KB: Known Experimental PDB Interface Residues
# ==============================================================================

def fetch_pdbekb_interface_residues(uniprot_id: str) -> Optional[dict]:
    """Fetch aggregated experimental PDB interface residues from PDBe-KB Graph API."""
    if not uniprot_id or uniprot_id.upper() in ["NONE", "NAN", ""]:
        return None

    cache_file = PDBEKB_CACHE_DIR / f"{uniprot_id}.json"
    if cache_file.is_file():
        try:
            with open(cache_file, "r", encoding="utf-8") as f:
                return json.load(f)
        except Exception:
            pass

    url = f"https://www.ebi.ac.uk/pdbe/graph-api/uniprot/interface_residues/{uniprot_id}"
    try:
        resp = requests.get(url, timeout=15)
        if resp.status_code == 200:
            data = resp.json()
            with open(cache_file, "w", encoding="utf-8") as f:
                json.dump(data, f)
            return data
    except Exception as e:
        pass
    return None


def evaluate_pdbekb_overlap(pdbekb_data: Optional[dict], partner_contact_residues: List[int], uniprot_id: str) -> dict:
    """Evaluate overlap percentage between predicted partner contact residues and known PDB interfaces."""
    if not pdbekb_data or uniprot_id not in pdbekb_data:
        return {
            "has_pdbekb_data": False,
            "known_pdb_interface_residues_count": 0,
            "known_pdb_ppi_overlap_count": 0,
            "known_pdb_ppi_overlap_pct": 0.0,
            "representative_pdb_ids": []
        }

    items = pdbekb_data[uniprot_id].get("data", [])
    known_interface_positions = set()
    pdb_id_matches = set()

    for entry in items:
        for res_item in entry.get("residues", []):
            start_idx = res_item.get("startIndex")
            end_idx = res_item.get("endIndex", start_idx)
            if start_idx is not None:
                for p in range(start_idx, end_idx + 1):
                    known_interface_positions.add(p)
                    if p in partner_contact_residues:
                        for pdb_entry in res_item.get("interactingPDBEntries", []):
                            pid = pdb_entry.get("pdbId")
                            if pid:
                                pdb_id_matches.add(pid.upper())

    partner_contacts = set(partner_contact_residues)
    overlap = partner_contacts.intersection(known_interface_positions)
    overlap_pct = (len(overlap) / len(partner_contacts)) * 100.0 if partner_contacts else 0.0

    return {
        "has_pdbekb_data": True,
        "known_pdb_interface_residues_count": len(known_interface_positions),
        "known_pdb_ppi_overlap_count": len(overlap),
        "known_pdb_ppi_overlap_pct": round(overlap_pct, 1),
        "representative_pdb_ids": sorted(list(pdb_id_matches))[:10]
    }


# ==============================================================================
# 5. Unified Molecular Mechanism Classifier
# ==============================================================================

def classify_binding_mechanism(stats: dict, uniprot_res: dict, pdbekb_res: dict) -> str:
    """Classify the likely molecular mechanism of the microprotein based on its binding epitope."""
    # 0. Check for non-interacting complex
    if stats.get("total_contacts", 0) == 0:
        return "Non-Interacting / Dissociated"

    # 1. Active site cleft direct contact (<= 5.0 A)
    act_dist = uniprot_res.get("active_site_min_dist_A")
    bind_dist = uniprot_res.get("binding_site_min_dist_A")

    if uniprot_res.get("is_active_site_cleft") or (act_dist is not None and act_dist <= 5.0):
        return "Active Site Blocker / Competitive Decoy"
    
    if uniprot_res.get("is_binding_site_overlap") or (bind_dist is not None and bind_dist <= 5.0):
        return "Cofactor / Substrate Pocket Competitor"

    # 2. Flanking the active site cleft (<= 8.0 A)
    if (act_dist is not None and act_dist <= 8.0) or (bind_dist is not None and bind_dist <= 8.0):
        return "Active Site Cleft Occlusion / Proximal Inhibitor"
    
    # 3. High overlap with known physiological PDB interface
    overlap_pct = pdbekb_res.get("known_pdb_ppi_overlap_pct", 0.0)
    if overlap_pct >= 60.0:
        return "Known Macromolecular Subunit Competitor / Adaptor"
    elif overlap_pct >= 30.0:
        return "Shared Physiological Complex Interface"
    
    # 4. Low overlap (< 20%) but robust buried SASA
    if stats.get("buried_sasa_A2", 0) >= 800.0:
        return "Allosteric Regulatory Surface / Novel Epitope"
    
    # 5. Low buried SASA or peripheral
    return "Peripheral Surface Interaction"


# ==============================================================================
# 6. Pipeline Orchestrator for a Single Complex
# ==============================================================================

def annotate_complex(pair_id: str, pdb_path: Path, partner_uniprot_id: str, pae_path: Optional[Path] = None) -> dict:
    """Run full end-to-end interface annotation for a single complex."""
    pdb_path = Path(pdb_path)
    if not pdb_path.is_file():
        raise FileNotFoundError(f"PDB does not exist: {pdb_path}")

    # 1. Coordinates and contacts
    chains = parse_pdb_coordinates(pdb_path)
    contacts_df, contact_stats = extract_contacts(chains, pae_path=pae_path, dist_cutoff=4.5)

    # Fast path for non-interacting complexes
    if contact_stats["total_contacts"] == 0:
        return {
            "pair_id": pair_id,
            "partner_uniprot_id": partner_uniprot_id,
            "total_contacts": 0,
            "microprotein_contact_residues_count": 0,
            "partner_contact_residues_count": 0,
            "microprotein_surface_involvement_pct": 0.0,
            "mean_interface_pae": None,
            "hbond_count": 0,
            "salt_bridge_count": 0,
            "hydrophobic_count": 0,
            "aromatic_count": 0,
            "vdw_count": 0,
            "prodigy_link_density": 0.0,
            "prodigy_predicted_dg_kcal": 0.0,
            "buried_sasa_A2": 0.0,
            "sasa_microprotein_A2": None,
            "sasa_partner_A2": None,
            "catalytic_sites_count": 0,
            "active_site_min_dist_A": None,
            "is_active_site_cleft": False,
            "active_site_details": "None",
            "binding_sites_count": 0,
            "binding_site_min_dist_A": None,
            "is_binding_site_overlap": False,
            "binding_site_details": "None",
            "partner_domains_contacted": "N/A",
            "known_pdb_interface_residues_count": 0,
            "known_pdb_ppi_overlap_count": 0,
            "known_pdb_ppi_overlap_pct": 0.0,
            "representative_pdb_ids": "None",
            "binding_mechanism_category": "Non-Interacting / Dissociated"
        }

    # 2. Buried Surface Area (SASA)
    sasa_stats = compute_buried_sasa(pdb_path)

    # 3. UniProt functional features
    uniprot_data = fetch_uniprot_features(partner_uniprot_id)
    uniprot_stats = evaluate_uniprot_proximity(uniprot_data, chains)

    # Domain mapping
    domains_contacted = []
    if uniprot_stats.get("has_uniprot_annotation") and "uniprot_domains_raw" in uniprot_stats:
        contact_b = set(contact_stats["partner_contact_residues"])
        for start, end, desc in uniprot_stats["uniprot_domains_raw"]:
            dom_span = set(range(start, end + 1))
            if contact_b.intersection(dom_span):
                domains_contacted.append(f"{desc} ({start}-{end})")
    domains_contacted_str = "; ".join(domains_contacted) if domains_contacted else "N/A"

    # 4. PDBe-KB known experimental PDB interface overlap
    pdbekb_data = fetch_pdbekb_interface_residues(partner_uniprot_id)
    pdbekb_stats = evaluate_pdbekb_overlap(pdbekb_data, contact_stats["partner_contact_residues"], partner_uniprot_id)

    # 5. Mechanism classification
    mechanism = classify_binding_mechanism(
        {**contact_stats, **sasa_stats},
        uniprot_stats,
        pdbekb_stats
    )

    # Save detailed contact TSV and PyMOL visualization script
    detail_tsv = INTERFACES_DIR / f"{pair_id}_contacts.tsv"
    if not contacts_df.empty:
        contacts_df.to_csv(detail_tsv, sep="\t", index=False)

    pml_file = INTERFACES_DIR / f"{pair_id}_interface.pml"
    res_a_str = "+".join(str(r) for r in contact_stats["microprotein_contact_residues"])
    res_b_str = "+".join(str(r) for r in contact_stats["partner_contact_residues"])
    pml_content = f"""# PyMOL Interface Visualization Script for {pair_id}
load {pdb_path.resolve()}, complex
hide everything, complex
show cartoon, complex
color marine, chain B
color forest, chain A

# Interface Residues
select if_uprot, chain A and resi {res_a_str}
select if_partner, chain B and resi {res_b_str}

show sticks, if_uprot
show sticks, if_partner
color tv_yellow, if_uprot
color tv_orange, if_partner

# Hydrogen bonds
distance hbonds, chain A, chain B, 3.5, mode=2
color cyan, hbonds

set cartoon_transparency, 0.2
zoom if_uprot or if_partner
"""
    pml_file.write_text(pml_content, encoding="utf-8")

    return {
        "pair_id": pair_id,
        "partner_uniprot_id": partner_uniprot_id,
        "total_contacts": contact_stats["total_contacts"],
        "microprotein_contact_residues_count": contact_stats["microprotein_contact_count"],
        "partner_contact_residues_count": contact_stats["partner_contact_count"],
        "microprotein_surface_involvement_pct": contact_stats["microprotein_surface_involvement_pct"],
        "mean_interface_pae": contact_stats["mean_interface_pae"],
        "hbond_count": contact_stats["hbond_count"],
        "salt_bridge_count": contact_stats["salt_bridge_count"],
        "hydrophobic_count": contact_stats["hydrophobic_count"],
        "aromatic_count": contact_stats["aromatic_count"],
        "vdw_count": contact_stats["vdw_count"],
        "prodigy_link_density": contact_stats["link_density"],
        "prodigy_predicted_dg_kcal": contact_stats["predicted_dg_kcal"],
        "buried_sasa_A2": sasa_stats["buried_sasa_A2"],
        "sasa_microprotein_A2": sasa_stats["sasa_microprotein_A2"],
        "sasa_partner_A2": sasa_stats["sasa_partner_A2"],
        "catalytic_sites_count": uniprot_stats["catalytic_sites_count"],
        "active_site_min_dist_A": uniprot_stats["active_site_min_dist_A"],
        "is_active_site_cleft": uniprot_stats["is_active_site_cleft"],
        "active_site_details": "; ".join(uniprot_stats["active_site_descriptions"]) if uniprot_stats["active_site_descriptions"] else "None",
        "binding_sites_count": uniprot_stats["binding_sites_count"],
        "binding_site_min_dist_A": uniprot_stats["binding_site_min_dist_A"],
        "is_binding_site_overlap": uniprot_stats["is_binding_site_overlap"],
        "binding_site_details": "; ".join(uniprot_stats["binding_site_descriptions"]) if uniprot_stats["binding_site_descriptions"] else "None",
        "partner_domains_contacted": domains_contacted_str,
        "known_pdb_interface_residues_count": pdbekb_stats["known_pdb_interface_residues_count"],
        "known_pdb_ppi_overlap_count": pdbekb_stats["known_pdb_ppi_overlap_count"],
        "known_pdb_ppi_overlap_pct": pdbekb_stats["known_pdb_ppi_overlap_pct"],
        "representative_pdb_ids": "; ".join(pdbekb_stats["representative_pdb_ids"]) if pdbekb_stats["representative_pdb_ids"] else "None",
        "binding_mechanism_category": mechanism
    }


# ==============================================================================
# 7. Batch Driver with Multithreading & Resilient Checkpointing
# ==============================================================================

def process_single_row(row: dict, checkpoint_file: Path, counter: dict, total_targets: int) -> Optional[dict]:
    pair_id = str(row["pair_id"])
    uprot = str(row["microprotein_gene"])
    partner = str(row["partner_gene_symbol"])
    uniprot_id = str(row["partner_uniprot_id"])
    iptm_val = float(row["iptm"]) if row.get("iptm") is not None else 0.0

    pdb_file = PDBS_DIR / f"{pair_id}.pdb"
    pae_file = PAES_DIR / f"{pair_id}_pae.json"

    if not pdb_file.is_file():
        return None

    try:
        res = annotate_complex(pair_id, pdb_file, uniprot_id, pae_path=pae_file)
        res["microprotein_gene"] = uprot
        res["microprotein_biotype"] = row.get("microprotein_biotype", "")
        res["microprotein_tier"] = row.get("microprotein_tier", "")
        res["partner_gene_symbol"] = partner
        res["iptm"] = iptm_val
        res["ptm"] = float(row["ptm"]) if row.get("ptm") is not None else 0.0
        res["partner_plddt"] = float(row["partner_plddt"]) if row.get("partner_plddt") is not None else 0.0
        res["microprotein_plddt"] = float(row["microprotein_plddt"]) if row.get("microprotein_plddt") is not None else 0.0

        # Append to checkpoint
        with CHECKPOINT_LOCK:
            with open(checkpoint_file, "a", encoding="utf-8") as f:
                f.write(json.dumps(res) + "\n")
            counter["done"] += 1
            cur_done = counter["done"]

        # Progress reporting
        if cur_done % 25 == 0 or cur_done == total_targets or res["total_contacts"] > 0:
            with PRINT_LOCK:
                if res["total_contacts"] > 0:
                    print(f"[{cur_done}/{total_targets}] {uprot} + {partner} ({pair_id}) | ipTM: {iptm_val:.3f} | Contacts: {res['total_contacts']} | Delta-SASA: {res['buried_sasa_A2']} A^2 | Overlap: {res['known_pdb_ppi_overlap_pct']}% | {res['binding_mechanism_category']}")
                elif cur_done % 100 == 0:
                    print(f"[{cur_done}/{total_targets}] ... {cur_done} complexes processed.")

        return res
    except Exception as e:
        with PRINT_LOCK:
            print(f"  [ERROR] {pair_id} failed: {e}")
        return None


def main():
    parser = argparse.ArgumentParser(description="Automated Protein Interface Site & Epitope Annotation Pipeline")
    parser.add_argument("--pair-id", default=None, help="Process a single pair_id")
    parser.add_argument("--min-iptm", type=float, default=None, help="Process all complexes with ipTM >= threshold")
    parser.add_argument("--limit", type=int, default=None, help="Limit number of complexes to process")
    parser.add_argument("--workers", type=int, default=4, help="Number of worker threads (default: 4)")
    parser.add_argument("--checkpoint", default=str(REPORTS_DIR / "interface_site_annotations_checkpoint.jsonl"), help="Checkpoint JSONL path")
    parser.add_argument("--out-tsv", default=str(REPORTS_DIR / "interface_site_annotations_summary.tsv"), help="Path to output TSV")
    parser.add_argument("--out-parquet", default=str(REPORTS_DIR / "interface_site_annotations_summary.parquet"), help="Path to output Parquet")
    args = parser.parse_args()

    atlas_path = REPORTS_DIR / "master_microprotein_ppi_structural_atlas.parquet"
    if not atlas_path.is_file():
        atlas_path = REPORTS_DIR / "master_microprotein_ppi_structural_atlas.tsv"
        if not atlas_path.is_file():
            print(f"Error: Master atlas not found at {atlas_path}")
            sys.exit(1)
        atlas_df = pd.read_csv(atlas_path, sep="\t")
    else:
        atlas_df = pd.read_parquet(atlas_path)

    # Filter targets
    if args.pair_id:
        targets = atlas_df[atlas_df["pair_id"] == args.pair_id].copy()
    elif args.min_iptm is not None:
        targets = atlas_df[atlas_df["iptm"] >= args.min_iptm].copy()
    else:
        # Default: Process ALL completed complexes!
        targets = atlas_df.copy()

    if args.limit:
        targets = targets.head(args.limit)

    targets = targets.sort_values("iptm", ascending=False).reset_index(drop=True)
    checkpoint_file = Path(args.checkpoint)

    # Load existing checkpoint
    completed_pairs = set()
    checkpoint_records = []
    if checkpoint_file.is_file():
        with open(checkpoint_file, "r", encoding="utf-8") as f:
            for line in f:
                line = line.strip()
                if line:
                    try:
                        record = json.loads(line)
                        completed_pairs.add(record["pair_id"])
                        checkpoint_records.append(record)
                    except Exception:
                        pass
        print(f"[CHECKPOINT] Loaded {len(completed_pairs)} already annotated complexes from {checkpoint_file.name}.")

    # Exclude already processed targets
    pending_targets = targets[~targets["pair_id"].isin(completed_pairs)].copy()

    print(f"================================================================================")
    print(f"  PROTEIN INTERFACE SITE & BINDING EPITOPE ANNOTATION PIPELINE")
    print(f"  Total Targets in Scope: {len(targets)} | Already Done: {len(completed_pairs)} | Pending: {len(pending_targets)}")
    print(f"  Worker Threads: {args.workers}")
    print(f"================================================================================")

    t_start = time.time()
    counter = {"done": len(completed_pairs)}

    if len(pending_targets) > 0:
        rows_to_process = pending_targets.to_dict("records")
        with concurrent.futures.ThreadPoolExecutor(max_workers=args.workers) as executor:
            futures = [
                executor.submit(process_single_row, r, checkpoint_file, counter, len(targets))
                for r in rows_to_process
            ]
            for f in concurrent.futures.as_completed(futures):
                pass

    # Read back all records from checkpoint file to build definitive summary
    final_records = []
    if checkpoint_file.is_file():
        with open(checkpoint_file, "r", encoding="utf-8") as f:
            for line in f:
                line = line.strip()
                if line:
                    try:
                        final_records.append(json.loads(line))
                    except Exception:
                        pass

    summary_df = pd.DataFrame(final_records)
    if not summary_df.empty:
        # Drop duplicates by pair_id if any
        summary_df = summary_df.drop_duplicates(subset=["pair_id"]).reset_index(drop=True)

        lead_cols = [
            "pair_id", "microprotein_gene", "microprotein_biotype", "microprotein_tier",
            "partner_gene_symbol", "partner_uniprot_id", "iptm", "ptm", "binding_mechanism_category",
            "buried_sasa_A2", "known_pdb_ppi_overlap_pct", "is_active_site_cleft", "active_site_min_dist_A",
            "active_site_details", "is_binding_site_overlap", "binding_site_min_dist_A",
            "partner_domains_contacted", "representative_pdb_ids", "total_contacts", "hbond_count",
            "salt_bridge_count", "hydrophobic_count", "prodigy_link_density", "prodigy_predicted_dg_kcal"
        ]
        present_leads = [c for c in lead_cols if c in summary_df.columns]
        other_cols = [c for c in summary_df.columns if c not in present_leads]
        summary_df = summary_df[present_leads + other_cols]

        out_tsv = Path(args.out_tsv)
        out_parquet = Path(args.out_parquet)
        summary_df.to_csv(out_tsv, sep="\t", index=False)
        summary_df.to_parquet(out_parquet, index=False)
        print(f"\n[DONE] Successfully consolidated {len(summary_df)} complexes in {time.time() - t_start:.1f}s.")
        print(f"TSV Saved:     {out_tsv.resolve()}")
        print(f"Parquet Saved: {out_parquet.resolve()}")
    else:
        print("\n[!] No results found.")


if __name__ == "__main__":
    main()
