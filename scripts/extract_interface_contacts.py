#!/usr/bin/env python3
"""
Residue-Level Protein-Protein Interface Contact Extractor

Analyzes 3D coordinate PDBs and PAE matrices from ESMFold2 cofolding:
1. Calculates all inter-chain atomic distances between Chain A (microprotein) and Chain B (partner).
2. Identifies specific residue-residue interactions:
   - Hydrogen bonds (donor/acceptor N-O <= 3.5 A)
   - Salt bridges (acidic Asp/Glu to basic Arg/Lys/His <= 4.0 A)
   - Hydrophobic packing (aliphatic/aromatic sidechains <= 4.5 A)
   - Van der Waals contacts (heavy atom distance <= 4.0 A)
3. Cross-references 3D contacts with inter-chain PAE confidence matrices.
4. Computes interface metrics (buried residues, interface footprint %, interface pLDDT).
5. Generates PyMOL (.pml) visualization script to inspect interface residues in 3D.
"""

import argparse
import json
import math
from pathlib import Path
import numpy as np
import pandas as pd

# Residue chemistry classifications
ACIDIC = {"ASP", "GLU"}
BASIC = {"ARG", "LYS", "HIS"}
HYDROPHOBIC = {"ALA", "VAL", "LEU", "ILE", "MET", "PHE", "TRP", "PRO"}
POLAR = {"SER", "THR", "CYS", "ASN", "GLN", "TYR"}
AROMATIC = {"PHE", "TYR", "TRP", "HIS"}

def parse_pdb_atoms(pdb_path: Path):
    """Parse ATOM/HETATM records from PDB file grouped by chain and residue."""
    chains = {"A": {}, "B": {}}
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

def load_pae_matrix(pae_path: Path):
    """Load PAE matrix if available."""
    if pae_path and pae_path.is_file():
        try:
            with open(pae_path, "r", encoding="utf-8") as f:
                data = json.load(f)
                return np.array(data)
        except Exception:
            pass
    return None

def classify_interaction(res_a_name, atom_a, res_b_name, atom_b, dist):
    """Classify physical nature of atomic contact."""
    elem_a = atom_a[0]
    elem_b = atom_b[0]
    
    # Salt bridge check
    if (res_a_name in ACIDIC and res_b_name in BASIC) or (res_a_name in BASIC and res_b_name in ACIDIC):
        acidic_atoms = {"OD1", "OD2", "OE1", "OE2", "OXT"}
        basic_atoms = {"NZ", "NH1", "NH2", "NE", "ND1", "NE2"}
        if (atom_a in acidic_atoms and atom_b in basic_atoms) or (atom_a in basic_atoms and atom_b in acidic_atoms):
            if dist <= 4.0:
                return "Salt Bridge"
                
    # Hydrogen bond check (donor-acceptor N-O / O-N / O-O <= 3.5 A)
    if (elem_a == "N" and elem_b == "O") or (elem_a == "O" and elem_b == "N") or (elem_a == "O" and elem_b == "O"):
        if dist <= 3.5:
            return "Hydrogen Bond"
            
    # Hydrophobic packing
    if res_a_name in HYDROPHOBIC and res_b_name in HYDROPHOBIC:
        sidechain_atoms_a = atom_a not in {"N", "CA", "C", "O"}
        sidechain_atoms_b = atom_b not in {"N", "CA", "C", "O"}
        if sidechain_atoms_a and sidechain_atoms_b and dist <= 4.5:
            return "Hydrophobic Packing"
            
    # Aromatic stacking
    if res_a_name in AROMATIC and res_b_name in AROMATIC:
        if atom_a not in {"N", "CA", "C", "O"} and atom_b not in {"N", "CA", "C", "O"} and dist <= 4.8:
            return "Aromatic Contact"
            
    # Van der Waals contact
    if dist <= 4.0:
        return "Van der Waals"
    elif dist <= 5.0:
        return "Peripheral Contact"
    else:
        return "Proximal"

def extract_interface(pdb_path: Path, pae_path: Path = None, dist_cutoff: float = 4.5):
    """Extract all residue-residue contacts between Chain A and Chain B."""
    chains = parse_pdb_atoms(pdb_path)
    res_a = chains["A"]
    res_b = chains["B"]
    
    if not res_a or not res_b:
        return None
        
    pae_matrix = load_pae_matrix(pae_path)
    len_a = len(res_a)
    
    contact_list = []
    interacting_a = set()
    interacting_b = set()
    
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
            
            # Compute distance matrix between all atoms of residue A and residue B
            dists = np.linalg.norm(coords_a[:, np.newaxis, :] - coords_b[np.newaxis, :, :], axis=2)
            min_idx = np.unravel_index(np.argmin(dists), dists.shape)
            min_dist = float(dists[min_idx])
            
            if min_dist <= dist_cutoff:
                atom_a_best = atoms_a[min_idx[0]]
                atom_b_best = atoms_b[min_idx[1]]
                
                # Inter-chain PAE lookup
                pae_val = None
                if pae_matrix is not None:
                    idx_a = r_a_num - 1
                    idx_b = len_a + (r_b_num - 1)
                    if idx_a < pae_matrix.shape[0] and idx_b < pae_matrix.shape[1]:
                        pae_ab = float(pae_matrix[idx_a, idx_b])
                        pae_ba = float(pae_matrix[idx_b, idx_a])
                        pae_val = round((pae_ab + pae_ba) / 2.0, 2)
                
                itype = classify_interaction(name_a, atom_a_best, name_b, atom_b_best, min_dist)
                interacting_a.add(r_a_num)
                interacting_b.add(r_b_num)
                
                contact_list.append({
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
                
    contacts_df = pd.DataFrame(contact_list).sort_values("min_distance_A").reset_index(drop=True)
    
    # Interface metrics
    uprot_contact_pct = (len(interacting_a) / len_a) * 100 if len_a > 0 else 0
    mean_if_plddt_a = float(np.mean([np.mean(res_a[num]["bfactors"]) for num in interacting_a])) if interacting_a else 0
    mean_if_plddt_b = float(np.mean([np.mean(res_b[num]["bfactors"]) for num in interacting_b])) if interacting_b else 0
    
    summary = {
        "pdb_file": pdb_path.name,
        "total_contacts": len(contacts_df),
        "microprotein_total_len": len_a,
        "partner_total_len": len(res_b),
        "microprotein_interface_residues": sorted(list(interacting_a)),
        "partner_interface_residues": sorted(list(interacting_b)),
        "microprotein_interface_count": len(interacting_a),
        "partner_interface_count": len(interacting_b),
        "microprotein_surface_involvement_pct": round(uprot_contact_pct, 1),
        "interface_plddt_microprotein": round(mean_if_plddt_a, 1),
        "interface_plddt_partner": round(mean_if_plddt_b, 1),
        "hydrogen_bonds_count": int((contacts_df["interaction_type"] == "Hydrogen Bond").sum()) if len(contacts_df) else 0,
        "salt_bridges_count": int((contacts_df["interaction_type"] == "Salt Bridge").sum()) if len(contacts_df) else 0,
        "hydrophobic_contacts_count": int((contacts_df["interaction_type"] == "Hydrophobic Packing").sum()) if len(contacts_df) else 0
    }
    
    return contacts_df, summary

def generate_pymol_script(pdb_path: Path, summary: dict, out_pml: Path):
    """Generate PyMOL visualization script highlighting interface residues."""
    res_a_str = "+".join(str(r) for r in summary["microprotein_interface_residues"])
    res_b_str = "+".join(str(r) for r in summary["partner_interface_residues"])
    
    pml = f"""# PyMOL Interface Visualization Script for {pdb_path.name}
load {pdb_path.resolve()}, complex
hide everything, complex
show cartoon, complex
color marine, chain B
color forest, chain A

# Select and show interface residues
select if_uprot, chain A and resi {res_a_str}
select if_partner, chain B and resi {res_b_str}

show sticks, if_uprot
show sticks, if_partner
color tv_yellow, if_uprot
color tv_orange, if_partner

# Highlight polar contacts / H-bonds
distance hbonds, chain A, chain B, 3.5, mode=2
color cyan, hbonds

set cartoon_transparency, 0.2
zoom if_uprot or if_partner
"""
    out_pml.write_text(pml, encoding="utf-8")

def main():
    parser = argparse.ArgumentParser(description="Extract residue-level protein interface contacts from PDB")
    parser.add_argument("pdb", help="Path to PDB coordinate file")
    parser.add_argument("--pae", default=None, help="Path to PAE matrix JSON")
    parser.add_argument("--cutoff", type=float, default=4.5, help="Distance cutoff in Angstroms (default: 4.5)")
    parser.add_argument("--out-tsv", default=None, help="Output TSV for residue contacts")
    parser.add_argument("--out-pml", default=None, help="Output PyMOL PML script")
    
    args = parser.parse_args()
    pdb_path = Path(args.pdb)
    pae_path = Path(args.pae) if args.pae else pdb_path.parent.parent / "pae" / f"{pdb_path.stem}_pae.json"
    
    contacts_df, summary = extract_interface(pdb_path, pae_path, dist_cutoff=args.cutoff)
    
    print(f"=== Interface Contact Analysis for {pdb_path.name} ===")
    print(f"Microprotein: {summary['microprotein_total_len']} aa | Partner: {summary['partner_total_len']} aa")
    print(f"Total Contacts (<= {args.cutoff} A): {summary['total_contacts']}")
    print(f"Microprotein Contact Residues: {summary['microprotein_interface_count']} ({summary['microprotein_surface_involvement_pct']}% of peptide)")
    print(f"Partner Contact Residues:      {summary['partner_interface_count']}")
    print(f"Interface Chemistry:           H-Bonds={summary['hydrogen_bonds_count']}, Salt Bridges={summary['salt_bridges_count']}, Hydrophobic={summary['hydrophobic_contacts_count']}")
    print(f"Interface pLDDT:               uProt={summary['interface_plddt_microprotein']}, Partner={summary['interface_plddt_partner']}")
    
    print("\n--- Top 20 Closest Residue-Residue Contacts ---")
    cols = ["uProt_res_num", "uProt_res_name", "uProt_atom", "partner_res_num", "partner_res_name", "partner_atom", "min_distance_A", "interaction_type", "pae_A"]
    print(contacts_df[cols].head(20).to_string(index=False))
    
    if args.out_tsv:
        contacts_df.to_csv(args.out_tsv, sep="\t", index=False)
        print(f"\nContacts TSV saved: {args.out_tsv}")
        
    if args.out_pml:
        generate_pymol_script(pdb_path, summary, Path(args.out_pml))
        print(f"PyMOL script saved: {args.out_pml}")

if __name__ == "__main__":
    main()
