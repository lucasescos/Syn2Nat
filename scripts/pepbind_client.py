#!/usr/bin/env python3
"""
PepBind & Protein-Peptide Complex API Client
============================================
Programmatic tools and APIs for querying peptide-binding protein complexes 
and computing inter-chain contact interfaces according to PepBind's PICI algorithm.

Reference:
Das AA, Sharma OP, Kumar MS, Krishna R, Mathur PP.
PepBind: A comprehensive database and computational tool for analysis of protein–peptide interactions.
Genomics, Proteomics & Bioinformatics (2013) 11(4):241–246.
DOI: 10.1016/j.gpb.2013.03.002
"""

from __future__ import annotations

import argparse
import io
import json
import math
import sys
import urllib.error
import urllib.parse
import urllib.request
from typing import Any, Dict, List, Optional, Set, Tuple

try:
    from Bio.PDB import PDBParser, MMCIFParser, NeighborSearch
    HAS_BIOPYTHON = True
except ImportError:
    HAS_BIOPYTHON = False

# Legacy PepBind Server Endpoints (Pondicherry University)
PEPBIND_BASE_URL = "http://pepbind.bicpu.edu.in"
RCSB_SEARCH_API = "https://search.rcsb.org/rcsbsearch/v2/query"
RCSB_DATA_API = "https://data.rcsb.org/rest/v1/core/entry"
PDBE_PISA_API = "https://www.ebi.ac.uk/pdbe/api/pisa/interactions"


class PepBindLegacyClient:
    """
    Client for interacting with the legacy PepBind web endpoints (Pondicherry University).
    Note: Institutional academic servers may experience connectivity timeouts or downtime.
    """

    def __init__(self, base_url: str = PEPBIND_BASE_URL, timeout: int = 8):
        self.base_url = base_url.rstrip("/")
        self.timeout = timeout
        self.headers = {
            "User-Agent": "Mozilla/5.0 (Windows NT 10.0; Win64; x64) PepBindClient/1.0",
            "Accept": "text/html,application/xhtml+xml,application/xml;q=0.9,*/*;q=0.8"
        }

    def fetch_entry(self, pdb_id: str) -> Dict[str, Any]:
        """Fetch details of a complex by 4-letter PDB ID from PepBind details.php."""
        url = f"{self.base_url}/details.php?id={pdb_id.upper()}"
        req = urllib.request.Request(url, headers=self.headers)
        try:
            with urllib.request.urlopen(req, timeout=self.timeout) as resp:
                content = resp.read().decode("utf-8", errors="ignore")
                return {
                    "status": "success",
                    "pdb_id": pdb_id.upper(),
                    "source": url,
                    "html_length": len(content),
                    "raw_content": content[:2000]
                }
        except Exception as err:
            return {
                "status": "error",
                "pdb_id": pdb_id.upper(),
                "error": str(err),
                "note": "Legacy PepBind server is currently unreachable. Use modern RCSB/PDBe APIs."
            }

    def predict_binding(self, peptide_sequence: str) -> Dict[str, Any]:
        """Query PepBind's beta domain prediction server."""
        url = f"{self.base_url}/PepBind_prediction_beta.php"
        data = urllib.parse.urlencode({"pep_seq": peptide_sequence.strip().upper()}).encode("utf-8")
        req = urllib.request.Request(url, data=data, headers=self.headers)
        try:
            with urllib.request.urlopen(req, timeout=self.timeout) as resp:
                content = resp.read().decode("utf-8", errors="ignore")
                return {
                    "status": "success",
                    "peptide_sequence": peptide_sequence,
                    "source": url,
                    "raw_content": content[:2000]
                }
        except Exception as err:
            return {
                "status": "error",
                "peptide_sequence": peptide_sequence,
                "error": str(err),
                "note": "Prediction server unreachable. Consider using Foldseek or ESMFold2 co-folding."
            }


class ModernPeptideComplexAPI:
    """
    Production-grade RESTful client to search and extract peptide-protein complexes 
    from the RCSB PDB and PDBe APIs matching PepBind's criteria (peptide length <= 35 residues).
    """

    @staticmethod
    def search_peptide_complexes(
        max_peptide_length: int = 35,
        min_peptide_length: int = 3,
        rows: int = 10
    ) -> List[str]:
        """
        Query RCSB PDB Search API for complexes containing at least one short peptide chain (<= max_length)
        and at least one receptor protein chain.
        """
        payload = {
            "query": {
                "type": "group",
                "logical_operator": "and",
                "nodes": [
                    {
                        "type": "terminal",
                        "service": "text",
                        "parameters": {
                            "attribute": "entity_poly.rcsb_sample_sequence_length",
                            "operator": "less_or_equal",
                            "value": max_peptide_length
                        }
                    },
                    {
                        "type": "terminal",
                        "service": "text",
                        "parameters": {
                            "attribute": "entity_poly.rcsb_sample_sequence_length",
                            "operator": "greater_or_equal",
                            "value": min_peptide_length
                        }
                    },
                    {
                        "type": "terminal",
                        "service": "text",
                        "parameters": {
                            "attribute": "rcsb_entry_info.polymer_entity_count_protein",
                            "operator": "greater_or_equal",
                            "value": 2
                        }
                    }
                ]
            },
            "return_type": "entry",
            "request_options": {
                "paginate": {"start": 0, "rows": rows}
            }
        }
        
        req = urllib.request.Request(
            RCSB_SEARCH_API,
            data=json.dumps(payload).encode("utf-8"),
            headers={"Content-Type": "application/json", "User-Agent": "PepBindClient/1.0"}
        )
        
        try:
            with urllib.request.urlopen(req, timeout=12) as resp:
                data = json.loads(resp.read().decode("utf-8"))
                results = [item["identifier"] for item in data.get("result_set", [])]
                return results
        except Exception as e:
            sys.stderr.write(f"[RCSB Search Error]: {e}\n")
            return []

    @staticmethod
    def get_entry_metadata(pdb_id: str) -> Dict[str, Any]:
        """Retrieve entry summary and polymer info via RCSB REST Data API."""
        url = f"{RCSB_DATA_API}/{pdb_id.upper()}"
        req = urllib.request.Request(url, headers={"User-Agent": "PepBindClient/1.0"})
        try:
            with urllib.request.urlopen(req, timeout=10) as resp:
                return json.loads(resp.read().decode("utf-8"))
        except Exception as e:
            return {"error": str(e)}


class PiciEngine:
    """
    Exact Python implementation of PepBind's PICI (Protein Inter-Chain Interaction) tool.
    Calculates interface contacts between a core protein chain and a short peptide chain:
      - Hydrogen bonds: O/N donor-acceptor distance <= 3.5 Å
      - Hydrophobic interactions: Nonpolar carbon-carbon contacts <= 5.0 Å
      - Ionic interactions / salt bridges: Opposite charge distance <= 6.0 Å
      - Disulfide bonds: SG-SG distance <= 2.2 Å
    """

    HYDROPHOBIC_RESIDUES: Set[str] = {
        "ALA", "VAL", "LEU", "ILE", "MET", "PHE", "TRP", "PRO", "TYR"
    }
    POS_CHARGED_RESIDUES: Set[str] = {"ARG", "LYS", "HIS"}
    NEG_CHARGED_RESIDUES: Set[str] = {"ASP", "GLU"}

    @classmethod
    def calculate_interface(
        cls,
        pdb_path: str,
        protein_chain_id: str,
        peptide_chain_id: str
    ) -> Dict[str, Any]:
        """
        Analyze the protein-peptide interface using PICI criteria.
        """
        if not HAS_BIOPYTHON:
            raise ImportError("Biopython is required for PiciEngine: run 'uv add biopython'")

        if pdb_path.endswith(".cif") or pdb_path.endswith(".mmcif"):
            parser = MMCIFParser(QUIET=True)
        else:
            parser = PDBParser(QUIET=True)

        structure = parser.get_structure("complex", pdb_path)
        # Select first model (as per PepBind specifications for NMR or crystal)
        model = next(iter(structure))

        if protein_chain_id not in model or peptide_chain_id not in model:
            available_chains = [c.id for c in model.get_chains()]
            raise KeyError(
                f"Specified chains ({protein_chain_id}, {peptide_chain_id}) not found. "
                f"Available chains: {available_chains}"
            )

        prot_chain = model[protein_chain_id]
        pep_chain = model[peptide_chain_id]

        prot_atoms = [a for a in prot_chain.get_atoms() if a.element != "H"]
        pep_atoms = [a for a in pep_chain.get_atoms() if a.element != "H"]

        ns = NeighborSearch(prot_atoms)

        h_bonds: List[Dict[str, Any]] = []
        hydrophobic: List[Dict[str, Any]] = []
        ionic: List[Dict[str, Any]] = []
        disulfides: List[Dict[str, Any]] = []

        seen_pairs: Set[Tuple[str, str]] = set()

        for atom_p in pep_atoms:
            res_p = atom_p.get_parent()
            res_p_name = res_p.get_resname().strip()
            res_p_id = res_p.id[1]

            # Neighbor query within 6.0 Å (max threshold for ionic)
            neighbors = ns.search(atom_p.coord, 6.0)

            for atom_r in neighbors:
                res_r = atom_r.get_parent()
                res_r_name = res_r.get_resname().strip()
                res_r_id = res_r.id[1]

                pair_key = (f"{res_p_name}{res_p_id}:{atom_p.name}", f"{res_r_name}{res_r_id}:{atom_r.name}")
                if pair_key in seen_pairs:
                    continue
                seen_pairs.add(pair_key)

                dist = float(atom_p - atom_r)

                # 1. Hydrogen Bond (O/N pairs <= 3.5 Å)
                if dist <= 3.5 and atom_p.element in ("N", "O") and atom_r.element in ("N", "O"):
                    h_bonds.append({
                        "peptide_res": f"{res_p_name}_{res_p_id}",
                        "protein_res": f"{res_r_name}_{res_r_id}",
                        "peptide_atom": atom_p.name,
                        "protein_atom": atom_r.name,
                        "distance_angstrom": round(dist, 2)
                    })

                # 2. Hydrophobic Interaction (<= 5.0 Å between nonpolar carbons)
                if dist <= 5.0 and res_p_name in cls.HYDROPHOBIC_RESIDUES and res_r_name in cls.HYDROPHOBIC_RESIDUES:
                    if atom_p.element == "C" and atom_r.element == "C":
                        hydrophobic.append({
                            "peptide_res": f"{res_p_name}_{res_p_id}",
                            "protein_res": f"{res_r_name}_{res_r_id}",
                            "peptide_atom": atom_p.name,
                            "protein_atom": atom_r.name,
                            "distance_angstrom": round(dist, 2)
                        })

                # 3. Ionic Interaction (<= 6.0 Å between opposite charges)
                if dist <= 6.0:
                    is_salt_bridge = (
                        (res_p_name in cls.POS_CHARGED_RESIDUES and res_r_name in cls.NEG_CHARGED_RESIDUES) or
                        (res_p_name in cls.NEG_CHARGED_RESIDUES and res_r_name in cls.POS_CHARGED_RESIDUES)
                    )
                    if is_salt_bridge:
                        ionic.append({
                            "peptide_res": f"{res_p_name}_{res_p_id}",
                            "protein_res": f"{res_r_name}_{res_r_id}",
                            "peptide_atom": atom_p.name,
                            "protein_atom": atom_r.name,
                            "distance_angstrom": round(dist, 2)
                        })

                # 4. Disulfide Bond (<= 2.2 Å between Cys SG atoms)
                if dist <= 2.2 and res_p_name == "CYS" and res_r_name == "CYS":
                    if atom_p.name == "SG" and atom_r.name == "SG":
                        disulfides.append({
                            "peptide_res": f"{res_p_name}_{res_p_id}",
                            "protein_res": f"{res_r_name}_{res_r_id}",
                            "distance_angstrom": round(dist, 2)
                        })

        # Aggregate by unique residue pairs
        res_h_bonds = {}
        for hb in h_bonds:
            key = (hb["peptide_res"], hb["protein_res"])
            if key not in res_h_bonds or hb["distance_angstrom"] < res_h_bonds[key]["min_distance"]:
                res_h_bonds[key] = {
                    "peptide_res": hb["peptide_res"],
                    "protein_res": hb["protein_res"],
                    "atoms": f"{hb['peptide_atom']}--{hb['protein_atom']}",
                    "min_distance": hb["distance_angstrom"]
                }

        res_hydrophobic = {}
        for hp in hydrophobic:
            key = (hp["peptide_res"], hp["protein_res"])
            if key not in res_hydrophobic or hp["distance_angstrom"] < res_hydrophobic[key]["min_distance"]:
                res_hydrophobic[key] = {
                    "peptide_res": hp["peptide_res"],
                    "protein_res": hp["protein_res"],
                    "min_distance": hp["distance_angstrom"]
                }

        res_ionic = {}
        for ion in ionic:
            key = (ion["peptide_res"], ion["protein_res"])
            if key not in res_ionic or ion["distance_angstrom"] < res_ionic[key]["min_distance"]:
                res_ionic[key] = {
                    "peptide_res": ion["peptide_res"],
                    "protein_res": ion["protein_res"],
                    "min_distance": ion["distance_angstrom"]
                }

        return {
            "summary": {
                "unique_hydrogen_bonded_residue_pairs": len(res_h_bonds),
                "unique_hydrophobic_residue_pairs": len(res_hydrophobic),
                "unique_ionic_residue_pairs": len(res_ionic),
                "disulfide_bonds_count": len(disulfides),
                "total_atom_contacts": len(h_bonds) + len(hydrophobic) + len(ionic) + len(disulfides)
            },
            "hydrogen_bonds_residue_pairs": list(res_h_bonds.values()),
            "hydrophobic_residue_pairs": list(res_hydrophobic.values()),
            "ionic_residue_pairs": list(res_ionic.values()),
            "disulfide_bonds": disulfides
        }


def main():
    parser = argparse.ArgumentParser(
        description="PepBind & Protein-Peptide Complex API Client and PICI Engine."
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    # Subcommand: query-legacy
    sub_leg = subparsers.add_parser("query-legacy", help="Query legacy PepBind web endpoint.")
    sub_leg.add_argument("pdb_id", help="4-letter PDB ID (e.g. 1CP3)")

    # Subcommand: search-rcsb
    sub_rcsb = subparsers.add_parser("search-rcsb", help="Search RCSB PDB for peptide-protein complexes.")
    sub_rcsb.add_argument("--max-length", type=int, default=35, help="Maximum peptide length (default: 35)")
    sub_rcsb.add_argument("--limit", type=int, default=10, help="Number of results to retrieve")

    # Subcommand: pici
    sub_pici = subparsers.add_parser("pici", help="Calculate interface contacts using PepBind PICI algorithm.")
    sub_pici.add_argument("pdb_file", help="Path to PDB or CIF structure file")
    sub_pici.add_argument("--protein-chain", required=True, help="Chain ID of receptor protein (e.g. A)")
    sub_pici.add_argument("--peptide-chain", required=True, help="Chain ID of bound peptide (e.g. B)")
    sub_pici.add_argument("--output", help="Optional output JSON file")

    args = parser.parse_args()

    if args.command == "query-legacy":
        client = PepBindLegacyClient()
        res = client.fetch_entry(args.pdb_id)
        print(json.dumps(res, indent=2))

    elif args.command == "search-rcsb":
        hits = ModernPeptideComplexAPI.search_peptide_complexes(
            max_peptide_length=args.max_length,
            rows=args.limit
        )
        print(json.dumps({"query": f"peptide_length <= {args.max_length}", "hits": hits}, indent=2))

    elif args.command == "pici":
        res = PiciEngine.calculate_interface(
            pdb_path=args.pdb_file,
            protein_chain_id=args.protein_chain,
            peptide_chain_id=args.peptide_chain
        )
        if args.output:
            with open(args.output, "w") as f:
                json.dump(res, f, indent=2)
            print(f"Interface contacts written to {args.output}")
        else:
            print(json.dumps(res, indent=2))


if __name__ == "__main__":
    main()
