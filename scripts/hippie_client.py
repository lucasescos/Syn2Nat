#!/usr/bin/env python3
"""
HIPPIE API Client
=================
Programmatic client for querying the Human Integrated Protein-Protein
Interaction rEference (HIPPIE v2.4, CBDM Group / JGU Mainz).

Web server: https://cbdm-01.zdv.uni-mainz.de/~mschaefer/hippie/
Documentation: docs/hippie.md
"""

from __future__ import annotations

import argparse
import io
import json
import sys
import urllib.parse
import urllib.request
from typing import Any, Dict, Iterable, List, Optional, Union

try:
    import pandas as pd
    HAS_PANDAS = True
except ImportError:
    HAS_PANDAS = False

HIPPIE_ENDPOINT = "https://cbdm-01.zdv.uni-mainz.de/~mschaefer/hippie/fast_query_tissue.php"

COLUMN_NAMES_CONCISE = [
    "uniprot_1",
    "entrez_1",
    "gene_1",
    "uniprot_2",
    "entrez_2",
    "gene_2",
    "score"
]

COLUMN_NAMES_MITAB = [
    "id_interactor_a",
    "id_interactor_b",
    "alt_ids_interactor_a",
    "alt_ids_interactor_b",
    "aliases_interactor_a",
    "aliases_interactor_b",
    "interaction_detection_methods",
    "first_author",
    "publication_ids",
    "taxid_interactor_a",
    "taxid_interactor_b",
    "interaction_types",
    "source_databases",
    "interaction_identifiers",
    "confidence_score"
]


class HippieClient:
    """Client for the HIPPIE REST interface."""

    def __init__(self, endpoint: str = HIPPIE_ENDPOINT, timeout: int = 30):
        self.endpoint = endpoint
        self.timeout = timeout

    def _prepare_gene_string(self, genes: Union[str, Iterable[str]]) -> str:
        """Standardize input genes into a newline-separated query string."""
        if isinstance(genes, str):
            # Split commas, semicolons, or pipes if pasted as a single string
            tokens = [g.strip() for g in genes.replace(";", ",").replace("|", ",").replace("\n", ",").split(",") if g.strip()]
            return "\n".join(tokens)
        return "\n".join(str(g).strip() for g in genes if str(g).strip())

    def query_raw(
        self,
        proteins: Union[str, Iterable[str]],
        layer: int = 1,
        conf_thres: float = 0.0,
        out_type: str = "conc_file",
        extra_params: Optional[Dict[str, Any]] = None
    ) -> str:
        """
        Send raw query to HIPPIE and return the plain response text.
        
        Parameters
        ----------
        proteins : str or iterable of str
            Gene symbols, UniProt IDs/accessions, or Entrez IDs (max 100).
        layer : int, optional
            0: interactions strictly within query set.
            1: interactions between query set and full HIPPIE interactome. Default is 1.
        conf_thres : float, optional
            Score threshold between 0.0 and 1.0 (e.g. 0.63 medium, 0.72 high). Default is 0.0.
        out_type : str, optional
            'conc_file' for TSV table, 'mitab' for PSI-MI TAB 2.5 format.
        extra_params : dict, optional
            Additional form parameters (e.g. tissue checkboxes, direction_type).
        """
        gene_str = self._prepare_gene_string(proteins)
        if not gene_str:
            raise ValueError("No valid protein identifiers provided.")

        payload = {
            "query_genes": gene_str,
            "layers": str(layer),
            "conf_thres": str(conf_thres),
            "out_type": out_type,
        }
        if extra_params:
            payload.update(extra_params)

        data = urllib.parse.urlencode(payload).encode("utf-8")
        req = urllib.request.Request(
            self.endpoint,
            data=data,
            headers={
                "User-Agent": "HIPPIE-Python-Client/2.4 (microproteinproject; research-use)",
                "Content-Type": "application/x-www-form-urlencoded"
            }
        )

        with urllib.request.urlopen(req, timeout=self.timeout) as resp:
            text = resp.read().decode("utf-8", errors="replace").strip()

        return text

    def query(
        self,
        proteins: Union[str, Iterable[str]],
        layer: int = 1,
        conf_thres: float = 0.0,
        out_type: str = "conc_file",
        extra_params: Optional[Dict[str, Any]] = None
    ) -> Any:
        """
        Query HIPPIE and parse the result into a pandas.DataFrame (if available)
        or a list of dictionaries.
        """
        raw_text = self.query_raw(
            proteins=proteins,
            layer=layer,
            conf_thres=conf_thres,
            out_type=out_type,
            extra_params=extra_params
        )

        if not raw_text or "No interactions passing the filter criteria" in raw_text:
            if HAS_PANDAS:
                cols = COLUMN_NAMES_CONCISE if out_type == "conc_file" else COLUMN_NAMES_MITAB
                return pd.DataFrame(columns=cols)
            return []

        if out_type == "conc_file":
            lines = raw_text.splitlines()
            if not lines:
                return pd.DataFrame(columns=COLUMN_NAMES_CONCISE) if HAS_PANDAS else []
            
            # First line is usually header: uniprot id 1 \t entrez gene id 1 ...
            data_lines = lines[1:] if "uniprot id 1" in lines[0].lower() else lines
            
            rows = []
            for line in data_lines:
                parts = line.split("\t")
                if len(parts) >= 7:
                    try:
                        score = float(parts[6])
                    except (ValueError, IndexError):
                        score = 0.0
                    rows.append({
                        "uniprot_1": parts[0].strip(),
                        "entrez_1": parts[1].strip(),
                        "gene_1": parts[2].strip(),
                        "uniprot_2": parts[3].strip(),
                        "entrez_2": parts[4].strip(),
                        "gene_2": parts[5].strip(),
                        "score": score
                    })

            if HAS_PANDAS:
                df = pd.DataFrame(rows)
                if not df.empty:
                    df["score"] = pd.to_numeric(df["score"], errors="coerce")
                return df
            return rows

        elif out_type == "mitab":
            lines = raw_text.splitlines()
            rows = []
            for line in lines:
                parts = line.split("\t")
                if len(parts) >= 15:
                    try:
                        score = float(parts[14])
                    except (ValueError, IndexError):
                        score = 0.0
                    rec = {COLUMN_NAMES_MITAB[i]: parts[i].strip() for i in range(min(len(parts), len(COLUMN_NAMES_MITAB)))}
                    rec["confidence_score"] = score
                    rows.append(rec)

            if HAS_PANDAS:
                df = pd.DataFrame(rows)
                if not df.empty:
                    df["confidence_score"] = pd.to_numeric(df["confidence_score"], errors="coerce")
                return df
            return rows

        return raw_text

    def get_interactors(
        self,
        protein: str,
        conf_thres: float = 0.63
    ) -> Any:
        """Retrieve all interactors for a single protein above confidence threshold."""
        return self.query(proteins=[protein], layer=1, conf_thres=conf_thres, out_type="conc_file")

    def query_subnetwork(
        self,
        proteins: Iterable[str],
        conf_thres: float = 0.0
    ) -> Any:
        """Retrieve interactions strictly internal to the provided protein list (layer 0)."""
        return self.query(proteins=proteins, layer=0, conf_thres=conf_thres, out_type="conc_file")

    def query_pair(
        self,
        protein_a: str,
        protein_b: str,
        conf_thres: float = 0.0
    ) -> Optional[Dict[str, Any]]:
        """
        Check if two proteins interact in HIPPIE and return the highest-scoring record.
        """
        df = self.query_subnetwork([protein_a, protein_b], conf_thres=conf_thres)
        if HAS_PANDAS:
            if df.empty:
                return None
            
            p_a_upper = protein_a.strip().upper()
            p_b_upper = protein_b.strip().upper()
            
            # Filter matches in either direction
            mask = (
                ((df["gene_1"].str.upper() == p_a_upper) | (df["uniprot_1"].str.upper() == p_a_upper) | (df["entrez_1"] == p_a_upper)) &
                ((df["gene_2"].str.upper() == p_b_upper) | (df["uniprot_2"].str.upper() == p_b_upper) | (df["entrez_2"] == p_b_upper))
            ) | (
                ((df["gene_1"].str.upper() == p_b_upper) | (df["uniprot_1"].str.upper() == p_b_upper) | (df["entrez_1"] == p_b_upper)) &
                ((df["gene_2"].str.upper() == p_a_upper) | (df["uniprot_2"].str.upper() == p_a_upper) | (df["entrez_2"] == p_a_upper))
            )
            matched = df[mask]
            if matched.empty:
                return None
            return matched.sort_values(by="score", ascending=False).iloc[0].to_dict()
        else:
            best = None
            for rec in df:
                # check match
                g1, g2 = rec["gene_1"].upper(), rec["gene_2"].upper()
                u1, u2 = rec["uniprot_1"].upper(), rec["uniprot_2"].upper()
                p1, p2 = protein_a.upper(), protein_b.upper()
                if (p1 in (g1, u1) and p2 in (g2, u2)) or (p2 in (g1, u1) and p1 in (g2, u2)):
                    if best is None or rec["score"] > best["score"]:
                        best = rec
            return best


# Convenient module-level functions
_default_client = HippieClient()

def query_hippie(
    proteins: Union[str, Iterable[str]],
    layer: int = 1,
    conf_thres: float = 0.0,
    out_type: str = "conc_file"
) -> Any:
    return _default_client.query(proteins, layer=layer, conf_thres=conf_thres, out_type=out_type)

def query_pair(protein_a: str, protein_b: str, conf_thres: float = 0.0) -> Optional[Dict[str, Any]]:
    return _default_client.query_pair(protein_a, protein_b, conf_thres=conf_thres)

def get_interactors(protein: str, conf_thres: float = 0.63) -> Any:
    return _default_client.get_interactors(protein, conf_thres=conf_thres)


def main():
    parser = argparse.ArgumentParser(description="Query HIPPIE human protein-protein interaction reference.")
    parser.add_argument(
        "-p", "--proteins",
        type=str,
        required=True,
        help="Comma-separated protein identifiers (Gene Symbol, UniProt ID, Entrez ID)."
    )
    parser.add_argument(
        "-l", "--layer",
        type=int,
        choices=[0, 1],
        default=1,
        help="0: interactions within set only; 1: interactions between set and all HIPPIE. (default: 1)"
    )
    parser.add_argument(
        "-c", "--conf", "--conf-thres",
        type=float,
        default=0.0,
        dest="conf_thres",
        help="Confidence threshold between 0.0 and 1.0 (e.g. 0.63 for medium, 0.72 for high)."
    )
    parser.add_argument(
        "-f", "--format",
        type=str,
        choices=["conc_file", "mitab", "raw"],
        default="conc_file",
        help="Output format from API (conc_file, mitab, raw)."
    )
    parser.add_argument(
        "-o", "--output",
        type=str,
        default=None,
        help="Path to save output table (supports .csv, .tsv, .json)."
    )
    parser.add_argument(
        "--pair",
        action="store_true",
        help="Evaluate as a pair query between exactly two proteins."
    )

    args = parser.parse_args()
    client = HippieClient()

    proteins = [p.strip() for p in args.proteins.split(",") if p.strip()]

    if args.pair:
        if len(proteins) != 2:
            print("Error: --pair requires exactly two proteins separated by a comma (e.g. --proteins TP53,MDM2)", file=sys.stderr)
            sys.exit(1)
        res = client.query_pair(proteins[0], proteins[1], conf_thres=args.conf_thres)
        if res:
            print(f"Interaction found between {proteins[0]} and {proteins[1]}:")
            print(json.dumps(res, indent=2))
        else:
            print(f"No interaction found between {proteins[0]} and {proteins[1]} at confidence threshold >= {args.conf_thres}.")
        return

    if args.format == "raw":
        result = client.query_raw(proteins, layer=args.layer, conf_thres=args.conf_thres, out_type="conc_file")
        print(result)
        return

    result = client.query(proteins, layer=args.layer, conf_thres=args.conf_thres, out_type=args.format)

    if HAS_PANDAS and isinstance(result, pd.DataFrame):
        print(f"Retrieved {len(result)} interactions (layer={args.layer}, conf_thres={args.conf_thres})")
        if not result.empty:
            print(result.head(20).to_string(index=False))
            if len(result) > 20:
                print(f"... and {len(result) - 20} more rows.")

        if args.output:
            if args.output.endswith(".csv"):
                result.to_csv(args.output, index=False)
            elif args.output.endswith(".json"):
                result.to_json(args.output, orient="records", indent=2)
            else:
                result.to_csv(args.output, sep="\t", index=False)
            print(f"Saved results to {args.output}")
    else:
        print(f"Retrieved {len(result)} interactions.")
        if args.output:
            with open(args.output, "w", encoding="utf-8") as f:
                json.dump(result, f, indent=2)
            print(f"Saved results to {args.output}")
        else:
            print(json.dumps(result[:10], indent=2))


if __name__ == "__main__":
    main()
