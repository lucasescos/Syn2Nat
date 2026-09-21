#!/usr/bin/env python3
"""
BioGRID API Client
==================
Programmatic Python client for querying the Biological General Repository for
Interaction Datasets (BioGRID REST Service API).

Web server: https://thebiogrid.org
REST API: https://webservice.thebiogrid.org
Documentation: docs/biogrid.md
"""

from __future__ import annotations

import argparse
import io
import json
import os
import sys
import urllib.parse
import urllib.request
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Union

try:
    import pandas as pd
    HAS_PANDAS = True
except ImportError:
    HAS_PANDAS = False

BIOGRID_BASE_URL = "https://webservice.thebiogrid.org"

TAB2_COLUMNS = [
    "BIOGRID_INTERACTION_ID",
    "ENTREZ_GENE_A",
    "ENTREZ_GENE_B",
    "BIOGRID_ID_A",
    "BIOGRID_ID_B",
    "SYSTEMATIC_NAME_A",
    "SYSTEMATIC_NAME_B",
    "OFFICIAL_SYMBOL_A",
    "OFFICIAL_SYMBOL_B",
    "SYNONYMS_A",
    "SYNONYMS_B",
    "EXPERIMENTAL_SYSTEM",
    "EXPERIMENTAL_SYSTEM_TYPE",
    "PUBMED_AUTHOR",
    "PUBMED_ID",
    "ORGANISM_A_ID",
    "ORGANISM_B_ID",
    "THROUGHPUT",
    "QUANTITATIVE_SCORE",
    "POST_TRANSLATIONAL_MODIFICATION",
    "PHENOTYPES",
    "QUALIFICATIONS",
    "TAGS",
    "SOURCE_DATABASE"
]


def load_env_key() -> Optional[str]:
    """Check environment variables or local .env file for BIOGRID_ACCESS_KEY."""
    if "BIOGRID_ACCESS_KEY" in os.environ:
        return os.environ["BIOGRID_ACCESS_KEY"]
    if "BIOGRID_API_KEY" in os.environ:
        return os.environ["BIOGRID_API_KEY"]
    
    # Search upwards for .env
    cur = Path.cwd()
    for _ in range(4):
        env_file = cur / ".env"
        if env_file.is_file():
            try:
                with open(env_file, "r", encoding="utf-8") as f:
                    for line in f:
                        line = line.strip()
                        if line.startswith("BIOGRID_ACCESS_KEY="):
                            return line.split("=", 1)[1].strip().strip('"').strip("'")
                        if line.startswith("BIOGRID_API_KEY="):
                            return line.split("=", 1)[1].strip().strip('"').strip("'")
            except Exception:
                pass
        if cur.parent == cur:
            break
        cur = cur.parent
    return None


class BioGridClient:
    """Production client for the BioGRID REST API."""

    def __init__(
        self,
        access_key: Optional[str] = None,
        base_url: str = BIOGRID_BASE_URL,
        timeout: int = 45
    ):
        self.access_key = access_key or load_env_key()
        self.base_url = base_url.rstrip("/")
        self.timeout = timeout

    def _ensure_access_key(self) -> str:
        """Ensure an access key is set before making requests."""
        if not self.access_key:
            raise ValueError(
                "BioGRID Access Key not found! "
                "Register for a free key at https://webservice.thebiogrid.org/ and "
                "either pass access_key to BioGridClient(), set os.environ['BIOGRID_ACCESS_KEY'], "
                "or add BIOGRID_ACCESS_KEY=your_key to your project .env file."
            )
        return self.access_key

    def _get(self, endpoint: str, params: Dict[str, Any]) -> str:
        """Perform HTTP GET request and return plain response text."""
        clean_params = {k: str(v) for k, v in params.items() if v is not None}
        encoded_query = urllib.parse.urlencode(clean_params)
        url = f"{self.base_url}/{endpoint.lstrip('/')}?{encoded_query}"

        req = urllib.request.Request(
            url,
            headers={
                "User-Agent": "BioGridPythonClient/1.0 (MicroproteinInteractomePipeline)"
            }
        )
        try:
            with urllib.request.urlopen(req, timeout=self.timeout) as resp:
                return resp.read().decode("utf-8", errors="replace")
        except urllib.error.HTTPError as err:
            err_body = err.read().decode("utf-8", errors="replace")
            raise RuntimeError(f"BioGRID API HTTP {err.code} error on {url}: {err_body}") from err
        except Exception as exc:
            raise RuntimeError(f"Failed to connect to BioGRID API ({url}): {exc}") from exc

    def get_version(self) -> str:
        """Retrieve active BioGRID database and REST version."""
        key = self._ensure_access_key()
        res = self._get("/version/", {"accesskey": key})
        return res.strip()

    def get_organisms(self, format: str = "json") -> Any:
        """Retrieve list of supported NCBI organisms."""
        key = self._ensure_access_key()
        res = self._get("/organisms/", {"accesskey": key, "format": format})
        return json.loads(res) if format == "json" else res

    def get_evidence_types(self, format: str = "json") -> Any:
        """Retrieve supported experimental evidence codes."""
        key = self._ensure_access_key()
        res = self._get("/evidence/", {"accesskey": key, "format": format})
        return json.loads(res) if format == "json" else res

    def get_identifier_types(self, format: str = "json") -> Any:
        """Retrieve supported additional identifier types (e.g. UNIPROT, ENSEMBL)."""
        key = self._ensure_access_key()
        res = self._get("/identifiers/", {"accesskey": key, "format": format})
        return json.loads(res) if format == "json" else res

    def get_interaction(self, interaction_id: Union[int, str], format: str = "json") -> Any:
        """Fetch a single interaction by its BioGRID Interaction ID."""
        key = self._ensure_access_key()
        res = self._get(f"/interaction/{interaction_id}", {"accesskey": key, "format": format})
        if format in ("json", "jsonExtended"):
            data = json.loads(res)
            return data.get(str(interaction_id), data)
        return res

    def count(
        self,
        genes: Optional[Union[str, Iterable[str]]] = None,
        tax_id: Optional[Union[int, str]] = None,
        search_names: bool = True,
        search_ids: bool = False,
        search_synonyms: bool = False,
        include_interactors: bool = True,
        inter_species_excluded: bool = False,
        self_interactions_excluded: bool = False,
        throughput: str = "any",
        exp_type: Optional[str] = None,
    ) -> int:
        """Return total count of interactions matching search criteria."""
        key = self._ensure_access_key()
        params: Dict[str, Any] = {
            "accesskey": key,
            "format": "count",
            "searchNames": "true" if search_names else "false",
            "searchIds": "true" if search_ids else "false",
            "searchSynonyms": "true" if search_synonyms else "false",
            "includeInteractors": "true" if include_interactors else "false",
            "interSpeciesExcluded": "true" if inter_species_excluded else "false",
            "selfInteractionsExcluded": "true" if self_interactions_excluded else "false",
            "throughputTag": throughput,
        }
        if genes:
            if isinstance(genes, str):
                params["geneList"] = genes.replace(",", "|").replace(";", "|")
            else:
                params["geneList"] = "|".join(str(g).strip() for g in genes if str(g).strip())
        if tax_id:
            params["taxId"] = tax_id

        res = self._get("/interactions/", params).strip()
        try:
            return int(res)
        except ValueError:
            # Check if response is JSON error
            try:
                err = json.loads(res)
                raise RuntimeError(f"BioGRID error: {err}")
            except Exception:
                raise RuntimeError(f"Unexpected response from BioGRID count: {res}")

    def query(
        self,
        genes: Optional[Union[str, Iterable[str]]] = None,
        tax_id: Optional[Union[int, str]] = 9606,
        search_names: bool = True,
        search_ids: bool = False,
        search_synonyms: bool = False,
        additional_identifier_types: Optional[str] = None,
        include_interactors: bool = True,
        include_interactor_interactions: bool = False,
        inter_species_excluded: bool = False,
        self_interactions_excluded: bool = False,
        evidence_list: Optional[Union[str, Iterable[str]]] = None,
        include_evidence: bool = False,
        throughput: str = "any",
        exp_type: Optional[str] = None,
        pubmed_list: Optional[Union[str, Iterable[str]]] = None,
        start: int = 0,
        max_results: int = 10000,
        paginate: bool = False,
        as_dataframe: bool = False
    ) -> Any:
        """
        Query BioGRID interactions with comprehensive multi-parameter filtering.

        Parameters
        ----------
        genes : str or iterable of str, optional
            Gene symbols, Entrez IDs, or identifiers to search.
        tax_id : int or str, optional
            NCBI Taxonomy ID (default 9606 for human). Pass "All" for all taxa.
        search_names : bool
            Search against official gene symbols (default True).
        search_ids : bool
            Search against Entrez IDs, systematic names, or ordered loci (default False).
        search_synonyms : bool
            Search against synonyms/aliases (default False).
        additional_identifier_types : str, optional
            e.g. 'UNIPROT' or 'ENSEMBL' when genes list contains UniProt/Ensembl IDs.
        include_interactors : bool
            If True (default), returns 1st-order interactors (edges where >=1 interactor is in genes).
            If False, returns edges strictly BETWEEN genes in the list (Layer 0 subnetwork).
        include_interactor_interactions : bool
            If True, returns interactions among the 1st-order interactors themselves.
        inter_species_excluded : bool
            If True, excludes interactions between different species.
        self_interactions_excluded : bool
            If True, excludes homomeric/self interactions (A-A).
        evidence_list : str or iterable of str, optional
            Experimental systems to filter by.
        include_evidence : bool
            If True, evidence_list is a whitelist. If False, evidence_list is a blacklist.
        throughput : str
            'any' (default), 'low', or 'high'.
        exp_type : str, optional
            'physical' or 'genetic'. Post-filtered if specified.
        pubmed_list : str or iterable of str, optional
            PubMed IDs to filter by.
        start : int
            Offset for pagination.
        max_results : int
            Max results per request (up to 10,000).
        paginate : bool
            If True, automatically retrieves all pages for queries with >10,000 hits.
        as_dataframe : bool
            If True, returns results as a pandas DataFrame (requires pandas).

        Returns
        -------
        list of dict, or pandas.DataFrame
        """
        key = self._ensure_access_key()
        gene_str = ""
        if genes:
            if isinstance(genes, str):
                gene_str = genes.replace(",", "|").replace(";", "|")
            else:
                gene_str = "|".join(str(g).strip() for g in genes if str(g).strip())

        params: Dict[str, Any] = {
            "accesskey": key,
            "format": "json",
            "searchNames": "true" if search_names else "false",
            "searchIds": "true" if search_ids else "false",
            "searchSynonyms": "true" if search_synonyms else "false",
            "includeInteractors": "true" if include_interactors else "false",
            "includeInteractorInteractions": "true" if include_interactor_interactions else "false",
            "interSpeciesExcluded": "true" if inter_species_excluded else "false",
            "selfInteractionsExcluded": "true" if self_interactions_excluded else "false",
            "throughputTag": throughput,
            "start": start,
            "max": min(max_results, 10000)
        }
        if gene_str:
            params["geneList"] = gene_str
        if tax_id is not None:
            params["taxId"] = tax_id
        if additional_identifier_types:
            params["additionalIdentifierTypes"] = additional_identifier_types
        if evidence_list:
            if isinstance(evidence_list, str):
                params["evidenceList"] = evidence_list.replace(",", "|")
            else:
                params["evidenceList"] = "|".join(evidence_list)
            params["includeEvidence"] = "true" if include_evidence else "false"
        if pubmed_list:
            if isinstance(pubmed_list, str):
                params["pubmedList"] = pubmed_list.replace(",", "|")
            else:
                params["pubmedList"] = "|".join(pubmed_list)

        # Handle pagination if requested
        collected_dict: Dict[str, Dict[str, Any]] = {}
        
        current_start = start
        while True:
            params["start"] = current_start
            raw = self._get("/interactions/", params)
            if not raw.strip():
                break
            try:
                page_data = json.loads(raw)
            except json.JSONDecodeError:
                raise RuntimeError(f"Invalid JSON returned from BioGRID: {raw[:300]}")

            if not isinstance(page_data, dict):
                break
            if not page_data:
                break

            collected_dict.update(page_data)

            if not paginate:
                break
            if len(page_data) < params["max"]:
                break
            current_start += len(page_data)
            if max_results and len(collected_dict) >= max_results:
                break

        # Convert dictionary to list of interaction records
        records = list(collected_dict.values())

        # Optional post-filter by physical/genetic experimental system type
        if exp_type:
            exp_type_clean = exp_type.strip().lower()
            records = [
                r for r in records
                if str(r.get("EXPERIMENTAL_SYSTEM_TYPE", "")).lower() == exp_type_clean
            ]

        if as_dataframe:
            if not HAS_PANDAS:
                raise ImportError("pandas is required for as_dataframe=True. Install with pip install pandas.")
            return pd.DataFrame(records)

        return records

    def query_pair(
        self,
        gene_a: str,
        gene_b: str,
        tax_id: Union[int, str] = 9606,
        exp_type: Optional[str] = None
    ) -> List[Dict[str, Any]]:
        """
        Check if two specific proteins physically or genetically interact.
        Returns all matching interaction records between gene_a and gene_b.
        """
        return self.query(
            genes=[gene_a, gene_b],
            tax_id=tax_id,
            include_interactors=False,  # Strictly interactions BETWEEN gene_a and gene_b
            exp_type=exp_type,
            as_dataframe=False
        )

    def query_subnetwork(
        self,
        genes: Iterable[str],
        tax_id: Union[int, str] = 9606,
        exp_type: Optional[str] = "physical",
        as_dataframe: bool = False
    ) -> Any:
        """
        Retrieve all interactions strictly internal to a candidate gene set (Layer 0 subnetwork).
        """
        return self.query(
            genes=genes,
            tax_id=tax_id,
            include_interactors=False,
            exp_type=exp_type,
            as_dataframe=as_dataframe
        )

    def query_interactors(
        self,
        gene: str,
        tax_id: Union[int, str] = 9606,
        throughput: str = "any",
        exp_type: Optional[str] = "physical",
        as_dataframe: bool = False
    ) -> Any:
        """
        Retrieve all 1st-order interactors for a given gene (Star subnetwork).
        """
        return self.query(
            genes=[gene],
            tax_id=tax_id,
            include_interactors=True,
            include_interactor_interactions=False,
            throughput=throughput,
            exp_type=exp_type,
            as_dataframe=as_dataframe
        )

    def query_by_identifier(
        self,
        identifiers: Union[str, Iterable[str]],
        id_type: str = "UNIPROT",
        tax_id: Union[int, str] = 9606,
        as_dataframe: bool = False
    ) -> Any:
        """
        Query using external identifiers (e.g. UniProt accessions or Ensembl IDs).
        """
        return self.query(
            genes=identifiers,
            tax_id=tax_id,
            search_names=False,
            additional_identifier_types=id_type,
            as_dataframe=as_dataframe
        )


def main():
    parser = argparse.ArgumentParser(
        description="BioGRID API Client: Query biological interactions programmatically."
    )
    parser.add_argument("--key", type=str, default=None, help="BioGRID 32-character Access Key")
    parser.add_argument("--genes", type=str, default=None, help="Comma-separated gene symbols")
    parser.add_argument("--taxid", type=str, default="9606", help="NCBI Tax ID (default: 9606 for human)")
    parser.add_argument("--pair", nargs=2, metavar=("GENE_A", "GENE_B"), help="Check interaction between two genes")
    parser.add_argument("--subnetwork", action="store_true", help="Extract internal subnetwork (includeInteractors=false)")
    parser.add_argument("--throughput", choices=["any", "low", "high"], default="any", help="Filter by throughput")
    parser.add_argument("--exp-type", choices=["physical", "genetic"], default=None, help="Filter physical or genetic")
    parser.add_argument("--count", action="store_true", help="Print total interaction count matching filters")
    parser.add_argument("--version", action="store_true", help="Print current BioGRID REST version")
    parser.add_argument("--interaction-id", type=int, default=None, help="Fetch single interaction by ID")
    parser.add_argument("--output", type=str, default=None, help="Path to save output table (.tsv, .csv, or .json)")
    parser.add_argument("--max", type=int, default=100, help="Max results to display/return")

    args = parser.parse_args()
    client = BioGridClient(access_key=args.key)

    if args.version:
        print(f"BioGRID REST Version: {client.get_version()}")
        return

    if args.interaction_id:
        data = client.get_interaction(args.interaction_id)
        print(json.dumps(data, indent=2))
        return

    if args.pair:
        gene_a, gene_b = args.pair
        hits = client.query_pair(gene_a, gene_b, tax_id=args.taxid, exp_type=args.exp_type)
        print(f"Found {len(hits)} interactions between {gene_a} and {gene_b}:")
        for h in hits:
            print(f"  - [{h.get('EXPERIMENTAL_SYSTEM_TYPE')}] {h.get('EXPERIMENTAL_SYSTEM')} | "
                  f"{h.get('PUBMED_AUTHOR')} (PMID: {h.get('PUBMED_ID')})")
        return

    if args.count:
        cnt = client.count(
            genes=args.genes,
            tax_id=args.taxid,
            include_interactors=not args.subnetwork,
            throughput=args.throughput,
            exp_type=args.exp_type
        )
        print(f"Total matching interactions: {cnt}")
        return

    if args.genes:
        genes = [g.strip() for g in args.genes.split(",") if g.strip()]
        if args.subnetwork:
            res = client.query_subnetwork(genes, tax_id=args.taxid, exp_type=args.exp_type, as_dataframe=HAS_PANDAS)
        else:
            res = client.query(
                genes=genes,
                tax_id=args.taxid,
                throughput=args.throughput,
                exp_type=args.exp_type,
                max_results=args.max,
                as_dataframe=HAS_PANDAS
            )

        if HAS_PANDAS and isinstance(res, pd.DataFrame):
            print(f"Retrieved {len(res)} interactions.")
            if not res.empty:
                cols = [c for c in ["OFFICIAL_SYMBOL_A", "OFFICIAL_SYMBOL_B", "EXPERIMENTAL_SYSTEM", "PUBMED_AUTHOR", "THROUGHPUT"] if c in res.columns]
                print(res[cols].head(args.max).to_string(index=False))
            if args.output:
                if args.output.endswith(".csv"):
                    res.to_csv(args.output, index=False)
                else:
                    res.to_csv(args.output, sep="\t", index=False)
                print(f"Saved to {args.output}")
        else:
            print(f"Retrieved {len(res)} interactions.")
            for r in res[:args.max]:
                print(f"  {r.get('OFFICIAL_SYMBOL_A')} <-> {r.get('OFFICIAL_SYMBOL_B')} "
                      f"({r.get('EXPERIMENTAL_SYSTEM')}, {r.get('PUBMED_AUTHOR')})")
            if args.output:
                with open(args.output, "w", encoding="utf-8") as f:
                    json.dump(res, f, indent=2)
                print(f"Saved to {args.output}")
        return

    parser.print_help()


if __name__ == "__main__":
    main()
