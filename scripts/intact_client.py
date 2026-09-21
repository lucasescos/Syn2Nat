#!/usr/bin/env python3
"""
IntAct API Client
==================
Programmatic client for querying the IntAct Molecular Interaction Database
(EMBL-EBI / IMEx Consortium).

Web Portal: https://www.ebi.ac.uk/intact/
REST API: https://www.ebi.ac.uk/intact/ws/
Documentation: docs/intact.md
"""

from __future__ import annotations

import argparse
import json
import sys
import urllib.error
import urllib.parse
import urllib.request
from typing import Any, Dict, List, Optional, Union

INTACT_BASE_URL = "https://www.ebi.ac.uk/intact/ws"


class IntActClient:
    """Python client for IntAct REST microservices."""

    def __init__(self, base_url: str = INTACT_BASE_URL, timeout: int = 30):
        self.base_url = base_url.rstrip("/")
        self.timeout = timeout
        self.user_agent = "IntAct-PythonClient/1.0 (microproteinproject)"

    def _request(
        self,
        endpoint: str,
        method: str = "GET",
        params: Optional[Dict[str, Any]] = None,
        json_data: Optional[Dict[str, Any]] = None,
        headers: Optional[Dict[str, str]] = None,
        raw_text: bool = False,
    ) -> Any:
        """Executes an HTTP request to the IntAct API."""
        url = f"{self.base_url}/{endpoint.lstrip('/')}"
        if params:
            encoded_params = urllib.parse.urlencode(
                {k: v for k, v in params.items() if v is not None}
            )
            url = f"{url}?{encoded_params}"

        req_headers = {
            "User-Agent": self.user_agent,
            "Accept": "application/json" if not raw_text else "*/*",
        }
        if headers:
            req_headers.update(headers)

        body_bytes = None
        if json_data is not None:
            body_bytes = json.dumps(json_data).encode("utf-8")
            req_headers["Content-Type"] = "application/json"
        elif method == "POST":
            # Some POST endpoints in IntAct expect POST with query parameters and empty body
            body_bytes = b""

        req = urllib.request.Request(
            url, data=body_bytes, headers=req_headers, method=method
        )

        try:
            with urllib.request.urlopen(req, timeout=self.timeout) as response:
                content = response.read().decode("utf-8")
                if raw_text:
                    return content
                try:
                    return json.loads(content)
                except json.JSONDecodeError:
                    return content
        except urllib.error.HTTPError as err:
            err_msg = err.read().decode("utf-8", errors="replace")
            raise RuntimeError(f"IntAct API HTTP error {err.code}: {err.reason} - {err_msg}")
        except urllib.error.URLError as err:
            raise RuntimeError(f"IntAct API connection error: {err.reason}")

    # -------------------------------------------------------------------------
    # Interactor Service
    # -------------------------------------------------------------------------
    def find_interactor(
        self, query: str, page: int = 0, page_size: int = 10
    ) -> Dict[str, Any]:
        """Search interactors by name, symbol, or accession."""
        encoded_query = urllib.parse.quote(query)
        endpoint = f"interactor/findInteractor/{encoded_query}"
        return self._request(endpoint, method="GET", params={"page": page, "pageSize": page_size})

    def resolve_interactors(
        self, identifiers: Union[str, List[str]], fuzzy: bool = False
    ) -> Dict[str, Any]:
        """Batch resolve identifier strings to canonical IntAct interactor entries."""
        if isinstance(identifiers, list):
            query_str = ",".join(identifiers)
        else:
            query_str = identifiers

        endpoint = "interactor/list/resolve"
        params = {"query": query_str, "fuzzySearch": str(fuzzy).lower(), "page": 0, "pageSize": 50}
        return self._request(endpoint, method="POST", params=params)

    # -------------------------------------------------------------------------
    # Interaction Service
    # -------------------------------------------------------------------------
    def search_interactions(
        self,
        query: str,
        min_miscore: float = 0.0,
        max_miscore: float = 1.0,
        species: Optional[str] = None,
        expansion_filter: bool = False,
        negative_filter: str = "POSITIVE_ONLY",
        mutation_filter: bool = False,
        advanced_search: bool = True,
        page: int = 0,
        page_size: int = 20,
    ) -> Dict[str, Any]:
        """Search binary interactions with full filtering options.
        
        Args:
            query: MIQL expression or keyword query (e.g. 'geneName:TP53')
            min_miscore: Minimum MIscore (0.0 to 1.0)
            max_miscore: Maximum MIscore (0.0 to 1.0)
            species: Filter by species name or taxid (e.g. 'Homo sapiens' or '9606')
            expansion_filter: If True, exclude spoke-expanded complexes (true binary only)
            negative_filter: 'POSITIVE_ONLY', 'NEGATIVE_ONLY', or 'POSITIVE_AND_NEGATIVE'
            mutation_filter: If True, restrict to interactions affected by mutations
            advanced_search: If True, query is parsed as MIQL syntax
            page: 0-indexed page number
            page_size: Results per page
        """
        endpoint = "interaction/list"
        params: Dict[str, Any] = {
            "query": query,
            "advancedSearch": str(advanced_search).lower(),
            "minMIScore": min_miscore,
            "maxMIScore": max_miscore,
            "expansionFilter": str(expansion_filter).lower(),
            "negativeFilter": negative_filter,
            "mutationFilter": str(mutation_filter).lower(),
            "draw": 1,
            "page": page,
            "pageSize": page_size,
        }
        if species:
            params["interactorSpeciesFilter"] = species

        return self._request(endpoint, method="POST", params=params)

    def get_facets(self, query: str) -> Dict[str, Any]:
        """Fetch statistical facet breakdown for an interaction query."""
        endpoint = "interaction/findInteractionFacets"
        return self._request(endpoint, method="POST", params={"query": query})

    # -------------------------------------------------------------------------
    # Graph Details & Export Services
    # -------------------------------------------------------------------------
    def get_interaction_details(self, interaction_ac: str) -> Dict[str, Any]:
        """Fetch complete curation details for an interaction accession."""
        endpoint = f"graph/interaction/details/{interaction_ac}"
        return self._request(endpoint, method="GET")

    def get_participants(self, interaction_ac: str) -> Dict[str, Any]:
        """Fetch participant entities and experimental roles for an interaction."""
        endpoint = f"graph/participants/details/{interaction_ac}"
        return self._request(endpoint, method="GET")

    def get_features(self, interaction_ac: str) -> Dict[str, Any]:
        """Fetch participant features (mutations, binding regions, PTMs)."""
        endpoint = f"graph/features/details/{interaction_ac}"
        return self._request(endpoint, method="GET")

    def export_interaction(
        self, interaction_ac: str, export_format: str = "miTab27"
    ) -> str:
        """Export a single interaction in standard format (miTab27, miJSON, miXML25, etc.)."""
        endpoint = f"graph/export/interaction/{interaction_ac}"
        return self._request(endpoint, method="GET", params={"format": export_format}, raw_text=True)

    def export_interaction_list(
        self, query: str, export_format: str = "miTab27"
    ) -> str:
        """Batch export interactions matching query in standard format."""
        endpoint = "graph/export/interaction/list"
        return self._request(
            endpoint, method="POST", params={"query": query, "format": export_format}, raw_text=True
        )

    # -------------------------------------------------------------------------
    # Network Graph Service
    # -------------------------------------------------------------------------
    def get_network_graph(self, query: str, page_size: int = 50) -> Dict[str, Any]:
        """Fetch Cytoscape.js compatible network JSON."""
        endpoint = "network/getInteractions"
        params = {"query": query, "page": 0, "pageSize": page_size}
        return self._request(endpoint, method="POST", params=params)


# -----------------------------------------------------------------------------
# CLI Interface
# -----------------------------------------------------------------------------
def main():
    parser = argparse.ArgumentParser(
        description="IntAct Molecular Interaction Database CLI Client"
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    # Subcommand: interactor
    interactor_p = subparsers.add_parser("interactor", help="Search interactors")
    interactor_p.add_argument("query", help="Gene symbol, protein name, or accession")
    interactor_p.add_argument("--page-size", type=int, default=5, help="Number of results")

    # Subcommand: search
    search_p = subparsers.add_parser("search", help="Search binary interactions")
    search_p.add_argument("query", help="MIQL query or gene identifier (e.g. 'geneName:TP53')")
    search_p.add_argument("--min-score", type=float, default=0.0, help="Minimum MIscore (0.0-1.0)")
    search_p.add_argument("--species", help="Filter by species (e.g. 'Homo sapiens')")
    search_p.add_argument("--no-expansion", action="store_true", help="Exclude spoke expanded complexes")
    search_p.add_argument("--mutations-only", action="store_true", help="Only affected by mutations")
    search_p.add_argument("--page-size", type=int, default=10, help="Results per page")

    # Subcommand: export
    export_p = subparsers.add_parser("export", help="Export interactions in standard format")
    export_p.add_argument("query", help="Interaction AC (e.g. EBI-21978305) or query string")
    export_p.add_argument("--format", default="miTab27", choices=["miTab25", "miTab27", "miJSON", "featureTab", "miXML25", "miXML30"])
    export_p.add_argument("--output", help="Write export to destination file")

    # Subcommand: details
    details_p = subparsers.add_parser("details", help="Fetch detailed curation data for an interaction AC")
    details_p.add_argument("ac", help="Interaction AC (e.g. EBI-21978305)")

    args = parser.parse_args()
    client = IntActClient()

    try:
        if args.command == "interactor":
            res = client.find_interactor(args.query, page_size=args.page_size)
            items = res.get("content", [])
            print(f"Found {res.get('totalElements', 0)} interactors for '{args.query}':\n")
            for it in items:
                print(f"  * AC: {it.get('interactorAc')} | ID: {it.get('interactorPreferredIdentifier')} | Name: {it.get('interactorName')}")
                print(f"    Species: {it.get('interactorSpecies')} (taxid: {it.get('interactorTaxId')}) | Type: {it.get('interactorType')}")
                print(f"    Description: {it.get('interactorDescription')}\n")

        elif args.command == "search":
            res = client.search_interactions(
                query=args.query,
                min_miscore=args.min_score,
                species=args.species,
                expansion_filter=args.no_expansion,
                mutation_filter=args.mutations_only,
                page_size=args.page_size,
            )
            items = res.get("data", res.get("content", []))
            total = res.get("recordsFiltered", res.get("totalElements", len(items)))
            print(f"Found {total} interactions matching '{args.query}':\n")
            for it in items:
                print(f"  * {it.get('moleculeA')} ({it.get('idA')}) <-> {it.get('moleculeB')} ({it.get('idB')})")
                print(f"    AC: {it.get('ac')} | MIscore: {it.get('intactMiscore')} | Type: {it.get('type')}")
                print(f"    Detection: {it.get('detectionMethod')} | PMID: {it.get('publicationPubmedIdentifier')}")
                if it.get("expansionMethod"):
                    print(f"    Expansion: {it.get('expansionMethod')}")
                print()

        elif args.command == "export":
            if args.query.startswith("EBI-"):
                text = client.export_interaction(args.query, export_format=args.format)
            else:
                text = client.export_interaction_list(args.query, export_format=args.format)

            if args.output:
                with open(args.output, "w", encoding="utf-8") as f:
                    f.write(text)
                print(f"Export written to {args.output}")
            else:
                print(text[:2000])
                if len(text) > 2000:
                    print(f"\n... [truncated, {len(text)} bytes total] ...")

        elif args.command == "details":
            details = client.get_interaction_details(args.ac)
            print(json.dumps(details, indent=2))

    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
