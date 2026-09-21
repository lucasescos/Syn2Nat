#!/usr/bin/env python3
"""STRING Database API Client (STRING v12.0).

===================================
Programmatic Python client and CLI tool for querying the STRING database:
- Identifier resolution and disambiguation (`get_string_ids`)
- Network retrieval and topological neighborhoods (`network`)
- Proteome-wide interaction partner extraction (`interaction_partners`)
- Pairwise and cross-species sequence homology (`homology`, `homology_best`)
- Over-representation functional enrichment (`enrichment`, `ppi_enrichment`)
- Reverse term searching and automated gene set descriptions (`functional_terms`, `geneset_description`)
- Network and enrichment figure generation

Web portal: https://string-db.org
API docs: docs/string.md
"""

from __future__ import annotations

import argparse
import io
import json
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
from typing import Any, Dict, Iterable, List, Optional, Union

try:
    import pandas as pd
    HAS_PANDAS = True
except ImportError:
    HAS_PANDAS = False

STRING_API_VERSION = "12-0"
DEFAULT_BASE_URL = f"https://version-{STRING_API_VERSION}.string-db.org/api"


class StringClient:
    """Client for the STRING Database REST API (v12.0)."""

    def __init__(
        self,
        base_url: str = DEFAULT_BASE_URL,
        caller_identity: str = "microprotein_project",
        min_delay: float = 1.0,
        timeout: int = 60,
    ):
        self.base_url = base_url.rstrip("/")
        self.caller_identity = caller_identity
        self.min_delay = min_delay
        self.timeout = timeout
        self._last_request_time = 0.0

    def _rate_limit(self) -> None:
        """Enforce polite rate limiting between consecutive calls."""
        elapsed = time.time() - self._last_request_time
        if elapsed < self.min_delay:
            time.sleep(self.min_delay - elapsed)

    def _request(
        self,
        endpoint: str,
        params: Dict[str, Any],
        format_type: str = "tsv",
    ) -> bytes:
        """Execute a POST request to the STRING API."""
        self._rate_limit()
        url = f"{self.base_url}/{format_type}/{endpoint}"

        # Standard caller identification
        params["caller_identity"] = self.caller_identity

        # Encode body
        data = urllib.parse.urlencode(params).encode("utf-8")
        req = urllib.request.Request(
            url,
            data=data,
            headers={
                "User-Agent": f"StringClient/12.0 ({self.caller_identity})",
                "Content-Type": "application/x-www-form-urlencoded",
            },
        )

        try:
            with urllib.request.urlopen(req, timeout=self.timeout) as response:
                content = response.read()
                self._last_request_time = time.time()
                return content
        except urllib.error.HTTPError as err:
            self._last_request_time = time.time()
            if err.code == 404:
                raise ValueError(
                    f"STRING API returned 404: Identifiers could not be resolved in species {params.get('species')}."
                ) from err
            error_body = err.read().decode("utf-8", errors="replace")
            raise RuntimeError(f"STRING API error (HTTP {err.code}): {error_body}") from err
        except urllib.error.URLError as err:
            self._last_request_time = time.time()
            raise RuntimeError(f"Network error connecting to STRING: {err.reason}") from err

    @staticmethod
    def _parse_tsv(raw_bytes: bytes) -> Union["pd.DataFrame", List[Dict[str, str]]]:
        """Convert TSV bytes to pandas DataFrame if available, else list of dicts."""
        text = raw_bytes.decode("utf-8")
        if HAS_PANDAS:
            return pd.read_csv(io.StringIO(text), sep="\t")
        lines = [line.strip().split("\t") for line in text.splitlines() if line.strip()]
        if not lines:
            return []
        headers = lines[0]
        return [dict(zip(headers, row)) for row in lines[1:]]

    def map_identifiers(
        self,
        identifiers: Iterable[str],
        species: int = 9606,
        echo_query: bool = True,
    ) -> Union["pd.DataFrame", List[Dict[str, str]]]:
        """Map common gene symbols, UniProt IDs, or synonyms to STRING IDs."""
        id_list = list(identifiers)
        params = {
            "identifiers": "\r".join(id_list),
            "species": species,
            "echo_query": 1 if echo_query else 0,
        }
        res = self._request("get_string_ids", params, format_type="tsv")
        return self._parse_tsv(res)

    def get_network(
        self,
        identifiers: Iterable[str],
        species: int = 9606,
        required_score: int = 400,
        network_type: str = "functional",
        add_nodes: int = 0,
    ) -> Union["pd.DataFrame", List[Dict[str, str]]]:
        """Retrieve interaction network between input proteins."""
        params = {
            "identifiers": "\r".join(identifiers),
            "species": species,
            "required_score": required_score,
            "network_type": network_type,
            "add_nodes": add_nodes,
        }
        res = self._request("network", params, format_type="tsv")
        return self._parse_tsv(res)

    def get_interaction_partners(
        self,
        identifiers: Iterable[str],
        species: int = 9606,
        limit: int = 10,
        required_score: int = 400,
        network_type: str = "functional",
    ) -> Union["pd.DataFrame", List[Dict[str, str]]]:
        """Retrieve top interaction partners across the proteome for input proteins."""
        params = {
            "identifiers": "\r".join(identifiers),
            "species": species,
            "limit": limit,
            "required_score": required_score,
            "network_type": network_type,
        }
        res = self._request("interaction_partners", params, format_type="tsv")
        return self._parse_tsv(res)

    def get_homology(
        self,
        identifiers: Iterable[str],
        species: int = 9606,
    ) -> Union["pd.DataFrame", List[Dict[str, str]]]:
        """Retrieve pairwise Smith-Waterman bit scores between input proteins."""
        params = {
            "identifiers": "\r".join(identifiers),
            "species": species,
        }
        res = self._request("homology", params, format_type="tsv")
        return self._parse_tsv(res)

    def get_homology_best(
        self,
        identifiers: Iterable[str],
        species: int = 9606,
        species_b: Optional[Iterable[int]] = None,
    ) -> Union["pd.DataFrame", List[Dict[str, str]]]:
        """Retrieve best homology hits between input proteins and target species."""
        params: Dict[str, Any] = {
            "identifiers": "\r".join(identifiers),
            "species": species,
        }
        if species_b:
            params["species_b"] = "%0d".join(str(s) for s in species_b)
        res = self._request("homology_best", params, format_type="tsv")
        return self._parse_tsv(res)

    def get_enrichment(
        self,
        identifiers: Iterable[str],
        species: int = 9606,
        background_ids: Optional[Iterable[str]] = None,
    ) -> Union["pd.DataFrame", List[Dict[str, str]]]:
        """Perform functional enrichment analysis (GO, KEGG, Pfam, etc.)."""
        params: Dict[str, Any] = {
            "identifiers": "\r".join(identifiers),
            "species": species,
        }
        if background_ids:
            params["background_string_identifiers"] = "%0d".join(background_ids)
        res = self._request("enrichment", params, format_type="tsv")
        return self._parse_tsv(res)

    def get_functional_annotation(
        self,
        identifiers: Iterable[str],
        species: int = 9606,
        allow_pubmed: bool = False,
    ) -> Union["pd.DataFrame", List[Dict[str, str]]]:
        """Retrieve all mapped functional annotations for the given proteins."""
        params = {
            "identifiers": "\r".join(identifiers),
            "species": species,
            "allow_pubmed": 1 if allow_pubmed else 0,
        }
        res = self._request("functional_annotation", params, format_type="tsv")
        return self._parse_tsv(res)

    def get_functional_terms(
        self,
        term_text: str,
        species: int = 9606,
    ) -> Union["pd.DataFrame", List[Dict[str, str]]]:
        """Retrieve proteins annotated with a functional term or disease name."""
        params = {
            "term_text": term_text,
            "species": species,
        }
        res = self._request("functional_terms", params, format_type="tsv")
        return self._parse_tsv(res)

    def get_geneset_description(
        self,
        identifiers: Iterable[str],
        species: int = 9606,
    ) -> Union["pd.DataFrame", List[Dict[str, str]]]:
        """Retrieve high-level biological theme descriptions for a gene set."""
        params = {
            "identifiers": "\r".join(identifiers),
            "species": species,
        }
        res = self._request("geneset_description", params, format_type="tsv")
        return self._parse_tsv(res)

    def get_ppi_enrichment(
        self,
        identifiers: Iterable[str],
        species: int = 9606,
        required_score: int = 400,
    ) -> Union["pd.DataFrame", List[Dict[str, str]]]:
        """Perform PPI enrichment statistical testing."""
        params = {
            "identifiers": "\r".join(identifiers),
            "species": species,
            "required_score": required_score,
        }
        res = self._request("ppi_enrichment", params, format_type="tsv")
        return self._parse_tsv(res)

    def download_network_image(
        self,
        identifiers: Iterable[str],
        output_path: str,
        species: int = 9606,
        format_type: str = "highres_image",
        network_flavor: str = "confidence",
        add_color_nodes: int = 0,
        add_white_nodes: int = 0,
    ) -> str:
        """Download network figure (PNG or SVG)."""
        params = {
            "identifiers": "\r".join(identifiers),
            "species": species,
            "network_flavor": network_flavor,
        }
        if add_color_nodes > 0:
            params["add_color_nodes"] = add_color_nodes
        if add_white_nodes > 0:
            params["add_white_nodes"] = add_white_nodes

        content = self._request("network", params, format_type=format_type)
        with open(output_path, "wb") as f:
            f.write(content)
        return output_path


def main() -> int:
    """Command-line interface entry point."""
    parent = argparse.ArgumentParser(add_help=False)
    parent.add_argument("--species", type=int, default=9606, help="NCBI Taxon ID (default: 9606 human)")
    parent.add_argument("--caller", default="string_cli", help="Caller identity token")

    parser = argparse.ArgumentParser(description="STRING Database API Client", parents=[parent])
    subparsers = parser.add_subparsers(dest="command", required=True)

    # map
    p_map = subparsers.add_parser("map", parents=[parent], help="Map identifiers to STRING IDs")
    p_map.add_argument("identifiers", nargs="+", help="Protein names or accessions")

    # network
    p_net = subparsers.add_parser("network", parents=[parent], help="Get interaction network")
    p_net.add_argument("identifiers", nargs="+", help="Protein identifiers")
    p_net.add_argument("--score", type=int, default=400, help="Score threshold (0-1000)")
    p_net.add_argument("--type", choices=["functional", "physical"], default="functional")
    p_net.add_argument("--add-nodes", type=int, default=0)

    # partners
    p_part = subparsers.add_parser("partners", parents=[parent], help="Get top interaction partners")
    p_part.add_argument("identifiers", nargs="+", help="Protein identifiers")
    p_part.add_argument("--limit", type=int, default=10, help="Max interactors per protein")
    p_part.add_argument("--score", type=int, default=400)

    # enrichment
    p_enr = subparsers.add_parser("enrichment", parents=[parent], help="Functional pathway enrichment")
    p_enr.add_argument("identifiers", nargs="+", help="Protein identifiers")

    # terms
    p_term = subparsers.add_parser("terms", parents=[parent], help="Find proteins matching a term or disease")
    p_term.add_argument("term", help="Search string (e.g. 'Melanoma')")

    args = parser.parse_args()
    client = StringClient(caller_identity=args.caller)

    try:
        if args.command == "map":
            res = client.map_identifiers(args.identifiers, species=args.species)
        elif args.command == "network":
            res = client.get_network(
                args.identifiers,
                species=args.species,
                required_score=args.score,
                network_type=args.type,
                add_nodes=args.add_nodes,
            )
        elif args.command == "partners":
            res = client.get_interaction_partners(
                args.identifiers,
                species=args.species,
                limit=args.limit,
                required_score=args.score,
            )
        elif args.command == "enrichment":
            res = client.get_enrichment(args.identifiers, species=args.species)
        elif args.command == "terms":
            res = client.get_functional_terms(args.term, species=args.species)
        else:
            parser.print_help()
            return 1

        if HAS_PANDAS and isinstance(res, pd.DataFrame):
            print(res.to_string(index=False))
        else:
            print(json.dumps(res, indent=2))
        return 0

    except Exception as exc:
        print(f"Error: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    sys.exit(main())
