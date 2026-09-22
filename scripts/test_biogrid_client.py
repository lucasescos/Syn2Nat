#!/usr/bin/env python3
"""
Unit and smoke tests for BioGRID Client.
"""

import json
import unittest
from unittest.mock import MagicMock, patch
from scripts.biogrid_client import BioGridClient, TAB2_COLUMNS


class TestBioGridClient(unittest.TestCase):

    def test_missing_access_key_raises(self):
        """Verify that an explicit error is raised if no access key is configured."""
        client = BioGridClient(access_key=None)
        # Ensure environment doesn't have a key during this test
        client.access_key = None
        with self.assertRaises(ValueError) as ctx:
            client.get_version()
        self.assertIn("BioGRID Access Key not found", str(ctx.exception))

    def test_tab2_columns_completeness(self):
        """Verify the 24 standard Tab 2.0 columns are defined."""
        self.assertEqual(len(TAB2_COLUMNS), 24)
        self.assertIn("BIOGRID_INTERACTION_ID", TAB2_COLUMNS)
        self.assertIn("OFFICIAL_SYMBOL_A", TAB2_COLUMNS)
        self.assertIn("OFFICIAL_SYMBOL_B", TAB2_COLUMNS)
        self.assertIn("EXPERIMENTAL_SYSTEM", TAB2_COLUMNS)

    @patch.object(BioGridClient, "_get")
    def test_count_query(self, mock_get):
        """Verify count endpoint query formulation and integer parsing."""
        mock_get.return_value = " 42 \n"
        client = BioGridClient(access_key="dummy_key_0123456789abcdef0123")
        count = client.count(genes=["TP53"], tax_id=9606)
        self.assertEqual(count, 42)
        mock_get.assert_called_once()
        args, kwargs = mock_get.call_args
        endpoint, params = args
        self.assertEqual(endpoint, "/interactions/")
        self.assertEqual(params["format"], "count")
        self.assertEqual(params["geneList"], "TP53")
        self.assertEqual(params["taxId"], 9606)

    @patch.object(BioGridClient, "_get")
    def test_query_pair_parsing(self, mock_get):
        """Verify pairwise interaction extraction."""
        mock_response = {
            "101": {
                "BIOGRID_INTERACTION_ID": 101,
                "OFFICIAL_SYMBOL_A": "TP53",
                "OFFICIAL_SYMBOL_B": "MDM2",
                "EXPERIMENTAL_SYSTEM": "Two-hybrid",
                "EXPERIMENTAL_SYSTEM_TYPE": "physical",
                "PUBMED_AUTHOR": "Momand J (1992)",
                "PUBMED_ID": "1533783"
            }
        }
        mock_get.return_value = json.dumps(mock_response)
        client = BioGridClient(access_key="dummy_key_0123456789abcdef0123")
        hits = client.query_pair("TP53", "MDM2", tax_id=9606, exp_type="physical")
        self.assertEqual(len(hits), 1)
        self.assertEqual(hits[0]["OFFICIAL_SYMBOL_A"], "TP53")
        self.assertEqual(hits[0]["OFFICIAL_SYMBOL_B"], "MDM2")
        self.assertEqual(hits[0]["EXPERIMENTAL_SYSTEM"], "Two-hybrid")

    @patch.object(BioGridClient, "_get")
    def test_query_subnetwork_layer0(self, mock_get):
        """Verify subnetwork requests set includeInteractors=false."""
        mock_get.return_value = json.dumps({})
        client = BioGridClient(access_key="dummy_key_0123456789abcdef0123")
        client.query_subnetwork(["TP53", "MDM2", "CDKN1A"], tax_id=9606)
        args, kwargs = mock_get.call_args
        endpoint, params = args
        self.assertEqual(params["includeInteractors"], "false")
        self.assertEqual(params["geneList"], "TP53|MDM2|CDKN1A")

    def test_live_api_queries(self):
        """Live smoke test against BioGRID webservice if key is available in environment."""
        client = BioGridClient()
        if not client.access_key:
            self.skipTest("No BIOGRID_ACCESS_KEY found; skipping live query test.")
        
        # Test version endpoint
        version = client.get_version()
        self.assertTrue(len(version) > 0)
        print(f"\n[LIVE TEST] BioGRID REST version: {version}")

        # Test count endpoint
        count = client.count(genes="TP53", tax_id=9606)
        self.assertGreater(count, 100)
        print(f"[LIVE TEST] TP53 total interactions: {count}")

        # Test pairwise query
        hits = client.query_pair("TP53", "MDM2", tax_id=9606, exp_type="physical")
        self.assertGreater(len(hits), 0)
        print(f"[LIVE TEST] TP53 - MDM2 physical interactions count: {len(hits)}")


if __name__ == "__main__":
    unittest.main()
