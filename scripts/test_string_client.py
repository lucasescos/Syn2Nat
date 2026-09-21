#!/usr/bin/env python3
"""Integration tests for the STRING database client."""

import unittest
from string_client import StringClient


class TestStringClient(unittest.TestCase):
    def setUp(self):
        self.client = StringClient(caller_identity="string_client_test", min_delay=0.5)

    def test_map_identifiers(self):
        df = self.client.map_identifiers(["TP53", "CDK2"], species=9606)
        self.assertGreaterEqual(len(df), 2)
        mapped_ids = set(df["preferredName"].tolist())
        self.assertIn("TP53", mapped_ids)
        self.assertIn("CDK2", mapped_ids)
        print(f"[PASS] Mapped IDs: {mapped_ids}")

    def test_network_tp53_mdm2(self):
        df = self.client.get_network(["TP53", "MDM2"], species=9606, required_score=700)
        self.assertGreater(len(df), 0)
        score = float(df.iloc[0]["score"])
        self.assertGreaterEqual(score, 0.9)
        print(f"[PASS] Network edge TP53 - MDM2 combined score: {score}")

    def test_interaction_partners(self):
        df = self.client.get_interaction_partners(["BRCA1"], species=9606, limit=5, required_score=700)
        self.assertGreaterEqual(len(df), 3)
        print(f"[PASS] Top BRCA1 interactors returned: {len(df)}")

    def test_enrichment(self):
        df = self.client.get_enrichment(["TP53", "MDM2", "CDKN1A"], species=9606)
        self.assertGreater(len(df), 0)
        self.assertIn("fdr", df.columns)
        print(f"[PASS] Enrichment terms returned: {len(df)}")


if __name__ == "__main__":
    unittest.main()
