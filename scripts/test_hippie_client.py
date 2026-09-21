#!/usr/bin/env python3
"""
Unit tests for HIPPIE client.
"""

import unittest
from hippie_client import HippieClient, query_pair, get_interactors

class TestHippieClient(unittest.TestCase):
    def setUp(self):
        self.client = HippieClient()

    def test_query_pair_tp53_mdm2(self):
        result = query_pair("TP53", "MDM2")
        self.assertIsNotNone(result)
        self.assertIn("score", result)
        self.assertGreaterEqual(result["score"], 0.9)
        print(f"[PASS] TP53 - MDM2 score: {result['score']}")

    def test_query_single_protein_layer1(self):
        df = self.client.query("TP53", layer=1, conf_thres=0.9)
        self.assertGreater(len(df), 5)
        print(f"[PASS] TP53 high-confidence (>=0.9) interactors count: {len(df)}")

    def test_query_layer0_subnetwork(self):
        genes = ["TP53", "MDM2", "CDKN1A", "ATM"]
        df = self.client.query_subnetwork(genes, conf_thres=0.7)
        self.assertGreater(len(df), 0)
        print(f"[PASS] Layer 0 subnetwork edges for {genes}: {len(df)}")

    def test_query_mitab(self):
        df = self.client.query("TP53", layer=1, conf_thres=0.95, out_type="mitab")
        self.assertGreater(len(df), 0)
        self.assertIn("confidence_score", df.columns)
        print(f"[PASS] MITAB format returned {len(df)} rows with columns: {list(df.columns[:4])}")

if __name__ == "__main__":
    unittest.main()
