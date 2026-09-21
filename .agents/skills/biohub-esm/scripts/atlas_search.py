#!/usr/bin/env python3
"""
ESM Atlas search & lookup script (public, keyless API).
Queries https://biohub.ai/esm/protein/api/v1alpha1 for similarity searches, protein hashes,
UniProt mapping, clusters, thumbnails, batch queries, and feature dossiers.
"""

import argparse
import hashlib
import json
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path

ATLAS_BASE = "https://biohub.ai/esm/protein/api/v1alpha1"


def query_atlas(endpoint: str, params: dict = None, method: str = "GET", data: bytes = None, is_json: bool = True) -> dict | bytes:
    url = f"{ATLAS_BASE}/{endpoint.lstrip('/')}"
    if params:
        query_string = urllib.parse.urlencode({k: v for k, v in params.items() if v is not None})
        url = f"{url}?{query_string}"

    headers = {"User-Agent": "Antigravity-BiohubESM-Skill/1.0"}
    if is_json:
        headers["Accept"] = "application/json"
    if data:
        headers["Content-Type"] = "application/json"

    req = urllib.request.Request(url, data=data, headers=headers, method=method)
    try:
        with urllib.request.urlopen(req, timeout=45) as resp:
            content = resp.read()
            if is_json:
                return json.loads(content.decode("utf-8"))
            return content
    except urllib.error.HTTPError as e:
        err = e.read().decode("utf-8", errors="replace")
        print(f"[ERROR] HTTP {e.code} for {url}: {err}", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"[ERROR] Failed to connect to Atlas API: {e}", file=sys.stderr)
        sys.exit(1)


def read_sequence(input_str: str) -> str:
    p = Path(input_str)
    if p.is_file():
        seq_lines = []
        with open(p, "r", encoding="utf-8") as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith(">"):
                    continue
                seq_lines.append(line)
        return "".join(seq_lines).upper()
    return "".join(input_str.split()).upper()


def main():
    parser = argparse.ArgumentParser(description="Query the public Biohub ESM Atlas API")
    subparsers = parser.add_subparsers(dest="command", required=True)

    # Subcommand: similarity
    sim_p = subparsers.add_parser("similarity", help="Search structurally and functionally similar proteins in the Atlas")
    sim_p.add_argument("sequence", help="Amino acid sequence or FASTA file")
    sim_p.add_argument("--topk", type=int, default=5, help="Number of similar proteins to return (default: 5)")
    sim_p.add_argument("--include-cluster", action="store_true", default=True, help="Include cluster metadata")

    # Subcommand: protein
    prot_p = subparsers.add_parser("protein", help="Retrieve protein metadata and SAE features by sequence hash")
    prot_p.add_argument("--hash", help="MD5 hex digest of the sequence")
    prot_p.add_argument("--sequence", help="Protein sequence (computes hash automatically)")
    prot_p.add_argument("--topk-features", type=int, default=10, help="Number of top features to return")
    prot_p.add_argument("--fold-on-miss", action="store_true", default=True, help="Fold on-the-fly with ESMFold2 if missing (<700 aa)")

    # Subcommand: uniprot
    uni_p = subparsers.add_parser("uniprot", help="Map UniProt Accession to Atlas metadata & protein hash")
    uni_p.add_argument("accession", help="UniProt Accession ID (e.g. 'P00520')")

    # Subcommand: feature
    feat_p = subparsers.add_parser("feature", help="Retrieve feature dossier from the 16,384 SAE dictionary")
    feat_p.add_argument("index", type=int, help="Feature index (0 to 16383)")

    # Subcommand: bulk-features
    bulk_p = subparsers.add_parser("bulk-features", help="Fetch all 16,384 SAE feature dossiers in one unpaginated call")
    bulk_p.add_argument("-o", "--output", default="all_sae_features.json", help="Output JSON filename")

    # Subcommand: cluster
    clust_p = subparsers.add_parser("cluster", help="Inspect an Atlas cluster by representative protein hash")
    clust_p.add_argument("hash", help="Cluster representative protein hash")

    # Subcommand: thumbnail
    thumb_p = subparsers.add_parser("thumbnail", help="Download PNG structure thumbnail by protein hash")
    thumb_p.add_argument("hash", help="Protein MD5 hash")
    thumb_p.add_argument("-o", "--output", default="thumbnail.png", help="Output PNG filename")
    thumb_p.add_argument("--type", default="plddt", choices=["plddt", "pct-characterized"], help="Thumbnail coloration type")

    # Subcommand: batch
    batch_p = subparsers.add_parser("batch", help="Submit batch queries for multiple protein hashes (up to 500)")
    batch_p.add_argument("hashes_file", help="Text file with one protein hash per line, or comma-separated hashes")
    batch_p.add_argument("-o", "--output", default="batch_results.zip", help="Output zip filename")
    batch_p.add_argument("--include-structure", action="store_true", help="Include predicted structure CIF/PDB files in batch")

    args = parser.parse_args()

    if args.command == "similarity":
        seq = read_sequence(args.sequence)
        print(f"Searching ESM Atlas for sequence of {len(seq)} aa (top {args.topk} hits)...")
        res = query_atlas("similarity-search", {
            "sequence": seq,
            "topk_results": args.topk,
            "include_cluster_info": args.include_cluster
        })

        similar = res.get("similar_proteins", [])
        print(f"\nFound {len(similar)} similar proteins:")
        print("-" * 75)
        print(f"{'Accession':<16} {'Similarity':<12} {'Length':<8} {'Mean pLDDT':<12} {'Name'}")
        print("-" * 75)
        for p in similar:
            acc = p.get("protein_accession") or p.get("protein_hash")[:12]
            sim = p.get("similarity_score", 0.0)
            length = p.get("sequence_length", "-")
            plddt = p.get("mean_plddt")
            plddt_str = f"{plddt:.1f}" if plddt is not None else "N/A"
            name = (p.get("protein_name") or "Uncharacterized")[:30]
            print(f"{acc:<16} {sim:<12.4f} {length:<8} {plddt_str:<12} {name}")

        top_feats = res.get("top_features_across_results", [])
        if top_feats:
            print("\nTop Firing SAE Features Across Results:")
            for f in top_feats[:5]:
                f_idx = f.get("feature_index")
                count = f.get("occurrence_count")
                mean_act = f.get("mean_activation", 0.0)
                print(f"  Feature #{f_idx}: fired in {count} hits (mean act: {mean_act:.2f})")

    elif args.command == "uniprot":
        acc = args.accession.strip()
        print(f"Resolving UniProt ID '{acc}' in ESM Atlas...")
        res = query_atlas(f"uniprot/{acc}")
        print("\n--- UniProt Mapping Result ---")
        print(f"Accession:    {res.get('accession')}")
        print(f"Protein Hash: {res.get('protein_hash')}")
        print(f"Gene:         {res.get('gene', 'N/A')}")
        print(f"Organism:     {res.get('organism', 'N/A')}")
        print(f"Name:         {res.get('protein_name', 'N/A')}")
        print(f"Length:       {res.get('sequence_length')} aa")
        if res.get("function"):
            print(f"\nFunction:     {res.get('function')[:300]}...")

    elif args.command == "protein":
        h = args.hash
        if not h and args.sequence:
            seq = read_sequence(args.sequence)
            h = hashlib.md5(seq.encode("utf-8")).hexdigest()
        if not h:
            print("[ERROR] Specify either --hash or --sequence.", file=sys.stderr)
            sys.exit(1)

        print(f"Looking up protein hash: {h}...")
        res = query_atlas(f"proteins/{h}", {
            "topk_features": args.topk_features,
            "fold_on_miss": str(args.fold_on_miss).lower()
        })
        print("\n--- Protein Atlas Details ---")
        print(f"Accession:    {res.get('protein_accession', 'N/A')}")
        print(f"Cluster Rep:  {res.get('cluster_rep_protein_hash', 'N/A')}")
        print(f"Mean pLDDT:   {res.get('mean_plddt', 'N/A')}")
        print(f"Length:       {res.get('sequence_length', 'N/A')}")

        feats = res.get("features", [])
        if feats:
            print("\nTop Features:")
            for f in feats[:args.topk_features]:
                idx = f.get("feature_index")
                score = f.get("normalized_score") or f.get("activation", 0.0)
                print(f"  Feature #{idx}: {score:.3f}")

    elif args.command == "feature":
        idx = args.index
        print(f"Fetching dossier for SAE Feature #{idx}...")
        res = query_atlas(f"features/{idx}")
        print("\n--- SAE Feature Dossier ---")
        print(f"Index:       #{idx}")
        print(f"Label:       {res.get('label', 'N/A')}")
        print(f"Category:    {res.get('category', 'N/A')}")
        print(f"Summary:     {res.get('summary', 'N/A')}")
        print(f"Threshold:   {res.get('threshold', 'N/A')}")
        if "top_swissprot_activations" in res:
            print("\nTop SwissProt Activations:")
            for sp in res["top_swissprot_activations"][:5]:
                print(f"  {sp.get('uniprot_id')}: act={sp.get('activation', 0.0):.2f} - {sp.get('protein_name')}")

    elif args.command == "bulk-features":
        print("Fetching all 16,384 SAE features from Atlas (single unpaginated call)...")
        res = query_atlas("features")
        out_path = Path(args.output).resolve()
        with open(out_path, "w", encoding="utf-8") as f:
            json.dump(res, f, indent=2)
        print(f"Saved {len(res) if isinstance(res, list) else 'feature dictionary'} to: {out_path}")

    elif args.command == "cluster":
        h = args.hash
        print(f"Fetching cluster for representative hash: {h}...")
        res = query_atlas(f"clusters/{h}")
        print("\n--- Cluster Summary ---")
        print(f"Representative: {h}")
        print(f"Size:           {res.get('cluster_size', 'N/A')}")
        print(f"Characterized:  {res.get('characterized_ratio', 'N/A')}")

    elif args.command == "thumbnail":
        h = args.hash
        print(f"Fetching {args.type} thumbnail for {h}...")
        content = query_atlas(f"proteins/{h}/thumbnail/{args.type}", is_json=False)
        out_path = Path(args.output).resolve()
        with open(out_path, "wb") as f:
            f.write(content)
        print(f"Saved thumbnail ({len(content)} bytes) to: {out_path}")

    elif args.command == "batch":
        hf = Path(args.hashes_file)
        if hf.is_file():
            hashes = [line.strip() for line in open(hf) if line.strip()]
        else:
            hashes = [x.strip() for x in args.hashes_file.split(",") if x.strip()]

        if len(hashes) > 500:
            print("[WARN] Truncating batch to maximum allowed 500 hashes.", file=sys.stderr)
            hashes = hashes[:500]

        print(f"Submitting batch request for {len(hashes)} protein hashes...")
        body = {
            "protein_hashes": hashes,
            "topk_features": 10,
            "include_structure": args.include_structure,
            "include_cluster_info": True,
            "include_sequence": True,
            "include_features": {"protein_level": True, "per_residue": False}
        }
        raw_res = query_atlas("proteins/batch", method="POST", data=json.dumps(body).encode("utf-8"), is_json=False)

        # Check if response is immediate zip (200) or job status (202)
        if raw_res.startswith(b"PK"):
            out_path = Path(args.output).resolve()
            with open(out_path, "wb") as f:
                f.write(raw_res)
            print(f"Saved immediate batch zip ({len(raw_res)} bytes) to: {out_path}")
        else:
            try:
                job_info = json.loads(raw_res.decode("utf-8"))
                job_id = job_info.get("job_id") or job_info.get("id")
                print(f"Batch job queued (ID: {job_id}). Polling status...")
                while True:
                    time.sleep(3)
                    poll_res = query_atlas(f"proteins/batch/jobs/{job_id}", is_json=False)
                    if poll_res.startswith(b"PK"):
                        out_path = Path(args.output).resolve()
                        with open(out_path, "wb") as f:
                            f.write(poll_res)
                        print(f"Batch completed! Saved archive to: {out_path}")
                        break
                    else:
                        st = json.loads(poll_res.decode("utf-8"))
                        print(f"  Status: {st.get('status', 'running')}...")
                        if st.get("status") in ("failed", "error", "expired"):
                            print(f"[ERROR] Batch job failed: {st}", file=sys.stderr)
                            sys.exit(1)
            except Exception as e:
                print(f"[ERROR] Failed handling batch response: {e}", file=sys.stderr)
                sys.exit(1)


if __name__ == "__main__":
    main()
