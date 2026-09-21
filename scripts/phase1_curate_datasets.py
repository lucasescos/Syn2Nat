#!/usr/bin/env python3
"""
Phase 1: Data Curation & Target Selection
Pipeline for the GENCODE Peptidein Deorphanization Framework.

Curates:
1. Microprotein/Peptidein Query Sets:
   - Tier 1 Peptideins (84 targets: 12 Tier 1A + 72 Tier 1B)
   - All Annotated Peptideins (121 targets: Tiers 1A, 1B, 2A, 2B)
   - Tier 1 Microprotein Catalog (711 targets: 18 Tier 1A + 693 Tier 1B)
2. Canonical Human Reference Proteome:
   - UniProt UP000005640 (canonical reviewed proteins, ~20,416 targets)
   - Generates indexed FASTA and Parquet tables for Phase 2 pre-filter screening.
"""

import argparse
import gzip
import os
import re
import sys
import time
from pathlib import Path
import pandas as pd
import requests


BASE_DIR = Path(__file__).resolve().parent.parent
DATA_DIR = BASE_DIR / "data"
PARSED_DIR = DATA_DIR / "parsed"
QUERIES_DIR = DATA_DIR / "queries"
TARGETS_DIR = DATA_DIR / "targets"
TARGETS_RAW_DIR = TARGETS_DIR / "raw"

UNIPROT_FTP_URL = (
    "https://ftp.uniprot.org/pub/databases/uniprot/current_release/"
    "knowledgebase/reference_proteomes/Eukaryota/UP000005640/UP000005640_9606.fasta.gz"
)


def write_fasta(records, output_path: Path, wrap: int = 60):
    """Write list of (header, sequence) tuples to FASTA with line wrapping."""
    with open(output_path, "w", encoding="utf-8") as f:
        for header, seq in records:
            f.write(f">{header}\n")
            for i in range(0, len(seq), wrap):
                f.write(f"{seq[i:i+wrap]}\n")


def curate_queries():
    """Extract and curate microprotein / peptidein query sets."""
    print("=" * 70)
    print("Phase 1: Curating Microprotein / Peptidein Query Datasets")
    print("=" * 70)

    orf_tiers_path = PARSED_DIR / "orf_tiers.parquet"
    orfs_path = PARSED_DIR / "orfs.parquet"

    if not orf_tiers_path.exists() or not orfs_path.exists():
        raise FileNotFoundError(
            f"Required parsed input files not found in {PARSED_DIR}. "
            f"Expected {orf_tiers_path.name} and {orfs_path.name}."
        )

    print(f"Loading {orf_tiers_path.name} and {orfs_path.name}...")
    df_tiers = pd.read_parquet(orf_tiers_path)
    df_orfs = pd.read_parquet(orfs_path)

    # Merge tier metadata with sequences and quality metrics
    merged = pd.merge(
        df_tiers,
        df_orfs[["orf_id", "sequence", "aa_valid", "looks_like_dna", "illegal_chars", "is_queryable"]],
        on="orf_id",
        how="inner",
    )
    print(f"Total merged ORFs in catalog: {len(merged):,}")

    QUERIES_DIR.mkdir(parents=True, exist_ok=True)

    # 1. Tier 1 Peptideins (12 Tier 1A + 72 Tier 1B = 84 targets)
    df_tier1_peps = merged[
        merged["is_peptidein"] & merged["final_tier"].isin(["1A", "1B"])
    ].copy().sort_values(["final_tier", "orf_id"])

    # 2. All Peptideins (121 targets)
    df_all_peps = merged[
        merged["is_peptidein"]
    ].copy().sort_values(["final_tier", "orf_id"])

    # 3. All Tier 1 Microproteins (18 Tier 1A + 693 Tier 1B = 711 targets)
    df_tier1_all = merged[
        merged["final_tier"].isin(["1A", "1B"])
    ].copy().sort_values(["final_tier", "orf_id"])

    # Build FASTA records
    def to_fasta_records(df):
        records = []
        for _, row in df.iterrows():
            hdr = (
                f"{row['orf_id']} gene={row['gene']} tier={row['final_tier']} "
                f"peptidein={'yes' if row['is_peptidein'] else 'no'} "
                f"length={row['length']} biotype={row['orf_biotype']}"
            )
            records.append((hdr, row["sequence"]))
        return records

    p1 = QUERIES_DIR / "tier1_peptideins.fasta"
    p2 = QUERIES_DIR / "all_peptideins.fasta"
    p3 = QUERIES_DIR / "tier1_all_microproteins.fasta"

    write_fasta(to_fasta_records(df_tier1_peps), p1)
    write_fasta(to_fasta_records(df_all_peps), p2)
    write_fasta(to_fasta_records(df_tier1_all), p3)

    print(f"Saved: {p1} ({len(df_tier1_peps)} sequences)")
    print(f"Saved: {p2} ({len(df_all_peps)} sequences)")
    print(f"Saved: {p3} ({len(df_tier1_all)} sequences)")

    # Save comprehensive metadata
    meta_df = merged[
        merged["is_peptidein"] | merged["final_tier"].isin(["1A", "1B"])
    ].copy().sort_values(["is_peptidein", "final_tier", "orf_id"], ascending=[False, True, True])

    meta_parquet = QUERIES_DIR / "query_metadata.parquet"
    meta_tsv = QUERIES_DIR / "query_metadata.tsv"
    meta_df.to_parquet(meta_parquet, index=False)
    meta_df.to_csv(meta_tsv, sep="\t", index=False)

    print(f"Saved metadata: {meta_parquet} ({len(meta_df)} rows)")
    print(f"Saved metadata: {meta_tsv}")

    print("\n--- Query Dataset Summary ---")
    print(f"• Tier 1A Peptideins (Tryptic MS + Ribo-seq): {(df_tier1_peps['final_tier'] == '1A').sum()}")
    print(f"• Tier 1B Peptideins (HLA + Ribo-seq):       {(df_tier1_peps['final_tier'] == '1B').sum()}")
    print(f"• Total Tier 1 Peptideins:                  {len(df_tier1_peps)}")
    print(f"• Total All Peptideins (Tiers 1A-2B):       {len(df_all_peps)}")
    print(f"• Total Tier 1 smORF Microproteins:         {len(df_tier1_all)}")
    print(f"• Peptidein length range: {df_tier1_peps['length'].min()} - {df_tier1_peps['length'].max()} aa (median {df_tier1_peps['length'].median():.0f} aa)")
    print("=" * 70 + "\n")


def download_canonical_proteome(url: str = UNIPROT_FTP_URL):
    """Download canonical human reference proteome archive from UniProt."""
    print("=" * 70)
    print("Phase 1: Downloading Canonical Human Reference Proteome")
    print("=" * 70)

    TARGETS_RAW_DIR.mkdir(parents=True, exist_ok=True)
    out_archive = TARGETS_RAW_DIR / "UP000005640_9606.fasta.gz"

    if out_archive.exists() and out_archive.stat().st_size > 5_000_000:
        print(f"Archive already exists ({out_archive.stat().st_size / 1e6:.2f} MB): {out_archive}")
        return out_archive

    print(f"Downloading from: {url}")
    print(f"Destination:      {out_archive}")

    start_time = time.time()
    response = requests.get(url, stream=True, timeout=60)
    response.raise_for_status()

    total_size = int(response.headers.get("content-length", 0))
    downloaded = 0
    chunk_size = 1024 * 512  # 512 KB

    with open(out_archive, "wb") as f:
        for chunk in response.iter_content(chunk_size=chunk_size):
            if chunk:
                f.write(chunk)
                downloaded += len(chunk)
                if total_size > 0:
                    pct = (downloaded / total_size) * 100
                    print(f"\rProgress: {downloaded / 1e6:.2f} / {total_size / 1e6:.2f} MB ({pct:.1f}%)", end="", flush=True)

    elapsed = time.time() - start_time
    print(f"\nDownload complete in {elapsed:.1f}s ({out_archive.stat().st_size / 1e6:.2f} MB).")
    print("=" * 70 + "\n")
    return out_archive


def index_canonical_proteome():
    """Extract and parse canonical FASTA into indexed metadata table."""
    print("=" * 70)
    print("Phase 1: Parsing and Indexing Canonical Human Proteome")
    print("=" * 70)

    raw_archive = TARGETS_RAW_DIR / "UP000005640_9606.fasta.gz"
    if not raw_archive.exists():
        download_canonical_proteome()

    TARGETS_DIR.mkdir(parents=True, exist_ok=True)
    out_fasta = TARGETS_DIR / "human_canonical_proteome.fasta"
    out_parquet = TARGETS_DIR / "human_canonical_proteome.parquet"
    out_tsv = TARGETS_DIR / "human_canonical_proteome_summary.tsv"

    print(f"Reading compressed archive: {raw_archive}")
    records = []
    current_header = None
    current_seq_chunks = []

    # Regex to parse UniProt FASTA headers
    # >db|Accession|EntryName ProteinDescription OS=Organism OX=TaxID [GN=Gene] PE=Evidence SV=Version
    header_pattern = re.compile(
        r"^>(?:sp|tr)\|(?P<accession>[A-Z0-9_-]+)\|(?P<entry_name>\S+)\s+(?P<desc>.*)$"
    )

    with gzip.open(raw_archive, "rt", encoding="utf-8") as gz_in, \
         open(out_fasta, "w", encoding="utf-8") as fasta_out:

        for line in gz_in:
            fasta_out.write(line)
            line = line.strip()
            if not line:
                continue

            if line.startswith(">"):
                if current_header is not None:
                    seq = "".join(current_seq_chunks)
                    records.append((current_header, seq))
                current_header = line
                current_seq_chunks = []
            else:
                current_seq_chunks.append(line)

        if current_header is not None:
            seq = "".join(current_seq_chunks)
            records.append((current_header, seq))

    print(f"Extracted {len(records):,} sequences to {out_fasta}")

    # Parse metadata fields into DataFrame
    parsed_entries = []
    for header, seq in records:
        match = header_pattern.match(header)
        if match:
            acc = match.group("accession")
            entry = match.group("entry_name")
            desc_full = match.group("desc")

            # Extract GN (gene symbol)
            gn_match = re.search(r"GN=([^\s]+)", desc_full)
            gene = gn_match.group(1) if gn_match else None

            # Extract OS (organism)
            os_match = re.search(r"OS=(.*?)\s+[A-Z]{2}=", desc_full)
            organism = os_match.group(1) if os_match else "Homo sapiens"

            # Protein name is description up to OS=
            pname_match = re.split(r"\s+OS=", desc_full, maxsplit=1)
            protein_name = pname_match[0].strip() if pname_match else desc_full

            # Protein evidence PE
            pe_match = re.search(r"PE=(\d+)", desc_full)
            pe = int(pe_match.group(1)) if pe_match else None

            parsed_entries.append({
                "uniprot_id": acc,
                "entry_name": entry,
                "gene_symbol": gene,
                "protein_name": protein_name,
                "pe_level": pe,
                "length": len(seq),
                "sequence": seq,
            })
        else:
            # Fallback parsing
            parts = header[1:].split(None, 1)
            acc = parts[0]
            parsed_entries.append({
                "uniprot_id": acc,
                "entry_name": acc,
                "gene_symbol": None,
                "protein_name": parts[1] if len(parts) > 1 else "",
                "pe_level": None,
                "length": len(seq),
                "sequence": seq,
            })

    df_targets = pd.DataFrame(parsed_entries)
    df_targets.to_parquet(out_parquet, index=False)

    # Save lightweight summary table without full sequence string
    summary_cols = ["uniprot_id", "entry_name", "gene_symbol", "protein_name", "pe_level", "length"]
    df_targets[summary_cols].to_csv(out_tsv, sep="\t", index=False)

    print(f"Saved indexed target table: {out_parquet} ({len(df_targets):,} rows)")
    print(f"Saved target summary TSV:   {out_tsv}")

    print("\n--- Target Proteome Summary ---")
    print(f"• Total canonical targets:  {len(df_targets):,}")
    print(f"• Unique UniProt IDs:       {df_targets['uniprot_id'].nunique():,}")
    print(f"• Unique Gene Symbols:      {df_targets['gene_symbol'].nunique():,}")
    print(f"• Target length min/median/max: {df_targets['length'].min()} / {df_targets['length'].median():.0f} / {df_targets['length'].max()} aa")
    print("=" * 70 + "\n")


def main():
    parser = argparse.ArgumentParser(
        description="Phase 1: Curate Microprotein Queries & Human Canonical Proteome Targets"
    )
    parser.add_argument("--curate-queries", action="store_true", help="Curate microprotein/peptidein query sets")
    parser.add_argument("--download-proteome", action="store_true", help="Download UniProt reference proteome archive")
    parser.add_argument("--index-proteome", action="store_true", help="Parse and index canonical proteome into Parquet/FASTA")
    parser.add_argument("--all", action="store_true", help="Run full Phase 1 curation pipeline")

    args = parser.parse_args()

    if not any([args.curate_queries, args.download_proteome, args.index_proteome, args.all]):
        parser.print_help()
        sys.exit(1)

    if args.all or args.curate_queries:
        curate_queries()

    if args.all or args.download_proteome:
        download_canonical_proteome()

    if args.all or args.index_proteome:
        index_canonical_proteome()

    print("Phase 1 curation complete!")


if __name__ == "__main__":
    main()
