#!/usr/bin/env python3
import sys
from pathlib import Path
BASE_DIR = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(BASE_DIR))

from scripts.annotate_interface_sites import fetch_uniprot_features, fetch_pdbekb_interface_residues

for gene, acc in [('PSMA1', 'P25786'), ('RPRD1A', 'Q96P16'), ('CIB2', 'O75838')]:
    print('=' * 60)
    print(f"Target: {gene} ({acc})")
    u_data = fetch_uniprot_features(acc)
    p_data = fetch_pdbekb_interface_residues(acc)

    # UniProt features
    feats = u_data.get('features', []) if u_data else []
    print(f"Total UniProt features: {len(feats)}")
    for f in feats:
        ftype = f.get('type')
        if ftype in ['Active site', 'Binding site', 'Domain', 'Region']:
            start = f.get('location', {}).get('start', {}).get('value')
            end = f.get('location', {}).get('end', {}).get('value')
            desc = f.get('description', '')
            print(f"  [{ftype}] Res {start}–{end}: {desc}")

    # PDBe-KB residues
    if p_data and acc in p_data:
        entries = p_data[acc].get('data', [])
        known_res = set()
        for e in entries:
            for r in e.get('residues', []):
                s = r.get('startIndex')
                end = r.get('endIndex', s)
                if s is not None:
                    for pos in range(s, end + 1):
                        known_res.add(pos)
        sorted_res = sorted(list(known_res))
        print(f"PDBe-KB: {len(sorted_res)} experimental contact residues across {len(entries)} PDB interfaces!")
        print(f"Sample interface residues: {sorted_res[:25]}")
    else:
        print("PDBe-KB: No entries found.")
