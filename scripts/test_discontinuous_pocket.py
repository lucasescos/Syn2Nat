import json
import math
import sys
import urllib.request
from pathlib import Path

BASE_DIR = Path(__file__).resolve().parent.parent
if str(BASE_DIR) not in sys.path:
    sys.path.insert(0, str(BASE_DIR))

from esm.sdk.api import ESMProtein
from scripts.test_esm3_api import get_key
from scripts.annotate_interface_sites import fetch_pdbekb_interface_residues

def test_rprd1a_authentic_pocket():
    api_key = get_key()
    
    # 1. Target: RPRD1A (Q96P16)
    pdb_path = "structures/pilot_15_targets/RPRD1A_Q96P16.pdb"
    prot = ESMProtein.from_pdb(pdb_path)
    full_seq = prot.sequence
    full_coords = prot.coordinates.tolist()
    
    # 2. Get real PDBe-KB experimental interface residues
    pdbekb_data = fetch_pdbekb_interface_residues("Q96P16")
    entries = pdbekb_data["Q96P16"].get("data", [])
    interface_pos = set()
    for e in entries:
        for r in e.get("residues", []):
            s = r.get("startIndex")
            end = r.get("endIndex", s)
            if s is not None:
                for pos in range(s, end + 1):
                    interface_pos.add(pos)
                    
    sorted_pos = sorted(list(interface_pos))
    print(f"[RPRD1A] Experimental PDBe-KB interface residues ({len(sorted_pos)} total): {sorted_pos}")
    
    # 3. Build sub-structure of authentic pocket residues
    # Extract only the pocket residues and their 3D coordinates!
    pocket_seq_list = []
    pocket_coords_list = []
    
    for pos in sorted_pos:
        idx = pos - 1 # 0-indexed
        if idx < len(full_seq):
            pocket_seq_list.append(full_seq[idx])
            pocket_coords_list.append(full_coords[idx])
            
    pocket_seq = "".join(pocket_seq_list)
    print(f"Authentic pocket sequence string ({len(pocket_seq)} aa): {pocket_seq}")
    
    # Clean coordinates
    clean_pocket_coords = []
    for res in pocket_coords_list:
        clean_res = []
        for atom in res:
            if atom is None:
                clean_res.append([None, None, None])
            else:
                clean_res.append([
                    None if (v is None or (isinstance(v, float) and math.isnan(v))) else float(v)
                    for v in atom
                ])
        clean_pocket_coords.append(clean_res)
        
    # 4. Add 30-aa masked binder
    binder_len = 30
    prompt_seq = pocket_seq + ("_" * binder_len)
    empty_res = [[None, None, None]] * 37
    prompt_coords = clean_pocket_coords + ([empty_res] * binder_len)
    
    print(f"Prompting ESM3 generate on authentic pocket ({len(prompt_seq)} total residues)...")
    body = {
        'model': 'esm3-open-2024-03',
        'track': 'sequence',
        'inputs': {
            'sequence': prompt_seq,
            'coordinates': prompt_coords
        },
        'num_steps': 15,
        'temperature': 0.5,
        'temperature_annealing': True
    }
    req = urllib.request.Request(
        'https://biohub.ai/api/v1/generate',
        data=json.dumps(body).encode('utf-8'),
        headers={'Authorization': f'Bearer {api_key}', 'Content-Type': 'application/json'},
        method='POST'
    )
    with urllib.request.urlopen(req) as resp:
        res = json.loads(resp.read().decode('utf-8'))
        full_out_seq = res['outputs']['sequence']
        designed_binder = full_out_seq[len(pocket_seq):]
        print(f"[SUCCESS] Designed binder for authentic RPRD1A pocket: {designed_binder}")
        return designed_binder

if __name__ == "__main__":
    test_rprd1a_authentic_pocket()
