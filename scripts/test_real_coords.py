#!/usr/bin/env python3
import json
import urllib.request
from test_esm3_api import get_key

api_key = get_key()
# Fetch a real 37-atom coordinate tensor from structure track for a 4-aa peptide
body = {
    'model': 'esm3-open-2024-03',
    'track': 'structure',
    'inputs': {
        'sequence': 'MGLK'
    },
    'num_steps': 1,
    'temperature': 0.0
}
req = urllib.request.Request(
    'https://biohub.ai/api/v1/generate',
    data=json.dumps(body).encode('utf-8'),
    headers={'Authorization': f'Bearer {api_key}', 'Content-Type': 'application/json'},
    method='POST'
)
with urllib.request.urlopen(req) as resp:
    res = json.loads(resp.read().decode('utf-8'))
    coords_4 = res['outputs']['coordinates']
    print('Coords len:', len(coords_4))
    print('Residue 0 atom count:', len(coords_4[0]))
    print('Residue 0 non-null atom count:', sum(1 for a in coords_4[0] if a is not None))
    # print what an unassigned atom looks like in coords_4[0]
    unassigned = [a for a in coords_4[0] if a is None or a == [0,0,0] or any(x is None for x in a)]
    print('Unassigned atom sample:', unassigned[:2] if unassigned else 'None are None')

    # Now let's test feeding this into track='sequence' with a masked residue
    # 4 residues with coords, followed by 4 residues with unassigned/null coords
    # Let's see what unassigned residue coords look like
    empty_res = [[None, None, None]] * 37
    test_body = {
        'model': 'esm3-open-2024-03',
        'track': 'sequence',
        'inputs': {
            'sequence': 'MGLK____',
            'coordinates': coords_4 + [empty_res]*4
        },
        'num_steps': 4,
        'temperature': 0.5
    }
    req2 = urllib.request.Request(
        'https://biohub.ai/api/v1/generate',
        data=json.dumps(test_body).encode('utf-8'),
        headers={'Authorization': f'Bearer {api_key}', 'Content-Type': 'application/json'},
        method='POST'
    )
    try:
        with urllib.request.urlopen(req2) as resp2:
            res2 = json.loads(resp2.read().decode('utf-8'))
            print('[SUCCESS] Infill conditioned on real coords:', res2['outputs']['sequence'])
    except urllib.error.HTTPError as e:
        print('[FAIL] Infill conditioned on real coords (HTTP):', e.read().decode('utf-8', errors='replace'))
