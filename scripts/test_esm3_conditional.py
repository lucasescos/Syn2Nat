#!/usr/bin/env python3
import json
import urllib.request
from test_esm3_api import get_key

api_key = get_key()
sample_coord = [
    [-7.5, 0.8, 12.2], [-6.8, 1.6, 11.2], [-6.0, 0.8, 10.1]
]

for test_name, empty_coord in [
    ("empty_list", []),
    ("list_of_none", [[None, None, None]] * 3),
    ("list_of_zeros", [[0.0, 0.0, 0.0]] * 37),
    ("list_of_3_zeros", [[0.0, 0.0, 0.0]] * 3),
]:
    body = {
        'model': 'esm3-open-2024-03',
        'track': 'sequence',
        'inputs': {
            'sequence': 'ACDEF_____',
            'coordinates': [sample_coord]*5 + [empty_coord]*5
        },
        'num_steps': 3,
        'temperature': 0.5
    }
    req = urllib.request.Request(
        'https://biohub.ai/api/v1/generate',
        data=json.dumps(body).encode('utf-8'),
        headers={'Authorization': f'Bearer {api_key}', 'Content-Type': 'application/json'},
        method='POST'
    )
    try:
        with urllib.request.urlopen(req) as resp:
            res = json.loads(resp.read().decode('utf-8'))
            print(f'[OK] {test_name}:', res['outputs']['sequence'])
            break
    except urllib.error.HTTPError as e:
        print(f'[FAIL] {test_name}:', e.read().decode('utf-8', errors='replace')[:180])
    except Exception as e:
        print(f'[FAIL] {test_name}:', e)
