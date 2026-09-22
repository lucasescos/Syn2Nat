#!/usr/bin/env python3
import json
import urllib.request
from test_esm3_api import get_key

def test_structure():
    api_key = get_key()
    body = {
        'model': 'esm3-open-2024-03',
        'track': 'structure',
        'inputs': {
            'sequence': 'MSHHWGYGKH'
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
        coords = res['outputs']['coordinates']
        print('Returned coordinates length:', len(coords))
        print('Single residue coordinate item shape/type:', type(coords[0]), len(coords[0]))
        print('Sample atom (N, CA, C):', coords[0][:3])

if __name__ == '__main__':
    test_structure()
