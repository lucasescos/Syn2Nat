#!/usr/bin/env python3
import json
import os
import sys
import urllib.request
from pathlib import Path

def get_key():
    for p in [Path.cwd() / '.env', Path.home() / '.env']:
        if p.exists():
            with open(p, 'r', encoding='utf-8') as f:
                for line in f:
                    if line.startswith('ESM_API_KEY='):
                        return line.strip().split('=', 1)[1].strip('"\'')
    return os.environ.get('ESM_API_KEY')

def main():
    api_key = get_key()
    body = {
        'model': 'esm3-open-2024-03',
        'track': 'sequence',
        'inputs': {
            'sequence': 'M___G__L___K'
        },
        'num_steps': 5,
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
            print('[SUCCESS] Output sequence:', res['outputs']['sequence'])
    except urllib.error.HTTPError as e:
        print(f'[ERROR] HTTP {e.code}:', e.read().decode('utf-8', errors='replace'))
    except Exception as e:
        print('[ERROR]:', e)

if __name__ == '__main__':
    main()
