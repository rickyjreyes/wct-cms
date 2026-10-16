#!/usr/bin/env python3
"""Recover exactly the four historical inputs without rewriting frozen receipts."""
from __future__ import annotations
import argparse
import concurrent.futures
import hashlib
import json
import subprocess
from pathlib import Path
from trace_cms_verified_root_stages import ROOT, PREFIX, published_identity, hash_input


def recover(label, directory, download=False):
    expected = published_identity(ROOT, label)
    path = directory / expected['filename']
    url = 'https://opendata.cern.ch' + expected['uri'].split('eospublic.cern.ch/')[1]
    row = {'label': label, 'expected': expected, 'url': url, 'path': str(path)}
    try:
        if download and not path.exists():
            partial = path.with_suffix('.root.partial')
            subprocess.run(['curl', '-L', '--fail', '--silent', '--show-error',
                            '--connect-timeout', '20', '--max-time', '1800',
                            '--retry', '2', '-C', '-', url, '-o', str(partial)], check=True)
            # Validate under the official filename before accepting the download.
            if partial.stat().st_size != expected['size']:
                raise ValueError('incomplete download; partial retained')
            partial.rename(path)
        row['receipt'] = hash_input(path, expected)
        row['status'] = 'VERIFIED'
    except (OSError, ValueError, subprocess.CalledProcessError) as exc:
        row['status'] = 'BLOCKED_EXTERNAL' if not path.exists() else 'INVALID_INPUT'
        row['error'] = str(exc)
    print(label, row['status'], flush=True)
    return row


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--directory', type=Path, default=ROOT/'data/audit_inputs')
    p.add_argument('--download', action='store_true')
    p.add_argument('--output', required=True, type=Path)
    args = p.parse_args()
    if args.output.exists():
        p.error('refusing to overwrite an earlier receipt')
    args.directory.mkdir(parents=True, exist_ok=True)
    with concurrent.futures.ThreadPoolExecutor(max_workers=4) as pool:
        rows = list(pool.map(lambda label: recover(label, args.directory, args.download), PREFIX))
    complete = all(r['status'] == 'VERIFIED' for r in rows)
    out = {'schema': 'CMS_EVENT_INPUT_RECOVERY_V1', 'complete': complete,
           'sampling': 'first 100000 ROOT entries BEFORE luminosity certification',
           'official_manifest_sha256': {
               record: hashlib.sha256((ROOT/f'audit/2026-09-06/baseline/cms{record}.json').read_bytes()).hexdigest()
               for record in ('30555','30522')}, 'files': rows}
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open('x') as f:
        json.dump(out, f, indent=2)
        f.write('\n')
    if not complete:
        raise SystemExit(2)

if __name__ == '__main__':
    main()
