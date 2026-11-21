"""Resume a large public archive with bounded parallel byte-range requests."""
from __future__ import annotations
import argparse
import concurrent.futures as cf
import gzip
import hashlib
import json
from pathlib import Path
import shutil
import zipfile

import requests
from download_encori_bulk_ranges import chunk


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--url', required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--workers', type=int, default=8)
    args = parser.parse_args()
    path = args.output.resolve()
    if path.exists():
        raise FileExistsError('Preserving existing final file: ' + str(path))
    path.parent.mkdir(parents=True, exist_ok=True)
    headers = {'Range': 'bytes=0-0', 'Accept-Encoding': 'identity'}
    with requests.get(args.url, headers=headers, timeout=(15, 30), stream=True) as response:
        response.raise_for_status()
        value = response.headers.get('Content-Range', '')
        if response.status_code != 206 or not value.startswith('bytes 0-0/'):
            raise ValueError('Provider does not support verified byte ranges')
        total = int(value.split('/')[1])
        metadata = {'url': args.url, 'etag': response.headers.get('ETag'),
                    'last_modified': response.headers.get('Last-Modified'), 'expected_bytes': total}
    prefix = path.with_name(path.name + '.part')
    start = prefix.stat().st_size if prefix.exists() else 0
    if start > total:
        raise ValueError('Existing partial file exceeds provider size')
    tasks = [(args.url, left, min(total-1, left+4*1024**2-1),
              path.with_name(path.name + f'.range_{left}_{min(total-1,left+4*1024**2-1)}'), total)
             for left in range(start, total, 4*1024**2)]
    print(json.dumps({'prefix_bytes': start, 'expected_bytes': total, 'range_jobs': len(tasks)}), flush=True)
    with cf.ThreadPoolExecutor(args.workers) as pool:
        for n, _ in enumerate(pool.map(chunk, tasks), 1):
            print(json.dumps({'ranges_complete': n, 'total': len(tasks)}), flush=True)
    assembled = path.with_name(path.name + '.assembling')
    with assembled.open('wb') as out:
        for part in ([prefix] if prefix.exists() else []) + [t[3] for t in tasks]:
            with part.open('rb') as source:
                shutil.copyfileobj(source, out, 1024**2)
    if assembled.stat().st_size != total:
        raise ValueError('Assembled length mismatch')
    if path.name.endswith('.zip'):
        with zipfile.ZipFile(assembled) as archive:
            if archive.testzip():
                raise ValueError('ZIP CRC check failed')
    elif path.name.endswith('.gz'):
        with gzip.open(assembled, 'rb') as stream:
            for _ in iter(lambda: stream.read(8*1024**2), b''):
                pass
    with assembled.open('rb') as stream:
        metadata['sha256'] = hashlib.file_digest(stream, 'sha256').hexdigest()
    assembled.replace(path)
    path.with_name(path.name + '.range_provenance.json').write_text(json.dumps(metadata, indent=2))
    for part in ([prefix] if prefix.exists() else []) + [t[3] for t in tasks]:
        part.unlink()
    print(json.dumps({'verified': str(path), **metadata}), flush=True)


if __name__ == '__main__':
    main()
