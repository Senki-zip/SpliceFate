"""Resume the 2026-11-19 resource restoration and publish only validated outputs.

Run after source discovery/downloads have been started. Logs, stage markers and
provider failures are retained; unavailable sources are never silently omitted.
"""
from __future__ import annotations

import concurrent.futures as cf
import fcntl
import gzip
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
import threading
import time
from urllib.parse import urljoin
import zipfile

import requests
from bs4 import BeautifulSoup

PROJECT = Path(__file__).resolve().parents[1]
RUN = PROJECT / 'runs/database_rebuild_20261119'
EXTERNAL = PROJECT / 'data/external'
ATLAS = EXTERNAL / 'rbp_atlas_20261115'
ORN = EXTERNAL / 'ornament_v2_20261115'
MOUSE = EXTERNAL / 'mouse_rbp_unified_20261116'
INTER = EXTERNAL / 'rbp_interactions_20261116'
GENES = MOUSE / 'reference/mm10_genes.tsv'
LOCK = threading.RLock()
STATE = {}
MOUSE_ORNAMENT_INPUTS = frozenset({
    'raw/oRNAment_motif_ppms.meme', 'raw/Motif.Info.tsv',
    'raw/gene_biotype.csv', 'raw/transcript_biotype.csv',
    'raw/data/tables/mus_musculus.csv.xz',
    'raw/data/tables/mus_musculus.ids.csv.xz',
    'raw/data/tables/mus_musculus.structure.csv.xz',
})


def stamp():
    return time.strftime('%Y-%m-%dT%H:%M:%S%z')


def save(**updates):
    with LOCK:
        STATE.update(updates, updated_at=stamp())
        tmp = RUN / 'status.json.tmp'
        tmp.write_text(json.dumps(STATE, indent=2, ensure_ascii=False))
        tmp.replace(RUN / 'status.json')


def stage(name, status, **details):
    with LOCK:
        STATE.setdefault('stages', {})[name] = {'status': status, 'updated_at': stamp(), **details}
        save()
    print(json.dumps({'stage': name, 'status': status, **details}), flush=True)


def active(pidfile):
    if not pidfile.exists():
        return False
    try:
        pid = int(pidfile.read_text().strip())
        stat = Path(f'/proc/{pid}/stat').read_text()
        return stat.split(') ', 1)[1].split()[0] != 'Z'
    except (OSError, ValueError):
        return False


def wait_for_downloads(files):
    while any(active(p) for p in files):
        save()
        time.sleep(20)


def command(name, script, *args):
    marker = RUN / (name + '.complete.json')
    if marker.exists():
        return
    stage(name, 'running')
    log = RUN / (name + '.log')
    pidfile = RUN / (name + '.pid')
    if name == 'ornament_download' and active(pidfile):
        # A restarted coordinator must not duplicate an ongoing large download.
        wait_for_downloads([pidfile])
        state_path = ORN / 'download_status.json'
        if not state_path.exists() or json.loads(state_path.read_text())['status'] == 'running':
            raise RuntimeError('Previous oRNAment downloader exited without a final status')
    else:
        with log.open('a') as stream:
            process = subprocess.Popen([sys.executable, '-u', str(PROJECT / 'scripts' / script), *map(str, args)],
                                       cwd=PROJECT, stdout=stream, stderr=subprocess.STDOUT)
            pidfile.write_text(str(process.pid) + '\n')
            code = process.wait()
            if code:
                raise subprocess.CalledProcessError(code, process.args)
    marker.write_text(json.dumps({'finished_at': stamp(), 'script': script, 'args': list(map(str, args))}, indent=2))
    stage(name, 'complete')


def retry_encori():
    from collect_rbp_resources import download
    for kind in ['peak', 'motif']:
        status_file = ATLAS / f'encori_{kind}_download_status.json'
        data = json.loads(status_file.read_text())
        rows = data['failed'] if kind == 'peak' else [r for r in data if r['status'] == 'failed']
        for row in rows:
            result = download('ENCORI', row['url'], row['path'], compress=kind == 'peak')
            if kind == 'motif' and '404' in result.get('error', ''):
                result['status'] = 'not_published'
            row.update(result)
        if kind == 'peak':
            data['failed'] = [r for r in rows if r['status'] == 'failed']
            failed = data['failed']
        else:
            failed = [r for r in data if r['status'] == 'failed']
        status_file.write_text(json.dumps(data, indent=2))
        if failed:
            raise RuntimeError(f'Unresolved ENCORI {kind} downloads: {len(failed)}')
    # Some provider HTTP-200 responses contain PHP failures rather than BED data.
    # Preserve that raw response and use the existing documented API repair.
    for meta in sorted((ATLAS / 'raw/encori').glob('*_experiments.json')):
        assembly = meta.name.split('_experiments')[0]
        for row in json.loads(meta.read_text())['data']:
            dataset = row['dataSetId']
            path = ATLAS / f'raw/encori/peaks/{assembly}/{dataset}.bed.gz'
            invalid = False
            with gzip.open(path, 'rt') as stream:
                for line in stream:
                    if not line.strip() or line.lstrip().startswith('#'):
                        continue
                    fields = line.rstrip().split('\t')
                    try:
                        assert len(fields) >= 6 and int(fields[1]) >= 0 and int(fields[2]) > int(fields[1])
                        assert fields[5] in {'+', '-', '.'}
                    except (ValueError, AssertionError):
                        invalid = True
                        break
            if invalid:
                command('repair_' + dataset, 'repair_encori_peak_payload.py', '--assembly', assembly, '--dataset', dataset)


def atlas_pipeline():
    stage('atlas_downloads', 'waiting_for_downloads')
    wait_for_downloads([ATLAS / (task + '.pid') for task in ['basic', 'cisbp', 'postar', 'encori-peaks', 'encori-motifs', 'integration_core']])
    if not (ATLAS / 'integrated/core_stats.json').exists():
        command('atlas_core', 'integrate_rbp_resources.py', '--phase', 'core', '--local-genes', GENES)
    urls = set(re.findall(r'https://rnasysu\.com/encori/download/[^\s<>]+\.RBP_[^\s<>]+\.tar\.gz',
                         (ATLAS / 'discovery/encori_api.html').read_text()))
    if any(not (ATLAS / 'raw/encori' / u.rsplit('/', 1)[-1]).exists() for u in urls):
        command('encori_bulk_resume', 'download_encori_bulk_ranges.py')
    retry_encori()
    for p in [ATLAS / 'raw/cisbp/All_species_RBP_information_PWMs.zip',
              ATLAS / 'raw/cisbp/SQLArchive_cisbp-rna_2_00.zip',
              ATLAS / 'raw/attract/ATtRACT.zip']:
        with zipfile.ZipFile(p) as archive:
            bad = archive.testzip()
            if bad:
                raise ValueError(f'ZIP CRC failure: {p}: {bad}')
    stage('atlas_downloads', 'complete')
    command('atlas_finish', 'integrate_rbp_resources.py', '--phase', 'finish', '--local-genes', GENES)
    command('atlas_validation', 'validate_rbp_collection.py', '--local-genes', GENES,
            '--report-path', PROJECT / 'docs/RBP_RESOURCE_RESTORATION_20261119.md')


def interaction_pipeline():
    stage('interaction_downloads', 'waiting_for_downloads')
    wait_for_downloads([INTER / 'download.pid'])
    if json.loads((INTER / 'download_status.json').read_text())['status'] != 'complete':
        command('interaction_download_retry', 'collect_mouse_interaction_resources.py')
    checks = []
    manifest = json.loads((INTER / 'download_manifest.json').read_text())
    plan = json.loads((INTER / 'download_plan.json').read_text())
    if {r['path'] for r in manifest} != {r[2] for r in plan}:
        raise ValueError('Interaction manifest does not cover the download plan')
    for row in manifest:
        p = INTER / row['path']
        with p.open('rb') as stream:
            digest = hashlib.file_digest(stream, 'sha256').hexdigest()
        passed = p.stat().st_size == row['bytes'] and digest == row['sha256']
        if p.name.endswith('.gz'):
            with gzip.open(p, 'rb') as stream:
                for _ in iter(lambda: stream.read(8 * 1024**2), b''):
                    pass
        elif p.name.endswith('.zip'):
            with zipfile.ZipFile(p) as archive:
                passed &= archive.testzip() is None
        checks.append({'path': row['path'], 'passed': passed})
    passed = all(r['passed'] for r in checks)
    (INTER / 'validation.json').write_text(json.dumps({'passed': passed, 'checked_at': stamp(), 'checks': checks}, indent=2))
    if not passed:
        raise ValueError('Interaction raw-source validation failed')
    stage('interaction_downloads', 'complete', verified_files=len(checks))


def discover_ornament():
    discovery = ORN / 'discovery'
    discovery.mkdir(parents=True, exist_ok=True)
    url = 'https://rnabiology.ircm.qc.ca/oRNAmentv2/downloads/'
    response = requests.get(url, timeout=(12, 30))
    response.raise_for_status()
    (discovery / 'downloads.html').write_text(response.text)
    soup = BeautifulSoup(response.text, 'html.parser')
    urls = set()
    for tag in soup.select('a[href]'):
        link = urljoin(response.url, tag['href'])
        if '/static/download/' in link and '/data/legacy/' not in link and not link.endswith('/'):
            urls.add(link)
    if not urls or not any('mus_musculus.csv.xz' in u for u in urls):
        raise ValueError('Official page has no expected V2 mouse exports; discovery requires review')
    def metadata(link):
        head = requests.head(link, timeout=(12, 30), allow_redirects=True)
        if head.status_code in (403, 404):
            # Some official track archives are listed but unavailable. Preserve
            # their status without blocking discovery of accessible table data.
            return {'url': link, 'bytes': None, 'http_status': head.status_code,
                    'discovery_error': f'Official endpoint returned HTTP {head.status_code}'}
        head.raise_for_status()
        row = {'url': link, 'bytes': int(head.headers['Content-Length']) if head.headers.get('Content-Length') else None,
               'etag': head.headers.get('ETag'), 'modified': head.headers.get('Last-Modified')}
        if link.endswith('.md5'):
            result = requests.get(link, timeout=(12, 30)); result.raise_for_status()
            (discovery / link.rsplit('/', 1)[-1]).write_text(result.text)
        return row
    with cf.ThreadPoolExecutor(4) as pool:
        rows = list(pool.map(metadata, sorted(urls)))
    (discovery / 'official_files.json').write_text(json.dumps(rows, indent=2))


def mouse_ornament_readiness(root=ORN):
    """Use atomically published status/manifest records, never partial files."""
    status_path = root / 'download_status.json'
    manifest_path = root / 'download_manifest.json'
    if not status_path.exists() or not manifest_path.exists():
        return sorted(MOUSE_ORNAMENT_INPUTS), []
    state = json.loads(status_path.read_text())
    verified = {r['path'] for r in state['files'] if r['status'] == 'verified'}
    records = {r['path']: r for r in json.loads(manifest_path.read_text())}
    pending, ready = [], []
    for relative in sorted(MOUSE_ORNAMENT_INPUTS):
        row = records.get(relative)
        path = root / relative
        if (relative not in verified or not row or not row.get('sha256')
                or not path.is_file() or path.stat().st_size != row['bytes']):
            pending.append(relative)
        else:
            ready.append(row)
    return pending, ready


def ornament_pipeline():
    # Bounded durable retries allow a transient provider outage to recover while
    # independent sources are integrated. Every attempt is visible in status.
    for attempt in range(1, 49):
        try:
            if not (ORN / 'discovery/official_files.json').exists():
                discover_ornament()
            pidfile = RUN / 'ornament_download.pid'
            pending, records = mouse_ornament_readiness()
            if pending and not active(pidfile):
                with (RUN / 'ornament_download.log').open('a') as stream:
                    process = subprocess.Popen(
                        [sys.executable, '-u', str(PROJECT / 'scripts/download_ornament_v2.py')],
                        cwd=PROJECT, stdout=stream, stderr=subprocess.STDOUT,
                        start_new_session=True)
                pidfile.write_text(str(process.pid) + '\n')
            stage('ornament_download', 'running' if active(pidfile) else 'complete',
                  pid=int(pidfile.read_text()) if pidfile.exists() else None,
                  scope='All species continue independently of mouse database build')
            while True:
                pending, records = mouse_ornament_readiness()
                if not pending:
                    stage('ornament', 'checking_mouse_inputs', files=len(records))
                    for row in records:
                        with (ORN / row['path']).open('rb') as stream:
                            digest = hashlib.file_digest(stream, 'sha256').hexdigest()
                        if digest != row['sha256']:
                            raise ValueError('Mouse input checksum mismatch: ' + row['path'])
                        if row['path'].endswith('.xz') and row.get('official_md5_match') is not True:
                            raise ValueError('Mouse table lacks official MD5 verification: ' + row['path'])
                    # Freeze only the exact inputs used by the mouse build; the
                    # full download manifest will keep changing in the background.
                    snapshot = ORN / 'mouse_build_input_manifest.json'
                    temp = snapshot.with_suffix('.json.tmp')
                    temp.write_text(json.dumps(records, indent=2))
                    temp.replace(snapshot)
                    stage('ornament', 'complete', verified_mouse_files=len(records),
                          other_species_downloading=active(pidfile))
                    return
                stage('ornament', 'waiting_for_mouse_inputs', pending=pending,
                      poll_seconds=20)
                if not active(pidfile):
                    raise RuntimeError('Downloader stopped before required mouse files were verified')
                time.sleep(20)
        except Exception as exc:
            stage('ornament', 'waiting_for_provider', attempt=attempt, error=str(exc), retry_seconds=300)
            if attempt == 48:
                raise
            time.sleep(300)


def backup():
    mount = Path('/media/huangtao/ae2c2db0-d1db-4f02-9cc6-79d9bd823577')
    if not os.path.ismount(mount):
        raise RuntimeError('10 TB backup drive is not mounted')
    dest = mount / 'home/hcq/SpliceFate_backups/mouse_rbp_rebuilt_20261119'
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists():
        raise FileExistsError('Preserving existing backup: ' + str(dest))
    temporary = dest.with_name(dest.name + '.copying')
    stage('backup', 'running', destination=str(dest))
    shutil.copytree(MOUSE, temporary, dirs_exist_ok=True, ignore=shutil.ignore_patterns('_tmp', '*.wal'))
    for row in (MOUSE / 'SHA256SUMS').read_text().splitlines():
        expected, relative = row.split(maxsplit=1)
        p = temporary / relative.lstrip('*')
        with p.open('rb') as stream:
            actual = hashlib.file_digest(stream, 'sha256').hexdigest()
        if expected != actual:
            raise ValueError('Backup checksum mismatch: ' + str(p))
    temporary.rename(dest)
    stage('backup', 'complete', destination=str(dest))


def main():
    RUN.mkdir(parents=True, exist_ok=True)
    with (RUN / '.restore.lock').open('w') as guard:
        fcntl.flock(guard, fcntl.LOCK_EX | fcntl.LOCK_NB)
        old = RUN / 'status.json'
        if old.exists():
            STATE.update(json.loads(old.read_text()))
        save(status='running', pid=os.getpid(), started_at=STATE.get('started_at', stamp()),
             scope='Restore original public RBP resources, rebuild mouse database, validate and run real-data smoke',
             database_published=False)
        save(mouse_first=True, monitor_interval_seconds=20)
        failures = []
        with cf.ThreadPoolExecutor(3) as pool:
            futures = {pool.submit(fun): name for name, fun in [
                ('atlas', atlas_pipeline), ('interactions', interaction_pipeline), ('ornament', ornament_pipeline)]}
            for future in cf.as_completed(futures):
                name = futures[future]
                try:
                    future.result()
                    stage(name, 'complete')
                except Exception as exc:
                    failures.append({'stage': name, 'error': str(exc)})
                    stage(name, 'failed', error=str(exc))
                save(failures=failures)
        if failures:
            save(status='incomplete', failures=failures)
            return 1
        wait_for_downloads([MOUSE / 'references.pid'])
        command('mouse_database', 'build_mouse_rbp_database.py', '--stage', 'all')
        # The extension scripts update existing checksums after committing.
        sums = MOUSE / 'SHA256SUMS'
        if not sums.exists():
            sums.write_text('')
        command('family_motifs', 'extend_mouse_family_motifs.py')
        command('interaction_audit', 'audit_mouse_interaction_resources.py')
        command('interaction_integration', 'integrate_mouse_interactions.py')
        command('mouse_validation', 'validate_mouse_rbp_database.py')
        if not json.loads((MOUSE / 'validation.json').read_text())['passed']:
            raise ValueError('Mouse database validation failed')
        save(database_published=True)
        backup()
        command('real_smoke', 'run_model_smoke.py', '--config', 'configs/nemo_p21_v0.3.yaml',
                '--output', 'runs/NEMO_P21_restored_database_smoke_20261119')
        save(status='complete', finished_at=stamp(),
             other_species_downloading=active(RUN / 'ornament_download.pid'))
    return 0


if __name__ == '__main__':
    try:
        raise SystemExit(main())
    except Exception as exc:
        save(status='failed', error=str(exc))
        raise
