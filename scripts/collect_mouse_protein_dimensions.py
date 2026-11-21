"""Download sequence-resolved mouse RBP annotations, retaining full evidence.

UniProt requests use the resource_http REST client. InterPro
member-signature and integrated matches remain distinct; source absences and
network failures are recorded separately. Restart to reuse verified payloads.
"""
from __future__ import annotations
import concurrent.futures as cf
import hashlib
import json
from pathlib import Path
import threading
import time

import duckdb
import requests
import resource_http

PROJECT = Path(__file__).resolve().parents[1]
ROOT = PROJECT / 'data/external/rbp_protein_dimensions_20261120'
DB = PROJECT / 'data/external/mouse_rbp_unified_20261116/mouse_rbp.duckdb'
LOCK = threading.RLock()
RECORDS = {}
STATE = {}


def stamp():
    return time.strftime('%Y-%m-%dT%H:%M:%S%z')


def save():
    with LOCK:
        for name, value in [('download_manifest.json', list(RECORDS.values())),
                            ('download_status.json', {**STATE, 'updated_at': stamp()})]:
            p = ROOT / name; tmp = p.with_suffix('.json.tmp')
            tmp.write_text(json.dumps(value, indent=2)); tmp.replace(p)


def record(path, source, query, **extra):
    p = ROOT / path
    item = {'source': source, 'path': path, 'request': query, 'retrieved_at': stamp(), **extra}
    if p.exists():
        with p.open('rb') as stream:
            item.update(bytes=p.stat().st_size, sha256=hashlib.file_digest(stream, 'sha256').hexdigest())
    with LOCK:
        RECORDS[path] = item; save()
    return item


def cached(path):
    p = ROOT / path; old = RECORDS.get(path)
    if old and old.get('status') in ('verified', 'not_in_source') and p.exists():
        with p.open('rb') as stream:
            return hashlib.file_digest(stream, 'sha256').hexdigest() == old.get('sha256')
    return False


def uniprot(payload, path, client):
    if cached(path) and RECORDS[path].get('request') == payload:
        return
    destination = ROOT / path
    destination.parent.mkdir(parents=True, exist_ok=True)
    last = None
    for attempt in range(4):
        request = {'base_url': 'https://rest.uniprot.org', **payload, 'save_raw': True,
                   'raw_output_path': str(destination.with_suffix(destination.suffix+'.part')),
                   'max_items': 1, 'max_depth': 1, 'timeout_sec': 180}
        result = client.execute(request)
        if result['ok']:
            temporary = Path(request['raw_output_path'])
            if payload.get('response_format') == 'text':
                if temporary.stat().st_size and not temporary.read_text().startswith('>'):
                    raise ValueError('Expected FASTA: ' + path)
            else:
                data = json.loads(temporary.read_text())
                if not isinstance(data.get('results'), list):
                    raise ValueError('Expected complete UniProt results: ' + path)
                if any(r.get('organism', {}).get('taxonId') != 10090 for r in data['results']):
                    raise ValueError('Non-mouse UniProt record')
            temporary.replace(destination)
            record(path, 'UniProtKB', payload, status='verified', client='resource_http')
            return
        last = result.get('error'); time.sleep(2**attempt)
    record(path, 'UniProtKB', payload, status='failed', error=last)


def interpro(uid):
    # Protein sequence is retained to check release consistency with UniProt.
    for kind, url in [
        ('protein', f'https://www.ebi.ac.uk/interpro/api/protein/uniprot/{uid}/?extra_fields=sequence'),
        ('entries', f'https://www.ebi.ac.uk/interpro/api/entry/all/protein/uniprot/{uid}/?page_size=200')]:
        relative = f'raw/interpro/{uid}.{kind}.json'
        if cached(relative):
            continue
        p = ROOT / relative; p.parent.mkdir(parents=True, exist_ok=True)
        last = None
        for attempt in range(3):
            try:
                pages = []; next_url = url; no_record = False
                with requests.Session() as session:
                    session.headers['User-Agent'] = 'SpliceFate-protein-annotation/1.0'
                    while next_url:
                        response = session.get(next_url, timeout=(10, 45))
                        if response.status_code in (204, 404):
                            pages.append({'source_http_status': response.status_code, 'results': []})
                            no_record = True; break
                        response.raise_for_status(); data = response.json(); pages.append(data)
                        next_url = data.get('next')
                        if next_url and not next_url.startswith('https://www.ebi.ac.uk/interpro/api/'):
                            raise ValueError('Unexpected pagination host')
                tmp = p.with_suffix('.json.tmp'); tmp.write_text(json.dumps({'pages': pages}))
                tmp.replace(p)
                record(relative, 'InterPro', {'url': url},
                       status='not_in_source' if no_record else 'verified')
                break
            except Exception as exc:
                last = str(exc); time.sleep(2**attempt)
        else:
            record(relative, 'InterPro', {'url': url}, status='failed', error=last)


def main():
    ROOT.mkdir(parents=True, exist_ok=True)
    if (ROOT/'download_manifest.json').exists():
        RECORDS.update({r['path']: r for r in json.loads((ROOT/'download_manifest.json').read_text())})
    db = duckdb.connect(str(DB), read_only=True)
    rows = db.execute('''SELECT DISTINCT p.uniprot_id, f.factor_id, f.gene_id, f.gene_symbol,
        f.broad_splicing_candidate FROM rbp2go_profiles p JOIN factor_catalog f USING(factor_id)
        WHERE f.rbp2go_listed OR f.other_rbp_source OR f.broad_splicing_candidate''').fetchall()
    ids = sorted({r[0] for r in rows}, key=lambda u: (not any(r[0] == u and r[4] for r in rows), u))
    old_plan=ROOT/'download_plan.json'
    if old_plan.exists():
        previous=json.loads(old_plan.read_text())['proteins']
        ids=previous+[u for u in ids if u not in set(previous)]
    # Query stable mouse genes for source-supported factors lacking a protein link.
    missing = db.execute('''SELECT DISTINCT f.gene_id FROM factor_catalog f
        WHERE (f.rbp2go_listed OR f.other_rbp_source OR f.broad_splicing_candidate) AND starts_with(f.gene_id,'ENSMUSG')
        AND NOT EXISTS(SELECT 1 FROM rbp2go_profiles p WHERE p.factor_id=f.factor_id)''').fetchall()
    db.close()
    plan = {'scope': 'Source-supported mouse RBPs plus all GO splicing candidates, including non-RNA-binding splice factors; unrelated background-only proteins excluded',
            'proteins': ids, 'factor_protein_links': rows, 'missing_stable_gene_ids': [r[0] for r in missing],
            'sources': ['UniProtKB full JSON and alternative-isoform FASTA', 'InterPro integrated and member-signature positions'],
            'coordinate_policy': 'protein source 1-based inclusive; normalized start0/end0 half-open; fuzzy positions retained'}
    (ROOT/'download_plan.json').write_text(json.dumps(plan, indent=2))
    client = resource_http
    jobs = []
    for n, start in enumerate(range(0, len(ids), 80)):
        query = 'organism_id:10090 AND ('+' OR '.join('accession:'+u for u in ids[start:start+80])+')'
        jobs += [({'path':'uniprotkb/stream','params':{'query':query,'format':'json'}},f'raw/uniprot/batch_{n:03d}.json'),
                 ({'path':'uniprotkb/stream','params':{'query':query,'format':'fasta','includeIsoform':'true'},'response_format':'text'}, f'raw/uniprot/batch_{n:03d}.fasta')]
    for n, start in enumerate(range(0, len(missing), 30)):
        query='organism_id:10090 AND ('+' OR '.join('xref:ensembl-'+r[0] for r in missing[start:start+30])+')'
        jobs.append(({'path':'uniprotkb/stream','params':{'query':query,'format':'json'}},f'raw/uniprot/missing_gene_{n:03d}.json'))
        jobs.append(({'path':'uniprotkb/stream','params':{'query':query,'format':'fasta','includeIsoform':'true'},'response_format':'text'},f'raw/uniprot/missing_gene_{n:03d}.fasta'))
    STATE.update(status='running', started_at=stamp(), uniprot_jobs=len(jobs), interpro_proteins=len(ids)); save()
    # Independent resource lanes; moderate request concurrency per provider.
    with cf.ThreadPoolExecutor(2) as lanes:
        def uni_lane():
            with cf.ThreadPoolExecutor(3) as pool:
                for i, _ in enumerate(pool.map(lambda job: uniprot(*job, client), jobs), 1):
                    with LOCK: STATE['uniprot_completed_jobs']=i; save()
        def ipr_lane():
            with cf.ThreadPoolExecutor(6) as pool:
                for i, _ in enumerate(pool.map(interpro, ids), 1):
                    if i%25 == 0 or i==len(ids):
                        with LOCK: STATE['interpro_completed_proteins']=i; save()
                        print(json.dumps({'interpro_completed_proteins':i,'total':len(ids),'time':stamp()}),flush=True)
        futures=[lanes.submit(uni_lane),lanes.submit(ipr_lane)]
        for future in futures: future.result()
    failures=[r for r in RECORDS.values() if r['status']=='failed']
    STATE.update(status='complete' if not failures else 'partial', failures=failures, finished_at=stamp()); save()
    print(json.dumps(STATE),flush=True)


if __name__ == '__main__':
    main()
