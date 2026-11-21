"""Download minimal P21 author processed data, resume ranges, verify gzip/TAR."""
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
import tarfile
import threading
import time
import requests

PROJECT=Path(__file__).resolve().parents[1]
ROOT=PROJECT/'test/processed'
DISCOVERY=PROJECT/'test/processed_discovery'
LOCK=threading.Lock();STATE={};CHUNK=8*1024**2;START=time.monotonic()


def stamp():return time.strftime('%Y-%m-%dT%H:%M:%S%z')


def atomic(name,value):
    p=ROOT/name;t=p.with_suffix(p.suffix+'.tmp');t.write_text(json.dumps(value,indent=2));t.replace(p)


def status():
    with LOCK:
        vals=list(STATE.values());complete=all(v['status']=='verified' for v in vals)
        atomic('download_status.json',{'status':'complete' if complete else 'running','updated_at':stamp(),
             'elapsed_seconds':round(time.monotonic()-START),'verified_files':sum(v['status']=='verified' for v in vals),
             'total_files':len(vals),'downloaded_bytes':sum(v.get('downloaded_bytes',0) for v in vals),
             'expected_bytes':sum(v.get('expected_bytes',0) for v in vals),'files':vals,
             'raw_FASTQ_download':'cancelled_and_deleted_by_user','free_disk_bytes':shutil.disk_usage(ROOT).free})


def build_jobs():
    jobs=[]
    for row in json.loads((DISCOVERY/'P21_processed_files.json').read_text()):
        if row['file'].endswith('_AllInfo.tsv.tar'):
            jobs.append({'file':row['file'],'url':row['url'],'expected_bytes':row['bytes'],'kind':'AllInfo'})
    for row in json.loads((DISCOVERY/'P21_shortread_processed_files.json').read_text()):
        if row.get('type')=='counts':
            jobs.append({'file':row['sample']+'.mex.tar.gz','url':row['url'],'expected_bytes':row['bytes'],'kind':'gene_counts'})
    # Pin metadata source to the author's current repository commit, not a moving raw URL.
    version_file=ROOT/'author_repository.json'
    if version_file.exists():meta=json.loads(version_file.read_text())
    else:
        r=requests.get('https://api.github.com/repos/noush-joglekar/biccn_tilgner_scisorseq/commits/main',timeout=(15,45));r.raise_for_status()
        meta={'repository':'https://github.com/noush-joglekar/biccn_tilgner_scisorseq','commit':r.json()['sha'],'retrieved_at':stamp()};atomic('author_repository.json',meta)
    for name,rel,kind in [('SR_fullData_combined_barcodeCTassignments.gz','data/SR_fullData_combined_barcodeCTassignments.gz','cell_annotation'),
                          ('AUTHOR_README.md','README.md','source_documentation')]:
        jobs.append({'file':name,'url':f"https://raw.githubusercontent.com/noush-joglekar/biccn_tilgner_scisorseq/{meta['commit']}/{rel}",'expected_bytes':0,'kind':kind})
    return jobs


def validate(path,job):
    name=job['file'];report={'validation':'local_size_SHA256_and_compression_CRC','official_checksum_available':False}
    if job['expected_bytes'] and path.stat().st_size!=job['expected_bytes']:raise ValueError('Provider size mismatch')
    with path.open('rb') as f:report['sha256']=hashlib.file_digest(f,'sha256').hexdigest()
    if name.endswith('.tar'):
        members=[];preview=[]
        with tarfile.open(path,'r:') as tar:
            for m in tar:
                if m.offset_data+m.size>path.stat().st_size:raise ValueError('Truncated TAR member')
                members.append({'name':m.name,'bytes':m.size,'is_file':m.isfile()})
                if m.isfile() and m.name.endswith('.tsv.gz'):
                    with tar.extractfile(m) as stream,gzip.GzipFile(fileobj=stream) as gz:
                        for line in gz:
                            if len(preview)<3:preview.append(line.decode().rstrip('\n'))
                        # Iteration to EOF validates embedded gzip CRC without extracting.
        if not preview:raise ValueError('No readable AllInfo gzip payload')
        widths={len(s.split('\t')) for s in preview}
        if widths!={9}:raise ValueError('Unexpected AllInfo preview field count: '+str(widths))
        report.update(tar_members=members,AllInfo_preview=preview,field_count=9,embedded_gzip_CRC_passed=True)
    elif name.endswith('.mex.tar.gz'):
        members=[];matrix_header=[]
        with tarfile.open(path,'r|gz') as tar:
            for m in tar:
                members.append(m.name)
                if m.isfile():
                    with tar.extractfile(m) as stream:
                        first=stream.read(65536)
                        if m.name.endswith('matrix.mtx'):matrix_header=first.decode(errors='replace').splitlines()[:4]
                        while stream.read(4*1024**2):pass
        # tar streaming can stop at end-of-archive before outer gzip trailer: check EOF explicitly.
        with gzip.open(path,'rb') as gz:
            while gz.read(4*1024**2):pass
        if not any('matrix.mtx' in m for m in members) or not any('barcodes.tsv' in m for m in members) or not any('features.tsv' in m or 'genes.tsv' in m for m in members):
            raise ValueError('Missing 10x MEX components')
        report.update(tar_members=members,matrix_header=matrix_header,outer_gzip_CRC_passed=True)
    elif name.endswith('.gz'):
        n=0;preview=[]
        with gzip.open(path,'rt') as f:
            for line in f:
                if len(preview)<3:preview.append(line.rstrip('\n'))
                n+=1
        if not n:raise ValueError('Empty cell annotation')
        report.update(annotation_rows=n,preview=preview,gzip_CRC_passed=True,
                      cell_key='sample identity parsed from label plus barcode; barcode alone not globally unique')
    elif not path.read_text().startswith('#'):raise ValueError('Invalid author README')
    return report


def download(job):
    name=job['file'];dest=ROOT/name;part=ROOT/(name+'.part');sourcefile=ROOT/(name+'.source.json')
    def update(**kw):
        with LOCK:STATE[name].update(kw)
    oldmanifest=ROOT/'download_manifest.json'
    records=json.loads(oldmanifest.read_text()) if oldmanifest.exists() else []
    previous=next((r for r in records if r['file']==name),None)
    if dest.exists() and previous and dest.stat().st_size==previous['bytes']:
        with dest.open('rb') as f:sha=hashlib.file_digest(f,'sha256').hexdigest()
        if sha==previous['sha256']:update(status='verified',downloaded_bytes=previous['bytes']);return
        raise ValueError('Existing final checksum mismatch; preserving '+name)
    source=json.loads(sourcefile.read_text()) if sourcefile.exists() else None
    attempt=0
    while not dest.exists():
        try:
            offset=part.stat().st_size if part.exists() else 0
            total=source['bytes'] if source else job['expected_bytes']
            if total and offset==total:break
            if total and offset>total:raise ValueError('Partial file larger than source')
            right=min(offset+CHUNK-1,total-1) if total else offset+CHUNK-1
            headers={'Range':f'bytes={offset}-{right}','Accept-Encoding':'identity','User-Agent':'SpliceFate-processed-data/1.0'}
            if source and source.get('etag'):headers['If-Range']=source['etag']
            update(status='downloading')
            with requests.get(job['url'],headers=headers,stream=True,timeout=(20,90)) as r:
                r.raise_for_status()
                if r.status_code==206:
                    match=re.fullmatch(r'bytes (\d+)-(\d+)/(\d+)',r.headers.get('Content-Range',''))
                    if not match or int(match[1])!=offset:raise ValueError('Invalid HTTP range response')
                    actual_total=int(match[3]);response_bytes=int(match[2])-offset+1
                elif r.status_code==200 and offset==0:
                    actual_total=int(r.headers.get('Content-Length') or 0);response_bytes=actual_total or None
                else:raise ValueError('Server ignored resume range; keeping partial file')
                if total and actual_total!=total:raise ValueError('Provider file size changed')
                etag=r.headers.get('ETag')
                if source and source.get('etag') and etag!=source['etag']:raise ValueError('Provider ETag changed')
                source={'url':job['url'],'bytes':actual_total,'etag':etag,'last_modified':r.headers.get('Last-Modified'),'final_url':r.url}
                sourcefile.write_text(json.dumps(source,indent=2));total=actual_total
                update(expected_bytes=total)
                n=0
                with part.open('ab') as f:
                    for block in r.iter_content(256*1024):
                        if shutil.disk_usage(ROOT).free<5*1024**3:raise OSError('Low disk: preserving partial downloads')
                        f.write(block);n+=len(block);update(downloaded_bytes=offset+n)
                if response_bytes is not None and n!=response_bytes:raise requests.ConnectionError('Truncated transfer')
            attempt=0
            if not total or part.stat().st_size==total:break
        except Exception as exc:
            transient=isinstance(exc,(requests.ConnectionError,requests.Timeout))
            if isinstance(exc,requests.HTTPError):transient=exc.response.status_code in (408,429) or exc.response.status_code>=500
            update(last_error=str(exc),retry_attempt=attempt+1)
            print(json.dumps({'file':name,'error':str(exc),'transient':transient,'time':stamp()}),flush=True)
            if not transient:raise
            attempt+=1;time.sleep(min(2**min(attempt,6),60))
    update(status='checking_compression_and_content')
    file=dest if dest.exists() else part
    report=validate(file,job)
    report.update(file=name,url=job['url'],kind=job['kind'],bytes=file.stat().st_size,verified_at=stamp(),source=source)
    if file==part:part.replace(dest)
    with LOCK:
        manifests=json.loads(oldmanifest.read_text()) if oldmanifest.exists() else []
        manifests=[r for r in manifests if r['file']!=name]+[report]
        atomic('download_manifest.json',manifests)
        (ROOT/'SHA256SUMS').write_text(''.join(r['sha256']+'  '+r['file']+'\n' for r in manifests))
        STATE[name].update(status='verified',downloaded_bytes=report['bytes'],sha256=report['sha256'])
        STATE[name].pop('last_error',None)
    print(json.dumps({'verified':name,'bytes':report['bytes'],'time':stamp()}),flush=True)
    status()


def main():
    ROOT.mkdir(exist_ok=True,parents=True)
    lock=(ROOT/'.download.lock').open('w');fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
    (ROOT/'download.pid').write_text(str(os.getpid())+'\n')
    jobs=build_jobs();atomic('download_plan.json',jobs)
    for job in jobs:STATE[job['file']]={**job,'status':'pending','downloaded_bytes':0}
    status()
    try:
        with cf.ThreadPoolExecutor(max_workers=2) as pool:
            pending={pool.submit(download,job) for job in jobs}
            while pending:
                done,pending=cf.wait(pending,timeout=20,return_when=cf.FIRST_COMPLETED)
                for f in done:f.result()
                status()
                print(json.dumps({'time':stamp(),'remaining_files':len(pending),'downloaded_bytes':sum(v.get('downloaded_bytes',0) for v in STATE.values())}),flush=True)
    except BaseException as exc:
        failure={'error':str(exc),'time':stamp()}
        atomic('failure.json',failure)
        old=json.loads((ROOT/'download_status.json').read_text())
        old.update(status='failed',error=str(exc),updated_at=stamp())
        atomic('download_status.json',old)
        raise
    status();print('All processed files downloaded and verified.',flush=True)


if __name__=='__main__':main()
