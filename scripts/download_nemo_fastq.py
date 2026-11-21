"""Resume validated HTTP ranges directly into TAR files without assembly copies."""
from __future__ import annotations
import concurrent.futures as cf
from datetime import datetime, timezone
import fcntl
import hashlib
import json
import os
from pathlib import Path
import shutil
import tarfile
import threading
import time
from urllib.parse import urlparse
import requests

ROOT=Path(__file__).resolve().parents[1]/'test'
CHUNK=32*1024**2
WORKERS=2
LOCK=threading.Lock()
PROGRESS={}
START=time.monotonic()


def now():return datetime.now(timezone.utc).isoformat()


def save(status='downloading',error=None):
    with LOCK:
        value={'status':status,'updated_at':now(),'elapsed_seconds':round(time.monotonic()-START),
            'files':list(PROGRESS.values()),'error':error,'workers':WORKERS,'chunk_bytes':CHUNK,
            'completed_range_bytes':sum(v.get('completed_range_bytes',0) for v in PROGRESS.values()),
            'expected_bytes':sum(v['expected_bytes'] for v in PROGRESS.values()),
            'free_disk_bytes':shutil.disk_usage(ROOT).free}
        temp=ROOT/'download_status.json.tmp';temp.write_text(json.dumps(value,indent=2));temp.replace(ROOT/'download_status.json')


def fetch_range(job):
    row,left,right=job
    name=Path(urlparse(row['url']).path).name
    key=f'{left}-{right}'
    part=ROOT/(name+'.part');journal=ROOT/(name+'.ranges.jsonl')
    attempt=0
    while True:
        try:
            sha=hashlib.sha256();written=0;started=time.monotonic()
            with requests.get(row['url'],headers={'Range':f'bytes={left}-{right}',
                 'If-Range':row['etag'],'Accept-Encoding':'identity','User-Agent':'SpliceFate-public-data-download/1.0'},
                 stream=True,timeout=(20,90)) as r:
                r.raise_for_status()
                if r.status_code!=206 or r.headers.get('Content-Range')!=f"bytes {left}-{right}/{row['bytes']}":
                    raise ValueError('Server did not honor exact range or file changed')
                fd=os.open(part,os.O_WRONLY|os.O_CREAT,0o664)
                try:
                    for block in r.iter_content(256*1024):
                        if not block:continue
                        if written+len(block)>right-left+1:raise ValueError('Range overflow')
                        if shutil.disk_usage(ROOT).free<1024**3:raise OSError('Less than 1 GiB free; stopping without deleting files')
                        view=memoryview(block);offset=left+written
                        while view:
                            count=os.pwrite(fd,view,offset)
                            if count<=0:raise OSError('Short disk write')
                            offset+=count;view=view[count:]
                        sha.update(block);written+=len(block)
                        if time.monotonic()-started>1800:raise TimeoutError('Range too slow; retrying')
                    if written!=right-left+1:raise ValueError('Truncated range')
                    os.fsync(fd)
                finally:os.close(fd)
            rec={'range':key,'start':left,'end':right,'sha256':sha.hexdigest()}
            with LOCK:
                with journal.open('a') as f:
                    f.write(json.dumps(rec)+'\n');f.flush();os.fsync(f.fileno())
                PROGRESS[name]['completed_range_bytes']+=written
                PROGRESS[name]['completed_ranges']+=1
                PROGRESS[name].pop('last_error',None)
                PROGRESS[name].pop('retry_attempt',None)
            return
        except Exception as exc:
            attempt+=1
            transient=isinstance(exc,(requests.ConnectionError,requests.Timeout,TimeoutError))
            if isinstance(exc,requests.HTTPError):
                transient=exc.response is not None and (exc.response.status_code in (408,429) or exc.response.status_code>=500)
            print(json.dumps({'time':now(),'file':name,'range':key,'attempt':attempt,'retryable':transient,'error':str(exc)}),flush=True)
            if not transient:raise
            with LOCK:
                PROGRESS[name]['last_error']=str(exc)
                PROGRESS[name]['retry_attempt']=attempt
                PROGRESS[name]['last_error_at']=now()
            # Persistent network outages wait and resume; disk/identity failures stop.
            time.sleep(min(2**min(attempt,6),60))


def main():
    ROOT.mkdir(exist_ok=True)
    scope_file=ROOT/'download_scope.json'
    if scope_file.exists() and json.loads(scope_file.read_text()).get('status')=='cancelled_by_user':
        raise SystemExit('Raw download cancelled by user; use the processed dataset downloader instead.')
    lock=(ROOT/'.download.lock').open('w')
    fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
    (ROOT/'download.pid').write_text(str(os.getpid())+'\n')
    rows=json.loads((ROOT/'download_discovery.json').read_text())
    (ROOT/'download_status.json').write_text(json.dumps({'status':'verifying_partial_ranges','updated_at':now(),
        'pid':os.getpid(),'message':'Checking saved per-range SHA256 before resuming; partial files are preserved.'},indent=2))
    state_path=ROOT/'download_manifest.json'
    manifest=json.loads(state_path.read_text()) if state_path.exists() else []
    expected_urls=[r['url'] for r in rows]
    (ROOT/'download_scope.json').write_text(json.dumps({'authorization':'User confirmed all 5 URLs','urls':expected_urls,'extract':False},indent=2))
    jobs=[]
    for row in rows:
        if row['status']!=200 or not row.get('etag'):raise ValueError('Missing provider size/ETag')
        row['bytes']=int(row['bytes']);name=Path(urlparse(row['url']).path).name
        final=ROOT/name;part=ROOT/(name+'.part');journal=ROOT/(name+'.ranges.jsonl')
        if final.exists():
            rec=next((r for r in manifest if r['file']==name and r['url']==row['url']),None)
            if not rec or final.stat().st_size!=row['bytes']:raise ValueError('Existing file lacks matching verified manifest: '+name)
            PROGRESS[name]={'file':name,'url':row['url'],'expected_bytes':row['bytes'],'status':'verified',
                             'completed_range_bytes':row['bytes'],'sha256':rec['sha256']};continue
        identity=ROOT/(name+'.source.json')
        if identity.exists() and json.loads(identity.read_text())!=row:raise ValueError('Provider identity changed; preserve existing partial: '+name)
        if not identity.exists():identity.write_text(json.dumps(row,indent=2))
        complete={}
        if journal.exists():
            for line in journal.read_text().splitlines():
                try:r=json.loads(line)
                except json.JSONDecodeError:continue
                complete[r['range']]=r
        if complete and not part.exists():raise ValueError('Range journal exists without partial file')
        valid={}
        if complete:
            with part.open('rb') as f:
                for key,r in complete.items():
                    f.seek(r['start']);b=f.read(r['end']-r['start']+1)
                    if hashlib.sha256(b).hexdigest()==r['sha256']:valid[key]=r
        remaining=0
        for left in range(0,row['bytes'],CHUNK):
            right=min(row['bytes']-1,left+CHUNK-1)
            if f'{left}-{right}' not in valid:jobs.append((row,left,right));remaining+=right-left+1
        PROGRESS[name]={'file':name,'url':row['url'],'expected_bytes':row['bytes'],'status':'downloading',
            'completed_range_bytes':row['bytes']-remaining,'completed_ranges':len(valid),
            'expected_ranges':(row['bytes']+CHUNK-1)//CHUNK}
    # Range writes are sparse; reserve only additional bytes, not file apparent sizes.
    remaining=sum(r['expected_bytes']-r['completed_range_bytes'] for r in PROGRESS.values())
    if shutil.disk_usage(ROOT).free<remaining+5*1024**3:raise OSError('Insufficient disk for requested archives plus 5 GiB margin')
    # Interleave files to keep useful progress on every requested archive.
    jobs.sort(key=lambda j:(j[1],j[0]['url']))
    save();print(json.dumps({'started':now(),'jobs':len(jobs),'remaining_bytes':remaining}),flush=True)
    try:
        with cf.ThreadPoolExecutor(max_workers=WORKERS) as pool:
            iterator=iter(jobs)
            futures={pool.submit(fetch_range,job) for job in [next(iterator,None) for _ in range(WORKERS)] if job is not None}
            finished_ranges=0
            while futures:
                done,futures=cf.wait(futures,timeout=20,return_when=cf.FIRST_COMPLETED)
                for f in done:
                    f.result();finished_ranges+=1
                    job=next(iterator,None)
                    if job is not None:futures.add(pool.submit(fetch_range,job))
                save()
                print(json.dumps({'time':now(),'completed_bytes':sum(v['completed_range_bytes'] for v in PROGRESS.values()),'remaining_ranges':len(jobs)-finished_ranges}),flush=True)
        for row in rows:
            name=Path(urlparse(row['url']).path).name
            if PROGRESS[name]['status']=='verified':continue
            part=ROOT/(name+'.part')
            if part.stat().st_size!=row['bytes']:raise ValueError('Final size mismatch: '+name)
            PROGRESS[name]['status']='checking';save('checking')
            with part.open('rb') as f:sha=hashlib.file_digest(f,'sha256').hexdigest()
            members=0;payload=0
            with tarfile.open(part,'r:') as tar:
                for member in tar:
                    if member.size<0 or member.offset_data+member.size>row['bytes']:raise ValueError('Invalid/truncated TAR member')
                    members+=1
                    if member.isfile():payload+=member.size
            if not members:raise ValueError('Empty TAR archive')
            part.replace(ROOT/name)
            rec={'file':name,'url':row['url'],'bytes':row['bytes'],'etag':row['etag'],'last_modified':row['last_modified'],
                'sha256':sha,'tar_members':members,'tar_payload_bytes':payload,'verified_at':now(),
                'checks':'exact HTTP ranges; provider size/ETag; per-range SHA256; local whole-file SHA256; TAR headers/bounds',
                'official_checksum_available':False}
            manifest=[r for r in manifest if r['file']!=name]+[rec]
            state_path.write_text(json.dumps(manifest,indent=2))
            (ROOT/'SHA256SUMS').write_text(''.join(r['sha256']+'  '+r['file']+'\n' for r in manifest))
            PROGRESS[name].update(status='verified',sha256=sha,tar_members=members)
            save('checking');print(json.dumps(rec),flush=True)
        save('complete')
    except BaseException as exc:
        save('failed',str(exc));raise


if __name__=='__main__':main()
