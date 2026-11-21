"""Download all GEO GSE295932 processed WT/CRISPR files with integrity audits."""
from __future__ import annotations
import concurrent.futures as cf
import csv
import fcntl
import gzip
import hashlib
import json
from pathlib import Path
import re
import threading
import time

import requests

ROOT=Path(__file__).resolve().parents[1]/'test/GSE295932'
LOCK=threading.RLock()
ROWS={}


def stamp():return time.strftime('%Y-%m-%dT%H:%M:%S%z')


def save():
    with LOCK:
        values=list(ROWS.values())
        state={'status':'complete' if values and all(r['status']=='verified' for r in values) else
               'partial' if any(r['status']=='failed' for r in values) else 'running',
               'updated_at':stamp(),'expected_files':len(values),
               'verified_files':sum(r['status']=='verified' for r in values),
               'expected_bytes':sum(r['expected_bytes'] for r in values),
               'downloaded_bytes':sum(r.get('downloaded_bytes',0) for r in values),'files':values}
        for name,data in [('download_status.json',state),('download_manifest.json',[r for r in values if r['status']=='verified'])]:
            p=ROOT/name;tmp=p.with_suffix('.json.tmp');tmp.write_text(json.dumps(data,indent=2));tmp.replace(p)


def download(row):
    name=row['file'];p=ROOT/'files'/name;p.parent.mkdir(parents=True,exist_ok=True)
    part=p.with_name(p.name+'.part')
    for attempt in range(5):
        try:
            if not p.exists():
                offset=part.stat().st_size if part.exists() else 0
                headers={'Accept-Encoding':'identity','User-Agent':'SpliceFate-GEO-download/1.0'}
                if offset:headers['Range']=f'bytes={offset}-'
                with requests.get(row['url'],headers=headers,stream=True,timeout=(20,60)) as response:
                    response.raise_for_status()
                    if offset and response.status_code==200:offset=0
                    if offset:
                        expected=f'bytes {offset}-{row["expected_bytes"]-1}/{row["expected_bytes"]}'
                        if response.status_code!=206 or response.headers.get('Content-Range')!=expected:
                            raise ValueError('Unexpected resume range')
                    with part.open('ab' if offset else 'wb') as out:
                        for block in response.iter_content(512*1024):
                            out.write(block)
                            with LOCK:row.update(status='downloading',downloaded_bytes=out.tell())
                    row.update(final_url=response.url,etag=response.headers.get('ETag'),last_modified=response.headers.get('Last-Modified'))
                if part.stat().st_size!=row['expected_bytes']:raise ValueError('File length differs from GEO filelist')
                target=part
            else:target=p
            if target.stat().st_size!=row['expected_bytes']:raise ValueError('Existing file length mismatch; preserved for review')
            with gzip.open(target,'rb') as stream:
                head=stream.read(1000)
                for _ in iter(lambda:stream.read(8*1024**2),b''):pass
            with target.open('rb') as stream:digest=hashlib.file_digest(stream,'sha256').hexdigest()
            if target==part:part.replace(p)
            with LOCK:
                row.update(status='verified',downloaded_bytes=p.stat().st_size,sha256=digest,
                           verified_at=stamp(),gzip_CRC_passed=True,
                           content_preview=head.decode('utf-8',errors='replace')[:700],
                           official_independent_checksum_available=False)
                row.pop('error',None);save()
            print(json.dumps({'file':name,'status':'verified','bytes':p.stat().st_size}),flush=True)
            return
        except Exception as exc:
            with LOCK:row.update(status='retrying',error=str(exc),attempt=attempt+1);save()
            time.sleep(min(2**attempt,15))
    with LOCK:row['status']='failed';save()


def main():
    ROOT.mkdir(parents=True,exist_ok=True)
    with (ROOT/'.download.lock').open('w') as guard:
        fcntl.flock(guard,fcntl.LOCK_EX|fcntl.LOCK_NB)
        inventory=list(csv.DictReader((ROOT/'metadata/filelist.txt').open(),delimiter='\t'))
        previous={r['file']:r for r in json.loads((ROOT/'download_manifest.json').read_text())} if (ROOT/'download_manifest.json').exists() else {}
        plan=[]
        for record in inventory:
            if record['#Archive/File']!='File':continue
            name=record['Name'];gsm=name.split('_')[0]
            url=f'https://ftp.ncbi.nlm.nih.gov/geo/samples/{gsm[:-3]}nnn/{gsm}/suppl/{name}'
            meta=(ROOT/f'metadata/{gsm}.soft.txt').read_text()
            urls=re.findall(r'!Sample_supplementary_file_\d+ = (\S+)',meta)
            if url.replace('https://','ftp://') not in urls:raise ValueError('File missing from sample metadata: '+name)
            condition='mixed_WT_CRISPR_KO' if int(gsm[3:])<=8961555 else 'WT'
            row={**previous.get(name,{}),'file':name,'sample':gsm,'condition':condition,
                 'stimulated':gsm in ['GSM8961554','GSM8961555'],'url':url,
                 'path':'files/'+name,'expected_bytes':int(record['Size']),'status':'pending'}
            ROWS[name]=row;plan.append(row)
        (ROOT/'download_plan.json').write_text(json.dumps(plan,indent=2));save()
        with cf.ThreadPoolExecutor(4) as pool:
            futures=[pool.submit(download,r) for r in plan]
            while futures:
                done,pending=cf.wait(futures,timeout=15,return_when=cf.FIRST_COMPLETED)
                for f in done:f.result()
                futures=list(pending);save()
        with (ROOT/'SHA256SUMS').open('w') as out:
            for row in sorted(ROWS.values(),key=lambda r:r['file']):
                if row['status']=='verified':out.write(row['sha256']+'  '+row['path']+'\n')
        if any(r['status']!='verified' for r in ROWS.values()):raise SystemExit('Some files failed; see download_status.json')
        print('All 11 processed GEO files verified.',flush=True)


if __name__=='__main__':main()
