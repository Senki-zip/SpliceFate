"""Download mouse PPI/complex resources and preserve source provenance."""
import concurrent.futures as cf
import gzip
import hashlib
import json
from pathlib import Path
import threading
import time
import zipfile
import requests

ROOT=Path(__file__).resolve().parents[1]/'data/external/rbp_interactions_20261116'
LOCK=threading.Lock()
JOBS=[]
for typ in ['protein.physical.links.full','protein.links.detailed','protein.info','protein.aliases']:
    name=f'10090.{typ}.v12.5.txt.gz'
    JOBS.append(('STRING',f'https://stringdb-downloads.org/download/{typ}.v12.5/{name}','raw/string/'+name))
for tax in ['10090','9606']:
    JOBS.append(('ComplexPortal',f'https://ftp.ebi.ac.uk/pub/databases/intact/complex/current/complextab/{tax}.tsv',f'raw/complexportal/{tax}.tsv'))
JOBS.append(('IntAct','https://ftp.ebi.ac.uk/pub/databases/intact/current/psimitab/species/mouse.zip','raw/intact/mouse.zip'))
JOBS.append(('IntAct_negative','https://ftp.ebi.ac.uk/pub/databases/intact/current/psimitab/species/mouse_negative.txt','raw/intact/mouse_negative.txt'))
for species,tax in [('mouse','10090'),('human','9606')]:
    JOBS.append(('SpliceosomeDB',f'https://spliceosomedb.ucsc.edu/proteins/export/?organism={tax}',f'raw/spliceosomedb/proteins_{species}.csv'))
JOBS.append(('SpliceosomeDB','https://spliceosomedb.ucsc.edu/experiments/export/?quick_search=','raw/spliceosomedb/experiments_all.csv'))


def download(job):
    source,url,rel=job;path=ROOT/rel;path.parent.mkdir(parents=True,exist_ok=True)
    part=path.with_name(path.name+'.part')
    error=''
    for attempt in range(5):
        try:
            meta={}
            if not path.exists():
                offset=part.stat().st_size if part.exists() else 0
                headers={'Accept-Encoding':'identity','User-Agent':'SpliceFate-resource-cache/1.0'}
                if offset:headers['Range']=f'bytes={offset}-'
                with requests.get(url,headers=headers,stream=True,timeout=(20,45)) as r:
                    r.raise_for_status()
                    if offset and (r.status_code!=206 or not r.headers.get('Content-Range','').startswith(f'bytes {offset}-')):
                        if r.status_code==200:offset=0
                        else:raise ValueError('Unexpected resume response')
                    expected=int(r.headers['Content-Length']) if r.headers.get('Content-Length') else None
                    n=0
                    with part.open('ab' if offset else 'wb') as f:
                        for block in r.iter_content(256*1024):
                            f.write(block);n+=len(block)
                    if expected is not None and n!=expected:raise ValueError('Truncated HTTP payload')
                    meta={'final_url':r.url,'etag':r.headers.get('ETag'),'last_modified':r.headers.get('Last-Modified')}
                with part.open('rb') as f:
                    if f.read(128).lstrip().lower().startswith((b'<html',b'<!doctype')):raise ValueError('HTML instead of data')
                if path.name.endswith('.gz'):
                    with gzip.open(part,'rb') as f:
                        for _ in iter(lambda:f.read(4*1024**2),b''):pass
                if path.name.endswith('.zip'):
                    with zipfile.ZipFile(part) as z:
                        bad=z.testzip()
                        if bad:raise ValueError('ZIP CRC failure '+bad)
                part.replace(path)
            with path.open('rb') as f:sha=hashlib.file_digest(f,'sha256').hexdigest()
            rec={'source':source,'url':url,'path':rel,'bytes':path.stat().st_size,'sha256':sha,
                 'status':'verified','compression_integrity_checked':path.name.endswith(('.gz','.zip')),
                 'timestamp':time.strftime('%Y-%m-%dT%H:%M:%S%z'),**meta}
            with LOCK:
                old=ROOT/'download_manifest.json'
                rows=json.loads(old.read_text()) if old.exists() else []
                prior=next((x for x in rows if x['path']==rel),None)
                if prior and not meta:rec={**prior,**rec}
                rows=[x for x in rows if x['path']!=rel]+[rec]
                old.write_text(json.dumps(rows,indent=2))
            print(json.dumps(rec),flush=True);return rec
        except Exception as exc:
            error=str(exc)
            print(json.dumps({'retry':attempt+1,'path':rel,'error':error}),flush=True)
            time.sleep(min(2**attempt,15))
    return {'source':source,'url':url,'path':rel,'status':'failed','error':error}


def main():
    ROOT.mkdir(parents=True,exist_ok=True)
    (ROOT/'download_plan.json').write_text(json.dumps(JOBS,indent=2))
    (ROOT/'download_status.json').write_text(json.dumps({'status':'running','expected':len(JOBS)}))
    with cf.ThreadPoolExecutor(max_workers=5) as p:rows=list(p.map(download,JOBS))
    state={'status':'complete' if all(x['status']=='verified' for x in rows) else 'partial',
           'verified':sum(x['status']=='verified' for x in rows),'expected':len(rows),'results':rows}
    (ROOT/'download_status.json').write_text(json.dumps(state,indent=2))
    print(json.dumps(state),flush=True)


if __name__=='__main__':main()
