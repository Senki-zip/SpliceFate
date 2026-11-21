"""Monitor source downloads and refresh the additive protein database extension."""
from __future__ import annotations
import collections
import fcntl
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import sys
import time

from collect_mouse_protein_dimensions import ROOT, PROJECT
from build_mouse_rbp_database import ROOT as MAIN, DBPATH, sha256


def stamp():
    return time.strftime('%Y-%m-%dT%H:%M:%S%z')


def active(pid):
    p=Path(f'/proc/{pid}/stat')
    return p.exists() and p.read_text().split(') ',1)[1].split()[0]!='Z'


def save(state):
    state['updated_at']=stamp();p=ROOT/'pipeline_status.json';tmp=p.with_suffix('.json.tmp')
    tmp.write_text(json.dumps(state,indent=2));tmp.replace(p)


def run(name,args):
    with (ROOT/(name+'.log')).open('a') as out:
        subprocess.run([sys.executable,'-u',*args],cwd=PROJECT,stdout=out,stderr=subprocess.STDOUT,check=True)


def export_coverage():
    import duckdb
    from integrate_mouse_protein_dimensions import coverage, validate_extension
    d=duckdb.connect(str(DBPATH),read_only=True);d.execute("SET file_search_path='"+str(MAIN)+"'")
    report={'checked_at':stamp(),'coverage':coverage(d),'checks':validate_extension(d)}
    report['passed']=all(v['passed'] for v in report['checks'].values())
    (ROOT/'validation.json').write_text(json.dumps(report,indent=2))
    cur=d.execute('SELECT * FROM protein_dimension_coverage ORDER BY gene_symbol')
    import csv
    with (ROOT/'factor_protein_coverage.tsv').open('w') as out:
        writer=csv.writer(out,delimiter='\t');writer.writerow([c[0] for c in cur.description]);writer.writerows(cur.fetchall())
    d.close()
    if not report['passed']:raise ValueError('Live protein extension validation failed')


def backup(label):
    mount=Path('/media/huangtao/ae2c2db0-d1db-4f02-9cc6-79d9bd823577')
    import os
    if not os.path.ismount(mount):raise ValueError('10 TB backup disk not mounted')
    # The validated baseline snapshot already contains the immutable sites store.
    # Save a new DB plus its annotations and link their checksums to that snapshot.
    dest=mount/'home/hcq/SpliceFate_backups'/('protein_dimensions_20261120_'+label)
    dest.mkdir(parents=True,exist_ok=True)
    target=dest/'mouse_rbp.duckdb'
    tmp=target.with_suffix('.duckdb.copying');shutil.copy2(DBPATH,tmp)
    expected=sha256(DBPATH)
    if sha256(tmp)!=expected:raise ValueError('Protein dimension backup mismatch')
    tmp.replace(target)
    for p in [MAIN/'SHA256SUMS',MAIN/'protein_dimensions.json',ROOT/'download_manifest.json',
              ROOT/'integration_status.json',ROOT/'validation.json',ROOT/'factor_protein_coverage.tsv']:
        if p.exists():shutil.copy2(p,dest/p.name)
    (dest/'backup.json').write_text(json.dumps({'verified_at':stamp(),'database_sha256':expected,
        'sites_snapshot':str(mount/'home/hcq/SpliceFate_backups/mouse_rbp_rebuilt_20261119/sites'),
        'restore_note':'Copy this DuckDB beside the unchanged sites/ directory from the validated 20261119 snapshot; set file_search_path to that restored root.'},indent=2))
    return str(dest)


def main():
    with (ROOT/'.pipeline.lock').open('w') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
        state={'status':'monitoring','started_at':stamp(),'poll_seconds':30,'source_scope':'all requested mouse RBP and GO splice-factor proteins'}
        initial=ROOT/'integration_initial.pid'
        while initial.exists() and active(int(initial.read_text())):
            save(state);time.sleep(10)
        # Ensure current code and the completed UniProt batch are present.
        run('integration_checkpoint',['scripts/integrate_mouse_protein_dimensions.py','--allow-partial'])
        export_coverage();state['initial_backup']=backup('checkpoint')
        integrated=int(json.loads((ROOT/'download_status.json').read_text()).get('interpro_completed_proteins',0))
        save(state)
        retries=0
        while True:
            source=json.loads((ROOT/'download_status.json').read_text())
            state.update(source_status=source['status'],interpro_completed_proteins=source.get('interpro_completed_proteins',0),
                         interpro_expected_proteins=source['interpro_proteins'],source_failures=len(source.get('failures',[])))
            save(state)
            if source['status']=='complete':
                run('integration_final',['scripts/integrate_mouse_protein_dimensions.py'])
                export_coverage()
                run('main_database_validation',['scripts/validate_mouse_rbp_database.py'])
                state['final_backup']=backup('complete')
                state.update(status='complete',finished_at=stamp());save(state);return
            if source['status']=='partial' and retries<3:
                retries+=1;state.update(status='retrying_sources',retry=retries);save(state)
                run('download_retry',['scripts/collect_mouse_protein_dimensions.py'])
                continue
            if source['status']=='partial':
                run('integration_latest',['scripts/integrate_mouse_protein_dimensions.py','--allow-partial'])
                export_coverage();state['partial_backup']=backup('latest_partial')
                state.update(status='incomplete_provider_failures',failures=source.get('failures',[]));save(state);return
            if not active(int((ROOT/'download.pid').read_text())):
                state.update(status='failed',error='Source downloader stopped before final status');save(state);return
            done=int(source.get('interpro_completed_proteins',0))
            if done-integrated>=1000:
                run('integration_checkpoint',['scripts/integrate_mouse_protein_dimensions.py','--allow-partial'])
                export_coverage();integrated=done;state['last_integrated_proteins']=integrated
            time.sleep(30)


if __name__=='__main__':
    try:main()
    except Exception as exc:
        state={'status':'failed','error':str(exc),'updated_at':stamp()}
        (ROOT/'pipeline_status.json').write_text(json.dumps(state,indent=2));raise
