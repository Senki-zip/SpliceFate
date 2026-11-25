"""Project held-out author KO structures onto the fixed WT common-anchor axis."""
from __future__ import annotations
import argparse
import collections
import csv
import json
from pathlib import Path
import sys

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
import duckdb
import numpy as np
import pandas as pd
from scipy import sparse

from convert_gse295932 import parse_blocks,import_source,condition_samples,condition_output
from splicefate.convert import observation_path_groups,count_event_observations
from splicefate.io import digest,dump

PROJECT=Path(__file__).resolve().parents[1]


def project(source,processed,output,condition='stimulated',min_cells=8,min_molecules=20,prior_dir=None):
    samples=condition_samples(condition)
    if len(samples)!=1:raise ValueError('Project each condition independently')
    gsm,_,_=samples[0];stem=gsm+('_FT_stim' if condition=='stimulated' else '_FT')
    raw=source/'files'/(stem+'_sc_iso_count.txt.gz');labels_path=source/'files'/(stem+'_CRISPR_barcode.tsv.gz')
    labels=pd.read_csv(labels_path,sep='\t',dtype=str)
    if labels.barcode.duplicated().any() or labels.target.isna().any():raise ValueError('Invalid official barcode labels')
    cells=[gsm+':'+b for b in labels.barcode];axis={cell:i for i,cell in enumerate(cells)}
    cell_targets=labels.target.to_numpy();targets=sorted(set(cell_targets)-{'Nontarget'})
    events=json.loads((processed/'events/groups.json').read_text());offsets=np.load(processed/'events/branch_offsets.npy')
    C=sparse.load_npz(processed/'events/branch_counts.npz').tocsc()
    wt_cells=pd.read_csv(processed/'events/cells.tsv',sep='\t',dtype=str).cell.tolist()
    if not set(wt_cells)<=set(cells) or any(cell_targets[axis[cell]]!='Nontarget' for cell in wt_cells):
        raise ValueError('Processed WT axis is not the official same-condition Nontarget group')
    observed_wt=set(wt_cells)
    by_gene=collections.defaultdict(list)
    for e,event in enumerate(events):by_gene[event['gene']].append((e,event))
    output.mkdir(parents=True,exist_ok=True);db=duckdb.connect(str(output/'projection_staging.duckdb'))
    db.execute('SET threads=4');db.execute("SET memory_limit='8GB'")
    db.execute('CREATE OR REPLACE TABLE observations(sample VARCHAR,cell VARCHAR,gene VARCHAR,isoform VARCHAR,count BIGINT)')
    source_audit=import_source(db,raw,gsm);db.execute('CREATE INDEX IF NOT EXISTS projection_gene ON observations(gene)')
    audit=collections.Counter();coverage=collections.Counter();rows_written=0
    columns=['targets','strength','event','gene_id','cell_stage','branch_index','WT_branch_molecules','KO_branch_molecules',
        'WT_event_molecules','KO_event_molecules','WT_observed_cells','KO_observed_cells',
        'probability_WT_truth','probability_KO_truth','delta_probability_truth','truth_status']
    with (output/'ko_truth.tsv').open('w',newline='') as handle:
        writer=csv.DictWriter(handle,fieldnames=columns,delimiter='\t');writer.writeheader()
        for gene,items in by_gene.items():
            observations=[];cache={}
            for cell,text,count in db.execute('SELECT cell,isoform,sum(count) FROM observations WHERE gene=? GROUP BY cell,isoform',[gene]).fetchall():
                count=int(count)
                if cell not in axis:audit['unassigned_clusters_excluded']+=count;continue
                if not count:continue
                if text not in cache:
                    try:cache[text]=parse_blocks(text)
                    except ValueError:cache[text]=None
                blocks=cache[text]
                if blocks is None:audit['invalid_structure_clusters_excluded']+=count;continue
                junctions=[(a[1],b[0]) for a,b in zip(blocks,blocks[1:])]
                observations.append((axis[cell],count,('', '',junctions,blocks)))
            groups=observation_path_groups(observations)
            for e,event in items:
                counts=count_event_observations(groups,event);k=len(event['branches'])
                dense=np.zeros((len(cells),k),np.int64)
                for (i,b),value in counts.items():dense[i,b]=value
                wt=dense[[axis[cell] for cell in wt_cells]]
                if not np.array_equal(wt,C[:,offsets[e]:offsets[e+1]].toarray()):
                    raise ValueError('Projected WT counts differ from fixed training axis at '+event['event'])
                # Cells labelled Nontarget but absent from the converted source axis
                # cannot contribute external WT reference counts.
                wt_sum=wt.sum(0);wt_total=int(wt_sum.sum());wt_n=int((wt.sum(1)>0).sum())
                for target in targets:
                    ko=dense[cell_targets==target];ko_sum=ko.sum(0);ko_total=int(ko_sum.sum());ko_n=int((ko.sum(1)>0).sum())
                    enough=wt_n>=min_cells and ko_n>=min_cells and min(wt_total,ko_total)>=min_molecules
                    status='ok' if enough else 'insufficient_observed_counts'
                    coverage[target,'events_observed']+=int(ko_total>0);coverage[target,'events_evaluable']+=int(enough)
                    for b in range(k):
                        pw=float(wt_sum[b]/wt_total) if wt_total else np.nan;pk=float(ko_sum[b]/ko_total) if ko_total else np.nan
                        writer.writerow(dict(zip(columns,[target,1.,event['event'],gene,condition,b,int(wt_sum[b]),int(ko_sum[b]),
                            wt_total,ko_total,wt_n,ko_n,pw,pk,pk-pw if enough else np.nan,status])))
                        rows_written+=1
    db.close()
    prior_counts={};universe=set()
    if prior_dir is not None:
        from splicefate.storage import load_priors
        prior=load_priors(prior_dir);universe=set(prior['symbols'])
        ids={s:f for f,s in zip(prior['factor_ids'],prior['symbols'])}
        with duckdb.connect(str(prior_dir/'priors.duckdb'),read_only=True) as prior_db:
            for target in ['BUB3','VIM']:
                prior_counts[target]=prior_db.execute('SELECT count(*) FROM event_candidates WHERE factor_id=?',[ids.get(target,'')]).fetchone()[0]
    rows=[]
    for target in targets:
        rows.append({'target':target,'official_label_cells':int((cell_targets==target).sum()),
            'events_observed':coverage[target,'events_observed'],'events_evaluable':coverage[target,'events_evaluable'],
            'audit_role':'negative_control' if target=='VIM' else 'external_coverage' if target=='BUB3' else 'splice_factor_KO',
            'in_prepared_factor_universe':target in universe if prior_dir else None,
            'prior_event_paths':prior_counts.get(target)})
    pd.DataFrame(rows).to_csv(output/'target_coverage.tsv',sep='\t',index=False)
    pd.DataFrame([r for r in rows if r['target'] in ('BUB3','VIM')]).to_csv(output/'BUB3_VIM_external_audit.tsv',sep='\t',index=False)
    report={'condition':condition,'sample':gsm,'WT_cells':len(wt_cells),'KO_targets':targets,'events':len(events),
        'rows':rows_written,'WT_integer_counts_exactly_reproduced':True,'filtering':dict(audit),
        'min_cells':min_cells,'min_molecules':min_molecules,'count_field':'cluster','source_audit':source_audit,
        'fixed_axis':'WT groups.json anchors and exact branches; no KO-driven event discovery',
        'training_or_selection_use':False,'strength':1.,
        'sources':{p.name:digest(p) for p in [raw,labels_path,processed/'events/groups.json',processed/'events/branch_counts.npz']}}
    dump(output/'projection_report.json',report)
    print(json.dumps(report,indent=2),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--source',type=Path,default=PROJECT/'test/GSE295932')
    p.add_argument('--condition',choices=['stimulated','unstimulated'],default='stimulated')
    p.add_argument('--processed-dir',type=Path);p.add_argument('--output',type=Path)
    p.add_argument('--min-cells',type=int,default=8);p.add_argument('--min-molecules',type=int,default=20)
    p.add_argument('--prior-dir',type=Path);a=p.parse_args()
    if a.min_cells<1 or a.min_molecules<1:p.error('Positive count support thresholds required')
    a.processed_dir=a.processed_dir or condition_output(a.condition)
    a.output=a.output or a.processed_dir/'external_KO_truth'
    project(a.source,a.processed_dir,a.output,a.condition,a.min_cells,a.min_molecules,a.prior_dir)
