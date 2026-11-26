"""Observed-chain reconstruction with event-specific eligibility and simplex constraints."""
from __future__ import annotations
import json
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import sparse
from scipy.optimize import minimize
from scipy.special import softmax,log_softmax


def chain_event_branch(chain,event):
    exons=json.loads(chain['exons']) if isinstance(chain['exons'],str) else chain['exons']
    # No unobserved flank is treated as exclusion.
    if not all(any(int(a)<=left and int(b)>=right for a,b in exons) for left,right in event['anchors']):return -1
    junctions=[]
    tokens=json.loads(chain['structure']) if isinstance(chain['structure'],str) else chain['structure']
    for token in tokens:
        if not token.startswith('junction:'):continue
        interval=token.split(':')[2];a,b=map(int,interval.split('-'))
        if a>=event['left'] and b<=event['right']:junctions.append((a,b))
    observed=sorted(junctions)
    matches=[i for i,branch in enumerate(event['branches']) if observed==sorted(map(tuple,branch))]
    return matches[0] if len(matches)==1 else -1


def load_structures(data,models,c):
    events_by_gene = {}
    for m in models.values():
        if m['status'] == 'ok':
            events_by_gene.setdefault(m['gene_id'], []).append(m['event_index'])
    genes = set(events_by_gene)
    root=Path(c['input']['genes']);counts=sparse.load_npz(root/'isoform_counts.npz').tocsr()
    records_by_gene={}
    for chunk in pd.read_csv(root/'isoforms.tsv',sep='\t',chunksize=20000):
        chosen=chunk[chunk.gene.isin(genes)]
        for index,r in zip(chosen.index,chosen.to_dict('records')):records_by_gene.setdefault(r['gene'],[]).append((index,r))
    totals=np.asarray(counts.sum(0)).ravel();result={}
    for gene in sorted(genes):
        group=records_by_gene.get(gene,[])
        group.sort(key=lambda x:(-totals[x[0]],x[1]['isoform']))
        all_total=sum(totals[i] for i,_ in group)
        selected=group[:c['simulation']['max_chain_states']]
        if not selected:continue
        idx=[i for i,_ in selected];chainrecords=[r for _,r in selected]
        eventids=events_by_gene[gene]
        maps={str(e):[chain_event_branch(r,data['events'][e]) for r in chainrecords] for e in eventids}
        C=counts[:,idx].astype(float)
        result[gene]={'chain_ids':[r['isoform'] for r in chainrecords],
            'counts':C,'event_branch_maps':maps,
            'retained_chain_molecule_fraction':float(totals[idx].sum()/all_total) if all_total else 0.,
            'structure_scope':'observed_splice_chain; not guaranteed full_length',
            'eligibility_basis':'representative_corrected_exon_blocks; molecule_specific_completeness_unresolved',
            'observed_chain_count':len(group),'retained_chain_count':len(selected)}
    return result


def reconstruct(prior,maps,targets,lam):
    prior=np.asarray(prior,float)
    if np.any(prior<=0) or not np.isclose(prior.sum(),1):raise ValueError('Strictly positive reference composition required')
    if not targets:return prior.copy(),{'status':'no_event_constraint','rank':0,'conflict':None}
    constraints=[];J=[]
    for e,target in targets.items():
        assignment=np.asarray(maps[e],int);target=np.asarray(target,float)
        H=np.array([assignment==k for k in range(len(target))],float);D=H.sum(0)
        den=D@prior
        if den<=1e-12:continue
        pred=H@prior/den
        if np.any(H.sum(1)==0):
            # Legal observed structures cannot realize all branch targets; report residual later.
            pass
        constraints.append((H,D,target))
        J.extend((H[:-1]*den-(H[:-1]@prior)[:,None]*D)/(den**2))
    if not constraints:return prior.copy(),{'status':'structure_unresolved','rank':0,'conflict':None}
    def objective(z):
        log_pi=log_softmax(z);pi=np.exp(log_pi)
        log_ratio=log_pi-np.log(prior)
        value=float(np.sum(pi*log_ratio));grad=log_ratio+1
        for H,D,target in constraints:
            num=np.maximum(H@pi,1e-12);den=max(D@pi,1e-12)
            value+=lam/len(constraints)*float(-np.sum(target*np.log(num/den)))
            grad+=lam/len(constraints)*(-H.T@(target/num)+D/den)
        return value,pi*(grad-grad@pi)
    fit=minimize(objective,np.log(prior),jac=True,method='L-BFGS-B',options={'maxiter':300,'ftol':1e-10})
    if not fit.success:return prior.copy(),{'status':'structure_optimizer_not_converged','rank':0,'conflict':None}
    pi = softmax(fit.x)
    residuals = []
    for H, D, target in constraints:
        den = max(float(D @ pi), 1e-12)
        residuals.append(np.sum(np.abs((H @ pi) / den - target)))
    conflict = float(np.mean(residuals))
    rank=int(np.linalg.matrix_rank(np.asarray(J),tol=1e-7)) if len(J) else 0
    return pi,{'status':'ok' if rank>=len(prior)-1 else 'structure_unresolved','rank':rank,'conflict':conflict}


def structure_predictions(structures,event_wt,event_new,c):
    result={};stats=[]
    for gene,structure in structures.items():
        C=structure['counts'].toarray() if sparse.issparse(structure['counts']) else np.asarray(structure['counts']);total=C.sum(1);k=C.shape[1]
        pooled=(C.sum(0)+.5)/(C.sum()+k*.5);s=c['simulation']['structure_prior_strength']
        ref=(C+s*pooled)/(total[:,None]+s);out=ref.copy();valid=total>0
        for i in range(len(C)):
            if not valid[i]:stats.append({'gene':gene,'cell_index':i,'status':'structure_unobserved','conflict':None});continue
            targets={e:event_new[e][i] for e in structure['event_branch_maps'] if e in event_new}
            wt={e:event_wt[e][i] for e in targets}
            if not targets or all(np.array_equal(targets[e],wt[e]) for e in targets):
                stats.append({'gene':gene,'cell_index':i,'status':'no_change','conflict':0.});continue
            # Reconstruct feasible WT and KO optima, then apply the relative composition
            # change to the observed reference. This guarantees q=0 is the identity even
            # if event marginals and truncated chain estimates disagree at WT.
            base,diag0=reconstruct(ref[i],structure['event_branch_maps'],wt,c['simulation']['structure_lambda'])
            pred,diag=reconstruct(ref[i],structure['event_branch_maps'],targets,c['simulation']['structure_lambda'])
            if diag['status']=='structure_optimizer_not_converged' or diag0['status']=='structure_optimizer_not_converged':
                stats.append({'gene':gene,'cell_index':i,'status':'structure_optimizer_not_converged','conflict':None});continue
            out[i]=softmax(np.log(ref[i])+np.log(np.maximum(pred,1e-300))-np.log(np.maximum(base,1e-300)))
            stats.append({'gene':gene,'cell_index':i,**diag})
        result[gene]={'reference':ref,'perturbed':out,'valid':valid,'chain_ids':structure['chain_ids'],
                      'scope':structure['structure_scope'],'retained_fraction':structure['retained_chain_molecule_fraction']}
    return result,stats


def _gene_prediction(job):
    gene,structure,wt,new,c=job
    result,audit=structure_predictions({gene:structure},wt,new,c)
    return gene,result[gene],audit


def _structure_worker_init():
    import os
    for name in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMBA_NUM_THREADS']:os.environ[name]='1'
    from threadpoolctl import threadpool_limits
    threadpool_limits(1)


def iter_structure_predictions(jobs,workers=1):
    if workers<=1:
        for job in jobs:yield _gene_prediction(job)
        return
    from concurrent.futures import ProcessPoolExecutor,wait,FIRST_COMPLETED
    import multiprocessing as mp
    with ProcessPoolExecutor(max_workers=workers,mp_context=mp.get_context('fork'),initializer=_structure_worker_init) as pool:
        iterator=iter(jobs);pending=set()
        for job in iterator:
            pending.add(pool.submit(_gene_prediction,job))
            if len(pending)>=2*workers:break
        while pending:
            done,_=wait(pending,return_when=FIRST_COMPLETED)
            for f in done:
                pending.remove(f);yield f.result()
                job=next(iterator,None)
                if job is not None:pending.add(pool.submit(_gene_prediction,job))
