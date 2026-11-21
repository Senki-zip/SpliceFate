"""Official CytoTRACE 0.3.1 scores on complete gene counts for DPT roots."""
from __future__ import annotations
import logging
import os
from pathlib import Path
import shutil
import subprocess
import tempfile

import numpy as np
import pandas as pd
from scipy import sparse
from scipy.io import mmwrite

from .io import digest

LOG=logging.getLogger('splicefate')
OFFICIAL_SOURCE_SHA256='7dafcfa42c792eaffe1048f4716cb54105b31f2ad01c4f0c26868c10593b155d'


def source_files(c):
    cfg=c['trajectory'].get('cytotrace',{})
    root=Path(cfg.get('source_dir',Path(__file__).resolve().parents[2]/'tools/CytoTRACE_official'))
    return [root/'R/CytoTRACE.R',root/'DESCRIPTION',root/'COPY_PROVENANCE.json']


def compute_scores(counts,cells,genes,cfg,output=None):
    cells=list(map(str,cells));genes=list(map(str,genes));matrix=sparse.csr_matrix(counts)
    if matrix.shape!=(len(cells),len(genes)) or len(set(cells))!=len(cells) or len(set(genes))!=len(genes):
        raise ValueError('CytoTRACE needs unique aligned full gene/cell axes')
    if any('\n' in v or '\r' in v for v in cells+genes):raise ValueError('Invalid cell/gene ID')
    if np.any(~np.isfinite(matrix.data)) or np.any(matrix.data<0) or np.any(matrix.data!=np.floor(matrix.data)):
        raise ValueError('CytoTRACE requires complete raw nonnegative integer gene counts')
    sums=np.asarray(matrix.sum(1)).ravel()
    if np.any(sums<=0):raise ValueError('CytoTRACE input contains empty cell libraries')
    root=Path(cfg['source_dir']);source=root/'R/CytoTRACE.R'
    if digest(source)!=OFFICIAL_SOURCE_SHA256 or 'Version: 0.3.1' not in (root/'DESCRIPTION').read_text():
        raise ValueError('Expected unmodified official CytoTRACE 0.3.1 source')
    executable=shutil.which(str(cfg.get('rscript_binary','Rscript')))
    if executable is None:raise RuntimeError('CytoTRACE requires a working trajectory.cytotrace.rscript_binary')
    prefix=Path(executable).absolute().parent.parent
    runner=Path(__file__).with_name('cytotrace_runner.R')
    env=os.environ.copy()
    for name in ['R_HOME','R_ARCH','R_LIBS']:env.pop(name,None)
    library=prefix/'lib/R/library'
    env.update(R_LIBS_USER=str(library),R_LIBS_SITE=str(library),LD_LIBRARY_PATH=str(prefix/'lib'),
               OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1')
    params={'seed':int(cfg.get('seed',42)), 'enable_fast':bool(cfg.get('enable_fast',True)),
            'ncores':int(cfg.get('ncores',1)),'subsample_size':int(cfg.get('subsample_size',1000))}
    with tempfile.TemporaryDirectory(prefix='splicefate_cytotrace_') as temp:
        folder=Path(temp)
        mmwrite(folder/'counts.mtx',matrix.T.astype(np.int64),field='integer')
        (folder/'cells.txt').write_text('\n'.join(cells)+'\n');(folder/'genes.txt').write_text('\n'.join(genes)+'\n')
        np.savetxt(folder/'cell_sums.txt',sums,fmt='%.17g')
        args=[executable,'--vanilla',str(runner),temp,str(source),str(library),str(params['seed']),
              'TRUE' if params['enable_fast'] else 'FALSE',str(params['ncores']),str(params['subsample_size'])]
        result=subprocess.run(args,env=env,capture_output=True,text=True)
        if output:
            output=Path(output);output.mkdir(parents=True,exist_ok=True)
            (output/'official.log').write_text(result.stdout+result.stderr)
        if result.returncode:
            raise RuntimeError('Official CytoTRACE failed:\n'+result.stdout[-1500:]+result.stderr[-4500:])
        scores=pd.read_csv(folder/'scores.tsv',sep='\t',dtype={'cell':str})
        if scores.cell.tolist()!=cells:raise ValueError('CytoTRACE changed cell axis')
        if not scores.cytotrace_qc.eq(scores.cytotrace_score.between(0,1)&np.isfinite(scores.cytotrace_score)).all():
            raise ValueError('CytoTRACE score/QC mismatch')
        runtime=(folder/'runtime.txt').read_text()
        if output:
            shutil.copy2(folder/'scores.tsv',output/'scores.tsv')
            (output/'runtime.txt').write_text(runtime)
    return scores,{'method':'official_CytoTRACE_0.3.1','source_sha256':digest(source),
                   'runner_sha256':digest(runner),'rscript':executable,'r_runtime':runtime,
                   'full_gene_count':len(genes),'cells':len(cells),'scored_cells':int(scores.cytotrace_qc.sum()),
                   'parameters':params,'batch_correction':False,'scope':'DPT_context_full_gene_counts'}


def score_contexts(data,frame,graph,c):
    from .trajectory import trajectory_contexts
    scores=np.full(len(frame),np.nan);counts=np.full(len(frame),np.nan);qc=np.zeros(len(frame),bool)
    audit=[];root=Path(c['output'])/'cytotrace'
    for number,(context,ix) in enumerate(trajectory_contexts(frame,graph,c)):
        if len(ix)<c['trajectory']['min_component_cells']:
            audit.append({'context':context,'cells':len(ix),'status':'insufficient_context_cells'});continue
        LOG.info('Official CytoTRACE %s: %d cells, %d full genes',context,len(ix),len(data['genes']))
        result,meta=compute_scores(data['G'][ix],frame.cell.iloc[ix],data['genes'].gene,
                                   c['trajectory']['cytotrace'],root/f'context_{number:03d}')
        scores[ix]=result.cytotrace_score;qc[ix]=result.cytotrace_qc;counts[ix]=result.cytotrace_detected_genes
        audit.append({'context':context,**meta,'status':'ok' if qc[ix].any() else 'no_finite_CytoTRACE_score'})
    frame['cytotrace_score']=scores;frame['cytotrace_qc']=qc;frame['cytotrace_detected_genes']=counts
    root.mkdir(parents=True,exist_ok=True)
    frame[['cell','cytotrace_score','cytotrace_qc','cytotrace_detected_genes']].to_csv(root/'scores.tsv',sep='\t',index=False)
    return {'method':'official_CytoTRACE_0.3.1','contexts':audit,
            'selection':'highest_finite_score_in_earliest_stage_per_DPT_context; ties_by_cell_ID',
            'all_cells_retained':True,'score_is_not_pseudotime':True}
