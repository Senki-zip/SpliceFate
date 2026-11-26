from __future__ import annotations
import hashlib
import importlib.metadata
import json
from pathlib import Path
import platform
import numpy as np
import pandas as pd
from scipy import sparse
import yaml


def dump(path, value):
    path = Path(path); path.parent.mkdir(parents=True, exist_ok=True)
    def convert(v):
        if isinstance(v, np.ndarray): return v.tolist()
        if isinstance(v, np.generic): return v.item()
        if isinstance(v, Path): return str(v)
        raise TypeError(type(v).__name__)
    tmp = path.with_suffix(path.suffix + '.tmp')
    tmp.write_text(json.dumps(value, indent=2, ensure_ascii=False, default=convert, allow_nan=False))
    tmp.replace(path)


def digest(path):
    with Path(path).open('rb') as f: return hashlib.file_digest(f, 'sha256').hexdigest()


def configuration(path):
    path = Path(path).resolve()
    config = yaml.safe_load(path.read_text())
    root = (path.parent / config.pop('project_root', '..')).resolve()
    for key, val in config['input'].items():
        if val is not None: config['input'][key] = str((root / val).resolve())
    if config['trajectory'].get('root_method')=='cytotrace':
        cfg=config['trajectory'].setdefault('cytotrace',{})
        cfg['source_dir']=str((root / cfg.get('source_dir','tools/CytoTRACE_official')).resolve())
        if '/' in cfg.get('rscript_binary','Rscript'):
            cfg['rscript_binary']=str((root / cfg['rscript_binary']).resolve())
    config['output'] = str((root / config['output']).resolve())
    config['_config_path'] = str(path)
    validate_config(config)
    return config


def validate_config(c):
    if not isinstance(c['clustering'].get('enabled',True),bool):
        raise ValueError('clustering.enabled must be true or false')
    if c['input'].get('functional_feedback'):
        raise ValueError('Manual functional_feedback is retired; isoform composition responses are learned from WT counts')
    if c['trajectory']['method']!='dpt' or c['clustering']['fallback_pseudotime_method']!='dpt':
        raise ValueError('WT ordering uses DPT, including the shared/global small-cluster context')
    if c['trajectory'].get('root_method','stage') not in ('stage','cytotrace'):
        raise ValueError('trajectory.root_method must be stage or cytotrace')
    if c['trajectory'].get('root_method')=='cytotrace':
        cfg=c['trajectory']['cytotrace']
        if cfg.get('ncores',1)<1 or cfg.get('subsample_size',1000)<2:
            raise ValueError('Invalid official CytoTRACE worker/subsample settings')
    if 'assembly' in c['input']:
        raise ValueError('assembly belongs at the top level, not in input paths')
    if c['assembly'] not in ('mm10','hg38'): raise ValueError('Event adapter requires native mm10 or hg38; no implicit liftover')
    if c['model'].get('count_unit','integer_umi') not in ('integer_umi','author_fractional_umi_estimates'):
        raise ValueError('Unsupported model.count_unit')
    for section, key in [('trajectory','splice_weight'), ('priors','motif_pvalue')]:
        v=c[section][key]
        if not 0 <= v <= 1: raise ValueError(f'{section}.{key} must be within [0,1]')
    if c['model']['response'] not in ('linear','spline'): raise ValueError('response must be linear or spline')
    for key in ['ridge','interaction_ridge','state_ridge','dm_concentration']:
        if c['model'][key] <= 0: raise ValueError(f'{key} must be positive')
    if (isinstance(c['model']['max_factors'],bool) or not isinstance(c['model']['max_factors'],int)
            or c['model']['max_factors'] < 1 or c['model']['max_interactions'] < 0):
        raise ValueError('Invalid model budget')
    from .binning import DEFAULTS as bin_defaults
    from .model import GROUP_DEFAULTS
    bins=c.setdefault('binning',{})
    model=c['model']
    gates=model.setdefault('group_gates',{})
    if not isinstance(bins,dict) or not isinstance(gates,dict):
        raise ValueError('binning and model.group_gates must be mappings')
    for key,value in bin_defaults.items():
        bins.setdefault(key,value)
    for key,value in GROUP_DEFAULTS.items():
        gates.setdefault(key,value)
    model.setdefault('state_dependent',False)
    if not isinstance(bins['enabled'],bool) or bins['mode']!='fixed_gaussian':
        raise ValueError('binning requires a boolean enabled and mode fixed_gaussian')
    if not isinstance(model['state_dependent'],bool):
        raise ValueError('model.state_dependent must be boolean')
    if model.setdefault('selection_mode','all_passing_singles_hybrid') not in ('all_passing_singles_hybrid','top_k'):
        raise ValueError('Invalid model.selection_mode')
    limit=model.setdefault('max_single_factors',None)
    if limit is not None and (type(limit) is not int or limit<1):
        raise ValueError('model.max_single_factors must be null or a positive integer')
    model.setdefault('spline_knots',0)
    model.setdefault('output_cap_lambda',.25)
    for section,values,minimums in [
        ('binning',bins,dict(n_bins=1,min_cells_per_bin=1,min_support_cells=1,min_oob_bins=1)),
        ('model',model,dict(spline_knots=0)),
        ('model.group_gates',gates,dict(min_group_rows=1)),
    ]:
        for key,minimum in minimums.items():
            if type(values[key]) is not int or values[key]<minimum:
                raise ValueError(f'{section}.{key} must be an integer >= {minimum}')
    for section,values,bounds in [
        ('binning',bins,{'overlap_factor':(0.,np.inf)}),
        ('model',model,{'output_cap_lambda':(0.,np.inf)}),
        ('model.group_gates',gates,{'min_group_cell_ess':(0.,np.inf),'min_group_effective_rows':(0.,np.inf),
            'max_geometry_condition':(1.,np.inf),'min_group_oob_gain':(0.,np.inf),
            'min_group_positive_gain_fraction':(0.,1.),'min_group_penalized_gain':(0.,1.)}),
    ]:
        for key,(low,high) in bounds.items():
            value=values[key]
            if (type(value) not in (int,float) or not np.isfinite(value) or not low<=value<=high
                    or (key=='overlap_factor' and value==0) or (key=='min_group_penalized_gain' and value==1)):
                raise ValueError(f'Invalid {section}.{key}')
    from .network import selection_settings,POLICY
    if not isinstance(model.get('network_selection',{}),dict):
        raise ValueError('model.network_selection must be a mapping')
    selection=selection_settings(c);c['model']['network_selection']=selection
    if selection['method']!=POLICY:
        raise ValueError('model.network_selection.method must be '+POLICY)
    if selection['screening_mode'] not in ('all_wt_penalized','cell_oob'):
        raise ValueError('Invalid model.network_selection.screening_mode')
    all_wt=selection['screening_mode']=='all_wt_penalized'
    for key in ['min_response_score','min_oob_fit_gain','min_penalized_fit_gain']:
        minimum=selection[key]
        if isinstance(minimum,bool) or not isinstance(minimum,(int,float)) or not np.isfinite(minimum) or minimum<0:
            raise ValueError('model.network_selection.'+key+' must be finite and nonnegative')
    for key,minimum in [('bootstrap_repeats',0 if all_wt else 2),('bootstrap_seed',0)]:
        value=selection[key]
        if isinstance(value,bool) or not isinstance(value,int) or value<minimum:
            raise ValueError('Invalid model.network_selection.'+key)
    weights=selection['weights']
    if not isinstance(weights,dict) or set(weights)!={'fit_gain','stability','response'}:
        raise ValueError('network selection weights require fit_gain, stability and response')
    if any(isinstance(w,bool) or not isinstance(w,(int,float)) or not np.isfinite(w) or w<0 for w in weights.values()) or sum(weights.values())<=0:
        raise ValueError('network selection weights must be finite, nonnegative with positive sum')
    if all_wt and weights['fit_gain']+weights['response']<=0:
        raise ValueError('all_wt_penalized network selection weights require positive fit_gain + response')
    for key in ['min_positive_gain_fraction','min_stability','min_partial_r2']:
        v=selection[key]
        if isinstance(v,bool) or not isinstance(v,(int,float)) or not np.isfinite(v) or not 0<=v<=1:
            raise ValueError('model.network_selection.'+key+' must be within [0,1]')
    for section,key,minimum,default in [('priors','partner_hops',0,1),('model','max_event_branches',2,8),('model','support_neighbors',1,10)]:
        value=c[section].setdefault(key,default)
        if isinstance(value,bool) or not isinstance(value,int) or value<minimum:
            raise ValueError(section+'.'+key+' must be an integer >= '+str(minimum))
    extras=c['priors'].setdefault('extra_targets',[])
    if not isinstance(extras,list) or any(not isinstance(t,str) or not t.strip() for t in extras) or len(set(extras))!=len(extras):
        raise ValueError('priors.extra_targets must contain unique identifiers')
    if not isinstance(c['simulation'].setdefault('write_unchanged_cells',False),bool):
        raise ValueError('simulation.write_unchanged_cells must be boolean')
    threshold=c['simulation'].setdefault('min_delta_threshold',1e-4)
    if isinstance(threshold,bool) or not isinstance(threshold,(int,float)) or not np.isfinite(threshold) or threshold<0:
        raise ValueError('simulation.min_delta_threshold must be finite and nonnegative')
    if 'gene_limit' in c['simulation']:raise ValueError('simulation.gene_limit is retired; remove it')
    for key in ['divergence_limit','convergence_tolerance']:
        value=c['simulation'][key]
        if isinstance(value,bool) or not isinstance(value,(int,float)) or not np.isfinite(value) or value<=0:
            raise ValueError('simulation.'+key+' must be finite and positive')
    if not 0<c['model']['support_distance_quantile']<=1:raise ValueError('Invalid support_distance_quantile')
    from .composition import settings
    iso=settings(c)
    for key in ['max_events_per_factor','min_branch_molecules','min_branch_cells','min_joint_cells','ridge','count_shrinkage']:
        if iso[key]<=0:raise ValueError('Positive model.isoform_state setting required: '+key)
    if type(iso['min_joint_bins']) is not int or iso['min_joint_bins']<1:
        raise ValueError('model.isoform_state.min_joint_bins must be a positive integer')
    if iso['max_branches']<2:raise ValueError('Isoform groups require at least two branches')
    if c['simulation']['max_rounds']<1 or c['simulation']['max_chain_states']<2:raise ValueError('Invalid cascade/structure budget')
    if c['simulation']['structure_lambda']<=0 or c['simulation']['structure_prior_strength']<=0:raise ValueError('Positive structure penalties required')
    if not 0<=c['trajectory'].get('velocity_weight',0)<=1:raise ValueError('velocity_weight must be within [0,1]')
    if c['trajectory'].get('S_calibration','composition') not in ('composition','intron_phase','hybrid'):
        raise ValueError('S_calibration must be composition, intron_phase or hybrid')
    if c['trajectory'].get('S_calibration') in ('intron_phase','hybrid'):
        ic=c['trajectory']['intron_phase']
        if not 0<=ic['weight']<=1 or not 0<=ic['reference_quantile']<1:raise ValueError('Invalid intron calibration weights/quantile')
        for key in ['min_gene_cells','min_region_counts','count_shrinkage','min_shared_genes','reliability_genes']:
            if ic[key]<=0:raise ValueError('Positive intron support settings required: '+key)


def provenance(c, inputs):
    packages={}
    for name in ['numpy','scipy','pandas','scanpy','scikit-learn','palantir','duckdb','numba']:
        try: packages[name]=importlib.metadata.version(name)
        except importlib.metadata.PackageNotFoundError: packages[name]='missing'
    source=Path(__file__).parent
    files={str(p):digest(p) for p in sorted(source.glob('*')) if p.is_file() and p.suffix in ('.py','.R')}
    return {'schema':'splicefate-v0.3','training_regime':'strict_WT_only',
        'config':c,'inputs':{str(p):digest(p) for p in inputs},'source_sha256':files,
        'python':platform.python_version(),'packages':packages,
        'cell_split':'all_WT_cells; no heldout cells','effect_interpretation':'conditional_expression_proxy_scenario',
        'S_completion':False,'IR_training':False}


def aligned_cells(table, reference):
    if table['cell'].duplicated().any() or len(set(reference)) != len(reference):
        raise ValueError('Duplicate cell IDs')
    if set(table.cell) != set(reference): raise ValueError('Cell sets differ')
    return table.set_index('cell').loc[reference].reset_index()


def load_data(c):
    gp=Path(c['input']['genes']);ep=Path(c['input']['events'])
    cells=pd.read_csv(gp/'cells.tsv',sep='\t',dtype=str)
    if 'cell' not in cells:raise ValueError('cells.tsv requires a cell column')
    for name in ['stage','batch']:
        if name not in cells:cells[name]='unknown'
    if 'stage_order' not in cells:
        cells['stage_order']='0'
        if not c['trajectory'].get('root_cell'):cells['trajectory_direction_available']='false'
    if 'trajectory_direction_available' not in cells:cells['trajectory_direction_available']='true'
    genes=pd.read_csv(gp/'genes.tsv',sep='\t',dtype=str)
    if genes.gene.duplicated().any(): raise ValueError('Duplicate stable gene IDs')
    ec=pd.read_csv(ep/'cells.tsv',sep='\t',dtype=str)
    aligned_cells(ec,list(cells.cell))
    order=ec.set_index('cell').loc[cells.cell].index
    rowmap=pd.Index(ec.cell).get_indexer(order)
    G=sparse.load_npz(gp/'gene_counts.npz').tocsr()
    # Event operations slice columns; CSC avoids rescanning every cell library.
    C=sparse.load_npz(ep/'branch_counts.npz').tocsr()[rowmap].tocsc()
    offsets=np.load(ep/'branch_offsets.npy')
    events=json.loads((ep/'groups.json').read_text())
    if G.shape != (len(cells),len(genes)) or C.shape[0]!=len(cells): raise ValueError('Count axes differ')
    if len(offsets)!=len(events)+1 or offsets[0]!=0 or offsets[-1]!=C.shape[1] or np.any(np.diff(offsets)<2):
        raise ValueError('Invalid event branch offsets')
    fractional=c['model'].get('count_unit','integer_umi')=='author_fractional_umi_estimates'
    if fractional:
        meta=json.loads((gp/'quantification.json').read_text())
        if meta.get('S_count_unit')!='author_fractional_umi_estimates' or meta.get('assembly')!=c['assembly'] or meta.get('validated') is not True:
            raise ValueError('Fractional source counts require validated matching quantification provenance')
    for label,X in [('G',G),('S',C)]:
        if np.any(X.data<0) or not np.all(np.isfinite(X.data)):
            raise ValueError('Finite nonnegative source counts required')
        if (label=='G' or not fractional) and not np.all(X.data==np.floor(X.data)):
            raise ValueError('Raw integer molecule counts required')
    library=np.asarray(G.sum(1)).ravel()
    if np.any(library<=0): raise ValueError('Zero complete library')
    scales=c['normalization']['target_sum']/library
    norms=ep/'normalization_factors.tsv'
    if norms.exists():
        old=aligned_cells(pd.read_csv(norms,sep='\t'),list(cells.cell))
        if not np.array_equal(library,old['complete_gene_UMI_library'].values):
            raise ValueError('Complete library differs from established normalization')
    hvg=pd.read_csv(ep/'genes.tsv',sep='\t',dtype=str).gene.tolist()
    if len(hvg)!=c['normalization']['hvg_count'] or not set(hvg)<=set(genes.gene):
        raise ValueError('Existing Seurat-v3 HVG axis mismatch; do not silently reselect')
    data={'cells':cells,'genes':genes,'G':G,'C':C,'offsets':offsets,'events':events,'scales':scales,'hvg':hvg}
    data['S_count_unit']='author_fractional_umi_estimates' if fractional else 'integer_umi'
    return data


def normalized(G,scales):
    x=G.multiply(scales[:,None]).tocsr().astype(float)
    x.data=np.log1p(x.data)
    return x

