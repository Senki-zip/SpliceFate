"""Configuration-driven launcher; set native thread limits before model imports."""
from __future__ import annotations
import argparse
import copy
import json
import logging
import os
from pathlib import Path
import shutil
import subprocess
import sys
import time


def launch_settings(path):
    import yaml
    path=Path(path).expanduser().resolve()
    raw=yaml.safe_load(path.read_text())
    root=(path.parent/raw.get('project_root','..')).resolve()
    runtime=raw.get('runtime',{})
    threads=runtime.get('threads',4)
    if isinstance(threads,bool) or not isinstance(threads,int) or threads<1:
        raise ValueError('runtime.threads must be a positive integer')
    binary=runtime.get('python_binary',sys.executable)
    if '/' in binary:binary=str((root/Path(binary).expanduser()).resolve())
    executable=shutil.which(binary)
    if not executable:raise ValueError('Python interpreter not found: '+binary)
    env=os.environ.copy()
    for name in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS',
                 'NUMEXPR_NUM_THREADS','NUMBA_NUM_THREADS','VECLIB_MAXIMUM_THREADS',
                 'BLIS_NUM_THREADS']:
        env[name]=str(threads)
    env['MPLBACKEND']='Agg'
    # Permit execution from any directory using this checkout's source.
    source=str(Path(__file__).resolve().parents[1])
    env['PYTHONPATH']=source+(os.pathsep+env['PYTHONPATH'] if env.get('PYTHONPATH') else '')
    return path,executable,env


def execution_plan(c):
    execution=c.get('execution',{})
    stages=execution.get('stages',['prepare','fit','simulate'])
    order=['prepare','fit','simulate']
    if not isinstance(stages,list) or not stages or any(s not in order for s in stages):
        raise ValueError('execution.stages must contain prepare, fit and/or simulate')
    if stages!=[s for s in order if s in stages]:
        raise ValueError('execution.stages must be unique and ordered: prepare, fit, simulate')
    limit=execution.get('event_limit')
    if limit is not None and (isinstance(limit,bool) or not isinstance(limit,int) or limit<1):
        raise ValueError('execution.event_limit must be null or a positive integer')
    bootstrap=execution.get('bootstrap',0)
    if isinstance(bootstrap,bool) or not isinstance(bootstrap,int) or bootstrap<0:
        raise ValueError('execution.bootstrap must be a nonnegative integer')
    jobs=[]
    if 'simulate' in stages:
        targets=c['simulation'].get('target_sets',[])
        strengths=c['simulation'].get('strengths',[])
        if not isinstance(targets,list) or not targets or not isinstance(strengths,list) or not strengths:
            raise ValueError('simulation.target_sets and strengths must be nonempty lists')
        known=set();normalized=[]
        for group in targets:
            if not isinstance(group,list) or not group or any(not isinstance(t,str) or not t.strip() or t!=t.strip() or any(k in t for k in ['/', '\\', '+', ';']) for t in group):
                raise ValueError('Each simulation.target_sets entry must be a list of nonempty target identifiers')
            if len(set(group))!=len(group):raise ValueError('Duplicate target in a target set')
            key=tuple(sorted(group))
            if key in known:raise ValueError('Duplicate target set')
            known.add(key);normalized.append(list(key))
        doses=[]
        for q in strengths:
            if isinstance(q,bool) or not isinstance(q,(int,float)) or not 0<=q<=1:
                raise ValueError('Each simulation strength must be a number in [0,1]')
            if q in doses:raise ValueError('Duplicate simulation strength')
            doses.append(float(q))
        jobs=[{'targets':group,'strength':q} for group in normalized for q in doses]
    return {'stages':stages,'event_limit':limit,'full_event_fit':limit is None,
            'bootstrap':bootstrap,'simulations':jobs,'output':c['output']}


def model_configuration(c):
    """Execution choices do not change the frozen model's identity."""
    result=copy.deepcopy(c)
    result.pop('execution',None);result.pop('runtime',None)
    result['simulation'].pop('target_sets',None)
    result['simulation'].pop('strengths',None)
    return result


def preflight(c,plan):
    from .cli import input_files
    missing=[str(p) for p in input_files(c) if not p.is_file()]
    if missing:raise FileNotFoundError('Missing input files:\n'+'\n'.join(missing))
    root=Path(c['output'])
    if 'prepare' in plan['stages'] and (root/'provenance.json').exists():
        raise FileExistsError('Prepared output already exists. Use a new output directory, or execution.stages: [simulate] / [fit, simulate]')
    if 'prepare' not in plan['stages'] and not (root/'provenance.json').is_file():
        raise FileNotFoundError('No prepared model at '+str(root))
    if 'simulate' in plan['stages'] and 'fit' not in plan['stages'] and not (root/'fit_summary.json').is_file():
        raise FileNotFoundError('No fitted model at '+str(root))
    if 'fit' not in plan['stages'] and 'prepare' not in plan['stages']:
        from .cli import check_prepared
        check_prepared(c)


def run_pipeline(c,plan,launch_config=None,runtime_info=None):
    import fcntl
    from .io import dump
    from .cli import prepare,fit,simulate_run
    root=Path(c['output']);root.mkdir(parents=True,exist_ok=True)
    with (root/'.pipeline.lock').open('w') as guard:
        fcntl.flock(guard,fcntl.LOCK_EX|fcntl.LOCK_NB)
        preflight(c,plan)
        if launch_config is not None:
            import yaml
            (root/'launch_config.resolved.yaml').write_text(yaml.safe_dump(launch_config,sort_keys=False))
        if runtime_info is not None:dump(root/'runtime.json',runtime_info)
        log=logging.getLogger('splicefate')
        handler=logging.FileHandler(root/'pipeline.log')
        handler.setFormatter(logging.Formatter('%(asctime)s %(levelname)s %(message)s'))
        log.addHandler(handler);log.setLevel(logging.INFO)
        state={'status':'running','started_at':time.strftime('%Y-%m-%dT%H:%M:%S%z'),
               'pid':os.getpid(),'python_binary':sys.executable,'plan':plan,
               'completed_stages':[],'completed_simulations':[]}
        dump(root/'run_status.json',state)
        dump(root/'execution_plan.json',plan)
        try:
            for stage in plan['stages']:
                state['current_stage']=stage;dump(root/'run_status.json',state)
                log.info('Starting stage: %s',stage)
                if stage=='prepare':prepare(c)
                elif stage=='fit':fit(c,plan['event_limit'])
                else:
                    for job in plan['simulations']:
                        state['current_simulation']=job;dump(root/'run_status.json',state)
                        simulate_run(c,job['targets'],job['strength'],bootstrap=plan['bootstrap'])
                        state['completed_simulations'].append(job)
                        dump(root/'run_status.json',state)
                state['completed_stages'].append(stage)
            state.update(status='complete',current_stage=None,current_simulation=None,
                         finished_at=time.strftime('%Y-%m-%dT%H:%M:%S%z'))
        except BaseException as exc:
            state.update(status='failed',error=str(exc),finished_at=time.strftime('%Y-%m-%dT%H:%M:%S%z'))
            log.exception('Pipeline failed');raise
        finally:
            dump(root/'run_status.json',state)
            log.removeHandler(handler);handler.close()
        return state


def main():
    parser=argparse.ArgumentParser(description='Run SpliceFate using paths, threads and KO scenarios from YAML')
    parser.add_argument('--config',default=str(Path(__file__).resolve().parents[2]/'configs/gse295932_nontarget_hg38.yaml'))
    parser.add_argument('--check',action='store_true',help='Validate configuration and inputs, print plan; do not run or write outputs')
    parser.add_argument('--worker',action='store_true',help=argparse.SUPPRESS)
    args=parser.parse_args()
    path,executable,env=launch_settings(args.config)
    if not args.worker:
        command=[executable,'-m','splicefate.runner','--worker','--config',str(path)]
        if args.check:command.append('--check')
        return subprocess.call(command,env=env)
    # All numerical imports occur after the child received configured thread limits.
    from .io import configuration
    from numba import get_num_threads
    from threadpoolctl import threadpool_info
    import scanpy as sc
    c=configuration(path);plan=execution_plan(c)
    threads=int(c.get('runtime',{}).get('threads',4));sc.settings.n_jobs=threads
    model=model_configuration(c);preflight(model,plan)
    info={'configured_threads':threads,'numba_threads':get_num_threads(),
          'native_threadpools':threadpool_info(),'python_binary':sys.executable}
    print(json.dumps({'config':str(path),'runtime':info,'plan':plan},indent=2),flush=True)
    if args.check:return 0
    logging.basicConfig(level=logging.INFO,format='%(asctime)s %(levelname)s %(message)s')
    run_pipeline(model,plan,launch_config=c,runtime_info=info)
    return 0


if __name__=='__main__':
    raise SystemExit(main())
