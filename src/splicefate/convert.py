"""Convert author molecule AllInfo + paired 10x MEX into model input schema.

No realignment, author significant-event selection, S imputation, or precursor
inference. Exact conflicting UMI structures are excluded and audited.
"""
from __future__ import annotations
import argparse
import collections
import csv
import fcntl
import gzip
import hashlib
import io
import json
import os
from pathlib import Path
import re
import subprocess
import tarfile
import time

import duckdb
import numpy as np
import pandas as pd
from scipy import sparse
from scipy.io import mmread

from .io import dump,digest
from .structure import chain_event_branch

SAMPLES={'P21_M1_HIPP':('P21_SC_M1_HIPP','M1','Hippocampus'),
         'P21_M1_VIS':('P21_SC_M1_VIS','M1','VisCortex'),
         'P21_M5_HIPP':('P21_SC_M5_HIPP','M5','Hippocampus'),
         'P21_M5_VIS':('P21_SC_M5_VIS','M5','VisCortex')}


def log(stage,**kwargs):print(json.dumps({'time':time.strftime('%Y-%m-%dT%H:%M:%S%z'),'stage':stage,**kwargs}),flush=True)


def stable_id(value):return str(value).split('.')[0]


def parse_regions(text):
    regions=[]
    for token in text.split(';%;'):
        if not token:continue
        chrom,left,right,strand=token.rsplit('_',3)
        start,end=int(left)-1,int(right)
        if strand not in ('+','-') or start<0 or end<=start:raise ValueError('Invalid AllInfo inclusive region')
        regions.append((chrom,start,end,strand))
    return regions


def parse_structure(introns,exons):
    """AllInfo 1-based closed -> model 0-based half-open, validate adjacency."""
    I=sorted(parse_regions(introns),key=lambda x:x[1]);E=sorted(parse_regions(exons),key=lambda x:x[1])
    if not E or len({(r[0],r[3]) for r in I+E})!=1:raise ValueError('Mixed chromosome/strand or no exon')
    chrom,_,_,strand=E[0]
    if any(a[2]>b[1] for a,b in zip(E,E[1:])):raise ValueError('Overlapping exon blocks')
    expected=[(chrom,a[2],b[1],strand) for a,b in zip(E,E[1:])]
    if I!=expected:raise ValueError('Intron chain and exon-block gaps disagree')
    return chrom,strand,[(i[1],i[2]) for i in I],[(e[1],e[2]) for e in E]


def chain_record(gene,chrom,strand,junctions,exons):
    ordered=junctions if strand=='+' else junctions[::-1]
    tokens=[f'junction:{chrom}:{a}-{b}:{strand}' for a,b in ordered]
    key=gene+'|'+chrom+'|'+strand+'|'+json.dumps(junctions,separators=(',',':'))
    return {'isoform':'SC_'+hashlib.sha256(key.encode()).hexdigest()[:24],'gene':gene,
            'chromosome':chrom,'strand':strand,'structure':tokens,
            'exons':exons if strand=='+' else exons[::-1],
            'structure_scope':'observed_splice_chain','full_length_validated':False}


def discover_events(gene,chrom,strand,paths,anchor=10):
    """Observed local competing junction paths; no significance/author HVG gate."""
    junctions={j for path in paths for j in path};result={}
    def add(kind,left,right,branches):
        branches=sorted(set(tuple(tuple(j) for j in b) for b in branches))
        if len(branches)<2:return
        key=json.dumps([gene,chrom,strand,branches],separators=(',',':'))
        eid='EV_'+hashlib.sha256(key.encode()).hexdigest()[:24]
        result[eid]={'event':eid,'gene':gene,'chromosome':chrom,'strand':strand,
          'kind':kind,'left':left,'right':right,'anchors':[[left-anchor,left],[right,right+anchor]],
          'anchor_bp':anchor,'branches':[[list(j) for j in b] for b in branches],
          'interpretation':'observed common-anchor conditional local branch choice; not unbiased transcript PSI'}
    for path in paths:
        for a,b in zip(path,path[1:]):
            if (a[0],b[1]) in junctions and a[1]<b[0]:add('cassette',a[0],b[1],[[a,b],[(a[0],b[1])]])
    byleft=collections.defaultdict(set);byright=collections.defaultdict(set)
    for a,b in junctions:byleft[a].add(b);byright[b].add(a)
    for left,ends in byleft.items():
        if len(ends)>1:add('shared_left',left,max(ends),[[(left,e)] for e in ends])
    for right,starts in byright.items():
        if len(starts)>1:add('shared_right',min(starts),right,[[(s,right)] for s in starts])
    return list(result.values())


def observation_path_groups(observations):
    groups=collections.defaultdict(list)
    for cell,count,info in observations:groups[tuple(info[2])].append((cell,count,info[3]))
    result={}
    for path,rows in groups.items():
        blocks=np.asarray([r[2] for r in rows],np.int64)
        result[path]={'cells':np.asarray([r[0] for r in rows],np.int64),
                      'counts':np.asarray([r[1] for r in rows]),
                      'starts':blocks[:,:,0],'ends':blocks[:,:,1]}
    return result


def count_event_observations(groups,event):
    """Same eligibility as chain_event_branch, vectorized over exon-end variants."""
    branches={tuple(map(tuple,branch)):k for k,branch in enumerate(event['branches'])}
    counts=collections.Counter()
    for path,data in groups.items():
        local=tuple(j for j in path if j[0]>=event['left'] and j[1]<=event['right'])
        branch=branches.get(local)
        if branch is None:continue
        mask=np.ones(len(data['cells']),bool)
        for left,right in event['anchors']:
            mask&=np.any((data['starts']<=left)&(data['ends']>=right),axis=1)
        if not mask.any():continue
        cells=data['cells'][mask];weights=data['counts'][mask]
        unique,inverse=np.unique(cells,return_inverse=True);sums=np.zeros(len(unique),np.result_type(weights.dtype,np.int64))
        np.add.at(sums,inverse,weights)
        for cell,count in zip(unique,sums):counts[int(cell),branch]+=count.item()
    return counts


def read_mex(path):
    with tarfile.open(path,'r:gz') as tar:
        members={Path(m.name).name:m for m in tar if m.isfile()}
        def stream(name):
            m=members.get(name) or members.get(name+'.gz')
            if m is None:raise ValueError('Missing MEX member '+name)
            s=tar.extractfile(m)
            return gzip.GzipFile(fileobj=s) if m.name.endswith('.gz') else s
        genes=pd.read_csv(stream('features.tsv'),sep='\t',header=None,names=['gene','gene_name','feature_type'])
        barcodes=pd.read_csv(stream('barcodes.tsv'),header=None,names=['barcode']).barcode.str.replace(r'-1$','',regex=True)
        matrix=mmread(stream('matrix.mtx')).T.tocsr()
    genes['gene']=genes.gene.map(stable_id)
    if genes.gene.duplicated().any() or barcodes.duplicated().any():raise ValueError('Duplicate MEX gene/barcode axis')
    if matrix.shape!=(len(barcodes),len(genes)):raise ValueError('MEX dimensions inconsistent')
    return matrix,genes,barcodes


def reference_genes(gtf):
    records={}
    opener=gzip.open if str(gtf).endswith('.gz') else open
    with opener(gtf,'rt') as f:
        for line in f:
            if line.startswith('#'):continue
            r=line.rstrip().split('\t')
            if r[2]!='gene':continue
            attr=dict(re.findall(r'(\w+) "([^"]*)"',r[8]));gid=stable_id(attr['gene_id'])
            records[gid]={'chromosome':r[0],'strand':r[6],'gene_type':attr.get('gene_type',''),
                          'gene_name':attr.get('gene_name',gid)}
    return records


def import_allinfo(db,input_root,output_root):
    names=['read_id','gene_versioned','label','barcode','umi','introns','tss','polya','exons']
    db.execute('CREATE TABLE IF NOT EXISTS raw_imports(sample VARCHAR PRIMARY KEY,rows BIGINT)')
    for sample,(lr,donor,region) in SAMPLES.items():
        if db.execute('SELECT 1 FROM raw_imports WHERE sample=?',[sample]).fetchone():continue
        path=input_root/(lr+'_AllInfo.tsv.tar');fifo=output_root/'allinfo_stream.tsv'
        fifo.unlink(missing_ok=True);os.mkfifo(fifo)
        process=subprocess.Popen(['bash','-o','pipefail','-c',
            'tar -xOf "$1" --wildcards "*.tsv.gz" | gzip -dc > "$2"','bash',str(path),str(fifo)],stderr=subprocess.PIPE)
        columns='{'+','.join("'"+n+"':'VARCHAR'" for n in names)+'}'
        try:
            db.execute('BEGIN')
            db.execute(f"CREATE OR REPLACE TABLE incoming_allinfo AS SELECT * FROM read_csv(?,header=false,delim='\t',quote='',escape='',parallel=false,columns={columns},strict_mode=true)",[str(fifo)])
            _,error=process.communicate(timeout=120)
            if process.returncode:raise RuntimeError(error.decode())
            db.execute('CREATE TABLE IF NOT EXISTS allinfo AS SELECT ?::VARCHAR sample,* FROM incoming_allinfo LIMIT 0',[sample])
            db.execute('INSERT INTO allinfo SELECT ?,* FROM incoming_allinfo',[sample])
            n=db.execute('SELECT count(*) FROM incoming_allinfo').fetchone()[0]
            db.execute('INSERT INTO raw_imports VALUES (?,?)',[sample,n]);db.execute('DROP TABLE incoming_allinfo');db.execute('COMMIT')
            db.execute('CHECKPOINT');log('AllInfo_imported',sample=sample,rows=n)
        except BaseException:
            process.terminate();process.wait()
            try:db.execute('ROLLBACK')
            except Exception:pass
            raise
        finally:fifo.unlink(missing_ok=True)


def read_annotations(path):
    annotations={};conflict=set()
    with gzip.open(path,'rt') as f:
        for line in f:
            barcode,label=line.rstrip('\n').split('\t');parts=label.split('::')
            if len(parts)!=6 or parts[0]!='P21':continue
            age,donor,region,broad,typ,subtype=parts
            sample='P21_'+donor+('_HIPP' if region=='Hippocampus' else '_VIS' if region=='VisCortex' else '_unsupported')
            if sample not in SAMPLES:continue
            key=(sample,barcode.removesuffix('-1'))
            value={'author_label':label,'cell_type_broad':broad,'cell_type':typ,'cell_subtype':subtype,
                   'donor':donor,'region':region}
            if key in annotations and annotations[key]!=value:conflict.add(key)
            annotations[key]=value
    if conflict:raise ValueError('Conflicting author cell labels: '+str(len(conflict)))
    return annotations


def build_molecules(db):
    # UMI uniqueness is per sample, cell and gene. Contradictory gene assignment
    # for the same read or contradictory chains for the same gene-UMI are excluded.
    db.execute("""CREATE OR REPLACE TABLE read_conflicts AS SELECT sample,read_id FROM allinfo
        GROUP BY sample,read_id HAVING count(DISTINCT (gene_versioned,barcode,umi,introns))>1""")
    db.execute("""CREATE OR REPLACE TABLE molecule_groups AS SELECT a.sample,a.barcode,
        split_part(gene_versioned,'.',1) gene,umi,count(*) supporting_read_rows,
        count(DISTINCT introns) distinct_chains,
        arg_max(struct_pack(introns:=introns,exons:=exons,label:=label,read_id:=read_id,tss:=tss,polya:=polya),
                (try_cast(regexp_extract(exons,'_([0-9]+)_[+-]$',1) AS BIGINT)
                 -try_cast(regexp_extract(exons,'^;%;[^;]+_([0-9]+)_[0-9]+_[+-]',1) AS BIGINT),read_id)) representative
        FROM allinfo a ANTI JOIN read_conflicts r USING(sample,read_id)
        WHERE regexp_full_match(barcode,'[ACGT]{16}') AND regexp_full_match(umi,'[ACGT]+')
          AND starts_with(gene_versioned,'ENSMUSG')
        GROUP BY a.sample,a.barcode,split_part(gene_versioned,'.',1),umi""")
    db.execute('CREATE OR REPLACE TABLE molecule_variants AS SELECT sample,barcode,gene,representative.introns introns,representative.exons exons,count(*) molecules FROM molecule_groups WHERE distinct_chains=1 GROUP BY ALL')
    db.execute('CHECKPOINT')


def convert(input_root,output_root,gtf,anchor=10,min_molecules=10,min_cells=3,reuse_molecules=False):
    input_root=Path(input_root).resolve();output_root=Path(output_root).resolve();output_root.mkdir(parents=True,exist_ok=True)
    conversion_lock=(output_root/'.conversion.lock').open('w')
    fcntl.flock(conversion_lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
    identity={'source_manifest_sha256':digest(input_root/'download_manifest.json'),'gtf_sha256':digest(gtf),
              'source_path':str(input_root),'schema':'allinfo_9cols_inclusive_to_halfopen_v1'}
    identity_file=output_root/'source_identity.json'
    if identity_file.exists() and json.loads(identity_file.read_text())!=identity:
        raise ValueError('Source manifest or annotation changed; use a new output directory')
    dump(identity_file,identity)
    state={'status':'running','stage':'import','started_at':time.strftime('%Y-%m-%dT%H:%M:%S%z')}
    dump(output_root/'conversion_status.json',state)
    db=duckdb.connect(str(output_root/'conversion_staging.duckdb'));db.execute("SET memory_limit='14GB'");db.execute('SET threads=4')
    db.execute('SET temp_directory=?',[str(output_root/'_tmp')])
    import_allinfo(db,input_root,output_root)
    log('deduplicating');state['stage']='deduplicate';dump(output_root/'conversion_status.json',state)
    if reuse_molecules:
        expected=db.execute('SELECT count(*) FROM molecule_groups WHERE distinct_chains=1').fetchone()[0]
        actual=db.execute('SELECT sum(molecules) FROM molecule_variants').fetchone()[0]
        if actual!=expected:raise ValueError('Existing molecule stage is inconsistent')
        log('reusing_validated_molecule_stage',molecules=actual)
    else:build_molecules(db)
    annotations=read_annotations(input_root/'SR_fullData_combined_barcodeCTassignments.gz');reference=reference_genes(gtf)
    good_cell_rows=db.execute("SELECT sample,barcode,count(*) FROM molecule_groups WHERE distinct_chains=1 GROUP BY sample,barcode").fetchall()
    lr_cells={(s,b):n for s,b,n in good_cell_rows};frames=[];matrices=[];mex_audit=[];feature_ref=None
    for sample,(_,donor,region) in SAMPLES.items():
        G,features,barcodes=read_mex(input_root/(sample+'.mex.tar.gz'))
        if feature_ref is None:feature_ref=features
        elif features[['gene','gene_name']].values.tolist()!=feature_ref[['gene','gene_name']].values.tolist():raise ValueError('MEX feature order differs across samples')
        selected=[j for j,b in enumerate(barcodes) if (sample,b) in lr_cells and (sample,b) in annotations]
        matrices.append(G[selected])
        for j in selected:
            barcode=barcodes.iloc[j];label=annotations[sample,barcode]
            frames.append({'cell':sample+':'+barcode,'sample_id':sample,'barcode':barcode,'batch':sample,
                'stage':'P21','stage_order':0,'trajectory_direction_available':'false',
                'LR_nonconflicting_gene_UMIs':lr_cells[sample,barcode],**label})
        mex_audit.append({'sample':sample,'shortread_cells':len(barcodes),'author_annotation_cells':sum(s==sample for s,b in annotations),
            'LR_UMI_cells':sum(s==sample for s,b in lr_cells),'retained_paired_annotated_cells':len(selected)})
        log('paired_cells',**mex_audit[-1])
    cells=pd.DataFrame(frames);G=sparse.vstack(matrices,format='csr');genes=feature_ref.drop(columns=['feature_type']).copy()
    for key in ['chromosome','strand','gene_type']:genes[key]=[reference.get(g,{}).get(key,'') for g in genes.gene]
    library=np.asarray(G.sum(1)).ravel();positive=library>0
    if not positive.all():G=G[positive];cells=cells.loc[positive].reset_index(drop=True);library=library[positive]
    gmap={g:i for i,g in enumerate(genes.gene)};cellmap={(r.sample_id,r.barcode):i for i,r in enumerate(cells.itertuples())}
    cells['complete_gene_UMI_library']=library
    db.execute('CREATE OR REPLACE TABLE paired_cells(sample VARCHAR,barcode VARCHAR,cell_index BIGINT)')
    db.execute('BEGIN')
    db.executemany('INSERT INTO paired_cells VALUES (?,?,?)',[(s,b,i) for (s,b),i in cellmap.items()])
    db.execute('COMMIT')
    db.execute('CREATE OR REPLACE TABLE paired_variants AS SELECT v.*,c.cell_index FROM molecule_variants v JOIN paired_cells c USING(sample,barcode) ORDER BY gene')
    db.execute('CREATE INDEX IF NOT EXISTS paired_variant_gene ON paired_variants(gene)')
    gene_totals=db.execute('SELECT gene,sum(molecules) FROM paired_variants GROUP BY gene ORDER BY gene').fetchall()
    state['stage']='structures_and_events';dump(output_root/'conversion_status.json',state)
    isoforms=[];chain_rows=[];chain_cols=[];chain_vals=[];events=[];event_rows=[];event_cols=[];event_vals=[];offsets=[0];qc=collections.Counter()
    for n,(gene,total) in enumerate(gene_totals):
        if gene not in gmap:qc['LR_gene_absent_from_G_molecules']+=int(total);continue
        variants=db.execute('SELECT introns,exons,cell_index,sum(molecules) FROM paired_variants WHERE gene=? GROUP BY ALL',[gene]).fetchall()
        parsed={};valid=[];paths=set();chain_groups=collections.defaultdict(list)
        for introntext,exontext,cell,count in variants:
            key=(introntext,exontext)
            if key not in parsed:
                try:parsed[key]=parse_structure(*key)
                except (ValueError,TypeError):parsed[key]=None
            info=parsed[key]
            if info is None:qc['invalid_structure_molecules']+=int(count);continue
            chrom,strand,junctions,exons=info
            if not junctions:qc['unspliced_structure_not_used_as_AS_molecules']+=int(count);continue
            ref=reference.get(gene,{})
            if ref and (ref['chromosome']!=chrom or ref['strand']!=strand):qc['gene_coordinate_conflict_molecules']+=int(count);continue
            qc['valid_structure_molecules']+=int(count);valid.append((int(cell),int(count),info));paths.add(tuple(junctions))
            chain_groups[(chrom,strand,tuple(junctions))].append((int(cell),int(count),exons))
        if not valid:continue
        for (chrom,strand,path),members in sorted(chain_groups.items()):
            # Conservative common exon coverage prevents representative-chain endpoints
            # from asserting coverage absent in another counted molecule.
            E=[list(v) for v in members[0][2]]
            for _,_,ex in members[1:]:
                E[0][0]=max(E[0][0],ex[0][0]);E[-1][1]=min(E[-1][1],ex[-1][1])
            if any(a>=b for a,b in E):qc['chain_common_coverage_unresolved']+=1;continue
            rec=chain_record(gene,chrom,strand,list(path),E);rec['gene_name']=genes.iloc[gmap[gene]].gene_name
            col=len(isoforms);isoforms.append(rec)
            bycell=collections.Counter()
            for cell,count,_ in members:bycell[cell]+=count
            for cell,count in bycell.items():chain_rows.append(cell);chain_cols.append(col);chain_vals.append(count)
        coordinate_groups=collections.defaultdict(list)
        for cell,count,info in valid:coordinate_groups[(info[0],info[1])].append((cell,count,info))
        for (chrom,strand),observations in coordinate_groups.items():
            these_paths={tuple(info[2]) for _,_,info in observations}
            path_groups=observation_path_groups(observations)
            for event in discover_events(gene,chrom,strand,these_paths,anchor):
                counts=count_event_observations(path_groups,event)
                branch_totals=collections.Counter();observed_cells=set()
                for (cell,branch),count in counts.items():branch_totals[branch]+=count;observed_cells.add(cell)
                if sum(counts.values())<min_molecules or len(observed_cells)<min_cells or len(branch_totals)<2:continue
                event.update(gene_name=genes.iloc[gmap[gene]].gene_name,branch_count=len(event['branches']),
                    total_UMIs=sum(counts.values()),informative_cells=len(observed_cells))
                events.append(event);base=offsets[-1]
                for (cell,branch),count in counts.items():event_rows.append(cell);event_cols.append(base+branch);event_vals.append(count)
                offsets.append(base+len(event['branches']))
        if (n+1)%250==0:log('gene_conversion',genes=n+1,total_genes=len(gene_totals),chains=len(isoforms),events=len(events))
    if not events or not isoforms:raise ValueError('No usable events or chains after conversion')
    I=sparse.coo_matrix((chain_vals,(chain_rows,chain_cols)),shape=(len(cells),len(isoforms)),dtype=np.int64).tocsr()
    C=sparse.coo_matrix((event_vals,(event_rows,event_cols)),shape=(len(cells),offsets[-1]),dtype=np.int64).tocsr()
    gp=output_root/'genes';ep=output_root/'events';gp.mkdir(exist_ok=True);ep.mkdir(exist_ok=True)
    for folder in [gp,ep]:cells.to_csv(folder/'cells.tsv',sep='\t',index=False)
    genes.to_csv(gp/'genes.tsv',sep='\t',index=False);sparse.save_npz(gp/'gene_counts.npz',G);sparse.save_npz(gp/'isoform_counts.npz',I)
    records=[]
    for r in isoforms:records.append({**r,'structure':json.dumps(r['structure']),'exons':json.dumps(r['exons'])})
    pd.DataFrame(records).to_csv(gp/'isoforms.tsv',sep='\t',index=False)
    sparse.save_npz(ep/'branch_counts.npz',C);np.save(ep/'branch_offsets.npy',np.asarray(offsets,np.int64));dump(ep/'groups.json',events)
    # Select HVGs from the complete raw paired gene counts; never the author significant list.
    import anndata as ad
    import scanpy as sc
    a=ad.AnnData(G.astype(float),var=genes.set_index('gene'))
    sc.pp.highly_variable_genes(a,flavor='seurat_v3',n_top_genes=3000,batch_key=None,subset=False)
    hv=a.var['highly_variable'].to_numpy();genes.loc[hv].to_csv(ep/'genes.tsv',sep='\t',index=False)
    pd.DataFrame({'cell':cells.cell,'complete_gene_UMI_library':library,'scale_to_100000':1e5/library}).to_csv(ep/'normalization_factors.tsv',sep='\t',index=False)
    qc.update(raw_AllInfo_rows=db.execute('SELECT count(*) FROM allinfo').fetchone()[0],
        ambiguous_read_ids=db.execute('SELECT count(*) FROM read_conflicts').fetchone()[0],
        unique_gene_UMI_groups=db.execute('SELECT count(*) FROM molecule_groups').fetchone()[0],
        conflicting_chain_UMI_groups=db.execute('SELECT count(*) FROM molecule_groups WHERE distinct_chains>1').fetchone()[0])
    report={'status':'complete','schema':'splicefate_author_allinfo_mex_v1','assembly':'mm10','source_coordinates':'AllInfo_1based_inclusive',
        'output_coordinates':'0based_halfopen','samples':mex_audit,'cells':len(cells),'genes':len(genes),'HVGs':int(hv.sum()),
        'chains':len(isoforms),'events':len(events),'branches':offsets[-1],'count_QC':dict(qc),
        'gene_counts':'paired_shortread_integer_UMIs','S_counts':'exact_gene_UMI_deduplicated_observed_common_anchor_branches',
        'UMI_rule':'sample+barcode+stable_gene+exact_UMI; conflicting intron chains excluded; no new edit-distance collapse',
        'duplicate_read_rule':'contradictory assignment for same read ID within sample excluded',
        'representative_rule':'same-chain UMI representative with widest observed genomic exon span (tie read ID); not full-length truth',
        'chain_exons':'intersection of observed terminal coverage within same junction chain, no terminal imputation',
        'event_definition':'observed cassette and shared donor/acceptor local competitions; shared 10bp anchors required',
        'event_filter':{'min_total_molecules':min_molecules,'min_cells':min_cells,'min_observed_branches':2,'author_significance_filter':False,'HVG_restricted_S':False},
        'HVG_method':'seurat_v3_3000_from_raw_paired_counts; stored axis order follows genes table',
        'normalization':'complete_gene_library_to_100000_then_log1p_in_model; matrices_saved_raw',
        'intron_calibration':'unavailable; intron chain lists excised gaps, not retained/nascent RNA',
        'trajectory':'P21_only; no root inferred from cell subtype or file ordering; specify a biological root for directed DPT',
        'no_cell_holdout':True,'no_S_imputation':True,'IR_training':False,
        'structure_full_length_claim':False,'overlapping_events':'may share molecules; composite likelihood',
        'source_manifest_sha256':digest(input_root/'download_manifest.json'),'annotation':str(gtf),'annotation_sha256':digest(gtf)}
    dump(output_root/'conversion_report.json',report)
    # Required files are audited with dimensions, count conservation and unique keys.
    column_totals=np.asarray(C.sum(0)).ravel()
    checks={'G_dimensions':G.shape==(len(cells),len(genes)),'I_dimensions':I.shape==(len(cells),len(isoforms)),
        'C_dimensions':C.shape==(len(cells),offsets[-1]),'unique_cells':not cells.cell.duplicated().any(),
        'HVG3000':int(hv.sum())==3000,'positive_G_libraries':bool(np.all(library>0)),
        'event_ids_unique':len({e['event'] for e in events})==len(events),
        'chain_ids_unique':len({r['isoform'] for r in isoforms})==len(isoforms),
        'event_totals_match':all(int(column_totals[offsets[j]:offsets[j+1]].sum())==e['total_UMIs'] for j,e in enumerate(events)),
        'integer_counts':all(np.all(x.data>=0) and np.all(x.data==np.floor(x.data)) for x in [G,I,C])}
    dump(output_root/'validation.json',{'passed':all(checks.values()),'checks':checks})
    state.update(status='complete' if all(checks.values()) else 'validation_failed',stage='complete',finished_at=time.strftime('%Y-%m-%dT%H:%M:%S%z'))
    dump(output_root/'conversion_status.json',state);db.close();log('complete',cells=len(cells),events=len(events),chains=len(isoforms))
    if not all(checks.values()):raise ValueError('Conversion validation failed')


def validate_output(output_root):
    """Linear sparse-column reduction, no repeated whole-CSR scans per event."""
    root=Path(output_root);gp=root/'genes';ep=root/'events'
    cells=pd.read_csv(gp/'cells.tsv',sep='\t',dtype=str);genes=pd.read_csv(gp/'genes.tsv',sep='\t',dtype=str)
    ec=pd.read_csv(ep/'cells.tsv',sep='\t',dtype=str);hvg=pd.read_csv(ep/'genes.tsv',sep='\t',dtype=str)
    chains=pd.read_csv(gp/'isoforms.tsv',sep='\t',usecols=['isoform','gene']);events=json.loads((ep/'groups.json').read_text())
    G=sparse.load_npz(gp/'gene_counts.npz');I=sparse.load_npz(gp/'isoform_counts.npz');C=sparse.load_npz(ep/'branch_counts.npz')
    offsets=np.load(ep/'branch_offsets.npy');norm=pd.read_csv(ep/'normalization_factors.tsv',sep='\t')
    column_totals=np.asarray(C.sum(0)).ravel();libraries=np.asarray(G.sum(1)).ravel()
    source_report=json.loads((root/'conversion_report.json').read_text())
    checks={'G_dimensions':G.shape==(len(cells),len(genes)),'I_dimensions':I.shape==(len(cells),len(chains)),
       'C_dimensions':C.shape==(len(cells),int(offsets[-1])),
       'identical_cell_axes':cells.cell.tolist()==ec.cell.tolist()==norm.cell.tolist(),
       'unique_cells':not cells.cell.duplicated().any(),'unique_genes':not genes.gene.duplicated().any(),
       'unique_chain_ids':not chains.isoform.duplicated().any(),'HVG3000':len(hvg)==3000 and hvg.gene.is_unique,
       'HVG_subset':set(hvg.gene)<=set(genes.gene),'chain_genes_available':set(chains.gene)<=set(genes.gene),
       'event_genes_available':{e['gene'] for e in events}<=set(genes.gene),
       'offset_shape':len(offsets)==len(events)+1 and offsets[0]==0,
       'branch_counts_match_offsets':all(offsets[j+1]-offsets[j]==len(e['branches']) for j,e in enumerate(events)),
       'event_ids_unique':len({e['event'] for e in events})==len(events),
       'event_totals_match':all(int(column_totals[offsets[j]:offsets[j+1]].sum())==e['total_UMIs'] for j,e in enumerate(events)),
       'positive_libraries':bool(np.all(libraries>0)),
       'normalization_denominators_match':np.array_equal(libraries,norm.complete_gene_UMI_library.to_numpy()),
       'normalization_scale_match':bool(np.allclose(1e5/libraries,norm.scale_to_100000.to_numpy())),
       'observed_chain_count_conservation':int(I.sum())==source_report['count_QC']['valid_structure_molecules'],
       'integer_nonnegative_counts':all(np.isfinite(x.data).all() and np.all(x.data>=0) and np.all(x.data==np.floor(x.data)) for x in [G,I,C])}
    report={'passed':all(checks.values()),'checks':checks,'cells':len(cells),'genes':len(genes),'chains':len(chains),
            'events':len(events),'branches':int(offsets[-1]),'G_UMIs':int(G.sum()),'observed_chain_UMIs':int(I.sum()),
            'event_support_sum_not_independent_UMIs':int(C.sum()),'validation_time':time.strftime('%Y-%m-%dT%H:%M:%S%z')}
    dump(root/'validation.json',report)
    status={'status':'complete' if report['passed'] else 'validation_failed','stage':'complete',
            'finished_at':report['validation_time']};dump(root/'conversion_status.json',status)
    if not report['passed']:raise ValueError('Converted dataset validation failed')
    log('validated',**{k:v for k,v in report.items() if k!='checks'})
    return report


def main():
    p=argparse.ArgumentParser(description='Author ScISOr-Seq2 AllInfo + matched MEX to SpliceFate inputs')
    p.add_argument('--input',required=True);p.add_argument('--output',required=True);p.add_argument('--gtf',required=True)
    p.add_argument('--anchor-bp',type=int,default=10);p.add_argument('--min-event-molecules',type=int,default=10);p.add_argument('--min-event-cells',type=int,default=3)
    p.add_argument('--reuse-molecule-stage',action='store_true',help='Resume validated existing source-identical molecule aggregation')
    p.add_argument('--validate-only',action='store_true',help='Validate already written model inputs without reprocessing molecules')
    a=p.parse_args()
    try:
        if a.validate_only:validate_output(a.output)
        else:convert(a.input,a.output,a.gtf,a.anchor_bp,a.min_event_molecules,a.min_event_cells,a.reuse_molecule_stage)
    except BaseException as exc:
        out=Path(a.output);out.mkdir(parents=True,exist_ok=True);dump(out/'conversion_status.json',{'status':'failed','error':str(exc),'time':time.strftime('%Y-%m-%dT%H:%M:%S%z')});raise


if __name__=='__main__':main()
