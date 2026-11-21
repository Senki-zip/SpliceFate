"""Build native human hg38 model priors from retained human source records.

Human identity/motifs/CLIP reuse the multispecies atlas; protein relations use
native human STRING physical and human curated complexes. No mouse projection.
"""
from pathlib import Path
import collections
import csv
import gzip
import json
import re
import sqlite3
import time

import duckdb

from build_mouse_rbp_database import create_rows,sha256,sqlstr,source_rows
from integrate_rbp_resources import rbp2go_rows
from splicefate.convert import reference_genes

PROJECT=Path(__file__).resolve().parents[1]
ROOT=PROJECT/'data/external/human_rbp_unified'
ATLAS=PROJECT/'data/external/rbp_atlas_20261115'
GTF=PROJECT.parent/'reference_documents/refdata-cellranger-arc-GRCh38-2024-A/genes/genes.gtf.gz'


def main():
    ROOT.mkdir(parents=True,exist_ok=True)
    status={'status':'running','started_at':time.strftime('%Y-%m-%dT%H:%M:%S%z'),'species':'Homo_sapiens','assembly':'hg38'}
    (ROOT/'build_status.json').write_text(json.dumps(status,indent=2))
    db=duckdb.connect(str(ROOT/'human_rbp.building.duckdb'));db.execute("SET memory_limit='12GB'");db.execute('SET threads=4')
    old=sqlite3.connect((ATLAS/'integrated/rbp_atlas.sqlite').as_uri()+'?mode=ro',uri=True)
    ref=reference_genes(GTF);symbols=collections.defaultdict(set)
    for gid,r in ref.items():symbols[r['gene_name'].upper()].add(gid)
    raw_factors=old.execute("SELECT * FROM factors WHERE species='Homo_sapiens'").fetchall()
    identities=[];factors={}
    for fid,species,gid,symbol in raw_factors:
        candidates=symbols.get(symbol.upper(),set())
        resolved=gid.split('.')[0] if gid.startswith('ENSG') else next(iter(candidates)) if len(candidates)==1 else gid
        target='Homo_sapiens:'+resolved;identities.append((fid,target))
        label=ref.get(resolved,{}).get('gene_name',symbol)
        factors[target]=(target,species,resolved,label)
    mapping=dict(identities)
    create_rows(db,'factors','factor_id VARCHAR,species VARCHAR,gene_id VARCHAR,gene_symbol VARCHAR',factors.values())
    create_rows(db,'family_sharing_factor_identity','original_factor_id VARCHAR,canonical_factor_id VARCHAR',
                sorted(set(identities+[(f,f) for f in factors])))
    names=['factor_sources','factor_motifs','functions','rbp2go_scores','experiments']
    for name in names:
        schema=old.execute(f'PRAGMA table_info({name})').fetchall();cols=[r[1] for r in schema]
        types={'TEXT':'VARCHAR','INTEGER':'BIGINT','REAL':'DOUBLE'}
        columns=','.join('"'+r[1]+'" '+types[r[2]] for r in schema)
        def rows(name=name,cols=cols):
            query=f"SELECT t.* FROM {name} t JOIN factors f ON t.factor_id=f.factor_id WHERE f.species='Homo_sapiens'"
            for row in old.execute(query):
                row=list(row);row[cols.index('factor_id')]=mapping[row[cols.index('factor_id')]];yield row
        n=create_rows(db,name,columns,rows());print(json.dumps({'table':name,'rows':n}),flush=True)
    for name,query in [
        ('motif_sources',"SELECT DISTINCT ms.* FROM motif_sources ms JOIN factor_motifs fm USING(source,source_motif_id) JOIN factors f USING(factor_id) WHERE f.species='Homo_sapiens'"),
        ('motifs',"SELECT DISTINCT m.* FROM motifs m JOIN motif_sources ms USING(motif_id) JOIN factor_motifs fm USING(source,source_motif_id) JOIN factors f USING(factor_id) WHERE f.species='Homo_sapiens'")]:
        schema=old.execute(f'PRAGMA table_info({name})').fetchall()
        columns=','.join('"'+r[1]+'" '+{'TEXT':'VARCHAR','INTEGER':'BIGINT','REAL':'DOUBLE'}[r[2]] for r in schema)
        create_rows(db,name,columns,old.execute(query))
    old.close()
    aliases=set()
    for f,_,g,s in factors.values():aliases.update([(f,g),(f,s)])
    for f,sid,symbol in db.execute('SELECT factor_id,source_id,gene_symbol FROM factor_sources').fetchall():
        if sid:aliases.add((f,sid))
        if symbol:aliases.add((f,symbol))
    create_rows(db,'factor_aliases','factor_id VARCHAR,alias VARCHAR',sorted(aliases))
    db.execute("""CREATE OR REPLACE VIEW candidate_splicing_factors AS SELECT DISTINCT f.* FROM factors f
       JOIN functions fn USING(factor_id) WHERE lower(fn.label) LIKE '%splic%' OR fn.term_id='GO:0008380'""")
    db.execute("""CREATE OR REPLACE VIEW valid_factor_motifs AS SELECT fm.*,ms.motif_id
        FROM factor_motifs fm JOIN motif_sources ms USING(source,source_motif_id) WHERE ms.status='valid'""")
    # Reuse the exact InterPro Type=Family rule with native human donors/recipients.
    ipr={r['ID']:r for r in source_rows(ATLAS/'raw/rbp2go/download_interpro_data.export')}
    uid_to_factor={u:f for f,u in db.execute('SELECT factor_id,uniprot_id FROM rbp2go_scores').fetchall()}
    members=collections.defaultdict(set);membership=[]
    for file in ['all-download_species_rbp.export','all-download_species_non_listed.export']:
        for row in rbp2go_rows(ATLAS/'raw/rbp2go'/file):
            if row['Species']!='Homo sapiens' or row['Uniprot_ID'] not in uid_to_factor:continue
            f=uid_to_factor[row['Uniprot_ID']]
            for ip in set(re.findall(r'IPR\d+',row.get('InterPro',''))):
                if ipr.get(ip,{}).get('Type')=='Family':
                    members[ip].add(f);membership.append((f,ip,ipr[ip]['Name'],row['Uniprot_ID']))
    create_rows(db,'motif_family_membership','factor_id VARCHAR,family_id VARCHAR,family_name VARCHAR,uniprot_id VARCHAR',membership)
    native=collections.defaultdict(set)
    for f,mid in db.execute("SELECT DISTINCT factor_id,motif_id FROM valid_factor_motifs WHERE source IN ('CISBP-RNA','ATtRACT') AND evidence!='mutant_binding_exclude_from_WT_default'").fetchall():native[f].add(mid)
    supported={f for (f,) in db.execute("SELECT DISTINCT factor_id FROM rbp2go_scores WHERE source_role='RBP_catalog' UNION SELECT DISTINCT factor_id FROM factor_sources WHERE source!='RBP2GO'").fetchall()}
    shared=set();shared_ev=[]
    for ip,group in members.items():
        for recipient in (group&supported)-set(native):
            for donor in group&set(native):
                for mid in native[donor]:
                    shared.add((recipient,mid));shared_ev.append((recipient,mid,donor,ip,'family_shared_candidate'))
    create_rows(db,'family_shared_motifs','factor_id VARCHAR,motif_id VARCHAR',sorted(shared))
    create_rows(db,'family_shared_motif_evidence','factor_id VARCHAR,motif_id VARCHAR,donor_factor_id VARCHAR,family_id VARCHAR,extension_label VARCHAR',shared_ev)
    db.execute("""CREATE OR REPLACE VIEW effective_factor_motifs AS SELECT DISTINCT factor_id,source,source_motif_id,motif_id,evidence FROM valid_factor_motifs
        UNION ALL SELECT factor_id,'FamilySharing',motif_id,motif_id,'family_shared_candidate' FROM family_shared_motifs""")
    # Source error replacements remain attached; no interval liftover required.
    overrides=json.loads((ATLAS/'peak_overrides.json').read_text()) if (ATLAS/'peak_overrides.json').exists() else {}
    experiments=db.execute('SELECT experiment_id,peak_file,peak_rows FROM experiments ORDER BY experiment_id').fetchall()
    def peaks():
        for n,(eid,rel,expected) in enumerate(experiments):
            path=ATLAS/(overrides.get(eid,{}).get('usable_peak_file',rel));actual=0
            with gzip.open(path,'rt') as stream:
                for line in stream:
                    if not line.strip() or line.lstrip().startswith(('#','track','browser')):continue
                    c=line.rstrip().split('\t')
                    if len(c)<6 or c[5] not in ['+','-','.']:raise ValueError('Invalid native hg38 peak '+eid)
                    start,end=int(c[1]),int(c[2])
                    if not 0<=start<end:raise ValueError('Invalid CLIP interval')
                    yield ('hg38',c[0],start,end,c[5],eid,eid+':'+str(actual));actual+=1
            if actual!=expected:raise ValueError('Peak count mismatch '+eid)
            if (n+1)%200==0:print(json.dumps({'CLIP_experiments':n+1,'total':len(experiments)}),flush=True)
    create_rows(db,'clip_peaks_original','assembly VARCHAR,chromosome VARCHAR,start0 BIGINT,end0 BIGINT,strand VARCHAR,experiment_id VARCHAR,peak_id VARCHAR',peaks())
    db.execute("""CREATE OR REPLACE VIEW candidate_nonperturbed_binding_experiments AS SELECT e.* FROM experiments e
        WHERE perturbation_flag=0 AND content_duplicate_of='' AND NOT EXISTS
        (SELECT 1 FROM experiments b WHERE b.independent_evidence_key=e.independent_evidence_key
         AND b.factor_id=e.factor_id AND b.assay=e.assay AND b.content_duplicate_of='' AND b.experiment_id<e.experiment_id)""")
    symbol_map=collections.defaultdict(set)
    for f,alias in aliases:symbol_map[alias.upper()].add(f)
    info={}
    with gzip.open(ROOT/'raw/string/9606.protein.info.v12.5.txt.gz','rt') as stream:
        for row in csv.DictReader(stream,delimiter='\t'):info[row['#string_protein_id']]=row['preferred_name']
    candidates=collections.defaultdict(lambda:collections.defaultdict(set))
    genes={g:f for f,_,g,s in factors.values() if g.startswith('ENSG')}
    with gzip.open(ROOT/'raw/string/9606.protein.aliases.v12.5.txt.gz','rt') as stream:
        for row in csv.DictReader(stream,delimiter='\t'):
            alias=row['alias'];pid=row['#string_protein_id']
            if alias.split('.')[0] in genes:candidates[pid][1].add(genes[alias.split('.')[0]])
            elif alias in uid_to_factor:candidates[pid][2].add(uid_to_factor[alias])
    for pid,name in info.items():
        if name.upper() in symbol_map:candidates[pid][3].update(symbol_map[name.upper()])
    smap={}
    for pid,ranked in candidates.items():
        hits=ranked[min(ranked)]
        if len(hits)==1:smap[pid]=next(iter(hits))
    create_rows(db,'string_identity_mapping','protein_id VARCHAR,factor_id VARCHAR',sorted(smap.items()))
    string_path=ROOT/'raw/string/9606.protein.physical.links.full.v12.5.txt.gz'
    db.execute('CREATE OR REPLACE TABLE string_physical AS SELECT * FROM read_csv(?,header=true,delim=\' \',parallel=false)',[str(string_path)])
    db.execute("""CREATE OR REPLACE TABLE candidate_experiment_supported_protein_interactions AS SELECT
        a.factor_id factor_a,b.factor_id factor_b,'STRING_physical12.5'::VARCHAR source,
        'physical_association'::VARCHAR relation_type,s.protein1||':'||s.protein2 evidence_id,
        experiments::BIGINT nontransferred_experiment_score,false construct_change_flag,false perturbation_keyword_flag
        FROM string_physical s JOIN string_identity_mapping a ON a.protein_id=s.protein1
        JOIN string_identity_mapping b ON b.protein_id=s.protein2 WHERE s.protein1<s.protein2 AND experiments>0""")
    comp=[]
    cp=ROOT.parent/'rbp_interactions_20261116/raw/complexportal/9606.tsv'
    with cp.open() as stream:
        for row in csv.DictReader(stream,delimiter='\t'):
            # Match the original curated physical subset; inferred classes stay out.
            if 'ECO:0000353' not in row['Evidence Code']:continue
            tokens=row.get('Expanded participant list') or row['Identifiers (and stoichiometry) of molecules in complex']
            for token in tokens.split('|'):
                uid=token.split('(')[0].strip()
                if uid in uid_to_factor:comp.append(('ComplexPortal',row['#Complex ac'],uid_to_factor[uid]))
    for row in source_rows(ATLAS/'raw/rbp2go/download_corum_complexes.export'):
        if row.get('Organism')!='Human':continue
        for uid in re.split('[; ]+',row.get('subunits(UniProt IDs)','')):
            if uid in uid_to_factor:comp.append(('CORUM',row['ComplexID'],uid_to_factor[uid]))
    create_rows(db,'curated_complex_memberships','source VARCHAR,complex_id VARCHAR,factor_id VARCHAR',sorted(set(comp)))
    db.execute('CREATE OR REPLACE VIEW candidate_curated_complex_memberships AS SELECT * FROM curated_complex_memberships')
    create_rows(db,'build_metadata','key VARCHAR,value VARCHAR',[
        ('species','Homo_sapiens'),('primary_assembly','hg38'),('schema_version','model_prior_human_1.0'),
        ('scope','Human CISBP/ATtRACT/RBP2GO/ENCORI; native STRING physical, CORUM and ComplexPortal'),
        ('family_sharing','human InterPro Type=Family one-hop motif candidates; explicit labels'),
        ('limitations','No human IntAct imported in this model prior build; oRNAment sites are not needed for native sequence rescanning'),
        ('reference_gtf',str(GTF)),('reference_gtf_sha256',sha256(GTF)),('created_at',time.strftime('%Y-%m-%dT%H:%M:%S%z'))])
    checks={
        'human_only':db.execute("SELECT count(*) FROM factors WHERE species!='Homo_sapiens'").fetchone()[0]==0,
        'factor_unique':db.execute('SELECT count(*)-count(DISTINCT factor_id) FROM factors').fetchone()[0]==0,
        'native_CLIP_rows':db.execute('SELECT count(*) FROM clip_peaks_original').fetchone()[0]==sum(x[2] for x in experiments),
        'native_hg38_only':db.execute("SELECT count(*) FROM clip_peaks_original WHERE assembly!='hg38'").fetchone()[0]==0,
        'motif_factor_fk':db.execute('SELECT count(*) FROM effective_factor_motifs e LEFT JOIN factors f USING(factor_id) WHERE f.factor_id IS NULL').fetchone()[0]==0,
        'motif_matrix_fk':db.execute('SELECT count(*) FROM effective_factor_motifs e LEFT JOIN motifs m USING(motif_id) WHERE m.motif_id IS NULL').fetchone()[0]==0,
        'PPI_factor_fk':db.execute('SELECT count(*) FROM candidate_experiment_supported_protein_interactions e LEFT JOIN factors a ON a.factor_id=e.factor_a LEFT JOIN factors b ON b.factor_id=e.factor_b WHERE a.factor_id IS NULL OR b.factor_id IS NULL').fetchone()[0]==0,
        'RBFOX2_resolves':db.execute("SELECT count(*) FROM factors WHERE gene_symbol='RBFOX2'").fetchone()[0]==1}
    tables=['factors','candidate_splicing_factors','motifs','effective_factor_motifs','family_shared_motifs','experiments',
            'clip_peaks_original','candidate_experiment_supported_protein_interactions','curated_complex_memberships']
    report={'passed':all(checks.values()),'checks':checks,'tables':{t:db.execute('SELECT count(*) FROM '+t).fetchone()[0] for t in tables},
            'species':'Homo_sapiens','assembly':'hg38','source_atlas_sha256':sha256(ATLAS/'integrated/rbp_atlas.sqlite')}
    db.execute('CHECKPOINT');db.close();(ROOT/'validation.json').write_text(json.dumps(report,indent=2))
    if not report['passed']:raise ValueError(report)
    temporary=ROOT/'human_rbp.building.duckdb';final=ROOT/'human_rbp.duckdb'
    if final.exists():raise FileExistsError('Preserving existing human model prior')
    temporary.rename(final);(ROOT/'SHA256SUMS').write_text(sha256(final)+'  human_rbp.duckdb\n')
    status.update(status='complete',validation_passed=True,finished_at=time.strftime('%Y-%m-%dT%H:%M:%S%z'))
    (ROOT/'build_status.json').write_text(json.dumps(status,indent=2));print(json.dumps(report,indent=2),flush=True)


if __name__=='__main__':main()
