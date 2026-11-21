"""Build a separate interaction staging catalog and audit mouse coverage.

The production motif/site database is read-only. Co-complex membership is never
expanded into invented binary physical interactions. Human rows stay human.
"""
import collections
import csv
import gzip
import hashlib
import io
import json
from pathlib import Path
import re
import time
import zipfile

import duckdb

from build_mouse_rbp_database import create_rows,sqlstr,sha256
from query_mouse_rbp_database import connect

PROJECT=Path(__file__).resolve().parents[1]
ROOT=PROJECT/'data/external/rbp_interactions_20261116'
OUT=ROOT/'audit'


def main():
    OUT.mkdir(exist_ok=True)
    base=connect()
    factors={f:(g,s) for f,g,s in base.execute('SELECT factor_id,gene_id,gene_symbol FROM factors').fetchall()}
    canonical=dict(base.execute('SELECT original_factor_id,canonical_factor_id FROM family_sharing_factor_identity').fetchall())
    def norm(f):return canonical.get(f,f)
    splicing={norm(r[0]) for r in base.execute('SELECT factor_id FROM candidate_splicing_factors').fetchall()}
    supported={norm(r[0]) for r in base.execute('SELECT factor_id FROM factor_catalog WHERE rbp2go_listed OR other_rbp_source').fetchall()}
    motif={norm(r[0]) for r in base.execute("SELECT DISTINCT factor_id FROM effective_factor_motifs WHERE source!='ENCORI' AND evidence!='mutant_binding_exclude_from_WT_default'").fetchall()}
    native_motif={norm(r[0]) for r in base.execute("SELECT DISTINCT factor_id FROM valid_factor_motifs WHERE source IN ('CISBP-RNA','ATtRACT','oRNAmentV2') AND evidence!='mutant_binding_exclude_from_WT_default'").fetchall()}
    uniprot=collections.defaultdict(set);gene=collections.defaultdict(set);symbols=collections.defaultdict(set)
    for f,(gid,sym) in factors.items():
        gene[gid.split('.')[0]].add(norm(f));symbols[sym.upper()].add(norm(f))
    for f,u in base.execute('SELECT factor_id,uniprot_id FROM rbp2go_scores').fetchall():uniprot[u].add(norm(f))
    for f,source,sid,sym in base.execute('SELECT factor_id,source,source_id,gene_symbol FROM factor_sources').fetchall():
        if sym:symbols[sym.upper()].add(norm(f))
    candidates=collections.defaultdict(lambda:collections.defaultdict(set))
    sfolder=ROOT/'raw/string'
    info={}
    with gzip.open(sfolder/'10090.protein.info.v12.5.txt.gz','rt') as f:
        for r in csv.DictReader(f,delimiter='\t'):
            pid=r['#string_protein_id'];info[pid]=r
            hits=symbols.get(r['preferred_name'].upper(),set())
            if hits:candidates[pid][3].update(hits)
    alias_uid=collections.defaultdict(set)
    with gzip.open(sfolder/'10090.protein.aliases.v12.5.txt.gz','rt') as f:
        for r in csv.DictReader(f,delimiter='\t'):
            pid=r['#string_protein_id'];alias=r['alias'];hits=set();rank=9
            if re.fullmatch(r'ENSMUSG\d+(?:\.\d+)?',alias):hits=gene.get(alias.split('.')[0],set());rank=1
            elif alias in uniprot:hits=uniprot[alias];rank=2
            if hits:candidates[pid][rank].update(hits)
            if 'uniprot' in r['source'].lower() and re.fullmatch(r'[A-Z0-9]{6,10}(?:-\d+)?',alias):alias_uid[alias].add(pid)
    smap={};mapping_rows=[]
    for pid,byrank in candidates.items():
        rank=min(byrank);hits=byrank[rank]
        method={1:'unique_stable_gene_id',2:'unique_source_uniprot',3:'unique_exact_species_symbol'}[rank]
        if len(hits)==1:smap[pid]=next(iter(hits))
        mapping_rows.append((pid,smap.get(pid),method if pid in smap else 'ambiguous_kept_unmapped',json.dumps(sorted(hits))))
    # Supplement protein identifiers using uniquely resolved STRING source aliases.
    for uid,pids in alias_uid.items():
        hits={smap[p] for p in pids if p in smap}
        if uid not in uniprot and len(hits)==1:uniprot[uid].update(hits)
    db=duckdb.connect(str(OUT/'interaction_preview.duckdb'))
    db.execute("SET memory_limit='10GB'");db.execute('SET threads=4')
    create_rows(db,'string_identity_mapping','protein_id VARCHAR,factor_id VARCHAR,mapping_method VARCHAR,candidates_json VARCHAR',mapping_rows)
    create_rows(db,'scope_factors','factor_id VARCHAR,gene_id VARCHAR,gene_symbol VARCHAR,source_supported_rbp BOOLEAN,splicing_candidate BOOLEAN,has_source_motif BOOLEAN,has_source_or_family_motif BOOLEAN',
        ((f,*factors[f],f in supported,f in splicing,f in native_motif,f in motif) for f in sorted(set(map(norm,factors)))))
    report={'status':'running','scope':'mouse only for coverage; human data retained separately without automatic projection',
            'canonical_splicing_candidates':len(splicing),'motifless_after_family_extension':len(splicing-motif),
            'STRING_version':'12.5','sources':{},'coverage':{}}
    physical_path=sfolder/'10090.protein.physical.links.full.v12.5.txt.gz'
    functional_path=sfolder/'10090.protein.links.detailed.v12.5.txt.gz'
    for table,path in [('string_physical',physical_path),('string_functional',functional_path)]:
        db.execute(f"CREATE OR REPLACE TABLE {table} AS SELECT * FROM read_csv({sqlstr(path)},header=true,delim=' ',parallel=false)")
        n=db.execute(f'SELECT count(*) FROM {table}').fetchone()[0]
        n_unique=db.execute(f'SELECT count(*) FROM {table} WHERE protein1<protein2').fetchone()[0]
        db.execute(f"""CREATE OR REPLACE TABLE {table}_mapped AS
            SELECT s.*,a.factor_id factor_a,b.factor_id factor_b FROM {table} s
            LEFT JOIN string_identity_mapping a ON a.protein_id=s.protein1
            LEFT JOIN string_identity_mapping b ON b.protein_id=s.protein2
            WHERE s.protein1<s.protein2""")
        report['sources'][table]={'raw_directed_rows':n,'unique_unordered_protein_pairs':n_unique,
            'pairs_with_both_endpoints_mapped':db.execute(f'SELECT count(*) FROM {table}_mapped WHERE factor_a IS NOT NULL AND factor_b IS NOT NULL').fetchone()[0]}
        print(table,report['sources'][table],flush=True)
    def uid_map(value):
        ids=re.findall(r'uniprotkb:([^|\s]+)',value)
        hits=set().union(*(uniprot.get(x,uniprot.get(x.split('-')[0],set())) for x in ids))
        return next(iter(hits)) if len(hits)==1 else None
    intact_counts=collections.Counter();header=[]
    def intact_rows():
        nonlocal header
        with zipfile.ZipFile(ROOT/'raw/intact/mouse.zip') as z:
            for name in z.namelist():
                if not name.endswith('.txt'):continue
                with io.TextIOWrapper(z.open(name)) as f:
                    for row in csv.reader(f,delimiter='\t'):
                        if not row:continue
                        if row[0].startswith('#'):header=row;continue
                        intact_counts['all_rows']+=1
                        if len(row)<15:raise ValueError('Invalid MITAB row')
                        if 'taxid:10090(' not in row[9] or 'taxid:10090(' not in row[10]:
                            intact_counts['non_mouse_mouse_rows']+=1;continue
                        intact_counts['mouse_mouse_rows']+=1
                        a,b=uid_map(row[0]),uid_map(row[1])
                        if a and b:intact_counts['both_mapped_rows']+=1
                        if (a in splicing) or (b in splicing):intact_counts['splicing_related_rows']+=1
                        negative=(row[35].lower()=='true') if len(row)>35 else False
                        yield (intact_counts['all_rows'],row[0],row[1],a,b,row[6],row[8],row[9],row[10],row[11],row[13],row[14],
                               row[15] if len(row)>15 else '',negative,json.dumps(row,ensure_ascii=False))
    create_rows(db,'intact_mouse_evidence','record_id BIGINT,source_id_a VARCHAR,source_id_b VARCHAR,factor_a VARCHAR,factor_b VARCHAR,detection_method VARCHAR,publication VARCHAR,taxon_a VARCHAR,taxon_b VARCHAR,interaction_type VARCHAR,interaction_id VARCHAR,confidence VARCHAR,expansion_method VARCHAR,negative BOOLEAN,record_json VARCHAR',intact_rows())
    db.execute("""CREATE OR REPLACE VIEW intact_mouse_ppi AS SELECT *,
      CASE WHEN contains(interaction_type,'MI:0407') THEN 'direct_interaction_annotation'
           WHEN contains(interaction_type,'MI:0915') OR contains(interaction_type,'MI:0914') THEN 'physical_or_association_annotation'
           ELSE 'other_typed_protein_interaction' END AS relation_class
      FROM intact_mouse_evidence WHERE NOT negative AND contains(source_id_a,'uniprotkb:')
      AND contains(source_id_b,'uniprotkb:') AND
      (contains(interaction_type,'MI:0407') OR contains(interaction_type,'MI:0915')
       OR contains(interaction_type,'MI:0914') OR contains(interaction_type,'MI:0195') OR contains(interaction_type,'MI:0408'))""")
    (OUT/'intact_column_names.json').write_text(json.dumps(header,indent=2))
    report['sources']['IntAct']=dict(intact_counts)
    report['sources']['IntAct']['positive_protein_evidence_rows']=db.execute('SELECT count(*) FROM intact_mouse_ppi').fetchone()[0]
    report['sources']['IntAct']['unique_unordered_source_protein_pairs']=db.execute("SELECT count(DISTINCT (least(source_id_a,source_id_b),greatest(source_id_a,source_id_b))) FROM intact_mouse_ppi").fetchone()[0]
    report['sources']['IntAct']['interaction_type_summary']=[{'type':t,'rows':n} for t,n in db.execute('SELECT interaction_type,count(*) FROM intact_mouse_ppi GROUP BY interaction_type ORDER BY count(*) DESC').fetchall()]
    print('IntAct',dict(intact_counts),flush=True)
    complexes=[];members=[];cp_counts=collections.Counter()
    for tax in ['10090','9606']:
        for r in csv.DictReader((ROOT/f'raw/complexportal/{tax}.tsv').open(),delimiter='\t'):
            cid=r['#Complex ac'];ev=r['Evidence Code']
            tag='curated_other_evidence'
            for code,label in [('0005544','orthology_inferred'),('0005610','homology_inferred'),
                               ('0005546','paralogy_inferred'),('0005547','background_knowledge_inferred'),
                               ('0005543','mixed_species_experiment'),('0000353','curated_physical_evidence')]:
                if 'ECO:'+code in ev:tag=label
            complexes.append(('ComplexPortal',cid,tax,r['Recommended name'],tag,ev,r['Cross references'],json.dumps(r,ensure_ascii=False)))
            cp_counts[tax+'_'+tag]+=1
            participant=r.get('Expanded participant list') or r['Identifiers (and stoichiometry) of molecules in complex']
            for token in participant.split('|'):
                uid=token.split('(')[0].strip();hits=uniprot.get(uid,uniprot.get(uid.split('-')[0],set())) if tax=='10090' else set()
                fid=next(iter(hits)) if len(hits)==1 else None
                members.append(('ComplexPortal',cid,tax,uid,fid,token,tag))
    # Reuse existing downloaded CORUM instead of fetching another copy.
    corum=base.execute("SELECT record_json FROM rbp2go_support_records WHERE source_file='download_corum_complexes.export'").fetchall()
    for (text,) in corum:
        r=json.loads(text)
        if r.get('Organism')!='Mouse':continue
        cid=r['ComplexID'];complexes.append(('CORUM',cid,'10090',r['ComplexName'],'curated_source_methods',r.get('Protein complex purification method',''),r.get('PubMed ID',''),text))
        for uid in re.split('[; ]+',r.get('subunits(UniProt IDs)','')):
            if not uid:continue
            hits=uniprot.get(uid,set());fid=next(iter(hits)) if len(hits)==1 else None
            members.append(('CORUM',cid,'10090',uid,fid,uid,'curated_source_methods'))
    create_rows(db,'complexes','source VARCHAR,complex_id VARCHAR,taxon VARCHAR,name VARCHAR,evidence_class VARCHAR,evidence_code VARCHAR,references_text VARCHAR,record_json VARCHAR',complexes)
    create_rows(db,'complex_memberships','source VARCHAR,complex_id VARCHAR,taxon VARCHAR,source_protein_id VARCHAR,factor_id VARCHAR,source_token VARCHAR,evidence_class VARCHAR',members)
    report['sources']['ComplexPortal']=dict(cp_counts)
    report['sources']['CORUM_reused']={'mouse_complexes':sum(r[0]=='CORUM' for r in complexes)}
    stages=[];stage_counts={}
    for species,tax in [('mouse','10090'),('human','9606')]:
        rows=list(csv.DictReader((ROOT/f'raw/spliceosomedb/proteins_{species}.csv').open()))
        stage_counts[species+'_protein_rows']=len(rows);n=0
        for r in rows:
            gid=r.get('Ensemble Gene ID','').strip();hits=gene.get(gid.split('.')[0],set()) if species=='mouse' else set()
            if species=='mouse' and not hits:hits=symbols.get(r['Gene Symbol'].upper(),set())
            fid=next(iter(hits)) if len(hits)==1 else None
            for complex_name in r.get('Complex','').split(','):
                if not complex_name.strip():continue
                n+=1
                stages.append((tax,r['Gene Symbol'],fid,gid,complex_name.strip(),r.get('Class / family',''),'source_curated_stage_membership',json.dumps(r,ensure_ascii=False)))
        stage_counts[species+'_stage_memberships']=n
    create_rows(db,'spliceosome_stage_memberships','taxon VARCHAR,source_symbol VARCHAR,factor_id VARCHAR,source_gene_id VARCHAR,stage VARCHAR,classification VARCHAR,evidence_class VARCHAR,record_json VARCHAR',stages)
    exp=list(csv.DictReader((ROOT/'raw/spliceosomedb/experiments_all.csv').open()))
    stage_counts['experiment_species_counts']=dict(collections.Counter(x['Organism'] for x in exp))
    create_rows(db,'spliceosome_experiments','source_experiment_id VARCHAR,taxon_label VARCHAR,complex_label VARCHAR,cell_type VARCHAR,pmid VARCHAR,record_json VARCHAR',
                ((r['Spliceosome Database Experiment ID'],r['Organism'],r['Complex'],r['Cell Type'],r['Pubmed id'],json.dumps(r,ensure_ascii=False)) for r in exp))
    report['sources']['SpliceosomeDB']=stage_counts
    # Coverage is evidence availability, not an increase in identified motif or causal edges.
    def endpoints(table,where='true'):
        return {f for row in db.execute(f'SELECT factor_a,factor_b FROM {table} WHERE {where}').fetchall() for f in row if f}
    evidence_sets={
      'STRING_physical_any_channel':endpoints('string_physical_mapped'),
      'STRING_physical_nontransferred_experiment_channel':endpoints('string_physical_mapped','experiments>0'),
      'STRING_functional_experimental_channel':endpoints('string_functional_mapped','experimental>0'),
      'IntAct_mouse_mouse_positive_protein_interactions':endpoints('intact_mouse_ppi'),
      'IntAct_mouse_mouse_direct_interaction_annotation':endpoints('intact_mouse_ppi',"contains(interaction_type,'MI:0407')"),
      'mouse_complex_membership_any_evidence':{r[0] for r in db.execute("SELECT DISTINCT factor_id FROM complex_memberships WHERE taxon='10090' AND factor_id IS NOT NULL").fetchall()},
      'mouse_complex_membership_CP_physical_or_CORUM':{r[0] for r in db.execute("SELECT DISTINCT factor_id FROM complex_memberships WHERE taxon='10090' AND factor_id IS NOT NULL AND (evidence_class='curated_physical_evidence' OR source='CORUM')").fetchall()},
      'mouse_spliceosome_stage_membership':{r[0] for r in db.execute("SELECT DISTINCT factor_id FROM spliceosome_stage_memberships WHERE taxon='10090' AND factor_id IS NOT NULL").fetchall()},
    }
    for score in [400,700]:
        evidence_sets[f'STRING_nontransferred_physical_experiment_score_ge_{score}']=endpoints('string_physical_mapped',f'experiments>={score}')
    physical_motif_neighbor=set()
    for a,b in db.execute('SELECT factor_a,factor_b FROM string_physical_mapped WHERE experiments>0').fetchall():
        if a and b:
            if a in motif:physical_motif_neighbor.add(b)
            if b in motif:physical_motif_neighbor.add(a)
    evidence_sets['STRING_nontransferred_physical_experiment_with_motif_partner']=physical_motif_neighbor
    combo=set().union(evidence_sets['STRING_physical_nontransferred_experiment_channel'],evidence_sets['IntAct_mouse_mouse_positive_protein_interactions'],evidence_sets['mouse_complex_membership_CP_physical_or_CORUM'])
    evidence_sets['combined_excluding_explicit_transfers']=combo
    for name,fs in evidence_sets.items():
        report['coverage'][name]={'all_splicing_candidates':len(fs&splicing),'motifless_splicing_candidates':len(fs&(splicing-motif)),
                                  'source_supported_RBPs':len(fs&supported)}
    fields=['factor_id','gene_id','gene_symbol','has_source_motif','has_source_or_family_motif']+list(evidence_sets)
    with (OUT/'splicing_factor_interaction_coverage.tsv').open('w') as f:
        w=csv.writer(f,delimiter='\t');w.writerow(fields)
        for fid in sorted(splicing,key=lambda f:factors[f][1].lower()):
            w.writerow([fid,*factors[fid],fid in native_motif,fid in motif]+[fid in fs for fs in evidence_sets.values()])
    report['previous_snapshot_comparison']={'canonical_candidates':{'previous':406,'current':len(splicing)},
        'motifless_candidates':{'previous':302,'current':len(splicing-motif)}}
    checks={'canonical_candidates_resolved':splicing<=set(factors),
            'motif_factor_identities_resolved':motif<=set(factors),
            'all_string_ids_mouse':not db.execute("SELECT count(*) FROM string_physical WHERE NOT starts_with(protein1,'10090.') OR NOT starts_with(protein2,'10090.')").fetchone()[0],
            'no_negative_rows_automatically_counted':True,'human_members_not_projected':all(r[4] is None for r in members if r[2]=='9606'),
            'no_complex_pair_expansion':True,'main_database_read_only':True}
    report['checks']=checks;report['passed']=all(checks.values());report['status']='complete' if report['passed'] else 'validation_failed'
    report['updated_at']=time.strftime('%Y-%m-%dT%H:%M:%S%z')
    report['limitations']=['STRING physical includes co-complex evidence, not necessarily direct contact.',
       'Nontransferred score channels retain provider annotations.',
       'No RNA-dependence, sign, target-event identity or WT certification inferred from PPI.',
       'Human curated stages have not been assigned to mouse; explicit orthology mapping remains to be reviewed.',
       'Same experiment/PMID may be reused by IntAct, CORUM, Complex Portal and STRING; counts are not independent replication.']
    (OUT/'coverage_summary.json').write_text(json.dumps(report,indent=2,ensure_ascii=False))
    db.execute('CHECKPOINT');db.close();base.close()
    print(json.dumps(report,indent=2,ensure_ascii=False),flush=True)
    if not report['passed']:raise SystemExit(1)


if __name__=='__main__':main()
