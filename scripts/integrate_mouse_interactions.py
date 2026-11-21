"""Import audited protein interactions into the production mouse evidence DB.

Transactional, idempotent extension; original factor/motif/site tables unchanged.
"""
from __future__ import annotations
import collections
import csv
import hashlib
import json
from pathlib import Path
import re
import shutil
import time

import duckdb

from build_mouse_rbp_database import ROOT,DBPATH,create_rows,sqlstr,sha256

INPUT=ROOT.parent/'rbp_interactions_20261116'
PREVIEW=INPUT/'audit/interaction_preview.duckdb'


def main():
    report_path=ROOT/'interaction_integration.json'
    report={'status':'running','started_at':time.strftime('%Y-%m-%dT%H:%M:%S%z'),
            'main_database':str(DBPATH),'preview_database':str(PREVIEW)}
    for p in [INPUT/'validation.json',INPUT/'audit/coverage_summary.json']:
        if not json.loads(p.read_text())['passed']:raise ValueError('Upstream validation failed: '+str(p))
    for r in json.loads((INPUT/'download_manifest.json').read_text()):
        p=INPUT/r['path']
        if p.stat().st_size!=r['bytes'] or sha256(p)!=r['sha256']:raise ValueError('Changed input: '+str(p))
    backup=ROOT/'backups/mouse_rbp.before_interactions_20261116.duckdb'
    backup.parent.mkdir(exist_ok=True)
    before_hash=sha256(DBPATH)
    if not backup.exists():
        shutil.copy2(DBPATH,backup)
        if sha256(backup)!=before_hash:raise ValueError('Backup checksum mismatch')
    report.update(backup=str(backup),backup_sha256=sha256(backup),input_database_sha256=before_hash,
                  preview_sha256=sha256(PREVIEW))
    report_path.write_text(json.dumps(report,indent=2))
    db=duckdb.connect(str(DBPATH));db.execute("SET memory_limit='10GB'");db.execute('SET threads=4')
    db.execute('SET file_search_path='+sqlstr(ROOT))
    db.execute('ATTACH '+sqlstr(PREVIEW)+' AS incoming (READ_ONLY)')
    protected=['factors','factor_sources','motifs','motif_sources','factor_motifs','motif_evidence',
               'experiments','gene_binding','family_shared_motifs','clip_peaks_original','clip_peaks_mm39']
    before={t:db.execute(f'SELECT count(*) FROM {t}').fetchone()[0] for t in protected}
    try:
        db.execute('BEGIN TRANSACTION')
        imported={}
        for src,dest in [
            ('string_identity_mapping','interaction_string_identity_mapping'),
            ('string_physical_mapped','interaction_string_physical_scores'),
            ('string_functional_mapped','interaction_string_functional_scores'),
            ('intact_mouse_evidence','interaction_intact_raw'),
            ('complexes','complexes'),('complex_memberships','complex_memberships'),
            ('spliceosome_stage_memberships','spliceosome_stage_memberships'),
            ('spliceosome_experiments','spliceosome_experiments')]:
            db.execute(f'CREATE OR REPLACE TABLE {dest} AS SELECT * FROM incoming.{src}')
            imported[dest]=db.execute(f'SELECT count(*) FROM {dest}').fetchone()[0]
            print(json.dumps({'imported':dest,'rows':imported[dest]}),flush=True)
        # Preserve explicit negatives in their own source table. They never become positive edges.
        def negatives():
            with (INPUT/'raw/intact/mouse_negative.txt').open() as f:
                for row in csv.reader(f,delimiter='\t'):
                    if row and not row[0].startswith('#'):
                        if len(row)!=42:raise ValueError('Unexpected negative MITAB width')
                        yield (hashlib.sha256(json.dumps(row,separators=(',',':')).encode()).hexdigest(),
                               row[0],row[1],row[9],row[10],row[13],row[8],row[35].lower()=='true',json.dumps(row))
        create_rows(db,'interaction_negative_records','record_key VARCHAR,source_id_a VARCHAR,source_id_b VARCHAR,taxon_a VARCHAR,taxon_b VARCHAR,interaction_id VARCHAR,publication VARCHAR,negative BOOLEAN,record_json VARCHAR',negatives())
        # Complete MITAB metadata stay accessible and mutation/condition flags are surfaced.
        db.execute("""CREATE OR REPLACE TABLE interaction_intact_evidence AS
            WITH metadata AS (SELECT *,json_extract_string(record_json,'$[28]') host_organism,
                json_extract_string(record_json,'$[12]') original_database,
                json_extract_string(record_json,'$[18]') experimental_role_a,
                json_extract_string(record_json,'$[19]') experimental_role_b,
                json_extract_string(record_json,'$[20]') interactor_type_a,
                json_extract_string(record_json,'$[21]') interactor_type_b,
                json_extract_string(record_json,'$[25]') annotation_a,
                json_extract_string(record_json,'$[26]') annotation_b,
                json_extract_string(record_json,'$[27]') interaction_annotation,
                json_extract_string(record_json,'$[36]') features_a,
                json_extract_string(record_json,'$[37]') features_b,
                json_extract_string(record_json,'$[30]') source_created_date,
                json_extract_string(record_json,'$[31]') source_updated_date FROM interaction_intact_raw)
            SELECT *,sha256(record_json) AS source_record_key,
              CASE WHEN contains(interaction_type,'MI:0407') THEN 'direct_physical'
                   WHEN contains(interaction_type,'MI:0915') THEN 'physical_association'
                   WHEN contains(interaction_type,'MI:0914') THEN 'association_unspecified'
                   WHEN contains(interaction_type,'MI:0403') THEN 'colocalization'
                   WHEN contains(interaction_type,'MI:2364') THEN 'proximity'
                   WHEN contains(interaction_type,'MI:0195') OR contains(interaction_type,'MI:0408') THEN 'covalent_binding'
                   ELSE 'other_typed_relation' END relation_type,
              contains(interactor_type_a,'MI:0326') AND contains(interactor_type_b,'MI:0326') both_annotated_proteins,
              regexp_matches(lower(coalesce(features_a,'')||' '||coalesce(features_b,'')),
                 'mutation|mutagenesis|deletion|substitution') construct_change_flag,
              regexp_matches(lower(coalesce(annotation_a,'')||' '||coalesce(annotation_b,'')||' '||
                 coalesce(interaction_annotation,'')),'knock.?out|knock.?down|crispr|sirna|shrna|over.?express|deplet|treat') perturbation_keyword_flag,
              regexp_matches(lower(coalesce(annotation_a,'')||' '||coalesce(annotation_b,'')||' '||
                 coalesce(interaction_annotation,'')),'rnase|rna.depend|rna.independ') rna_dependence_text_flag,
              'unknown'::VARCHAR rna_dependence,'not_reviewed'::VARCHAR wt_review_status,
              NULL::VARCHAR regulatory_sign,NULL::DOUBLE effect_size,
              'source_metadata_not_cell_type_matched'::VARCHAR context_status
            FROM metadata""")
        # Source evidence rows are not independent replication counts. Exact source duplicates
        # share a key; IntAct accessions and publications remain available for cross-source audits.
        db.execute("""CREATE OR REPLACE TABLE protein_interaction_evidence AS
            SELECT 'IntAct:'||source_record_key evidence_id,'IntAct'::VARCHAR source,
              factor_a,factor_b,source_id_a,source_id_b,relation_type,
              detection_method,interaction_id source_interaction_id,publication,confidence,
              host_organism,expansion_method,'10090'::VARCHAR taxon_a,'10090'::VARCHAR taxon_b,
              negative,construct_change_flag,perturbation_keyword_flag,rna_dependence_text_flag,
              rna_dependence,wt_review_status,regulatory_sign,effect_size,
              relation_type IN ('direct_physical','physical_association','association_unspecified','covalent_binding')
                AND contains(source_id_a,'uniprotkb:') AND contains(source_id_b,'uniprotkb:')
                AND NOT negative AS eligible_positive_protein_relation,
              NULL::BIGINT nontransferred_experiment_score,NULL::BIGINT transferred_experiment_score,
              NULL::BIGINT combined_score
            FROM interaction_intact_evidence QUALIFY row_number() OVER(PARTITION BY source_record_key ORDER BY record_id)=1
            UNION ALL
            SELECT 'STRING12.5:'||sha256(protein1||'|'||protein2),'STRING_physical12.5',
              factor_a,factor_b,protein1,protein2,'physical_network_association',
              NULL,NULL,NULL,NULL,NULL,NULL,'10090','10090',false,false,false,false,
              'unknown','not_reviewed',NULL,NULL,true,experiments,experiments_transferred,combined_score
            FROM interaction_string_physical_scores""")
        db.execute("""CREATE OR REPLACE VIEW positive_protein_interactions AS
            SELECT * FROM protein_interaction_evidence WHERE eligible_positive_protein_relation
            AND factor_a IS NOT NULL AND factor_b IS NOT NULL""")
        db.execute("""CREATE OR REPLACE VIEW candidate_experiment_supported_protein_interactions AS
            SELECT * FROM positive_protein_interactions WHERE source='IntAct' OR nontransferred_experiment_score>0""")
        db.execute("""CREATE OR REPLACE VIEW candidate_unmodified_protein_interactions AS
            SELECT * FROM candidate_experiment_supported_protein_interactions
            WHERE NOT construct_change_flag AND NOT perturbation_keyword_flag""")
        db.execute("""CREATE OR REPLACE TABLE protein_interaction_pairs AS
            SELECT least(factor_a,factor_b) factor_a,greatest(factor_a,factor_b) factor_b,
                count(*) source_record_support_count,
                list(DISTINCT source ORDER BY source) sources,
                bool_or(source='IntAct') has_IntAct_association,
                bool_or(relation_type='direct_physical') has_direct_physical_annotation,
                bool_or(source='STRING_physical12.5' AND nontransferred_experiment_score>0) has_STRING_nontransferred_experiment_channel,
                bool_or(source='STRING_physical12.5' AND transferred_experiment_score>0) has_STRING_transferred_experiment_channel,
                bool_or(construct_change_flag OR perturbation_keyword_flag) has_flagged_construct_or_condition,
                factor_a=factor_b is_same_canonical_factor,
                'support_count_not_independent_replications'::VARCHAR support_count_interpretation
            FROM positive_protein_interactions GROUP BY least(factor_a,factor_b),greatest(factor_a,factor_b),factor_a=factor_b""")
        db.execute("""CREATE OR REPLACE TABLE interaction_source_accessions AS
            SELECT e.evidence_id,regexp_extract(token,'^([^:]+):',1) namespace,
              regexp_replace(token,'^[^:]+:','') accession
            FROM protein_interaction_evidence e,unnest(string_split(e.source_interaction_id,'|')) t(token)
            WHERE e.source='IntAct' AND token IS NOT NULL AND token!='-'""")
        db.execute("""CREATE OR REPLACE TABLE interaction_publications AS
            SELECT DISTINCT e.evidence_id,unnest(regexp_extract_all(publication,'pubmed:([0-9]+)',1)) pmid
            FROM protein_interaction_evidence e WHERE publication IS NOT NULL""")
        db.execute("""CREATE OR REPLACE VIEW interaction_evidence_reuse_candidates AS
            SELECT 'shared_source_accession'::VARCHAR reason,namespace||':'||accession source_reference,
              count(DISTINCT evidence_id) source_rows FROM interaction_source_accessions
              GROUP BY namespace,accession HAVING count(DISTINCT evidence_id)>1""")
        # No pair expansion from complexes and no human-to-mouse stage projection.
        db.execute("""CREATE OR REPLACE VIEW mouse_complex_memberships AS
          SELECT m.*,c.name,c.evidence_code,c.references_text FROM complex_memberships m
          JOIN complexes c USING(source,complex_id,taxon) WHERE m.taxon='10090'""")
        db.execute("""CREATE OR REPLACE VIEW mouse_curated_complex_memberships AS
          SELECT * FROM mouse_complex_memberships
          WHERE evidence_class='curated_physical_evidence' OR source='CORUM'""")
        db.execute("""CREATE OR REPLACE VIEW mouse_spliceosome_stage_memberships AS
          SELECT * FROM spliceosome_stage_memberships WHERE taxon='10090' AND factor_id IS NOT NULL""")
        # Preserve raw source protein/genomic mappings; no need to invent factors for unknowns.
        db.execute("""CREATE OR REPLACE TABLE interaction_identity_map AS
          SELECT 'STRING12.5'::VARCHAR source,protein_id source_id,factor_id,mapping_method,
            candidates_json FROM interaction_string_identity_mapping""")
        provenance=[]
        for r in json.loads((INPUT/'download_manifest.json').read_text()):
            provenance.append((r['source'],r['path'],r['url'],r.get('timestamp'),r['bytes'],r['sha256'],
                               r.get('etag'),r.get('last_modified')))
        extra_manifest=INPUT/'spliceosome_download_manifest.json'
        for r in json.loads(extra_manifest.read_text()) if extra_manifest.exists() else []:
            if r.get('path') not in {x[1] for x in provenance} and 'sha256' in r:
                provenance.append((r['source'],r['path'],r['url'],None,r['bytes'],r['sha256'],None,None))
        create_rows(db,'interaction_input_assets','source VARCHAR,path VARCHAR,url VARCHAR,retrieved_at VARCHAR,bytes BIGINT,sha256 VARCHAR,etag VARCHAR,last_modified VARCHAR',provenance)
        snapshots={r['source']:r.get('timestamp','')[:10].replace('-','') for r in json.loads((INPUT/'download_manifest.json').read_text())}
        corum_dates=[json.loads(line)['recorded_at'][:10].replace('-','') for line in (ROOT.parent/'rbp_atlas_20261115/downloads.jsonl').read_text().splitlines()
                     if json.loads(line).get('path')=='raw/rbp2go/download_corum_complexes.export' and json.loads(line).get('status')=='downloaded']
        create_rows(db,'interaction_sources','source VARCHAR,version VARCHAR,scope VARCHAR,limitations VARCHAR',[
            ('IntAct','current_snapshot_'+snapshots.get('IntAct','unknown'),'mouse-mouse molecular origin; host retained separately','Association does not establish direct contact; construct changes and RNase status require review.'),
            ('STRING_physical12.5','12.5','mouse physical network','Contains co-complex, transferred and text-derived evidence; scores are not effects or independent experiments.'),
            ('STRING_functional12.5','12.5','mouse functional association','Separate score table only; not automatically converted into physical evidence.'),
            ('ComplexPortal','current_snapshot_'+snapshots.get('ComplexPortal','unknown'),'mouse and human complexes','Evidence tags retained; human members not assigned to mouse; no clique expansion.'),
            ('CORUM','RBP2GO_snapshot_'+(corum_dates[-1] if corum_dates else 'unknown'),'existing Mouse complex records','Reused existing records; not an additional independent dataset.'),
            ('SpliceosomeDB','snapshot_'+snapshots.get('SpliceosomeDB','unknown'),'native mouse/human catalog and source experiments','Mouse stage information is not inferred; human stage information not projected to mouse.')])
        create_rows(db,'interaction_policy','key VARCHAR,value VARCHAR',[
            ('schema_version','1.0'),('snapshot_date',time.strftime('%Y-%m-%d')),
            ('identity','Use existing canonical mouse factors; ambiguous source IDs remain unmapped.'),
            ('physical_evidence','Direct physical, physical association and unspecified association remain separate.'),
            ('functional_network','STRING functional scores are separate from physical evidence.'),
            ('complex_policy','Membership only; never infer every pair contacts directly.'),
            ('human_projection','Not performed; human records keep taxon 9606 and null mouse factor ID.'),
            ('RNA_dependence','Unknown unless curated from a specific experimental comparison; text flag only prompts review.'),
            ('WT_policy','Metadata flags are exploratory screens; no dataset is automatically certified WT.'),
            ('causal_effect','No regulatory sign, KO effect or new factor-RNA target inferred.'),
            ('dedup','AB/BA collapsed for STRING; exact IntAct source rows share a key; canonical factor pairs aggregate references, never independent replicate counts.'),
            ('source_reuse','Shared accessions and PMIDs permit review; same PMID alone does not collapse different experiments.'),
        ])
        # A compact per-factor coverage table is generated from production relations.
        db.execute("""CREATE OR REPLACE TABLE factor_interaction_coverage AS
            WITH endpoints AS (SELECT factor_a factor_id,source,relation_type,nontransferred_experiment_score
                FROM candidate_experiment_supported_protein_interactions UNION ALL
                SELECT factor_b,source,relation_type,nontransferred_experiment_score FROM candidate_experiment_supported_protein_interactions),
            interactions AS (SELECT factor_id,count(*) source_endpoint_records,
                bool_or(source='IntAct') has_IntAct_association,
                bool_or(relation_type='direct_physical') has_direct_physical_annotation,
                bool_or(source='STRING_physical12.5' AND nontransferred_experiment_score>0) has_STRING_nontransferred_experiment
                FROM endpoints GROUP BY factor_id),
            comp AS (SELECT factor_id,count(DISTINCT (source,complex_id)) curated_complexes
                FROM mouse_curated_complex_memberships WHERE factor_id IS NOT NULL GROUP BY factor_id)
            SELECT f.*,coalesce(i.source_endpoint_records,0) source_endpoint_records,
                coalesce(i.has_IntAct_association,false) has_IntAct_association,
                coalesce(i.has_direct_physical_annotation,false) has_direct_physical_annotation,
                coalesce(i.has_STRING_nontransferred_experiment,false) has_STRING_nontransferred_experiment,
                coalesce(c.curated_complexes,0) curated_complexes,
                i.factor_id IS NOT NULL OR c.factor_id IS NOT NULL AS has_experiment_channel_or_curated_complex,
                EXISTS(SELECT 1 FROM effective_factor_motifs v WHERE v.factor_id=f.factor_id
                     AND v.source!='ENCORI' AND v.evidence!='mutant_binding_exclude_from_WT_default') has_source_or_family_motif,
                EXISTS(SELECT 1 FROM candidate_splicing_factors s JOIN family_sharing_factor_identity m
                     ON m.original_factor_id=s.factor_id WHERE m.canonical_factor_id=f.factor_id) broad_splicing_candidate
            FROM factors f LEFT JOIN interactions i USING(factor_id) LEFT JOIN comp c USING(factor_id)
            WHERE f.factor_id IN (SELECT DISTINCT canonical_factor_id FROM family_sharing_factor_identity)""")
        after={t:db.execute(f'SELECT count(*) FROM {t}').fetchone()[0] for t in protected}
        checks={'protected_base_counts_unchanged':before==after,
            'evidence_ids_unique':db.execute('SELECT count(*)-count(DISTINCT evidence_id) FROM protein_interaction_evidence').fetchone()[0]==0,
            'complex_keys_unique':db.execute('SELECT count(*)-count(DISTINCT (source,complex_id,taxon)) FROM complexes').fetchone()[0]==0,
            'complex_member_fk':db.execute('SELECT count(*) FROM complex_memberships m LEFT JOIN complexes c USING(source,complex_id,taxon) WHERE c.complex_id IS NULL').fetchone()[0]==0,
            'human_not_projected':db.execute("SELECT count(*) FROM complex_memberships WHERE taxon='9606' AND factor_id IS NOT NULL").fetchone()[0]==0,
            'no_negative_in_positive':db.execute('SELECT count(*) FROM positive_protein_interactions WHERE negative').fetchone()[0]==0,
            'no_colocalization_or_proximity_in_positive':db.execute("SELECT count(*) FROM positive_protein_interactions WHERE relation_type IN ('colocalization','proximity','other_typed_relation')").fetchone()[0]==0,
            'RNA_dependence_not_invented':db.execute("SELECT count(*) FROM protein_interaction_evidence WHERE rna_dependence!='unknown'").fetchone()[0]==0,
            'regulatory_sign_not_invented':db.execute('SELECT count(*) FROM protein_interaction_evidence WHERE regulatory_sign IS NOT NULL OR effect_size IS NOT NULL').fetchone()[0]==0,
            'mouse_stage_still_empty':db.execute('SELECT count(*) FROM mouse_spliceosome_stage_memberships').fetchone()[0]==0,
            'functional_not_inserted_as_physical':db.execute("SELECT count(*) FROM protein_interaction_evidence WHERE source='STRING_functional12.5'").fetchone()[0]==0}
        for t,cols in [('protein_interaction_evidence',['factor_a','factor_b']),('complex_memberships',['factor_id']),
                       ('spliceosome_stage_memberships',['factor_id']),('interaction_identity_map',['factor_id'])]:
            for col in cols:
                checks[t+'_'+col+'_fk']=db.execute(f'SELECT count(*) FROM {t} e LEFT JOIN factors f ON e.{col}=f.factor_id WHERE e.{col} IS NOT NULL AND f.factor_id IS NULL').fetchone()[0]==0
        report.update(checks=checks,passed=all(checks.values()),protected_counts=after,imported_tables=imported)
        if not report['passed']:raise ValueError('Failed integration checks: '+str(checks))
        report['statistics']={
            'source_evidence_records':db.execute('SELECT count(*) FROM protein_interaction_evidence').fetchone()[0],
            'positive_canonical_pairs':db.execute('SELECT count(*) FROM protein_interaction_pairs').fetchone()[0],
            'negative_source_records':db.execute('SELECT count(*) FROM interaction_negative_records').fetchone()[0],
            'IntAct_construct_flagged_rows':db.execute('SELECT count(*) FROM interaction_intact_evidence WHERE construct_change_flag').fetchone()[0],
            'IntAct_condition_flagged_rows':db.execute('SELECT count(*) FROM interaction_intact_evidence WHERE perturbation_keyword_flag').fetchone()[0],
            'canonical_splicing_candidates':db.execute('SELECT count(*) FROM factor_interaction_coverage WHERE broad_splicing_candidate').fetchone()[0],
            'motifless_splicing_candidates_with_relations':db.execute('SELECT count(*) FROM factor_interaction_coverage WHERE broad_splicing_candidate AND NOT has_source_or_family_motif AND has_experiment_channel_or_curated_complex').fetchone()[0]}
        db.execute('COMMIT');db.execute('CHECKPOINT')
    except BaseException as e:
        try:db.execute('ROLLBACK')
        except Exception:pass
        report.update(status='failed',error=str(e));report_path.write_text(json.dumps(report,indent=2));raise
    finally:db.close()
    report.update(status='complete',finished_at=time.strftime('%Y-%m-%dT%H:%M:%S%z'),output_database_sha256=sha256(DBPATH))
    report_path.write_text(json.dumps(report,indent=2))
    # Refresh the production checksum; no Parquet parts have changed.
    p=ROOT/'SHA256SUMS';lines=p.read_text().splitlines()
    lines=[report['output_database_sha256']+'  mouse_rbp.duckdb' if x.endswith('  mouse_rbp.duckdb') else x for x in lines]
    p.write_text('\n'.join(lines)+'\n')
    print(json.dumps(report,indent=2),flush=True)


if __name__=='__main__':main()
