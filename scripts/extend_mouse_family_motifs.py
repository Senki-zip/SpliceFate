"""Add explicitly labeled, one-hop mouse family-sharing motif candidates.

Original motif associations remain unchanged. InterPro Family annotations, not
generic RRM/KH domains or gene-name prefixes, define sharing groups. This is a
candidate-generating rule with explicit source labels.
"""
from __future__ import annotations
import collections
import csv
import json
from pathlib import Path
import re
import time

import duckdb

from build_mouse_rbp_database import ROOT, DBPATH, ATLAS, create_rows, source_rows, sha256


def main():
    db=duckdb.connect(str(DBPATH))
    db.execute("SET file_search_path='"+str(ROOT).replace("'","''")+"'")
    factors={f:(g,s) for f,g,s in db.execute('SELECT factor_id,gene_id,gene_symbol FROM factors').fetchall()}
    profiles=[(f,json.loads(p)) for f,p in db.execute('SELECT factor_id,record_json FROM rbp2go_profiles').fetchall()]
    symbols=collections.defaultdict(set)
    for fid,(gid,symbol) in factors.items():
        if gid.startswith('ENSMUSG'):
            symbols[symbol.upper()].add(fid)
    # Only the unresolved protein's own explicit alias can point to an exact stable primary symbol.
    # Do not merge merely because some other protein happens to list a similar alias.
    equivalence={};identity_rows=[]
    for fid,r in profiles:
        if not factors[fid][0].startswith('symbol:'):continue
        hits=set().union(*(symbols.get(x.upper(),set()) for x in re.split(r'[;,\s]+',r.get('Alias_Names',''))))
        if len(hits)==1:
            target=next(iter(hits))
            if fid in equivalence and equivalence[fid]!=target:raise ValueError('Conflicting protein aliases')
            equivalence[fid]=target
            identity_rows.append((fid,target,r['Uniprot_ID'],r.get('Alias_Names',''),'explicit_protein_alias_to_unique_stable_primary_symbol'))
    def canonical(fid):return equivalence.get(fid,fid)
    all_identity=[(f,canonical(f)) for f in factors]
    create_rows(db,'family_sharing_factor_identity','original_factor_id VARCHAR,canonical_factor_id VARCHAR',all_identity)
    create_rows(db,'family_sharing_identity_evidence','original_factor_id VARCHAR,canonical_factor_id VARCHAR,uniprot_id VARCHAR,source_aliases VARCHAR,mapping_method VARCHAR',sorted(set(identity_rows)))
    supported_raw={r[0] for r in db.execute('SELECT factor_id FROM factor_catalog WHERE rbp2go_listed OR other_rbp_source').fetchall()}
    supported={canonical(f) for f in supported_raw}
    splicing_raw={r[0] for r in db.execute('SELECT factor_id FROM candidate_splicing_factors').fetchall()}
    splicing={canonical(f) for f in splicing_raw}
    existing=db.execute("""SELECT DISTINCT v.factor_id,v.source,v.source_motif_id,v.motif_id,v.evidence,v.experiment_type,v.pmid
      FROM valid_factor_motifs v WHERE source IN ('CISBP-RNA','ATtRACT','oRNAmentV2')
      AND evidence!='mutant_binding_exclude_from_WT_default'""").fetchall()
    donor_links=collections.defaultdict(set)
    for f,*r in existing:donor_links[canonical(f)].add(tuple(r))
    covered=set(donor_links)
    metadata={r['ID']:r for r in source_rows(ATLAS/'raw/rbp2go/download_interpro_data.export')}
    member_info=collections.defaultdict(lambda:{'family_ids':set(),'rbd_ids':set(),'uniprot_ids':set()})
    membership=set()
    for original,r in profiles:
        fid=canonical(original)
        if fid not in supported:continue
        info=member_info[fid]
        info['rbd_ids'].update(re.findall(r'IPR\d+',r.get('RBD_list','')))
        info['uniprot_ids'].add(r['Uniprot_ID'])
        for ip in r.get('InterPro','').split(';'):
            if metadata.get(ip,{}).get('Type')=='Family':
                info['family_ids'].add(ip)
                membership.add((fid,ip,metadata[ip]['Name'],r['Uniprot_ID'],'RBP2GO_InterPro_Family'))
    create_rows(db,'motif_family_membership','factor_id VARCHAR,family_id VARCHAR,family_name VARCHAR,uniprot_id VARCHAR,annotation_source VARCHAR',sorted(membership))
    families=collections.defaultdict(set)
    for f,info in member_info.items():
        for ip in info['family_ids']:families[ip].add(f)
    pair_rows=[];evidence=[]
    for family,members in sorted(families.items()):
        donors=members & covered;recipients=members-covered
        for recipient in sorted(recipients):
            for donor in sorted(donors):
                if recipient==donor:raise ValueError('Identity mismatch would create self-transfer')
                target_rbd=member_info[recipient]['rbd_ids'];donor_rbd=member_info[donor]['rbd_ids']
                common=sorted(target_rbd & donor_rbd)
                context=[]
                if not target_rbd:context.append('recipient_RBD_annotation_missing')
                elif not common:context.append('no_shared_annotated_RBD')
                if target_rbd!=donor_rbd:context.append('RBD_annotation_sets_differ')
                if 'YTH' in metadata[family]['Name']:context.append('RNA_modification_context_not_encoded_by_PWM')
                pair_rows.append((recipient,donor,family,metadata[family]['Name'],json.dumps(common),bool(common),json.dumps(context)))
                for source,source_mid,mid,donor_evidence,assay,pmid in sorted(donor_links[donor]):
                    evidence.append((recipient,mid,donor,source,source_mid,donor_evidence,assay,pmid,family,
                        metadata[family]['Name'],'family_shared_candidate',bool(common),json.dumps(context)))
    create_rows(db,'family_sharing_candidate_pairs',
        'recipient_factor_id VARCHAR,donor_factor_id VARCHAR,family_id VARCHAR,family_name VARCHAR,shared_RBD_ids_json VARCHAR,has_shared_annotated_RBD BOOLEAN,context_flags_json VARCHAR',pair_rows)
    create_rows(db,'family_shared_motif_evidence',
        'factor_id VARCHAR,motif_id VARCHAR,donor_factor_id VARCHAR,donor_source VARCHAR,donor_source_motif_id VARCHAR,donor_evidence VARCHAR,donor_assay VARCHAR,donor_pmid VARCHAR,family_id VARCHAR,family_name VARCHAR,extension_label VARCHAR,has_shared_annotated_RBD BOOLEAN,context_flags_json VARCHAR',evidence)
    db.execute("""CREATE OR REPLACE TABLE family_shared_motifs AS SELECT factor_id,motif_id,
        'family_shared_candidate'::VARCHAR extension_label,
        bool_or(has_shared_annotated_RBD) has_shared_annotated_RBD,
        list(DISTINCT donor_factor_id ORDER BY donor_factor_id) donor_factor_ids,
        list(DISTINCT family_id ORDER BY family_id) family_ids,
        count(*) provenance_links_not_independent_replications
        FROM family_shared_motif_evidence GROUP BY factor_id,motif_id""")
    db.execute("""CREATE OR REPLACE VIEW effective_factor_motifs AS
      SELECT DISTINCT i.canonical_factor_id factor_id,v.source,v.source_motif_id,v.motif_id,v.evidence,
          false is_family_extension,'source_association'::VARCHAR extension_label
      FROM valid_factor_motifs v JOIN family_sharing_factor_identity i ON v.factor_id=i.original_factor_id
      UNION ALL SELECT factor_id,'FamilySharing'::VARCHAR source,motif_id source_motif_id,motif_id,
          extension_label evidence,true is_family_extension,extension_label FROM family_shared_motifs""")
    # Link compatible matrices to the existing oRNAment sites, without materializing copies.
    # Matrices without an exact oRNAment counterpart explicitly need a new sequence scan.
    db.execute("""CREATE OR REPLACE TABLE family_shared_site_motif_mapping AS
      SELECT DISTINCT f.factor_id,f.motif_id,m.source_motif_id ornament_source_motif_id
      FROM family_shared_motifs f JOIN motif_sources m USING(motif_id)
      WHERE m.source='oRNAmentV2' AND m.status='valid'""")
    db.execute("""CREATE OR REPLACE VIEW family_shared_motif_sites AS SELECT s.*,m.factor_id,
      m.motif_id AS shared_matrix_id,true AS is_family_extension,'family_shared_candidate'::VARCHAR extension_label
      FROM motif_sites s JOIN family_shared_site_motif_mapping m ON s.source_motif_id=m.ornament_source_motif_id""")
    create_rows(db,'family_sharing_policy','key VARCHAR,value VARCHAR',[
      ('version','1.0'),('scope','mouse recipient and mouse motif donor; no new cross-species transfer'),
      ('family_definition','exact same InterPro annotation with Type=Family; generic Domain / Homologous_superfamily annotations excluded'),
      ('recipient_definition','source-supported mouse RBP without any existing valid curated motif after explicit alias resolution'),
      ('donor_definition','existing CISBP-RNA / ATtRACT / oRNAmentV2 association, including source-inferred; exclude mutant-only and ENCORI de novo candidates'),
      ('extension_label','family_shared_candidate'),('propagation','one hop from pre-existing source associations only; no inference cascades'),
      ('interpretation','exploratory same-family motif prior; no proof of conserved sequence preference, occupancy, binding strength or regulatory sign'),
      ('RBD_flags','shared RBD annotations and modification context are diagnostic labels, not a validated transfer threshold'),
      ('original_data','base factor_motifs and base genomic site data unchanged; extension available via separate tables and effective view'),
    ])
    recipients={r[0] for r in pair_rows};with_rbd={r[0] for r in pair_rows if r[5]}
    identity_recovered={f for f in supported_raw if f not in {x[0] for x in existing} and canonical(f) in covered}
    stats={
      'status':'complete','updated_at':time.strftime('%Y-%m-%dT%H:%M:%S%z'),
      'source_supported_catalog_rows':len(supported_raw),'canonical_source_supported_RBPs':len(supported),
      'original_curated_motif_factor_rows':len({x[0] for x in existing}),
      'canonical_existing_motif_factors':len(covered),
      'canonical_existing_motif_source_supported_RBPs':len(covered & supported),
      'identity_resolved_catalog_rows':len(equivalence),
      'identity_recovered_motif_catalog_rows':len(identity_recovered),
      'new_family_shared_RBPs':len(recipients),'new_family_shared_RBPs_with_shared_RBD_annotation':len(with_rbd),
      'canonical_motif_RBPs_after_extension':len((covered|recipients)&supported),
      'raw_splicing_candidates':len(splicing_raw),'canonical_splicing_candidates':len(splicing),
      'splicing_original_covered_rows':len(splicing_raw & {x[0] for x in existing}),
      'splicing_covered_after_identity_resolution':len(splicing&covered),
      'splicing_new_family_shared':len(splicing&recipients),
      'splicing_new_family_shared_with_shared_RBD_annotation':len(splicing&with_rbd),
      'new_family_groups':len({r[2] for r in pair_rows}),
      'recipient_donor_family_pairs':len(pair_rows),
      'provenance_links_not_independent_replications':len(evidence),
      'distinct_recipient_matrix_links':db.execute('SELECT count(*) FROM family_shared_motifs').fetchone()[0],
      'new_recipients_with_existing_ornament_matrix_sites':db.execute('SELECT count(DISTINCT factor_id) FROM family_shared_site_motif_mapping').fetchone()[0],
      'extension_matrices_requiring_new_scan':db.execute('SELECT count(DISTINCT motif_id) FROM family_shared_motifs WHERE motif_id NOT IN (SELECT motif_id FROM motif_sources WHERE source=\'oRNAmentV2\' AND status=\'valid\')').fetchone()[0],
    }
    detail=[]
    for recipient in sorted(recipients,key=lambda f:factors[f][1].lower()):
        pairs=[p for p in pair_rows if p[0]==recipient]
        detail.append({'factor_id':recipient,'gene_symbol':factors[recipient][1],
            'splicing_candidate':recipient in splicing,'extension_label':'family_shared_candidate',
            'families':sorted({p[3] for p in pairs}),'family_ids':sorted({p[2] for p in pairs}),
            'donor_symbols':sorted({factors[p[1]][1] for p in pairs}),
            'has_shared_annotated_RBD':recipient in with_rbd,
            'context_flags':sorted({flag for p in pairs for flag in json.loads(p[6])}),
            'matrix_count':db.execute('SELECT count(*) FROM family_shared_motifs WHERE factor_id=?',[recipient]).fetchone()[0]})
    report={'statistics':stats,'recipients':detail,'identity_recovered':[
        {'original_factor_id':f,'canonical_factor_id':canonical(f)} for f in sorted(identity_recovered)]}
    checks={
      'no_self_transfer':all(r[0]!=r[2] for r in evidence),
      'only_previous_motifless_recipients':not bool(recipients & covered),
      'all_donors_preexist':all(r[2] in covered for r in evidence),
      'all_recipients_supported':recipients <= supported,
      'only_exact_annotated_families':all(r[8] in member_info[r[0]]['family_ids'] & member_info[r[2]]['family_ids'] for r in evidence),
      'all_extension_labels_present':all(r[10]=='family_shared_candidate' for r in evidence),
      'no_dangling_matrices':db.execute('SELECT count(*) FROM family_shared_motifs f LEFT JOIN motifs m USING(motif_id) WHERE m.motif_id IS NULL').fetchone()[0]==0,
      'no_duplicate_recipient_matrix_links':db.execute('SELECT count(*)-count(DISTINCT (factor_id,motif_id)) FROM family_shared_motifs').fetchone()[0]==0,
    }
    report['checks']=checks;report['validation_passed']=all(checks.values())
    if not report['validation_passed']:raise ValueError(checks)
    db.execute('CHECKPOINT');db.close()
    (ROOT/'family_sharing_report.json').write_text(json.dumps(report,indent=2,ensure_ascii=False))
    with (ROOT/'family_shared_recipients.tsv').open('w') as f:
        w=csv.writer(f,delimiter='\t');w.writerow(detail[0].keys() if detail else [])
        for r in detail:w.writerow(json.dumps(v,ensure_ascii=False) if isinstance(v,list) else v for v in r.values())
    # Refresh the changed DB checksum; unchanged Parquet hashes are retained.
    sums=ROOT/'SHA256SUMS';lines=sums.read_text().splitlines()
    lines=[sha256(DBPATH)+'  mouse_rbp.duckdb' if l.endswith('  mouse_rbp.duckdb') else l for l in lines]
    sums.write_text('\n'.join(lines)+'\n')
    print(json.dumps(report,indent=2,ensure_ascii=False),flush=True)


if __name__=='__main__':main()
