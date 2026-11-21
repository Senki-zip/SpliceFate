"""Read-only mouse evidence queries. Run with .venv-rbp-db/bin/python."""
from __future__ import annotations
import argparse
import json
from pathlib import Path
import re
import duckdb

ROOT=Path(__file__).resolve().parents[1]/'data/external/mouse_rbp_unified_20261116'


def connect():
    db=duckdb.connect(str(ROOT/'mouse_rbp.duckdb'),read_only=True)
    db.execute("SET file_search_path='"+str(ROOT).replace("'","''")+"'")
    db.execute('SET threads=4')
    db.execute("SET memory_limit='8GB'")
    return db


def records(db,sql,params=None):
    cur=db.execute(sql,params or [])
    names=[r[0] for r in cur.description]
    return [dict(zip(names,row)) for row in cur.fetchall()]


def interaction_records(db,factor_ids,limit):
    """Source labels and unknown fields stay attached to each returned relation."""
    if not db.execute("SELECT count(*) FROM information_schema.tables WHERE table_name='protein_interaction_evidence'").fetchone()[0]:
        return {'status':'interaction_extension_not_installed'}
    exists=db.execute("SELECT count(*) FROM information_schema.tables WHERE table_name='family_sharing_factor_identity'").fetchone()[0]
    canonical_ids=factor_ids
    if exists and factor_ids:
        canonical_ids=[r[0] for r in db.execute('SELECT DISTINCT canonical_factor_id FROM family_sharing_factor_identity WHERE original_factor_id IN (SELECT unnest(?))',[factor_ids]).fetchall()]
    where='(e.factor_a IN (SELECT unnest(?)) OR e.factor_b IN (SELECT unnest(?)))'
    params=[canonical_ids,canonical_ids]
    result={'canonical_factor_ids':canonical_ids,
        'interpretation':'Protein relations and complex membership; no inherited motif, RNA targets or KO sign. Support rows are not independent replication. Mouse molecular origin does not establish mouse experimental host or target cell context.',
        'coverage':records(db,'SELECT * FROM factor_interaction_coverage WHERE factor_id IN (SELECT unnest(?))',[canonical_ids])}
    result['evidence_summary']=records(db,'SELECT source,relation_type,count(*) source_records,count(*) FILTER(WHERE construct_change_flag OR perturbation_keyword_flag) flagged_construct_or_condition_records FROM protein_interaction_evidence e WHERE '+where+' GROUP BY source,relation_type',params)
    result['pair_summary']=records(db,'SELECT count(*) distinct_canonical_pairs,count(*) FILTER(WHERE has_direct_physical_annotation) pairs_with_direct_contact_annotation FROM protein_interaction_pairs e WHERE '+where,params)[0]
    result['pairs']=records(db,"""SELECT p.*,fa.gene_symbol symbol_a,fb.gene_symbol symbol_b FROM protein_interaction_pairs p
       JOIN factors fa ON fa.factor_id=p.factor_a JOIN factors fb ON fb.factor_id=p.factor_b
       WHERE (p.factor_a IN (SELECT unnest(?)) OR p.factor_b IN (SELECT unnest(?)))
       ORDER BY has_direct_physical_annotation DESC,has_IntAct_association DESC,has_STRING_nontransferred_experiment_channel DESC,p.factor_a,p.factor_b LIMIT ?""",params+[limit])
    result['evidence_examples']=records(db,"""SELECT e.evidence_id,e.source,e.factor_a,e.factor_b,
       a.gene_symbol symbol_a,b.gene_symbol symbol_b,e.source_id_a,e.source_id_b,e.relation_type,
       e.detection_method,e.source_interaction_id,e.publication,e.confidence,e.host_organism,
       e.expansion_method,e.construct_change_flag,e.perturbation_keyword_flag,e.rna_dependence_text_flag,
       e.rna_dependence,e.wt_review_status,e.regulatory_sign,e.nontransferred_experiment_score,e.transferred_experiment_score,e.combined_score
       FROM protein_interaction_evidence e LEFT JOIN factors a ON a.factor_id=e.factor_a
       LEFT JOIN factors b ON b.factor_id=e.factor_b WHERE """+where+""" ORDER BY e.source='IntAct' DESC,
       e.relation_type='direct_physical' DESC,e.evidence_id LIMIT ?""",params+[limit])
    result['complex_membership_count']=records(db,"SELECT count(*) membership_rows,count(DISTINCT (source,complex_id)) complexes FROM mouse_complex_memberships WHERE factor_id IN (SELECT unnest(?))",[canonical_ids])[0]
    result['complex_memberships']=records(db,"SELECT * FROM mouse_complex_memberships WHERE factor_id IN (SELECT unnest(?)) ORDER BY source,complex_id LIMIT ?",[canonical_ids,limit])
    result['native_mouse_stage_memberships']=records(db,"SELECT * FROM mouse_spliceosome_stage_memberships WHERE factor_id IN (SELECT unnest(?)) LIMIT ?",[canonical_ids,limit])
    result['stage_note']='Human stage memberships retained in source species; no automatic human-to-mouse projection.'
    result['returned_row_limit_per_list']=limit
    return result


def protein_records(db,factor_ids,limit):
    if not db.execute("SELECT count(*) FROM information_schema.tables WHERE table_name='protein_entries'").fetchone()[0]:
        return {'status':'protein_dimension_extension_not_installed'}
    ids=[r[0] for r in db.execute('SELECT DISTINCT canonical_factor_id FROM family_sharing_factor_identity WHERE original_factor_id IN (SELECT unnest(?))',[factor_ids]).fetchall()]
    params=[ids,limit]
    result={'coordinate_convention':'protein and genomic intervals are 0-based half-open',
        'interpretation':'Domain matches and reference isoform consequences are structural evidence. Overlapping source models are not independent domains; no numeric activity weights or observed-chain completion inferred.',
        'collection_policy':records(db,"SELECT * FROM protein_dimension_policy WHERE key IN ('collection_status','created_at','genomics','isoforms')"),
        'coverage':records(db,'SELECT * FROM protein_dimension_coverage WHERE factor_id IN (SELECT unnest(?))',[ids]),
        'proteins':records(db,"SELECT m.*,p.entry_name,p.entry_type,p.length_aa,p.sequence_version,p.last_annotation_update FROM protein_factor_mapping m JOIN protein_entries p USING(accession) WHERE factor_id IN (SELECT unnest(?)) ORDER BY p.accession LIMIT ?",params),
        'domain_memberships_without_positions':records(db,"SELECT a.*,c.name,c.entry_type FROM factor_domain_annotations a LEFT JOIN domain_catalog c ON c.domain_id=a.domain_id AND c.source='RBP2GO_InterPro' WHERE a.factor_id IN (SELECT unnest(?)) ORDER BY a.uniprot_id,a.domain_id LIMIT ?",params),
        'located_domains':records(db,'SELECT * FROM factor_protein_domains WHERE factor_id IN (SELECT unnest(?)) ORDER BY accession,sequence_id,start0,source LIMIT ?',params),
        'isoforms':records(db,'SELECT * EXCLUDE(record_json) FROM factor_protein_isoforms WHERE factor_id IN (SELECT unnest(?)) ORDER BY accession,isoform_id LIMIT ?',params),
        'functional_regions':records(db,"SELECT * EXCLUDE(record_json,source_location_json) FROM factor_protein_features WHERE factor_id IN (SELECT unnest(?)) AND feature_type IN ('Region','Motif','Active site','Binding site','Site','Transmembrane','Signal','Transit peptide') ORDER BY accession,start0 LIMIT ?",params),
        'isoform_domain_consequences':records(db,"SELECT m.factor_id,c.*,f.description,f.feature_type FROM protein_isoform_feature_consequences c JOIN protein_features f USING(feature_key) JOIN protein_factor_mapping m ON m.accession=c.accession WHERE m.factor_id IN (SELECT unnest(?)) AND m.status='resolved' ORDER BY c.isoform_id,f.start0 LIMIT ?",params),
        'reference_transcripts':records(db,"SELECT DISTINCT t.* FROM protein_local_transcript_matches t JOIN protein_factor_mapping m USING(accession) WHERE m.factor_id IN (SELECT unnest(?)) AND m.status='resolved' ORDER BY t.assembly,t.transcript_id LIMIT ?",params),
        'genomic_domain_blocks':records(db,"SELECT b.* FROM protein_domain_genomic_blocks b JOIN protein_sequences s USING(sequence_id) JOIN protein_factor_mapping m USING(accession) WHERE m.factor_id IN (SELECT unnest(?)) AND m.status='resolved' ORDER BY b.assembly,b.transcript_id,b.match_id,b.fragment_index,b.block_index LIMIT ?",params)}
    result['returned_limit_per_list']=limit
    return result


def main():
    p=argparse.ArgumentParser()
    p.add_argument('--factor',help='Mouse factor symbol, stable gene ID or source protein ID')
    p.add_argument('--target-gene',help='Target RNA gene symbol or stable Ensembl ID')
    p.add_argument('--region',help='0-based half-open genomic interval, e.g. chr3:100000-101000')
    p.add_argument('--assembly',choices=['mm10','mm39'],default='mm39')
    p.add_argument('--limit',type=int,default=20)
    p.add_argument('--include-family-shared',action='store_true',help='Include labeled same-family motif candidates in site queries')
    p.add_argument('--interactions-only',action='store_true',help='Query protein interactions/complexes for --factor without scanning motif/CLIP sites')
    p.add_argument('--protein-only',action='store_true',help='Query domains, functional sites, reference isoforms and CDS mappings for --factor')
    p.add_argument('--sql',help='Read-only SQL; table paths resolve relative to the database root')
    a=p.parse_args()
    if not (a.factor or a.target_gene or a.region or a.sql):p.error('Specify --factor, --target-gene, --region or --sql')
    if not 1<=a.limit<=10000:p.error('--limit must be between 1 and 10000')
    if a.interactions_only and (not a.factor or a.target_gene or a.region or a.sql):p.error('--interactions-only requires --factor and cannot be combined with --region, --target-gene or --sql')
    if a.protein_only and (not a.factor or a.target_gene or a.region or a.sql or a.interactions_only):p.error('--protein-only requires --factor and cannot be combined with other query modes')
    db=connect()
    if a.sql:
        print(json.dumps(records(db,a.sql),ensure_ascii=False,indent=2,default=str));return
    result={'assembly':a.assembly,'coordinate_convention':'0-based half-open; each split motif block tested independently',
            'interpretation':'Predicted motif sites and experimental binding evidence; no signed KO effect. Counts across transcripts, databases and shared motifs are not independent replicates.'}
    factor_ids=[]
    motif_ids=[]
    if a.factor:
        factors=records(db,"""SELECT DISTINCT f.* FROM factor_catalog f LEFT JOIN factor_sources s USING(factor_id)
             WHERE upper(f.gene_symbol)=upper(?) OR f.gene_id=? OR upper(s.gene_symbol)=upper(?) OR s.source_id=?
             OR EXISTS(SELECT 1 FROM factor_aliases al WHERE al.factor_id=f.factor_id AND upper(al.alias)=upper(?))""",
             [a.factor,a.factor,a.factor,a.factor,a.factor])
        result['factors']=factors
        factor_ids=[f['factor_id'] for f in factors]
        extension_available=db.execute("SELECT count(*) FROM information_schema.tables WHERE table_name='family_shared_motifs'").fetchone()[0]>0
        if extension_available and factor_ids:
            mapped=db.execute('SELECT DISTINCT canonical_factor_id FROM family_sharing_factor_identity WHERE original_factor_id IN (SELECT unnest(?))',[factor_ids]).fetchall()
            canonical_ids=[r[0] for r in mapped]
            result['family_sharing_canonical_factor_ids']=canonical_ids
            result['family_shared_motifs']=records(db,'SELECT * FROM family_shared_motifs WHERE factor_id IN (SELECT unnest(?))',[canonical_ids])
            result['family_shared_site_queries_enabled']=a.include_family_shared
            result['family_shared_note']='Exploratory same-family candidates; shared domains do not establish identical specificity. Source links and context flags are in family_shared_motif_evidence. Matrices without exact oRNAment counterparts require a new sequence scan.'
        else:
            canonical_ids=factor_ids
        if a.protein_only:
            result['protein_annotations']=protein_records(db,factor_ids,a.limit)
            print(json.dumps(result,ensure_ascii=False,indent=2,default=str));return
        result['protein_interactions']=interaction_records(db,factor_ids,a.limit)
        if a.interactions_only:
            print(json.dumps(result,ensure_ascii=False,indent=2,default=str));return
        if not factor_ids:
            result['message']='No unambiguous or unresolved-source entry matched this identifier.'
        else:
            result['motif_links']=records(db,"""SELECT f.gene_symbol,v.* FROM valid_factor_motifs v JOIN factors f USING(factor_id)
                WHERE v.factor_id IN (SELECT unnest(?)) AND v.source!='ENCORI' ORDER BY v.source,v.source_motif_id""",[factor_ids])
            result['motif_link_summary']=records(db,"""SELECT source,evidence,count(*) links,count(DISTINCT motif_id) exact_matrices
                FROM valid_factor_motifs WHERE factor_id IN (SELECT unnest(?)) GROUP BY source,evidence""",[factor_ids])
            result['clip_experiments']=records(db,"""SELECT experiment_id,rbp_label,assay,cell_tissue,treatment,source_accession,
                perturbation_flag,peak_rows,content_duplicate_of FROM experiments WHERE factor_id IN (SELECT unnest(?))""",[factor_ids])
            result['target_gene_summary']=records(db,"""SELECT source,count(DISTINCT target_gene_id) source_target_ids,
                count(DISTINCT canonical_gene_id) mapped_ensembl_targets FROM resolved_gene_binding
                WHERE factor_id IN (SELECT unnest(?)) GROUP BY source""",[factor_ids])
            motif_ids=[r[0] for r in db.execute("SELECT DISTINCT source_motif_id FROM factor_motifs WHERE source='oRNAmentV2' AND factor_id IN (SELECT unnest(?))",[factor_ids]).fetchall()]
            if extension_available and a.include_family_shared:
                motif_ids=sorted(set(motif_ids)|{r[0] for r in db.execute('SELECT DISTINCT ornament_source_motif_id FROM family_shared_site_motif_mapping WHERE factor_id IN (SELECT unnest(?))',[canonical_ids]).fetchall()})
    target_ids=[]
    if a.target_gene:
        result['target_identities']=records(db,"SELECT DISTINCT * FROM gene_identities WHERE upper(gene_symbol)=upper(?) OR gene_id=?",[a.target_gene,a.target_gene])
        target_ids=sorted({r['gene_id'] for r in result['target_identities']})
        result['gene_binding_summary']=records(db,"""SELECT source,count(*) source_annotations,count(DISTINCT factor_id) factors
            FROM resolved_gene_binding WHERE canonical_gene_id IN (SELECT unnest(?)) GROUP BY source""",[target_ids])
    motif_filters=[];motif_params=[];peak_filters=['p.assembly=?'];peak_params=[a.assembly]
    if a.factor:
        motif_filters+=['s.source_motif_id IN (SELECT unnest(?))','s.species_specific'];motif_params.append(motif_ids)
        peak_filters.append('e.factor_id IN (SELECT unnest(?))');peak_params.append(factor_ids)
    if a.target_gene:
        motif_filters.append('s.gene_id IN (SELECT unnest(?))');motif_params.append(target_ids)
        target_chroms=[r[0] for r in db.execute('SELECT DISTINCT chromosome FROM motif_gene_loci WHERE gene_id IN (SELECT unnest(?))',[target_ids]).fetchall()]
        motif_filters.append('s.chromosome IN (SELECT unnest(?))');motif_params.append(target_chroms)
    if a.region:
        match=re.fullmatch(r'([^:]+):(\d+)-(\d+)',a.region)
        if not match or int(match[2])>=int(match[3]):p.error('Invalid --region; use chromosome:start-end with start<end')
        chrom,start,end=match[1],int(match[2]),int(match[3])
        canonical=db.execute('SELECT DISTINCT chromosome FROM chromosome_aliases WHERE assembly=? AND source_chromosome=?',[a.assembly,chrom]).fetchall()
        if len(canonical)==1:chrom=canonical[0][0]
        motif_filters+=['s.chromosome=?','s.start0<? AND s.end0>?',
            'len(list_filter(list_zip(s.block_starts0,s.block_ends0),lambda b: b[1] < ? AND b[2] > ?))>0']
        motif_params += [chrom,end,start,end,start]
        peak_filters+=['p.chromosome=?','p.start0<? AND p.end0>?'];peak_params += [chrom,end,start]
        result['interval']={'chromosome':chrom,'start0':start,'end0':end}
    if a.assembly=='mm39':
        where=' AND '.join(motif_filters) or 'true'
        result['motif_site_count']=records(db,'SELECT count(*) transcript_feature_site_rows FROM motif_sites s WHERE '+where,motif_params)[0]
        result['motif_sites']=records(db,"""SELECT s.site_id,s.gene_id,s.transcript_id,s.feature_id,s.source_motif_id,
           s.chromosome,s.strand,s.block_starts0,s.block_ends0,s.region,s.score,s.percentile,s.species_specific,s.CLIP_support
           FROM motif_sites s WHERE """+where+' LIMIT ?',motif_params+[a.limit])
        if a.factor:
            original_mids={r[0] for r in db.execute("SELECT DISTINCT source_motif_id FROM factor_motifs WHERE source='oRNAmentV2' AND factor_id IN (SELECT unnest(?))",[factor_ids]).fetchall()}
            for row in result['motif_sites']:
                row['is_family_extension']=row['source_motif_id'] not in original_mids
                row['extension_label']='family_shared_candidate' if row['is_family_extension'] else 'source_association'
                if row['is_family_extension']:
                    row['CLIP_support_scope']='original_source_motif_annotation; not evidence for the recipient RBP'
                    row['recipient_specific_CLIP_support']=None
    else:
        result['motif_site_note']='oRNAment V2 sites are stored in native mm39 only; mm10 was not silently substituted.'
    if a.region or a.factor:
        result['clip_peaks']=records(db,"""SELECT p.*,e.factor_id,e.rbp_label,e.cell_tissue,e.assay,e.source_accession,e.perturbation_flag
          FROM clip_peaks p JOIN experiments e USING(experiment_id) WHERE """+' AND '.join(peak_filters)+' LIMIT ?',peak_params+[a.limit])
    result['returned_row_limit']=a.limit
    print(json.dumps(result,ensure_ascii=False,indent=2,default=str))


if __name__=='__main__':main()
