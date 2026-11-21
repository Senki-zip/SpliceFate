"""Add typed protein/isoform/domain/CDS evidence to the existing mouse database.

No inferred KO signs or numeric activity weights are added. Source sequence
identity is required before coordinates are transferred between resources.
"""
from __future__ import annotations
import argparse
import collections
import csv
import gzip
import hashlib
import json
from pathlib import Path
import re
import shutil
import time

import duckdb
import numpy as np
import pysam
from Bio.Seq import Seq

from build_mouse_rbp_database import ROOT as MAIN, DBPATH, REFERENCES, create_rows, sha256, sqlstr

ROOT = MAIN.parent / 'rbp_protein_dimensions_20261120'


def js(value):
    return json.dumps(value, ensure_ascii=False, separators=(',', ':'))


def key(*parts):
    return hashlib.sha256(js(parts).encode()).hexdigest()


def stable(value):
    return value.split('.')[0] if value else ''


def protein_interval(location):
    """Keep uncertain endpoints distinct from exact amino-acid coordinates."""
    start, end = location.get('start', {}), location.get('end', {})
    lo, hi = start.get('value'), end.get('value')
    exact = (start.get('modifier') == 'EXACT' and end.get('modifier') == 'EXACT'
             and isinstance(lo, int) and isinstance(hi, int) and 1 <= lo <= hi)
    return (lo-1 if isinstance(lo, int) else None, hi,
            'exact' if exact else 'uncertain_or_unknown', location.get('sequence'))


def isoform_edit_map(canonical, isoform, edits):
    """Apply only explicit VAR_SEQ edits; verify the complete resulting sequence."""
    cursor = outpos = 0; unchanged = []; pieces = []
    for edit in sorted(edits, key=lambda x: x['start0']):
        a, b = edit['start0'], edit['end0']
        if (edit['status'] != 'exact' or a is None or b is None or a < cursor or b > len(canonical)
                or edit.get('ambiguous_alternative')):
            return None, 'unresolved_edit_coordinates'
        if edit.get('original') and canonical[a:b] != edit['original']:
            return None, 'original_sequence_mismatch'
        if a > cursor:
            piece = canonical[cursor:a]; pieces.append(piece)
            unchanged.append((cursor, a, outpos)); outpos += len(piece)
        pieces.append(edit['replacement']); outpos += len(edit['replacement']); cursor = b
    if cursor < len(canonical):
        pieces.append(canonical[cursor:]); unchanged.append((cursor, len(canonical), outpos))
    if ''.join(pieces) != isoform:
        return None, 'isoform_sequence_mismatch'
    return unchanged, 'sequence_verified'


def feature_consequence(start, end, edits, segments):
    overlapping = [e for e in edits if e['start0'] < end and e['end0'] > start]
    if overlapping:
        if any(e['start0'] <= start and e['end0'] >= end and e['replacement'] == '' for e in overlapping):
            return 'deleted', None, None
        return 'overlaps_splice_alteration', None, None
    for a, b, offset in segments:
        if a <= start and end <= b:
            return 'retained_unchanged', offset + start-a, offset + end-a
    return 'unresolved', None, None


def genomic_blocks(cds_segments, start0, end0):
    """Project a verified protein interval through spliced CDS offsets."""
    left, right = start0*3, end0*3; blocks = []
    for segment in cds_segments:
        a, b, chrom, gs, ge, strand, exon = segment
        lo, hi = max(left, a), min(right, b)
        if lo >= hi:
            continue
        if strand == '+':
            x, y = gs+lo-a, gs+hi-a
        else:
            x, y = ge-(hi-a), ge-(lo-a)
        blocks.append((chrom, x, y, strand, exon, lo, hi))
    return blocks


def parse_fasta(path):
    name = None; chunks = []
    with path.open() as stream:
        for line in stream:
            if line.startswith('>'):
                if name is not None:
                    yield name, ''.join(chunks)
                token = line[1:].split()[0]; name = token.split('|')[1] if '|' in token else token
                chunks = []
            else:
                chunks.append(line.strip())
    if name is not None:
        yield name, ''.join(chunks)


def verified_sources(allow_partial):
    state = json.loads((ROOT/'download_status.json').read_text())
    if not allow_partial and state['status'] != 'complete':
        raise RuntimeError('Source collection is incomplete; use --allow-partial only for a labeled incremental import')
    manifest = json.loads((ROOT/'download_manifest.json').read_text())
    good = []
    for row in manifest:
        if row['status'] not in ('verified', 'not_in_source'):
            continue
        p = ROOT/row['path']
        if p.stat().st_size != row['bytes'] or sha256(p) != row['sha256']:
            raise ValueError('Source checksum mismatch: '+str(p))
        good.append(row)
    return state, good


def integrate(allow_partial=False):
    state, manifest = verified_sources(allow_partial)
    paths = {r['path']: r for r in manifest}
    entries = {}
    for row in manifest:
        if row['source'] == 'UniProtKB' and row['path'].endswith('.json'):
            for entry in json.loads((ROOT/row['path']).read_text()).get('results', []):
                if entry['organism']['taxonId'] != 10090:
                    raise ValueError('Only mouse protein entries may be imported')
                uid = entry['primaryAccession']
                if uid in entries and entries[uid]['sequence']['value'] != entry['sequence']['value']:
                    raise ValueError('Conflicting UniProt releases: '+uid)
                entries[uid] = entry
    if not entries:
        raise ValueError('No verified UniProt entries')
    fasta = {}
    for row in manifest:
        if row['source'] == 'UniProtKB' and row['path'].endswith('.fasta'):
            for uid, sequence in parse_fasta(ROOT/row['path']):
                if uid in fasta and fasta[uid] != sequence:
                    raise ValueError('Conflicting FASTA protein sequence: '+uid)
                fasta[uid] = sequence
    backup = MAIN/'backups/mouse_rbp.before_protein_dimensions_20261120.duckdb'
    report = {'started_at':time.strftime('%Y-%m-%dT%H:%M:%S%z'), 'source_status':state['status'],
              'scope':'source-supported mouse RBP proteins plus GO splicing candidates; existing background domain records also normalized',
              'partial_import':state['status'] != 'complete', 'status':'running', 'input_files':len(manifest)}
    if not backup.exists():
        shutil.copy2(DBPATH, backup)
        if sha256(backup) != sha256(DBPATH):
            raise ValueError('Backup checksum mismatch')
    report['backup'] = str(backup); report['backup_sha256'] = sha256(backup)
    (ROOT/'integration_status.json').write_text(json.dumps(report,indent=2))
    db = duckdb.connect(str(DBPATH)); db.execute('SET threads=4'); db.execute("SET memory_limit='12GB'")
    db.execute('SET file_search_path='+sqlstr(MAIN))
    protected = ['factors','factor_sources','motifs','factor_motifs','family_shared_motifs',
                 'experiments','clip_peaks_original','protein_interaction_evidence','complex_memberships']
    before = {t:db.execute(f'SELECT count(*),bit_xor(hash(t)) FROM {t} t').fetchone() for t in protected}
    canonical = dict(db.execute('SELECT * FROM family_sharing_factor_identity').fetchall())
    factors = {f:(g,s) for f,g,s in db.execute('SELECT factor_id,gene_id,gene_symbol FROM factors').fetchall()}
    gene_factors = collections.defaultdict(set); old_proteins = collections.defaultdict(set)
    profiles = [(canonical.get(f,f), u, role, json.loads(raw)) for f,u,role,raw in db.execute('SELECT * FROM rbp2go_profiles').fetchall()]
    for f,(g,s) in factors.items():
        if g.startswith('ENSMUSG'):
            gene_factors[stable(g)].add(canonical.get(f,f))
    for f,u,_,_ in profiles:
        old_proteins[u].add(f)
    # Normalize existing InterPro catalog and protein memberships, with no invented positions.
    catalog = {}
    for (raw,) in db.execute("SELECT record_json FROM rbp2go_support_records WHERE source_file='download_interpro_data.export'").fetchall():
        row=json.loads(raw); catalog[('RBP2GO_InterPro',row['ID'])] = (row['ID'],row['Name'],row['Type'].lower(),row.get('RNA_binding','').lower()=='true',raw)
    memberships = set(); protein_stats=[]
    for fid,u,role,row in profiles:
        rbd = {m[0]:int(m[1]) if m[1] else 1 for m in re.findall(r'(IPR\d+)(?:\s*\[(\d+)\])?',row.get('RBD_list',''))}
        for ip in sorted(set(re.findall(r'IPR\d+',row.get('InterPro','')))|set(rbd)):
            memberships.add((fid,u,ip,ip in rbd,rbd.get(ip), 'RBP2GO', 'membership_only_no_position'))
        def numeric(k):
            try:return float(row.get(k,''))
            except (ValueError,TypeError):return None
        protein_stats.append((fid,u,numeric('Length_AA'),numeric('RBDs_nb'),numeric('RBDs_content_fraction'),
                              numeric('IDR_content_fraction'),row.get('has_IDR'),row.get('LCD_status'),row.get('coiled_coil_status')))
    protein_rows=[]; mapping=[]; aliases={}; seq_rows={}; iso_rows=[]; edit_rows=[]; iso_checks={}
    alias_candidates=collections.defaultdict(set)
    for uid,entry in entries.items():
        for alias in [uid]+entry.get('secondaryAccessions',[]):alias_candidates[alias].add(uid)
    accession_resolution=[]
    for alias,hits in sorted(alias_candidates.items()):
        resolved=alias if alias in entries else next(iter(hits)) if len(hits)==1 else None
        if resolved:aliases[alias]=resolved
        accession_resolution.append((alias,resolved,'exact_primary' if alias in entries else 'unique_secondary' if resolved else 'ambiguous_secondary',js(sorted(hits))))
    features=[]; feature_ev=[]; comments=[]; xrefs=[]; transcript_links=[]
    uid_factors=collections.defaultdict(set); sequences_by_gene=collections.defaultdict(set)
    edits_by_id={}; iso_specs=[]
    for uid,entry in sorted(entries.items()):
        seq=entry['sequence']['value'];genes=entry.get('genes',[]);audit=entry.get('entryAudit',{})
        protein_rows.append((uid,entry.get('uniProtkbId'),entry['entryType'],10090,len(seq),sha256_text(seq),
                             entry.get('annotationScore'),entry.get('proteinExistence'),audit.get('sequenceVersion'),
                             audit.get('entryVersion'),audit.get('lastAnnotationUpdateDate'),js(entry.get('proteinDescription',{})),js(genes)))
        if uid in fasta and fasta[uid] != seq:raise ValueError('UniProt JSON/FASTA disagree: '+uid)
        seq_rows[uid]=(uid,uid,'canonical',seq,len(seq),sha256_text(seq),'UniProtKB_JSON',True)
        ensembl_genes=set()
        for x in entry.get('uniProtKBCrossReferences',[]):
            props={p['key']:p['value'] for p in x.get('properties',[])}
            xrefs.append((key(uid,x),uid,x['database'],x['id'],x.get('isoformId'),js(props)))
            if x['database']=='Ensembl' and props.get('GeneId','').startswith('ENSMUSG'):
                gid=stable(props['GeneId']);ensembl_genes.add(gid)
                transcript_links.append((key(uid,x),uid,x.get('isoformId'),x['id'],stable(x['id']),props.get('ProteinId'),stable(props.get('ProteinId','')),props['GeneId'],gid,'UniProtKB_Ensembl_xref'))
        old=set().union(*(old_proteins.get(a,set()) for a in [uid]+entry.get('secondaryAccessions',[]) if aliases.get(a)==uid))
        xhits=set().union(*(gene_factors.get(g,set()) for g in ensembl_genes))
        for fid in sorted(old|xhits):
            conflict = bool(ensembl_genes and factors[fid][0].startswith('ENSMUSG') and stable(factors[fid][0]) not in ensembl_genes)
            ambiguous = len(old|xhits)>1
            status='gene_id_conflict' if conflict else 'ambiguous_factor_mapping' if ambiguous else 'resolved'
            method='existing_UniProt_accession' if fid in old else 'exact_mouse_Ensembl_gene_id'
            mapping.append((uid,fid,factors[fid][0],method,status,js(sorted(ensembl_genes))))
            if status=='resolved':
                uid_factors[uid].add(fid)
                sequences_by_gene[stable(factors[fid][0])].add(uid)
        for n,comment in enumerate(entry.get('comments',[])):
            comments.append((key(uid,n,comment),uid,comment['commentType'],js(comment)))
            if comment['commentType']=='ALTERNATIVE PRODUCTS':
                for iso in comment.get('isoforms',[]):
                    for iid in iso.get('isoformIds',[]):
                        status=iso.get('isoformSequenceStatus','Unknown')
                        sequence=seq if status=='Displayed' else fasta.get(iid)
                        if sequence:
                            if iid in fasta and sequence!=fasta[iid]:raise ValueError('Displayed isoform sequence mismatch '+iid)
                            seq_rows[iid]=(iid,uid,'annotated_isoform',sequence,len(sequence),sha256_text(sequence),'UniProtKB_displayed' if status=='Displayed' else 'UniProtKB_isoform_FASTA',True)
                        iso_rows.append((iid,uid,iso.get('name',{}).get('value'),status,js(iso.get('sequenceIds',[])),iid in seq_rows,js(iso)))
                        iso_specs.append((uid,iid,status,iso.get('sequenceIds',[])))
        for n,feat in enumerate(entry.get('features',[])):
            fid=key(uid,n,feat);a,b,status,seqref=protein_interval(feat.get('location',{}))
            features.append((fid,uid,seqref or uid,feat['type'],feat.get('description',''),a,b,status,
                             feat.get('featureId'),js(feat.get('location',{})),js(feat)))
            for ev in feat.get('evidences',[]):
                feature_ev.append((fid,ev.get('evidenceCode'),ev.get('source'),ev.get('id'),js(ev)))
            if feat['type']=='Alternative sequence' and feat.get('featureId'):
                alt=feat.get('alternativeSequence',{});replacements=alt.get('alternativeSequences',[])
                edit={'start0':a,'end0':b,'status':status,'original':alt.get('originalSequence'),
                      'replacement':replacements[0] if len(replacements)==1 else '', 'ambiguous_alternative':len(replacements)>1,
                      'feature_id':fid,'source_feature_id':feat['featureId']}
                if seqref and seqref!=uid:edit['status']='noncanonical_edit_reference'
                edits_by_id[(uid,feat['featureId'])]=edit
    # Resolve isoform transformations only if complete source edits reproduce FASTA exactly.
    feature_impacts=[]; projected=[]; features_by_protein=collections.defaultdict(list)
    for feat in features:features_by_protein[feat[1]].append(feat)
    for uid,iid,status,ids in iso_specs:
        canonical_seq=seq_rows[uid][3]
        if status=='Displayed':iso_checks[iid]='canonical_displayed';continue
        if iid not in seq_rows:iso_checks[iid]='sequence_not_available';continue
        if not ids or any((uid,sid) not in edits_by_id for sid in ids):iso_checks[iid]='source_edits_unavailable';continue
        edits=[edits_by_id[(uid,sid)] for sid in ids]
        segments,check=isoform_edit_map(canonical_seq,seq_rows[iid][3],edits);iso_checks[iid]=check
        for edit in edits:
            edit_rows.append((iid,uid,edit['source_feature_id'],edit['feature_id'],edit['start0'],edit['end0'],edit['replacement'],check))
        if segments is None:continue
        for feat in features_by_protein[uid]:
            fid,puid,seqid,typ,description,a,b,certainty,*_ = feat
            if puid!=uid or seqid!=uid or certainty!='exact' or typ not in {'Domain','Region','Motif','Repeat','Zinc finger','Coiled coil','Active site','Binding site','Site','Signal','Transit peptide','Topological domain','Transmembrane'}:continue
            effect,ia,ib=feature_consequence(a,b,edits,segments)
            feature_impacts.append((iid,fid,uid,effect,ia,ib,'source_edit_and_full_sequence_verified','not_a_functional_activity_weight'))
            if effect=='retained_unchanged' and typ in {'Domain','Repeat','Zinc finger'}:
                projected.append((key('isoform_projection',iid,fid),iid,uid,'UniProtKB_isoform_projection',typ.lower(),
                                  feat[8] or fid,description,ia,ib,'exact_source_edit_sequence_verified',None,None,[(ia,ib,'CONTINUOUS')]))
    # InterPro includes families, superfamilies and domains. Never label all as domains.
    domain_matches=[]; domain_fragments=[]; interpro_audit=[]
    for relative,row in paths.items():
        if row['source']!='InterPro' or not relative.endswith('.entries.json'):continue
        requested=Path(relative).name.split('.')[0];uid=aliases.get(requested)
        if not uid:continue
        protein_path=f'raw/interpro/{requested}.protein.json';pmeta={}
        if protein_path in paths:
            pages=json.loads((ROOT/protein_path).read_text())['pages']
            pmeta=pages[0].get('metadata',{}) if pages else {}
        source_seq=pmeta.get('sequence')
        sequence_status='exact_sequence_match' if source_seq==seq_rows[uid][3] else 'sequence_mismatch' if source_seq else 'source_sequence_unavailable'
        if pmeta and str(pmeta.get('source_organism',{}).get('taxId'))!='10090':raise ValueError('Non-mouse InterPro sequence')
        interpro_audit.append((requested,uid,sequence_status,len(source_seq) if source_seq else None,
                               sha256_text(source_seq) if source_seq else None,pmeta.get('is_fragment'),js(pmeta)))
        for page in json.loads((ROOT/relative).read_text())['pages']:
            for ent in page.get('results',[]):
                md=ent['metadata'];source=md.get('source_database','unknown');acc=md['accession']
                name=md.get('name','');name=js(name) if not isinstance(name,str) else name
                rna_binding=any(g.get('identifier') in {'GO:0003723','GO:0003729','GO:0003727'} for g in md.get('go_terms') or [])
                catalog[('InterPro_API:'+source,acc)]=(acc,name,md.get('type','unknown'),rna_binding,js(md))
                for prot in ent.get('proteins',[]):
                    if prot['accession'].upper()!=requested.upper():raise ValueError('InterPro returned different protein')
                    if str(prot.get('organism'))!='10090':raise ValueError('Non-mouse InterPro matches')
                    for n,loc in enumerate(prot.get('entry_protein_locations',[])):
                        frags=[(int(f['start'])-1,int(f['end']),f.get('dc-status','unknown')) for f in loc.get('fragments',[])]
                        if not frags:continue
                        mid=key(requested,source,acc,n,loc);a=min(f[0] for f in frags);b=max(f[1] for f in frags)
                        record=(mid,uid,uid,'InterPro_API:'+source,md.get('type','unknown'),acc,name,a,b,
                                sequence_status,loc.get('model'),loc.get('score'),frags)
                        domain_matches.append(record)
    # UniProt positions share the canonical sequence and retain their own evidence.
    for feat in features:
        fid,uid,seqid,typ,name,a,b,status,source_id,*_=feat
        if typ not in {'Domain','Repeat','Zinc finger'}:continue
        if a is None or b is None:continue
        domain_matches.append((key('UniProtKB',fid),seqid,uid,'UniProtKB',typ.lower(),source_id or fid,name,a,b,
                               'exact_sequence_match' if seqid==uid and status=='exact' else 'uncertain_or_isoform_specific',None,None,[(a,b,'CONTINUOUS')]))
    domain_matches+=projected
    for record in domain_matches:
        domain_fragments.extend((record[0],n,a,b,dc) for n,(a,b,dc) in enumerate(record[12]))
    # Experimental binding peptides remain evidence, never declared whole domains.
    peptides=[]
    for (raw,) in db.execute("SELECT record_json FROM rbp2go_support_records WHERE source_file='download_peptide_data.export'").fetchall():
        x=json.loads(raw);requested=x.get('UniProt_ID','');uid=aliases.get(requested,requested)
        try:a,b=int(x['Peptide_Start'])-1,int(x['Peptide_Stop'])
        except (ValueError,KeyError):a=b=None
        sequence=x.get('Peptide_Sequence','');canon=seq_rows.get(uid)
        status='reference_sequence_unavailable' if canon is None else 'exact_sequence_match' if a is not None and 0<=a<b<=canon[4] and canon[3][a:b]==sequence else 'coordinate_or_sequence_mismatch'
        peptides.append((key('RBP2GO_peptide',raw),requested,uid,sequence,a,b,status,x.get('Original_Datasets'),js(x)))
    tables={
        'protein_entries':('accession VARCHAR,entry_name VARCHAR,entry_type VARCHAR,taxon INTEGER,length_aa INTEGER,sequence_sha256 VARCHAR,annotation_score DOUBLE,protein_existence VARCHAR,sequence_version INTEGER,entry_version INTEGER,last_annotation_update VARCHAR,description_json VARCHAR,genes_json VARCHAR',protein_rows),
        'protein_accession_resolution':('source_accession VARCHAR,accession VARCHAR,status VARCHAR,candidates_json VARCHAR',accession_resolution),
        'protein_factor_mapping':('accession VARCHAR,factor_id VARCHAR,gene_id VARCHAR,mapping_method VARCHAR,status VARCHAR,source_gene_ids_json VARCHAR',mapping),
        'protein_sequences':('sequence_id VARCHAR,accession VARCHAR,sequence_type VARCHAR,sequence VARCHAR,length_aa INTEGER,sequence_sha256 VARCHAR,source VARCHAR,source_sequence_verified BOOLEAN',seq_rows.values()),
        'protein_isoforms':('isoform_id VARCHAR,accession VARCHAR,name VARCHAR,sequence_status VARCHAR,source_edit_ids_json VARCHAR,sequence_available BOOLEAN,record_json VARCHAR',iso_rows),
        'protein_features':('feature_key VARCHAR,accession VARCHAR,sequence_id VARCHAR,feature_type VARCHAR,description VARCHAR,start0 INTEGER,end0 INTEGER,coordinate_status VARCHAR,source_feature_id VARCHAR,source_location_json VARCHAR,record_json VARCHAR',features),
        'protein_feature_evidence':('feature_key VARCHAR,evidence_code VARCHAR,evidence_source VARCHAR,evidence_id VARCHAR,record_json VARCHAR',feature_ev),
        'protein_function_comments':('comment_key VARCHAR,accession VARCHAR,comment_type VARCHAR,record_json VARCHAR',comments),
        'protein_database_xrefs':('xref_key VARCHAR,accession VARCHAR,database_name VARCHAR,source_id VARCHAR,isoform_id VARCHAR,properties_json VARCHAR',xrefs),
        'protein_transcript_links':('link_key VARCHAR,accession VARCHAR,isoform_id VARCHAR,transcript_id_versioned VARCHAR,transcript_id VARCHAR,protein_id_versioned VARCHAR,protein_id VARCHAR,gene_id_versioned VARCHAR,gene_id VARCHAR,source VARCHAR',transcript_links),
        'protein_isoform_edits':('isoform_id VARCHAR,accession VARCHAR,source_edit_id VARCHAR,feature_key VARCHAR,start0 INTEGER,end0 INTEGER,replacement_sequence VARCHAR,verification_status VARCHAR',edit_rows),
        'protein_isoform_feature_consequences':('isoform_id VARCHAR,feature_key VARCHAR,accession VARCHAR,structural_consequence VARCHAR,isoform_start0 INTEGER,isoform_end0 INTEGER,evidence_scope VARCHAR,functional_interpretation VARCHAR',feature_impacts),
        'domain_catalog':('source VARCHAR,domain_id VARCHAR,name VARCHAR,entry_type VARCHAR,rna_binding_annotation BOOLEAN,record_json VARCHAR',[(source,*v) for (source,_),v in catalog.items()]),
        'factor_domain_annotations':('factor_id VARCHAR,uniprot_id VARCHAR,domain_id VARCHAR,is_annotated_RBD BOOLEAN,source_RBD_multiplicity INTEGER,source VARCHAR,coordinate_scope VARCHAR',sorted(memberships)),
        'protein_rbp2go_properties':('factor_id VARCHAR,uniprot_id VARCHAR,length_aa DOUBLE,RBD_count DOUBLE,RBD_content_fraction DOUBLE,IDR_content_fraction DOUBLE,has_IDR VARCHAR,LCD_status VARCHAR,coiled_coil_status VARCHAR',protein_stats),
        'protein_domain_matches':('match_id VARCHAR,sequence_id VARCHAR,accession VARCHAR,source VARCHAR,entry_type VARCHAR,domain_id VARCHAR,name VARCHAR,start0 INTEGER,end0 INTEGER,sequence_mapping_status VARCHAR,model VARCHAR,score DOUBLE,fragments_json VARCHAR',[( *r[:12],js(r[12])) for r in domain_matches]),
        'protein_domain_fragments':('match_id VARCHAR,fragment_index INTEGER,start0 INTEGER,end0 INTEGER,discontinuity_status VARCHAR',domain_fragments),
        'protein_interpro_sequence_audit':('requested_accession VARCHAR,accession VARCHAR,sequence_status VARCHAR,source_length INTEGER,source_sequence_sha256 VARCHAR,source_is_fragment BOOLEAN,source_metadata_json VARCHAR',interpro_audit),
        'rna_binding_peptides':('peptide_key VARCHAR,source_uniprot_id VARCHAR,accession VARCHAR,sequence VARCHAR,start0 INTEGER,end0 INTEGER,coordinate_status VARCHAR,original_datasets VARCHAR,record_json VARCHAR',peptides),
        'protein_isoform_sequence_checks':('isoform_id VARCHAR,status VARCHAR',sorted(iso_checks.items())),
    }
    db.execute('BEGIN TRANSACTION')
    try:
        report['tables']={}
        for table,(schema,rows) in tables.items():
            n=create_rows(db,table,schema,rows);report['tables'][table]=n
            print(js({'table':table,'rows':n}),flush=True)
        build_coding_layer(db,factors,uid_factors,seq_rows,domain_matches,report)
        db.execute("""CREATE OR REPLACE VIEW sequence_verified_domain_matches AS
            SELECT * FROM protein_domain_matches WHERE sequence_mapping_status IN
            ('exact_sequence_match','exact_source_edit_sequence_verified')
            AND entry_type IN ('domain','repeat','zinc finger')""")
        db.execute("""CREATE OR REPLACE TABLE protein_peptide_domain_overlaps AS
            SELECT DISTINCT p.peptide_key,d.match_id,
                greatest(p.start0,f.start0) overlap_start0,least(p.end0,f.end0) overlap_end0,
                'coordinate_overlap_only_not_domain_specific_binding_validation'::VARCHAR interpretation
            FROM rna_binding_peptides p JOIN sequence_verified_domain_matches d
            ON p.accession=d.accession AND d.sequence_id=p.accession
            JOIN protein_domain_fragments f USING(match_id)
            WHERE p.coordinate_status='exact_sequence_match' AND p.start0<f.end0 AND p.end0>f.start0""")
        create_rows(db,'protein_dimension_input_assets',
                    'source VARCHAR,path VARCHAR,bytes BIGINT,sha256 VARCHAR,retrieved_at VARCHAR',
                    [(r['source'],str(ROOT/r['path']),r['bytes'],r['sha256'],r['retrieved_at']) for r in manifest])
        create_rows(db,'protein_dimension_policy','key VARCHAR,value VARCHAR',[
            ('version','1.0'),('species','Mus_musculus'),('protein_coordinate_convention','0-based half-open; source positions preserved in raw JSON'),
            ('mapping','Existing UniProt accessions and explicit mouse Ensembl gene IDs; ambiguous or conflicting mappings excluded from resolved views'),
            ('domain','InterPro family, superfamily and member-signature matches remain distinct; overlapping source models are not independent domains'),
            ('isoforms','Only source VAR_SEQ edits that reproduce complete isoform FASTA are projected'),
            ('activity','No 0-1 activity weights, KO signs, motif-domain attribution or functional loss inferred from a coordinate overlap'),
            ('sequences','InterPro coordinates transfer only on complete source sequence identity'),
            ('genomics','CDS projection requires exact translation equality to a same-gene protein sequence; genome builds and transcript/protein versions retained'),
            ('WT_evidence','Mutagenesis and natural variants retained as reference annotations; not added to WT training priors'),
            ('dataset_chains','Reference isoforms only; observed truncated splice chains not automatically treated as coding transcripts'),
            ('collection_status',state['status']),('created_at',report['started_at']),
        ])
        build_views(db)
        after={t:db.execute(f'SELECT count(*),bit_xor(hash(t)) FROM {t} t').fetchone() for t in protected}
        if before!=after:raise ValueError('Original evidence tables changed')
        report['protected_original_tables_unchanged']=True
        report['checks']=validate_extension(db)
        if not all(v['passed'] for v in report['checks'].values()):
            raise ValueError('Protein dimension validation failed: '+js(report['checks']))
        db.execute('COMMIT');db.execute('CHECKPOINT')
        report['coverage']=coverage(db)
        report['status']='complete' if not report['partial_import'] else 'partial_sources_integrated'
        report['passed']=True
    except BaseException as exc:
        try:db.execute('ROLLBACK')
        except Exception:pass
        report.update(status='failed',error=str(exc));raise
    finally:
        db.close();report['updated_at']=time.strftime('%Y-%m-%dT%H:%M:%S%z')
        (ROOT/'integration_status.json').write_text(json.dumps(report,indent=2,ensure_ascii=False))
    checksum=sha256(DBPATH);sums=MAIN/'SHA256SUMS'
    lines=sums.read_text().splitlines()
    lines=[checksum+'  mouse_rbp.duckdb' if l.endswith('  mouse_rbp.duckdb') else l for l in lines]
    sums.write_text('\n'.join(lines)+'\n')
    (MAIN/'protein_dimensions.json').write_text(json.dumps(report,indent=2,ensure_ascii=False))
    return report


def sha256_text(value):
    return hashlib.sha256(value.encode()).hexdigest()


def build_coding_layer(db,factors,uid_factors,sequences,domain_matches,report):
    gene_set={stable(factors[f][0]) for ids in uid_factors.values() for f in ids if factors[f][0].startswith('ENSMUSG')}
    source_rows=[];segments_by_transcript=collections.defaultdict(list);transcript_meta={}
    assets=[]
    for assembly,root in REFERENCES.items():
        path=root/'genes/genes.gtf.gz';assets.append((assembly,str(path),sha256(path)))
        with gzip.open(path,'rt') as stream:
            for line in stream:
                if line.startswith('#'):continue
                c=line.rstrip().split('\t')
                if c[2] not in {'transcript','CDS','start_codon','stop_codon','UTR','Selenocysteine'}:continue
                attrs=dict(re.findall(r'(\w+) "([^"]*)"',c[8]));gene=stable(attrs.get('gene_id',''))
                exon_match=re.search(r'(?:^|;\s*)exon_number\s+(?:"(\d+)"|(\d+))',c[8])
                if exon_match:attrs['exon_number']=exon_match[1] or exon_match[2]
                if gene not in gene_set:continue
                tx=stable(attrs.get('transcript_id',''));tkey=(assembly,tx)
                if c[2]=='transcript':
                    transcript_meta[tkey]=attrs
                    continue
                row=(assembly,gene,tx,attrs.get('transcript_version'),stable(attrs.get('protein_id','')),attrs.get('protein_id'),
                     c[2],c[0],int(c[3])-1,int(c[4]),c[6],int(c[7]) if c[7]!='.' else None,
                     int(attrs['exon_number']) if attrs.get('exon_number') else None,js(attrs),c[8])
                source_rows.append(row)
                if c[2]=='CDS':segments_by_transcript[tkey].append(row)
    n=create_rows(db,'reference_coding_segments',
                 'assembly VARCHAR,gene_id VARCHAR,transcript_id VARCHAR,transcript_version VARCHAR,protein_id VARCHAR,protein_id_versioned VARCHAR,feature_type VARCHAR,chromosome VARCHAR,start0 BIGINT,end0 BIGINT,strand VARCHAR,phase INTEGER,exon_number INTEGER,attributes_json VARCHAR,raw_attributes VARCHAR',source_rows)
    report['tables']['reference_coding_segments']=n
    fasta={a:pysam.FastaFile(str(r/'fasta/genome.fa')) for a,r in REFERENCES.items()}
    candidates=collections.defaultdict(list)
    for sid,row in sequences.items():
        uid=row[1]
        for factor in uid_factors.get(uid,set()):
            gene=stable(factors[factor][0]);candidates[(gene,row[5])].append(sid)
    source_xrefs=collections.defaultdict(set)
    for uid,iso,tx in db.execute('SELECT accession,isoform_id,transcript_id FROM protein_transcript_links').fetchall():
        source_xrefs[tx].add(iso or uid)
    transcript_rows=[];matches=[];coding_layout={}
    try:
        for (assembly,tx),rows in sorted(segments_by_transcript.items()):
            strand=rows[0][10];rows.sort(key=lambda r:r[8],reverse=strand=='-')
            pieces=[];layout=[];offset=0;phase_valid=True
            first_phase=rows[0][11] or 0
            for n,row in enumerate(rows):
                chrom,a,b=row[7:10];phase=row[11] or 0
                if n and phase!=(3-offset%3)%3:phase_valid=False
                seq=fasta[assembly].fetch(chrom,a,b).upper()
                if strand=='-':seq=str(Seq(seq).reverse_complement())
                if n==0 and first_phase:
                    seq=seq[first_phase:]
                    if strand=='+':a+=first_phase
                    else:b-=first_phase
                layout.append((offset,offset+len(seq),chrom,a,b,strand,row[12]));offset+=len(seq);pieces.append(seq)
            cds=''.join(pieces);aa=str(Seq(cds[:len(cds)//3*3]).translate())
            if aa.endswith('*'):aa=aa[:-1]
            status='valid_translation' if phase_valid and len(cds)%3==0 and '*' not in aa else 'incomplete_or_inconsistent_CDS'
            attrs=transcript_meta.get((assembly,tx),{});raw_tags=';'.join(r[14] for r in rows)
            partial='cds_start_NF' in raw_tags or 'cds_end_NF' in raw_tags or first_phase!=0
            if partial and status=='valid_translation':status='partial_CDS'
            gene=rows[0][1];digest=sha256_text(aa)
            transcript_rows.append((assembly,gene,tx,attrs.get('transcript_version',rows[0][3]),rows[0][4],rows[0][5],
                                    len(cds),len(aa),digest,status,phase_valid,partial,first_phase))
            coding_layout[(assembly,tx)]=layout
            if status!='valid_translation':continue
            for sid in sorted(set(candidates.get((gene,digest),[]))):
                if sequences[sid][3]!=aa:raise ValueError('Sequence hash collision')
                matches.append((assembly,gene,tx,rows[0][4],sid,sequences[sid][1],
                                'source_transcript_xref_and_exact_sequence' if sid in source_xrefs[tx] else 'same_gene_exact_sequence',True))
    finally:
        for handle in fasta.values():handle.close()
    create_rows(db,'reference_coding_transcripts',
                'assembly VARCHAR,gene_id VARCHAR,transcript_id VARCHAR,transcript_version VARCHAR,protein_id VARCHAR,protein_id_versioned VARCHAR,CDS_length_nt INTEGER,translated_length_aa INTEGER,translation_sha256 VARCHAR,translation_status VARCHAR,phase_consistent BOOLEAN,partial_CDS BOOLEAN,first_phase INTEGER',transcript_rows)
    create_rows(db,'protein_local_transcript_matches',
                'assembly VARCHAR,gene_id VARCHAR,transcript_id VARCHAR,protein_id VARCHAR,sequence_id VARCHAR,accession VARCHAR,mapping_method VARCHAR,sequence_verified BOOLEAN',matches)
    db.execute("""CREATE OR REPLACE TABLE protein_transcript_link_audit AS
        SELECT x.link_key,x.accession,x.isoform_id,x.transcript_id,x.gene_id,r.assembly,
          x.transcript_id_versioned,x.protein_id_versioned,r.transcript_version,
          r.protein_id_versioned local_protein_id_versioned,
          CASE WHEN r.transcript_id IS NULL THEN 'not_in_local_coding_reference'
               WHEN x.gene_id!=r.gene_id THEN 'gene_id_mismatch'
               WHEN EXISTS(SELECT 1 FROM protein_local_transcript_matches m WHERE m.assembly=r.assembly
                    AND m.transcript_id=x.transcript_id AND m.accession=x.accession
                    AND (x.isoform_id IS NULL OR m.sequence_id=x.isoform_id)) THEN 'exact_sequence_verified'
               WHEN r.translation_status!='valid_translation' THEN r.translation_status
               ELSE 'no_exact_source_isoform_sequence_match' END mapping_status
          FROM protein_transcript_links x LEFT JOIN reference_coding_transcripts r USING(transcript_id)""")
    verified_status={'exact_sequence_match','exact_source_edit_sequence_verified'}
    domains_by_seq=collections.defaultdict(list)
    for dom in domain_matches:
        if dom[9] in verified_status and dom[4] in {'domain','repeat','zinc finger'}:domains_by_seq[dom[1]].append(dom)
    def blocks():
        for assembly,gene,tx,pid,sid,uid,method,_ in matches:
            for dom in domains_by_seq[sid]:
                for fragment_index,(a,b,dc) in enumerate(dom[12]):
                    projected=genomic_blocks(coding_layout[(assembly,tx)],a,b)
                    if sum(p[2]-p[1] for p in projected)!=(b-a)*3:
                        raise ValueError('Domain projection does not preserve CDS length')
                    for block_index,(chrom,start,end,strand,exon,cl,cr) in enumerate(projected):
                        yield (assembly,gene,tx,sid,dom[0],fragment_index,block_index,chrom,start,end,strand,exon,cl,cr,'exact_CDS_translation_and_protein_sequence')
    create_rows(db,'protein_domain_genomic_blocks',
                'assembly VARCHAR,gene_id VARCHAR,transcript_id VARCHAR,sequence_id VARCHAR,match_id VARCHAR,fragment_index INTEGER,block_index INTEGER,chromosome VARCHAR,start0 BIGINT,end0 BIGINT,strand VARCHAR,exon_number INTEGER,CDS_offset_start0 INTEGER,CDS_offset_end0 INTEGER,mapping_evidence VARCHAR',blocks())
    create_rows(db,'protein_reference_assets','assembly VARCHAR,path VARCHAR,sha256 VARCHAR',assets)


def build_views(db):
    db.execute("""CREATE OR REPLACE VIEW factor_protein_domains AS SELECT m.factor_id,m.gene_id,d.*
        FROM protein_factor_mapping m JOIN sequence_verified_domain_matches d USING(accession)
        WHERE m.status='resolved'""")
    db.execute("""CREATE OR REPLACE VIEW factor_protein_features AS SELECT m.factor_id,m.gene_id,f.*
        FROM protein_factor_mapping m JOIN protein_features f USING(accession) WHERE m.status='resolved'""")
    db.execute("""CREATE OR REPLACE VIEW factor_protein_isoforms AS SELECT m.factor_id,m.gene_id,i.*,
        c.status edit_sequence_check FROM protein_factor_mapping m JOIN protein_isoforms i USING(accession)
        LEFT JOIN protein_isoform_sequence_checks c USING(isoform_id) WHERE m.status='resolved'""")
    db.execute("""CREATE OR REPLACE VIEW reference_transcript_function_context AS
        SELECT r.*,f.transcript_type transcript_biotype,
          coalesce(f.transcript_type='nonsense_mediated_decay',false) source_annotated_NMD,
          'reference_annotation_not_NMD_prediction_for_observed_chains'::VARCHAR interpretation
        FROM reference_coding_transcripts r LEFT JOIN reference_features f
        ON f.assembly=r.assembly AND f.transcript_id=r.transcript_id AND f.feature_type='transcript'""")
    db.execute("""CREATE OR REPLACE VIEW protein_dimension_coverage AS WITH
        domains AS (SELECT factor_id,count(DISTINCT match_id) matches FROM factor_protein_domains GROUP BY factor_id),
        iso AS (SELECT factor_id,count(DISTINCT isoform_id) isoforms FROM factor_protein_isoforms GROUP BY factor_id),
        seq AS (SELECT factor_id,count(DISTINCT accession) proteins FROM protein_factor_mapping WHERE status='resolved' GROUP BY factor_id),
        canonical_flags AS (SELECT i.canonical_factor_id factor_id,
          bool_or(f.rbp2go_listed OR f.other_rbp_source) canonical_source_supported_RBP,
          bool_or(f.broad_splicing_candidate) canonical_splicing_candidate
          FROM factor_catalog f JOIN family_sharing_factor_identity i ON i.original_factor_id=f.factor_id
          GROUP BY i.canonical_factor_id)
        SELECT f.*,coalesce(s.proteins,0) sequenced_proteins,coalesce(d.matches,0) located_domain_matches,
        coalesce(i.isoforms,0) annotated_protein_isoforms,
        coalesce(c.canonical_source_supported_RBP,false) canonical_source_supported_RBP,
        coalesce(c.canonical_splicing_candidate,false) canonical_splicing_candidate
        FROM factor_catalog f LEFT JOIN seq s USING(factor_id) LEFT JOIN domains d USING(factor_id)
        LEFT JOIN iso i USING(factor_id) LEFT JOIN canonical_flags c USING(factor_id)""")


def validate_extension(db):
    queries={
        'protein_entries_mouse_only':"SELECT count(*) FROM protein_entries WHERE taxon!=10090",
        'protein_entry_unique':"SELECT count(*)-count(DISTINCT accession) FROM protein_entries",
        'sequence_unique':"SELECT count(*)-count(DISTINCT sequence_id) FROM protein_sequences",
        'sequence_lengths':"SELECT count(*) FROM protein_sequences WHERE length(sequence)!=length_aa OR sha256(sequence)!=sequence_sha256",
        'mapping_factor_fk':"SELECT count(*) FROM protein_factor_mapping m LEFT JOIN factors f USING(factor_id) WHERE f.factor_id IS NULL",
        'resolved_mapping_unique':"SELECT count(*) FROM (SELECT accession FROM protein_factor_mapping WHERE status='resolved' GROUP BY accession HAVING count(*)!=1)",
        'feature_keys_unique':"SELECT count(*)-count(DISTINCT feature_key) FROM protein_features",
        'domain_keys_unique':"SELECT count(*)-count(DISTINCT match_id) FROM protein_domain_matches",
        'domain_fragment_fk':"SELECT count(*) FROM protein_domain_fragments f LEFT JOIN protein_domain_matches d USING(match_id) WHERE d.match_id IS NULL",
        'verified_domain_bounds':"SELECT count(*) FROM sequence_verified_domain_matches d LEFT JOIN protein_sequences s USING(sequence_id) WHERE s.sequence_id IS NULL OR d.start0<0 OR d.end0<=d.start0 OR d.end0>s.length_aa",
        'verified_fragment_bounds':"SELECT count(*) FROM protein_domain_fragments f JOIN sequence_verified_domain_matches d USING(match_id) WHERE f.start0<d.start0 OR f.end0>d.end0 OR f.start0>=f.end0",
        'isoform_feature_retention_bounds':"SELECT count(*) FROM protein_isoform_feature_consequences c JOIN protein_sequences s ON s.sequence_id=c.isoform_id WHERE c.structural_consequence='retained_unchanged' AND (c.isoform_start0<0 OR c.isoform_end0<=c.isoform_start0 OR c.isoform_end0>s.length_aa)",
        'genomic_intervals':"SELECT count(*) FROM protein_domain_genomic_blocks WHERE start0<0 OR end0<=start0 OR end0-start0!=CDS_offset_end0-CDS_offset_start0",
        'genomic_domain_fk':"SELECT count(*) FROM protein_domain_genomic_blocks b LEFT JOIN sequence_verified_domain_matches d USING(match_id) WHERE d.match_id IS NULL",
        'genomic_reference_fk':"SELECT count(*) FROM protein_domain_genomic_blocks b LEFT JOIN reference_coding_transcripts t USING(assembly,transcript_id) WHERE t.transcript_id IS NULL OR t.translation_status!='valid_translation'",
        'no_invented_functional_weights':"SELECT count(*) FROM information_schema.columns WHERE table_name IN ('protein_isoform_feature_consequences','protein_domain_genomic_blocks','protein_features') AND column_name IN ('functional_weight','activity_score','KO_effect','regulatory_sign')",
    }
    results={}
    for name,sql in queries.items():
        value=db.execute(sql).fetchone()[0];results[name]={'value':value,'expected':0,'passed':value==0}
    return results


def coverage(db):
    def records(sql):
        c=db.execute(sql);return [dict(zip([x[0] for x in c.description],r)) for r in c.fetchall()]
    return {'catalog':records("""SELECT count(*) catalog_rows,
        count(*) FILTER(WHERE rbp2go_listed OR other_rbp_source) source_supported_RBP_rows,
        count(*) FILTER(WHERE broad_splicing_candidate) broad_splicing_rows,
        count(*) FILTER(WHERE broad_splicing_candidate AND sequenced_proteins>0) splicing_with_sequence,
        count(*) FILTER(WHERE broad_splicing_candidate AND located_domain_matches>0) splicing_with_located_domains,
        count(*) FILTER(WHERE broad_splicing_candidate AND annotated_protein_isoforms>0) splicing_with_annotated_isoforms,
        count(*) FILTER(WHERE (rbp2go_listed OR other_rbp_source) AND located_domain_matches>0) RBP_with_located_domains
        FROM protein_dimension_coverage""")[0],
        'canonical_splicing_coverage':records("""SELECT count(*) candidates,
          count(*) FILTER(WHERE sequenced_proteins>0) with_sequence,
          count(*) FILTER(WHERE located_domain_matches>0) with_located_domains,
          count(*) FILTER(WHERE annotated_protein_isoforms>0) with_annotated_isoforms
          FROM protein_dimension_coverage WHERE canonical_splicing_candidate""")[0],
        'domain_sources':records('SELECT source,entry_type,count(*) records FROM protein_domain_matches GROUP BY source,entry_type'),
        'interpro_sequence_consistency':records('SELECT sequence_status,count(*) proteins FROM protein_interpro_sequence_audit GROUP BY sequence_status'),
        'isoform_consequences':records('SELECT structural_consequence,count(*) records FROM protein_isoform_feature_consequences GROUP BY structural_consequence'),
        'isoform_sequence_checks':records('SELECT status,count(*) isoforms FROM protein_isoform_sequence_checks GROUP BY status'),
        'transcript_mapping':records('SELECT assembly,count(DISTINCT gene_id) genes,count(DISTINCT transcript_id) transcripts,count(DISTINCT sequence_id) sequences FROM protein_local_transcript_matches GROUP BY assembly'),
        'genomic_domains':records('SELECT assembly,count(*) blocks,count(DISTINCT gene_id) genes,count(DISTINCT match_id) domain_matches FROM protein_domain_genomic_blocks GROUP BY assembly'),
        'peptide_consistency':records('SELECT coordinate_status,count(*) peptides FROM rna_binding_peptides GROUP BY coordinate_status')}


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--allow-partial',action='store_true')
    args=parser.parse_args();print(json.dumps(integrate(args.allow_partial),indent=2,ensure_ascii=False))
