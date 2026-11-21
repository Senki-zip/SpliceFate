"""Audit transcript-context redundancy without altering the unified database."""
import sys,json,time
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parent))
from query_mouse_rbp_database import connect,ROOT
c=connect();c.execute("SET memory_limit='10GB'");c.execute('SET threads=2');c.execute('SET temp_directory=\''+str(ROOT/'_dedup_audit_tmp')+'\'')
rows=[]
chroms=[r[0] for r in c.execute('select distinct chromosome from motif_gene_loci order by chromosome').fetchall()]
for chrom in chroms:
 t=time.time()
 r=c.execute('''SELECT count(*) raw_rows,
 count(DISTINCT (strand,source_motif_id,list_sort(list_zip(block_starts0,block_ends0)))) unique_motif_genomic_sites,
 count(DISTINCT (strand,list_sort(list_zip(block_starts0,block_ends0)))) unique_genomic_locations,
 count(*) FILTER(WHERE species_specific) mouse_matched_rows,
 count(DISTINCT (strand,source_motif_id,list_sort(list_zip(block_starts0,block_ends0)))) FILTER(WHERE species_specific) mouse_matched_unique_motif_sites
 FROM motif_sites WHERE chromosome=?''',[chrom]).fetchone()
 item=dict(zip(['raw_rows','unique_motif_genomic_sites','unique_genomic_locations','mouse_matched_rows','mouse_matched_unique_motif_sites'],r));item['chromosome']=chrom;rows.append(item)
 print(json.dumps(item),flush=True)
 result={'status':'running','definition':'Deduplicate transcript/feature context by mm39 chromosome + strand + motif ID + exact sorted genomic blocks; distinct physical locations additionally ignore motif ID. Split blocks retained; these are candidate matches, not occupancy. Counts of different motifs are not statistically independent.','chromosomes_completed':len(rows),'chromosomes_planned':len(chroms),'by_chromosome':rows,'totals':{k:sum(x[k] for x in rows) for k in item if k!='chromosome'}}
 (ROOT/'site_redundancy_audit.json').write_text(json.dumps(result,indent=2))
result['status']='complete';result['transcript_context_redundancy_fraction']=1-result['totals']['unique_motif_genomic_sites']/result['totals']['raw_rows'];(ROOT/'site_redundancy_audit.json').write_text(json.dumps(result,indent=2));print(json.dumps(result['totals']),flush=True)
