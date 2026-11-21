"""Check full-store relations/bounds, sampled sequence scores and liftover integrity."""
from __future__ import annotations
import collections
import json
from pathlib import Path
import time

from build_mouse_rbp_database import ROOT, DBPATH, Fasta, matrices, sha256
from query_mouse_rbp_database import connect, records


def main():
    db=connect()
    report={'started_at':time.strftime('%Y-%m-%dT%H:%M:%S%z'),'checks':{},'tables':{},'statistics':{}}
    checks=report['checks'];stats=report['statistics']
    def check(name,sql,expected=0):
        value=db.execute(sql).fetchone()[0]
        checks[name]={'value':value,'expected':expected,'passed':value==expected}
        print(name,value,flush=True)
    tables=[r[0] for r in db.execute("SELECT table_name FROM information_schema.tables WHERE table_type='BASE TABLE'").fetchall()]
    for name in tables:
        report['tables'][name]=db.execute(f'SELECT count(*) FROM "{name}"').fetchone()[0]
    check('mouse_only_factors',"SELECT count(*) FROM factors WHERE species!='Mus_musculus'")
    check('factor_ids_unique',"SELECT count(*)-count(DISTINCT factor_id) FROM factors")
    check('motif_ids_unique',"SELECT count(*)-count(DISTINCT motif_id) FROM motifs")
    for table in ['factor_sources','factor_motifs','functions','experiments','gene_binding','rbp2go_scores','motif_evidence']:
        check(table+'_factor_fk',f'SELECT count(*) FROM {table} a LEFT JOIN factors f USING(factor_id) WHERE f.factor_id IS NULL')
    check('motif_source_keys_unique','SELECT count(*)-count(DISTINCT (source,source_motif_id)) FROM motif_sources')
    check('valid_motif_fk',"SELECT count(*) FROM motif_sources s LEFT JOIN motifs m USING(motif_id) WHERE s.status='valid' AND m.motif_id IS NULL")
    check('factor_motif_source_fk','SELECT count(*) FROM factor_motifs f LEFT JOIN motif_sources s USING(source,source_motif_id) WHERE s.source IS NULL')
    check('factor_motif_exact_rows_unique','SELECT count(*)-count(DISTINCT (factor_id,source,source_motif_id,evidence,experiment_type,pmid)) FROM factor_motifs')
    check('clip_rows_match_imported_experiments',
          'SELECT (SELECT count(*) FROM clip_peaks_original)-(SELECT coalesce(sum(peak_rows),0) FROM experiments)')
    check('clip_rows_match_each_experiment',"""SELECT count(*) FROM experiments e LEFT JOIN
          (SELECT experiment_id,count(*) n FROM clip_peaks_original GROUP BY experiment_id) p USING(experiment_id)
          WHERE coalesce(p.n,0)!=e.peak_rows""")
    stats['previous_snapshot_clip_rows']={'previous':5467198,
        'current':db.execute('SELECT count(*) FROM clip_peaks_original').fetchone()[0]}
    check('peak_ids_unique','SELECT count(*)-count(DISTINCT peak_id) FROM clip_peaks_original')
    check('liftover_result_single','SELECT count(*)-count(DISTINCT peak_id) FROM clip_peaks_lifted_raw')
    check('liftover_status_complete','SELECT (SELECT count(*) FROM clip_peaks_original)-(SELECT count(*) FROM clip_peak_mapping)')
    check('liftover_lengths_preserved',"SELECT count(*) FROM clip_peaks_mm39 m JOIN clip_peaks_original o USING(peak_id) WHERE m.end0-m.start0 != o.end0-o.start0")
    check('peak_bounds',"SELECT count(*) FROM clip_peaks p JOIN chromosomes c USING(assembly,chromosome) WHERE p.start0<0 OR p.end0<=p.start0 OR p.end0>c.length")
    check('peak_unknown_chromosomes',"SELECT count(*) FROM clip_peaks p LEFT JOIN chromosomes c USING(assembly,chromosome) WHERE c.chromosome IS NULL")
    check('rna_structure_feature_unique','SELECT count(*)-count(DISTINCT feature_id) FROM rna_structure')
    check('rna_structure_length','SELECT count(*) FROM rna_structure WHERE structure IS NOT NULL AND length(sequence)!=length(structure)')
    check('ornament_feature_unique','SELECT count(*)-count(DISTINCT feature_id) FROM ornament_features')
    check('perturbation_effect_columns_absent',"""SELECT count(*) FROM information_schema.columns
         WHERE table_name IN ('gene_binding','factor_motifs','motif_sites') AND
         (lower(column_name) LIKE '%delta_psi%' OR lower(column_name) LIKE '%logfc%' OR lower(column_name) LIKE '%crispr%')""")
    if 'protein_interaction_evidence' in tables:
        check('interaction_evidence_ids_unique','SELECT count(*)-count(DISTINCT evidence_id) FROM protein_interaction_evidence')
        check('interaction_pairs_unique','SELECT count(*)-count(DISTINCT (factor_a,factor_b)) FROM protein_interaction_pairs')
        for col in ['factor_a','factor_b']:
            check('interaction_'+col+'_fk',f'SELECT count(*) FROM protein_interaction_evidence e LEFT JOIN factors f ON e.{col}=f.factor_id WHERE e.{col} IS NOT NULL AND f.factor_id IS NULL')
        check('complex_member_factor_fk','SELECT count(*) FROM complex_memberships m LEFT JOIN factors f USING(factor_id) WHERE m.factor_id IS NOT NULL AND f.factor_id IS NULL')
        check('complex_member_complex_fk','SELECT count(*) FROM complex_memberships m LEFT JOIN complexes c USING(source,complex_id,taxon) WHERE c.complex_id IS NULL')
        check('human_complex_not_mapped_to_mouse',"SELECT count(*) FROM complex_memberships WHERE taxon='9606' AND factor_id IS NOT NULL")
        check('interaction_negatives_excluded','SELECT count(*) FROM positive_protein_interactions WHERE negative')
        check('interaction_nonphysical_excluded',"SELECT count(*) FROM positive_protein_interactions WHERE relation_type IN ('colocalization','proximity','other_typed_relation')")
        check('interaction_sign_not_inferred','SELECT count(*) FROM protein_interaction_evidence WHERE regulatory_sign IS NOT NULL OR effect_size IS NOT NULL')
        check('functional_scores_separate',"SELECT count(*) FROM protein_interaction_evidence WHERE source='STRING_functional12.5'")
        check('no_invented_mouse_stages','SELECT count(*) FROM mouse_spliceosome_stage_memberships')
        stats['interaction_extension']={
            'source_evidence_records':report['tables']['protein_interaction_evidence'],
            'canonical_pairs':report['tables']['protein_interaction_pairs'],
            'complexes':report['tables']['complexes'],
            'complex_memberships':report['tables']['complex_memberships'],
            'motifless_splicing_candidates_with_relations':db.execute('SELECT count(*) FROM factor_interaction_coverage WHERE broad_splicing_candidate AND NOT has_source_or_family_motif AND has_experiment_channel_or_curated_complex').fetchone()[0]}
    stats['site_scan']=records(db,"""SELECT count(*) rows,min(site_id) min_site_id,max(site_id) max_site_id,
        count(DISTINCT source_motif_id) distinct_source_motifs,
        count(*) FILTER(WHERE species_specific) species_matched_rows,
        count(*) FILTER(WHERE len(block_starts0)>1) split_site_rows,
        count(*) FILTER(WHERE contains(region,'I')) intronic_site_rows,
        count(*) FILTER(WHERE CLIP_support IS NULL) clip_support_missing,
        count(*) FILTER(WHERE CLIP_support=0) clip_support_zero,
        count(*) FILTER(WHERE CLIP_support=1) clip_support_cluster,
        count(*) FILTER(WHERE CLIP_support=2) clip_support_same_rbp,
        count(*) FILTER(WHERE len(block_starts0)!=len(block_ends0)
          OR list_sum(list_transform(list_zip(block_starts0,block_ends0),lambda b: b[2]-b[1])) != 7
          OR len(list_filter(list_zip(block_starts0,block_ends0),lambda b: b[1]<0 OR b[2]<=b[1]))>0
          OR start0!=list_min(block_starts0) OR end0!=list_max(block_ends0)) invalid_block_rows,
        count(*) FILTER(WHERE NOT isfinite(score) OR score<0 OR score>1) invalid_scores
        FROM motif_sites""")[0]
    for name in ['invalid_block_rows','invalid_scores']:
        value=stats['site_scan'][name]
        checks[name]={'value':value,'expected':0,'passed':value==0}
    checks['site_id_range']={'passed':stats['site_scan']['min_site_id']==1 and stats['site_scan']['max_site_id']==stats['site_scan']['rows']}
    print('sites',stats['site_scan'],flush=True)
    check('site_chromosome_bounds',"""SELECT count(*) FROM motif_sites s LEFT JOIN chromosomes c
       ON c.assembly=s.assembly AND c.chromosome=s.chromosome WHERE c.chromosome IS NULL OR s.start0<0 OR s.end0>c.length""")
    check('site_motif_fk',"""SELECT count(*) FROM (SELECT DISTINCT source,source_motif_id FROM motif_sites) s
       LEFT JOIN motif_sources m USING(source,source_motif_id) WHERE m.source IS NULL""")
    check('site_feature_fk',"""SELECT count(*) FROM (SELECT DISTINCT feature_id FROM motif_sites) s
       LEFT JOIN ornament_features f USING(feature_id) WHERE f.feature_id IS NULL""")
    check('site_gene_transcript_consistency',"""SELECT count(*) FROM
        (SELECT DISTINCT feature_id,transcript_id,gene_id FROM motif_sites) s JOIN ornament_features f USING(feature_id)
        WHERE s.gene_id!=f.ensembl_gene_id OR s.transcript_id!=f.ensembl_transcript_id""")
    # Samples span all motifs, both strands and split/single-block sites; no first-motif-only audit.
    sample=records(db,"""SELECT source_motif_id,strand,len(block_starts0)>1 split,
       min_by(struct_pack(chromosome:=chromosome,starts:=block_starts0,ends:=block_ends0,score:=score,site_id:=site_id),hash(site_id)) sample
       FROM motif_sites GROUP BY source_motif_id,strand,len(block_starts0)>1""")
    sample+=records(db,"""SELECT min_by(source_motif_id,hash(site_id)) source_motif_id,strand,len(block_starts0)>1 split,
       min_by(struct_pack(chromosome:=chromosome,starts:=block_starts0,ends:=block_ends0,score:=score,site_id:=site_id),hash(site_id)) sample
       FROM motif_sites GROUP BY chromosome,strand,len(block_starts0)>1""")
    fa=Fasta('mm39');mm=matrices();seq_results=collections.Counter();failed=[]
    for row in sample:
        r=row['sample'];chrom=r['chromosome']
        if chrom not in fa.index:
            seq_results['unknown_chromosome']+=1;continue
        seq=''.join(fa.fetch(chrom,s,e) for s,e in sorted(zip(r['starts'],r['ends'])))
        if row['strand']==-1:seq=seq.translate(str.maketrans('ACGT','TGCA'))[::-1]
        matrix=mm[int(row['source_motif_id'])][2]
        import math
        lo,hi=math.prod(min(x) for x in matrix),math.prod(max(x) for x in matrix)
        value=(math.prod(p['ACGT'.index(b)] if b in 'ACGT' else 0 for p,b in zip(matrix,seq))-lo)/(hi-lo)
        if len(seq)==7 and abs(value-r['score'])<1e-7:
            seq_results['matched']+=1
        else:
            seq_results['mismatched']+=1
            failed.append({**row,'sequence':seq,'reconstructed_score':value})
    fa.file.close()
    stats['sequence_score_validation']={'sampling':'min hash(site_id) per motif x strand x split flag, plus per chromosome x strand x split flag',
        'sampled_groups':len(sample),'motifs':len({r['source_motif_id'] for r in sample}),'results':dict(seq_results),'mismatches':failed[:50]}
    checks['sequence_score_reconstruction']={'passed':not seq_results['mismatched'] and not seq_results['unknown_chromosome']}
    # Deliberately query the gap of a real split site: outer-span overlap is insufficient.
    split=next(r for r in sample if r['split'])['sample']
    parts=sorted(zip(split['starts'],split['ends']))
    gap=next(((parts[i][1],parts[i+1][0]) for i in range(len(parts)-1) if parts[i][1]<parts[i+1][0]),None)
    if gap:
        n=db.execute('SELECT len(list_filter(list_zip(?::BIGINT[],?::BIGINT[]),lambda b: b[1] < ? AND b[2] > ?))',
                     [split['starts'],split['ends'],gap[1],gap[0]]).fetchone()[0]
        checks['split_site_gap_not_a_hit']={'passed':n==0,'gap':gap,'site_id':split['site_id']}
    stats['liftover']=records(db,'SELECT status,count(*) rows FROM clip_peak_mapping GROUP BY status')
    stats['mapped_clip_without_local_sequence']=records(db,"""SELECT p.chromosome,count(*) rows FROM clip_peaks_mm39 p
       LEFT JOIN local_reference_chromosomes c ON c.assembly='mm39' AND c.chromosome=p.chromosome
       WHERE c.chromosome IS NULL GROUP BY p.chromosome""")
    # Mapped CLIP intervals can contain real reference substitutions; report, rather than pretend to validate binding.
    fa10=Fasta('mm10');fa39=Fasta('mm39');comparison=collections.Counter()
    for row in records(db,"""SELECT o.chromosome old_chrom,o.start0 old_start,o.end0 old_end,o.strand old_strand,
           m.chromosome,m.start0,m.end0,m.strand FROM clip_peaks_mm39 m JOIN clip_peaks_original o USING(peak_id)
           WHERE o.end0-o.start0<=2000 ORDER BY hash(o.peak_id) LIMIT 500"""):
        old=fa10.fetch(row['old_chrom'],row['old_start'],row['old_end']);new=fa39.fetch(row['chromosome'],row['start0'],row['end0'])
        if row['old_strand']!=row['strand']:new=new.translate(str.maketrans('ACGT','TGCA'))[::-1]
        comparison['sequence_identical' if old==new else 'sequence_differs_between_references']+=1
    fa10.file.close();fa39.file.close();stats['liftover_sequence_sample']=dict(comparison)
    stats['source_factor_links']=records(db,'SELECT source,count(*) rows,count(DISTINCT factor_id) factors FROM factor_sources GROUP BY source')
    stats['motif_links']=records(db,'SELECT source,count(*) rows,count(DISTINCT factor_id) factors,count(DISTINCT motif_id) exact_matrices FROM valid_factor_motifs GROUP BY source')
    stats['factor_catalog']=records(db,"""SELECT count(*) catalog_entries,count(*) FILTER(WHERE rbp2go_listed OR other_rbp_source) source_supported_rbp_entries,
       count(*) FILTER(WHERE NOT rbp2go_listed AND NOT other_rbp_source) background_only_entries,
       count(*) FILTER(WHERE broad_splicing_candidate) broad_splicing_candidates FROM factor_catalog""")[0]
    stats['target_mapping']=records(db,'SELECT source,mapping_method,count(*) source_target_identifiers FROM target_identity_mapping GROUP BY source,mapping_method')
    stats['reference_gene_coverage']=records(db,"""SELECT assembly,count(DISTINCT gene_id) genes,count(*) gene_rows FROM genes GROUP BY assembly""")
    stats['site_annotation']=records(db,"""SELECT count(DISTINCT gene_id) target_genes,count(DISTINCT transcript_id) transcripts,
       count(DISTINCT feature_id) rna_features FROM motif_sites""")[0]
    stats['site_files']=len(list((ROOT/'sites').rglob('*.parquet')))
    stats['parquet_bytes']=sum(p.stat().st_size for p in (ROOT/'sites').rglob('*.parquet'))
    stats['duckdb_bytes']=DBPATH.stat().st_size
    report['passed']=all(c['passed'] for c in checks.values())
    report['finished_at']=time.strftime('%Y-%m-%dT%H:%M:%S%z')
    (ROOT/'validation.json').write_text(json.dumps(report,indent=2,ensure_ascii=False))
    # Hash the physical database and all site parts so future copying/corruption can be checked.
    with (ROOT/'SHA256SUMS').open('w') as f:
        for path in [DBPATH]+sorted((ROOT/'sites').rglob('*.parquet')):
            f.write(sha256(path)+'  '+str(path.relative_to(ROOT))+'\n')
    status=json.loads((ROOT/'build_status.json').read_text())
    status.update(status='complete' if report['passed'] else 'validation_failed',validation_passed=report['passed'],updated_at=report['finished_at'])
    (ROOT/'build_status.json').write_text(json.dumps(status,indent=2))
    site=stats['site_scan'];factor=stats['factor_catalog']
    lines=['# 小鼠 RBP 统一数据库验收记录','',
        f"验收时间：{report['finished_at']}；通过：{report['passed']}。",'',
        "统一入口：`mouse_rbp.duckdb`；位点分区：同目录 `sites/`。",'',
        '| 项目 | 数量 |','|---|---:|',
        f"| 因子/蛋白目录行（含背景） | {factor['catalog_entries']:,} |",
        f"| 来源支持的 RBP 目录行 | {factor['source_supported_rbp_entries']:,} |",
        f"| 仅背景行 | {factor['background_only_entries']:,} |",
        f"| GO 宽泛剪接候选 | {factor['broad_splicing_candidates']:,} |",
        f"| oRNAment 小鼠转录本/feature 位点行 | {site['rows']:,} |",
        f"| 其中 species_specific 为真的位点行 | {site['species_matched_rows']:,} |",
        f"| 跨区块位点行 | {site['split_site_rows']:,} |",
        f"| 含内含子区域的位点行 | {site['intronic_site_rows']:,} |",
        f"| 原始 mm10 CLIP 区间 | {report['tables']['clip_peaks_original']:,} |",
        f"| 严格转换后的 mm39 CLIP 区间 | {report['tables']['clip_peaks_mm39']:,} |",
        f"| 来源级基因结合注释 | {report['tables']['gene_binding']:,} |",'',
        '计数单位为保留转录本和来源上下文的记录；位点、实验及因子分别统计。', '',
        '| motif 来源 | 因子数 | 关联行 | 精确数值矩阵组 |','|---|---:|---:|---:|']
    for r in stats['motif_links']:
        lines.append(f"| {r['source']} | {r['factors']:,} | {r['rows']:,} | {r['exact_matrices']:,} |")
    lines+=['','ENCORI 为 HOMER 初始候选；oRNAment 来源中含同源推断。相同矩阵及来源重叠保留来源标签。','',
        '## 坐标与验证','',
        '规范化坐标为 0-based half-open，oRNAment 原始 1-based inclusive 字符串保留。跨外显子位点逐区块处理，mm10 CLIP 原始记录完整保留。','',
        'CLIP 转换：`'+json.dumps(stats['liftover'],ensure_ascii=False)+'`。','',
        '序列评分抽样：`'+json.dumps({k:v for k,v in stats['sequence_score_validation'].items() if k!='mismatches'},ensure_ascii=False)+'`。','',
        f"数据目录文件校验见 `SHA256SUMS`，完整验证结果见 `validation.json`。数据库约 {stats['duckdb_bytes']/1024**3:.2f} GiB，位点 Parquet 约 {stats['parquet_bytes']/1024**3:.2f} GiB。",'',
        '## 限制','',
        '位点层使用已下载的 ENCORI CLIP 与 oRNAment 支持字段。POSTAR3 全量坐标和 oRNAment 官方 CLIP/浏览器压缩包的下载状态为未取得。', '',
        '本地参考为过滤版 GENCODE，oRNAment 为 Ensembl115。数据库补充保存后者的 feature/基因身份；两套注释按来源分别记录。', '',
        '本流程输出统一数据库与来源、坐标及序列验收结果。']
    (ROOT/'SUMMARY.md').write_text('\n'.join(lines)+'\n')
    print(json.dumps({'passed':report['passed'],'statistics':stats},ensure_ascii=False,indent=2),flush=True)
    if not report['passed']:raise SystemExit(1)


if __name__=='__main__':main()
