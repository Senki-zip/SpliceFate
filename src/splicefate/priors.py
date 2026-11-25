"""Native-assembly motif/CLIP anchors and typed protein modules."""
from __future__ import annotations
import collections
import gzip
import json
from pathlib import Path
import numpy as np
import duckdb
from numba import njit


class Fasta:
    def __init__(self,path):
        self.file=Path(path).open('rb')
        self.index={l.split()[0]:tuple(map(int,l.split()[1:5])) for l in Path(str(path)+'.fai').read_text().splitlines()}
    def fetch(self,chrom,start,end):
        n,offset,bases,width=self.index[chrom];start=max(0,start);end=min(n,end)
        if end<=start:return b''
        left=offset+start//bases*width+start%bases;right=offset+(end-1)//bases*width+(end-1)%bases+1
        self.file.seek(left);return self.file.read(right-left).replace(b'\n',b'').replace(b'\r',b'').upper()


def pwm_threshold(matrix,pvalue,resolution=.05):
    p=np.maximum(np.asarray(matrix,float),1e-6);p/=p.sum(1,keepdims=True)
    weights=np.rint(np.log2(p/.25)/resolution).astype(np.int32)
    dist=np.ones(1);lower=0
    for row in weights:
        low,high=int(row.min()),int(row.max());out=np.zeros(len(dist)+high-low)
        for s in row:out[int(s)-low:int(s)-low+len(dist)]+=dist*.25
        dist=out;lower+=low
    eligible=np.flatnonzero(np.cumsum(dist[::-1])[::-1]<=pvalue)
    cutoff=lower+int(eligible[0]) if len(eligible) else lower+len(dist)
    return weights,cutoff


@njit(cache=True)
def scan(sequence,starts,ends,weights):
    scores=np.full(len(starts),-2000000000,dtype=np.int32)
    positions=np.full(len(starts),-1,dtype=np.int64)
    m=len(weights)
    for frag in range(len(starts)):
        for pos in range(starts[frag],ends[frag]-m+1):
            score=0
            for k in range(m):
                b=sequence[pos+k]
                if b>3:score=-2000000000;break
                score+=weights[k,b]
            if score>scores[frag]:scores[frag]=score;positions[frag]=pos-starts[frag]
    return scores,positions


def event_windows(event,width):
    bounds=sorted({p for b in event['branches'] for j in b for p in j});result=[]
    for boundary in bounds:
        start,end=max(0,boundary-width),boundary+width
        if result and start<=result[-1][1]:result[-1][1]=max(end,result[-1][1])
        else:result.append([start,end])
    return result


def build_priors(data,c,output=None):
    db=duckdb.connect(c['input']['database'],read_only=True)
    db.execute("SET file_search_path='"+str(Path(c['input']['database']).parent).replace("'","''")+"'")
    assembly=c['assembly'];species='Homo_sapiens' if assembly=='hg38' else 'Mus_musculus'
    names={r[0] for r in db.execute('SHOW TABLES').fetchall()}
    if 'build_metadata' in names:
        meta=dict(db.execute('SELECT key,value FROM build_metadata').fetchall())
        if meta.get('species')!=species or meta.get('primary_assembly') not in (assembly,'mm39' if assembly=='mm10' else assembly):
            raise ValueError('Prior database species/assembly differs from input configuration')
    identities=dict(db.execute('SELECT original_factor_id,canonical_factor_id FROM family_sharing_factor_identity').fetchall())
    factors={f:(g,s) for f,g,s in db.execute('SELECT factor_id,gene_id,gene_symbol FROM factors').fetchall()}
    available=set(data['genes'].gene);expressed=set(data['genes'].gene[np.asarray(data['G'].sum(0)).ravel()>0])
    targets=c['priors']['targets'];all_targets=targets==['ALL']
    extra=c['priors'].get('extra_targets',[])
    if all_targets:
        seeds=sorted({identities.get(f,f) for (f,) in db.execute('SELECT factor_id FROM candidate_splicing_factors').fetchall()
                      if factors[identities.get(f,f)][0] in expressed})
    else:
        seeds=[]
        for target in targets:
            hits={identities.get(f,f) for f,g,s in db.execute('''SELECT DISTINCT f.factor_id,f.gene_id,f.gene_symbol
                FROM factors f LEFT JOIN factor_aliases a USING(factor_id)
                WHERE upper(f.gene_symbol)=upper(?) OR f.gene_id=? OR upper(a.alias)=upper(?)''',[target,target,target]).fetchall()}
            if len(hits)!=1:raise ValueError(f'Target {target} resolves to {len(hits)} identities')
            fid=next(iter(hits))
            if factors[fid][0] not in expressed:raise ValueError(f'Target {target} absent/unexpressed in full gene input')
            seeds.append(fid)
    for target in extra:
        hits={identities.get(f,f) for (f,) in db.execute('''SELECT DISTINCT f.factor_id
            FROM factors f LEFT JOIN factor_aliases a USING(factor_id)
            WHERE upper(f.gene_symbol)=upper(?) OR f.gene_id=? OR f.factor_id=? OR upper(a.alias)=upper(?)''',[target,target,target,target]).fetchall()}
        if len(hits)!=1:raise ValueError('Extra target does not resolve uniquely: '+target)
        fid=next(iter(hits))
        if factors[fid][0] not in expressed:raise ValueError('Extra target not expressed: '+target)
        seeds.append(fid)
    seeds=sorted(set(seeds))
    # Keep evidence per source. No arbitrary multi-hop diffusion or motif inheritance via PPI.
    pair_evidence=collections.defaultdict(list)
    for a,b,source,rel,eid,score,construct,condition in db.execute('''SELECT factor_a,factor_b,source,relation_type,evidence_id,
       nontransferred_experiment_score,construct_change_flag,perturbation_keyword_flag
       FROM candidate_experiment_supported_protein_interactions''').fetchall():
        if construct or condition or not a or not b or a==b:continue
        if factors[a][0] not in expressed or factors[b][0] not in expressed:continue
        if source=='IntAct':
            allowed=['direct_physical','physical_association','covalent_binding']
            if c['priors']['allow_unspecified_association']:allowed.append('association_unspecified')
            if rel not in allowed:continue
            rank=3 if rel=='direct_physical' else 2
        else:
            if score is None or score<c['priors']['min_experiment_score']:continue
            rank=1+float(score)/1000
        pair_evidence[tuple(sorted([a,b]))].append({'source':source,'relation_type':rel,'evidence_id':eid,'score':score,'rank':rank})
    # Complexes are retained as member groups; pair-specific joint effects are candidate groups,
    # explicitly co_complex rather than invented pairwise contact.
    complexes=collections.defaultdict(set)
    complex_view='candidate_curated_complex_memberships' if 'candidate_curated_complex_memberships' in names else 'mouse_curated_complex_memberships'
    for source,cid,fid in db.execute('SELECT source,complex_id,factor_id FROM '+complex_view+' WHERE factor_id IS NOT NULL').fetchall():
        if factors[fid][0] in expressed:complexes[(source,cid)].add(fid)
    for (source,cid),members in complexes.items():
        for a in members:
            for b in members-{a}:
                pair=tuple(sorted([a,b]));record={'source':source,'relation_type':'co_complex','complex_id':cid,'rank':1.5}
                if record not in pair_evidence[pair]:pair_evidence[pair].append(record)
    selected_pairs={};universe=set(seeds);frontier=set(seeds);depth={f:0 for f in seeds}
    # Candidate retrieval is independent of the learned composition propagation budget.
    adjacency=collections.defaultdict(list)
    for pair,ev in pair_evidence.items():
        item=(max(x['rank'] for x in ev),pair,ev)
        for f in pair:adjacency[f].append(item)
    for layer in range(c['priors'].get('partner_hops',1)):
        discovered=set()
        for factor in sorted(frontier):
            for _,pair,ev in sorted(adjacency[factor],key=lambda r:(-r[0],r[1]))[:c['priors']['max_partners_per_factor']]:
                selected_pairs[pair]=ev
                for f in pair:
                    if f not in universe:discovered.add(f);depth[f]=layer+1
        universe.update(discovered);frontier=discovered
        if not frontier:break
    gids=[factors[f][0] for f in sorted(universe)]
    factor_ids=sorted(universe);gene_indices=[list(data['genes'].gene).index(factors[f][0]) for f in factor_ids]
    N=data['G'][:,gene_indices].multiply(data['scales'][:,None]).toarray()
    X=np.log1p(N)
    relation=db.execute('''SELECT DISTINCT v.factor_id,v.source,v.source_motif_id,v.motif_id,v.evidence,m.matrix_json
       FROM effective_factor_motifs v JOIN motifs m USING(motif_id) WHERE v.factor_id IN (SELECT unnest(?))
       AND v.source!='ENCORI' AND v.evidence!='mutant_binding_exclude_from_WT_default' ''',[factor_ids]).fetchall()
    if not c['priors']['include_family_shared']:relation=[r for r in relation if r[1]!='FamilySharing']
    mat={r[3]:json.loads(r[5]) for r in relation};links=collections.defaultdict(list)
    for f,source,sid,mid,ev,_ in relation:links[mid].append((f,source,sid,ev))
    windows=[];seqs=[];starts=[];ends=[];pos=0
    fasta=Fasta(c['input']['fasta']);encoding=np.full(256,4,dtype=np.uint8)
    for k,b in enumerate(b'ACGT'):encoding[b]=k
    eligible=[]
    for e,event in enumerate(data['events']):
        counts=data['C'][:,data['offsets'][e]:data['offsets'][e+1]]
        nn=np.asarray(counts.sum(1)).ravel()
        if (nn>0).sum()<c['priors']['min_event_cells'] or nn.sum()<c['priors']['min_event_molecules']:continue
        eligible.append(e)
        for left,right in event_windows(event,c['priors']['window_bp']):
            sequence=encoding[np.frombuffer(fasta.fetch(event['chromosome'],left,right),dtype=np.uint8)]
            right=left+len(sequence)
            if event['strand']=='-':
                sequence=sequence[::-1].copy();good=sequence<4;sequence[good]=3-sequence[good]
            windows.append((e,left,right));seqs.append(sequence);starts.append(pos);pos+=len(sequence);ends.append(pos)
    fasta.file.close()
    seq=np.concatenate(seqs) if seqs else np.zeros(0,dtype=np.uint8)
    anchors=collections.defaultdict(list)
    for mid,matrix in mat.items():
        weights,cut=pwm_threshold(matrix,c['priors']['motif_pvalue'])
        scores,locations=scan(seq,np.asarray(starts,np.int64),np.asarray(ends,np.int64),weights)
        for j in np.flatnonzero(scores>=cut):
            e,left,right=windows[j];event=data['events'][e];loc=int(locations[j]);length=len(matrix)
            start=left+loc if event['strand']=='+' else right-loc-length
            for fid,source,sid,ev in links[mid]:
                anchors[(e,fid)].append({'type':'motif','matrix_id':mid,'source':source,'source_motif_id':sid,
                    'evidence':ev,'assembly':assembly,'start0':start,'end0':start+length,'score_discrete':int(scores[j]),
                    'per_position_pvalue_threshold':c['priors']['motif_pvalue'],'candidate_not_binding':True})
    # Experiment-specific peak evidence in the native event assembly; no coordinate mixing.
    old_to_new=identities
    exp=db.execute('''SELECT experiment_id,factor_id,cell_tissue,assay,source_accession FROM candidate_nonperturbed_binding_experiments''').fetchall()
    exp={r[0]:(old_to_new.get(r[1],r[1]),*r[2:]) for r in exp if old_to_new.get(r[1],r[1]) in universe}
    for chrom in sorted({e['chromosome'] for e in data['events']}):
        rows=db.execute('''SELECT start0,end0,strand,experiment_id,peak_id FROM clip_peaks_original
              WHERE chromosome=? AND experiment_id IN (SELECT unnest(?)) ORDER BY start0''',[chrom,list(exp)]).fetchall()
        if not rows:continue
        begin=np.array([r[0] for r in rows]);end=np.array([r[1] for r in rows])
        maxlen=int(np.max(end-begin))
        for e,left,right in windows:
            event=data['events'][e]
            if event['chromosome']!=chrom:continue
            lo=np.searchsorted(begin,left-maxlen,side='left');hi=np.searchsorted(begin,right,side='left')
            for j in range(lo,hi):
                start,stop,strand,eid,pid=rows[j]
                if stop<=left or strand!=event['strand']:continue
                f,tissue,assay,accession=exp[eid]
                anchors[(e,f)].append({'type':'CLIP','assembly':assembly,'start0':start,'end0':stop,'experiment_id':eid,
                    'peak_id':pid,'tissue':tissue,'assay':assay,'accession':accession,'WT_review':'keyword_screen_only'})
    # Normalize evidence once; generate candidate rows one event at a time.
    by_event=collections.defaultdict(dict)
    for (e,f),records in anchors.items():
        by_event[e][f]=list({json.dumps(r,sort_keys=True):r for r in records}.values())
    del anchors
    neighbors=collections.defaultdict(list)
    for (a,b),ev in selected_pairs.items():
        neighbors[a].append((b,ev));neighbors[b].append((a,ev))
    from collections.abc import Mapping
    class GeneratedEvents(Mapping):
        def __len__(self):return len(by_event)
        def __iter__(self):return (str(e) for e in sorted(by_event))
        def __getitem__(self,key):
            e=int(key)
            if e not in by_event:raise KeyError(key)
            candidates={}
            for f,records in by_event[e].items():
                candidates[f]={'factor_id':f,'gene_id':factors[f][0],'symbol':factors[f][1],
                    'route':'direct_anchor_candidate','anchors':records,'paths':[],
                    'priority':3 if any(r['type']=='CLIP' for r in records) else 2}
            for partner,records in by_event[e].items():
                for recipient,ev in neighbors[partner]:
                    if recipient not in candidates:
                        candidates[recipient]={'factor_id':recipient,'gene_id':factors[recipient][0],
                            'symbol':factors[recipient][1],'route':'indirect_partner_candidate','anchors':[],'paths':[],'priority':1}
                    candidates[recipient]['paths'].append({'partner_factor_id':partner,'relation_evidence':ev,
                        'anchor_factor_id':partner,'anchor_count':len(records),'not_recipient_motif':True})
            return list(candidates.values())
    generated=GeneratedEvents()
    result={'factor_ids':factor_ids,'gene_ids':gids,'symbols':[factors[f][1] for f in factor_ids],
            'seed_factor_ids':seeds,'pairs':[{'members':list(k),'evidence':v} for k,v in selected_pairs.items()],
            'events':generated if output is not None else dict(generated.items()),
            'eligible_events':eligible,'matrices_scanned':len(mat),'assembly':assembly,
            'motif_threshold_interpretation':'per-position per-motif under iid uniform background; not FDR',
            'source':'native '+assembly+' rescanning plus native CLIP',
            'RNA_expression_proxy_not_measured_activity':True,
            'factor_retrieval_depth':depth,'partner_hops':c['priors'].get('partner_hops',1)}
    db.close()
    if output is not None:
        from .storage import save_priors,load_priors
        save_priors(output,result)
        result=load_priors(output)
    return result,X
