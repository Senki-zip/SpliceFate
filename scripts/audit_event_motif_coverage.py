"""Descriptive event-proximal motif coverage with a shuffled sequence control.

This reports sequence matches, not in vivo binding, regulatory signs or KO effects.
"""
from pathlib import Path
import csv
import json
import sqlite3
import collections
import numpy as np
from numba import njit

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'data/external/rbp_atlas_20261115/integrated'
FASTA = ROOT.parent / 'reference_documents/refdata-cellranger-arc-mm10-2020-A-2.0.0/fasta/genome.fa'


class Fasta:
    def __init__(self, path):
        self.file = path.open('rb')
        self.index = {}
        for line in Path(str(path)+'.fai').open():
            chrom, length, offset, bases, width, *_ = line.split()
            self.index[chrom] = tuple(map(int, (length,offset,bases,width)))

    def fetch(self, chrom, start, end):
        length, offset, bases, width = self.index[chrom]
        start, end = max(0,start), min(length,end)
        left = offset + start//bases*width + start%bases
        right = offset + (end-1)//bases*width + (end-1)%bases+1
        self.file.seek(left)
        return self.file.read(right-left).replace(b'\n',b'').replace(b'\r',b'').upper()


def threshold(matrix, pvalue, resolution=.05):
    p = np.maximum(matrix, 1e-6); p /= p.sum(axis=1,keepdims=True)
    weights = np.rint(np.log2(p/.25)/resolution).astype(np.int32)
    dist, lower = np.ones(1), 0
    for row in weights:
        low, high = int(row.min()), int(row.max())
        updated = np.zeros(len(dist)+high-low)
        for score in row:
            shift = int(score)-low
            updated[shift:shift+len(dist)] += dist*.25
        dist, lower = updated, lower+low
    tails = np.cumsum(dist[::-1])[::-1]
    eligible = np.flatnonzero(tails <= pvalue)
    cutoff = lower+int(eligible[0]) if len(eligible) else lower+len(dist)
    return weights, cutoff


@njit
def best_scores(sequence, starts, ends, event_ids, weights, event_count):
    result = np.full(event_count, -2000000000, dtype=np.int32)
    length = len(weights)
    for fragment in range(len(starts)):
        event = event_ids[fragment]
        for start in range(starts[fragment], ends[fragment]-length+1):
            score = 0
            for k in range(length):
                base = sequence[start+k]
                if base>3:
                    score = -2000000000
                    break
                score += weights[k,base]
            if score>result[event]:result[event] = score
    return result


def main():
    db=sqlite3.connect((OUT/'rbp_atlas.sqlite').as_uri()+'?mode=ro',uri=True)
    links=db.execute("""SELECT DISTINCT v.motif_id,v.factor_id,f.gene_symbol,m.matrix_json
        FROM valid_factor_motifs v JOIN mouse_candidate_splicing_factors f ON f.factor_id=v.factor_id
        JOIN motifs m ON v.motif_id=m.motif_id
        WHERE v.source IN ('CISBP-RNA','ATtRACT')
        AND v.evidence!='mutant_binding_exclude_from_WT_default'""").fetchall()
    matrices={r[0]:json.loads(r[3]) for r in links}
    events=json.loads((ROOT/'data/processed/GSE250381_local_events_v8_20261108/groups.json').read_text())
    fasta=Fasta(FASTA);rng=np.random.default_rng(0)
    encoding=np.full(256,4,dtype=np.uint8)
    for k,base in enumerate(b'ACGT'):encoding[base]=k
    sequences=[];controls=[];starts=[];ends=[];eventids=[];lengths=[];pos=0
    for eventid,event in enumerate(events):
        boundaries=sorted({coordinate for branch in event['branches'] for junction in branch for coordinate in junction})
        intervals=[]
        for boundary in boundaries:
            left,right=max(0,boundary-100),boundary+100
            if intervals and left<=intervals[-1][1]:intervals[-1][1]=max(right,intervals[-1][1])
            else:intervals.append([left,right])
        event_length=0
        for left,right in intervals:
            raw=fasta.fetch(event['chromosome'],left,right)
            seq=encoding[np.frombuffer(raw,dtype=np.uint8)]
            if event['strand']=='-':
                seq=seq[::-1].copy();valid=seq<4;seq[valid]=3-seq[valid]
            shuffled=seq.copy();valid=shuffled<4;shuffled[valid]=rng.permutation(shuffled[valid])
            sequences.append(seq);controls.append(shuffled);starts.append(pos);ends.append(pos+len(seq));eventids.append(eventid)
            pos+=len(seq);event_length+=len(seq)
        lengths.append(event_length)
    sequences=np.concatenate(sequences);controls=np.concatenate(controls)
    starts,ends,eventids=map(lambda x:np.asarray(x,dtype=np.int64),(starts,ends,eventids))
    levels=[1e-4,1e-5,1e-6];motifids=sorted(matrices)
    observed=np.zeros((len(levels),len(motifids),len(events)),dtype=bool);shuffled=np.zeros_like(observed)
    for i,mid in enumerate(motifids):
        matrix=np.asarray(matrices[mid]);weights,_=threshold(matrix,levels[0])
        score=best_scores(sequences,starts,ends,eventids,weights,len(events))
        nullscore=best_scores(controls,starts,ends,eventids,weights,len(events))
        for k,pvalue in enumerate(levels):
            _,cut=threshold(matrix,pvalue);observed[k,i]=score>=cut;shuffled[k,i]=nullscore>=cut
        if (i+1)%50==0:print('Scanned',i+1,'/',len(motifids),flush=True)
    factor_links=collections.defaultdict(set);symbols={}
    motif_index={m:i for i,m in enumerate(motifids)}
    for mid,fid,symbol,_ in links:factor_links[fid].add(motif_index[mid]);symbols[fid]=symbol
    result={'event_count':len(events),'factor_count':len(factor_links),'distinct_matrices':len(motifids),
        'sequence_definition':'Union of +/-100 genomic bases around each observed splice-junction endpoint; RNA strand orientation; no artificial concatenation across fragments',
        'total_scanned_bases_including_repeated_regions_across_events':pos,'median_bases_per_event':float(np.median(lengths)),
        'null':'Independent uniform A/C/G/U at each position; 0.05-bit discretized log-odds; PWM probability floor 1e-6',
        'pvalue_definition':'Per-position per-motif upper-tail score probability under the null; no multiple-testing correction; not a binding probability',
        'control':'One fixed-seed mononucleotide-composition-preserving shuffle per sequence fragment; not dinucleotide matched and not a formal FDR estimate',
        'thresholds':[], 'interpretation':'Motif presence coverage only; sequence matching and its shuffle background do not establish factor-specific occupancy or regulatory effects.'}
    rows=[]
    for k,pvalue in enumerate(levels):
        factor_counts=np.zeros(len(events),dtype=int);null_counts=np.zeros(len(events),dtype=int)
        for fid,indexes in factor_links.items():
            hit=observed[k,list(indexes)].any(axis=0);nullhit=shuffled[k,list(indexes)].any(axis=0)
            factor_counts+=hit;null_counts+=nullhit
            rows.append([pvalue,fid,symbols[fid],int(hit.sum()),int(nullhit.sum())])
        result['thresholds'].append({'per_position_pvalue':pvalue,'events_with_any_motif':int((factor_counts>0).sum()),
            'shuffled_events_with_any_motif':int((null_counts>0).sum()),'median_factor_labels_per_event':float(np.median(factor_counts)),
            'shuffled_median_factor_labels_per_event':float(np.median(null_counts))})
    np.savez_compressed(OUT/'event_motif_scan_hits.npz',observed=observed,shuffled=shuffled,motif_ids=np.asarray(motifids),
        event_ids=np.asarray([e['event'] for e in events]),thresholds=np.asarray(levels))
    (OUT/'event_motif_sequence_coverage_20261115.json').write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    with (OUT/'event_motif_factor_coverage.tsv').open('w') as f:
        w=csv.writer(f,delimiter='\t');w.writerow(['per_position_pvalue','factor_id','symbol','events_with_match','shuffled_events_with_match']);w.writerows(rows)
    print(json.dumps(result,ensure_ascii=False,indent=2),flush=True)


if __name__=='__main__':main()
