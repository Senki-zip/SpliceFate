"""Event-addressable relational artifacts and memory-mapped numerical arrays."""
from __future__ import annotations

from collections.abc import Mapping
import json
from pathlib import Path
import struct
import zipfile

import duckdb
import numpy as np
from scipy import sparse

from .io import dump


def compact_candidate(row):
    return {key: row[key] for key in ('factor_id','gene_id','symbol','route','priority',
        'anchor_count','path_count','partner_ids','WT_response_score','selection_rank',
        'WT_network_weight','selection_score','fit_gain','oob_fit_gain','stability') if key in row}


class PriorEvents(Mapping):
    """One event at a time; connections are local to a call/process."""
    def __init__(self, path):
        self.path = str(path)
        with duckdb.connect(self.path, read_only=True,config={'threads':1}) as db:
            self.keys_ = [str(r[0]) for r in db.execute('SELECT DISTINCT event_index FROM event_candidates ORDER BY event_index').fetchall()]

        self.keyset=set(self.keys_)

    def __contains__(self,event):return str(event) in self.keyset
    def __iter__(self): return iter(self.keys_)
    def __len__(self): return len(self.keys_)
    def __getitem__(self, event):
        with duckdb.connect(self.path, read_only=True,config={'threads':1}) as db:
            rows = db.execute('''SELECT factor_id,gene_id,symbol,route,priority,anchor_count,partner_ids_json,path_count
                FROM event_candidates WHERE event_index=? ORDER BY factor_id''',[int(event)]).fetchall()
        if not rows: raise KeyError(str(event))
        return [dict(zip(['factor_id','gene_id','symbol','route','priority','anchor_count','partner_ids','path_count'],
                         (*r[:6],json.loads(r[6]),r[7]))) for r in rows]

    def evidence(self,event,factor):
        with duckdb.connect(self.path,read_only=True,config={'threads':1}) as db:
            anchor=db.execute('SELECT anchors_json FROM event_anchors WHERE event_index=? AND factor_id=?',[int(event),factor]).fetchone()
            candidate=next(r for r in self[event] if r['factor_id']==factor)
            paths=[]
            for partner in candidate['partner_ids']:
                a,b=sorted([factor,partner])
                evidence=db.execute('SELECT evidence_json FROM pair_evidence WHERE factor_a=? AND factor_b=?',[a,b]).fetchone()
                count=db.execute('SELECT anchor_count FROM event_anchors WHERE event_index=? AND factor_id=?',[int(event),partner]).fetchone()
                paths.append({'partner_factor_id':partner,'anchor_factor_id':partner,
                    'relation_evidence':json.loads(evidence[0]),'anchor_count':count[0], 'not_recipient_motif':True})
        return {**candidate,'anchors':json.loads(anchor[0]) if anchor else [],'paths':paths}


def save_priors(root, prior):
    root=Path(root);path=root/'priors.duckdb'
    db=duckdb.connect(str(path))
    try:
        db.execute('CREATE OR REPLACE TABLE pair_evidence(pair_id INTEGER,factor_a VARCHAR,factor_b VARCHAR,max_rank DOUBLE,evidence_json VARCHAR)')
        db.execute('CREATE OR REPLACE TABLE event_anchors(event_index INTEGER,factor_id VARCHAR,anchor_count INTEGER,has_clip BOOLEAN,anchors_json VARCHAR)')
        db.execute('CREATE OR REPLACE TABLE event_candidates(event_index INTEGER,factor_id VARCHAR,gene_id VARCHAR,symbol VARCHAR,route VARCHAR,priority INTEGER,anchor_count INTEGER,partner_ids_json VARCHAR,path_count INTEGER)')
        import pandas as pd
        def insert(table,rows):
            if not rows:return
            cols=[r[0] for r in db.execute('DESCRIBE '+table).fetchall()]
            batch=pd.DataFrame(rows,columns=cols);db.register('prior_batch',batch)
            db.execute('INSERT INTO '+table+' SELECT * FROM prior_batch');db.unregister('prior_batch')
        insert('pair_evidence',[(i,*p['members'],max((e['rank'] for e in p.get('evidence',[])),default=0),json.dumps(p.get('evidence',[]))) for i,p in enumerate(prior['pairs'])])
        rows=[];anchors=[]
        for event,candidates in prior['events'].items():
            for r in candidates:
                a=r.get('anchors',[]);partners=sorted({p['partner_factor_id'] for p in r.get('paths',[])})
                rows.append((int(event),r['factor_id'],r.get('gene_id'),r.get('symbol',r['factor_id']),r.get('route'),r.get('priority',0),len(a),json.dumps(partners),len(r.get('paths',[]))))
                if a:anchors.append((int(event),r['factor_id'],len(a),any(x['type']=='CLIP' for x in a),json.dumps(a)))
            if len(rows)>=50000:
                insert('event_candidates',rows);insert('event_anchors',anchors);rows=[];anchors=[]
        insert('event_candidates',rows);insert('event_anchors',anchors)
        db.execute('CREATE INDEX IF NOT EXISTS candidates_event ON event_candidates(event_index)')
        db.execute('CREATE INDEX IF NOT EXISTS anchors_event ON event_anchors(event_index,factor_id)')
        db.execute('CHECKPOINT')
    finally:db.close()
    dump(root/'priors_meta.json',{k:v for k,v in prior.items() if k not in ('events','pairs')})


def load_priors(root):
    root=Path(root);prior=json.loads((root/'priors_meta.json').read_text())
    prior['events']=PriorEvents(root/'priors.duckdb')
    with duckdb.connect(str(root/'priors.duckdb'),read_only=True,config={'threads':1}) as db:
        prior['pairs']=[{'members':[a,b],'max_rank':rank} for a,b,rank in db.execute('SELECT factor_a,factor_b,max_rank FROM pair_evidence').fetchall()]
    return prior


def mmap_npz(path, name):
    """Map a ZIP_STORED .npy member without decompressing/copying the archive."""
    with zipfile.ZipFile(path) as archive:
        info=archive.getinfo(name+'.npy')
        if info.compress_type!=zipfile.ZIP_STORED:raise ValueError('Model arrays must use uncompressed NPZ')
        with open(path,'rb') as f:
            f.seek(info.header_offset+26);n,e=struct.unpack('<HH',f.read(4))
            f.seek(info.header_offset+30+n+e)
            version=np.lib.format.read_magic(f)
            reader=np.lib.format.read_array_header_1_0 if version==(1,0) else np.lib.format.read_array_header_2_0
            shape,fortran,dtype=reader(f)
            offset=f.tell()
    if not np.prod(shape):return np.empty(shape,dtype=dtype)
    return np.memmap(path,mode='r',dtype=dtype,offset=offset,shape=shape,order='F' if fortran else 'C')


class ModelWriter:
    """Append each joint fit immediately; pack arrays on successful completion."""
    def __init__(self,root):
        self.root=Path(root);self.root.mkdir(parents=True,exist_ok=True)
        self.db=duckdb.connect(str(self.root/'models.duckdb'))
        self.db.execute('CREATE OR REPLACE TABLE models(event_index INTEGER PRIMARY KEY,metadata_json VARCHAR,wt_offset BIGINT,wt_rows INTEGER,wt_cols INTEGER,support_offset BIGINT,support_rows INTEGER,support_cols INTEGER)')
        self.wt=(self.root/'.wt.bin').open('wb');self.support=(self.root/'.support.bin').open('wb')
        self.wsize=0;self.ssize=0

    def append(self,e,model):
        m=dict(model);wt=np.asarray(m.pop('WT_prediction',[]),dtype=np.float64)
        sup=dict(m.get('support',{}));ref=np.asarray(sup.pop('reference_inputs',[]),dtype=np.float64)
        if 'support' in m:m['support']=sup
        m.pop('observed_denominators',None)
        m['candidates']=[compact_candidate(r) for r in m.get('candidates',[])]
        if wt.ndim!=2:wt=np.empty((0,0))
        if ref.ndim!=2:ref=np.empty((0,0))
        self.db.execute('INSERT INTO models VALUES (?,?,?,?,?,?,?,?)',[int(e),json.dumps(m,allow_nan=False),self.wsize,*wt.shape,self.ssize,*ref.shape])
        wt.tofile(self.wt);ref.tofile(self.support);self.wsize+=wt.size;self.ssize+=ref.size

    def finish(self):
        self.wt.close();self.support.close();self.db.execute('CHECKPOINT');self.db.close()
        wt=np.memmap(self.root/'.wt.bin',mode='r',dtype='float64',shape=(self.wsize,)) if self.wsize else np.empty(0)
        ref=np.memmap(self.root/'.support.bin',mode='r',dtype='float64',shape=(self.ssize,)) if self.ssize else np.empty(0)
        np.savez(self.root/'model_arrays.npz',wt=wt,support=ref)
        del wt,ref
        (self.root/'.wt.bin').unlink();(self.root/'.support.bin').unlink()
        return ModelStore(self.root)

    def close(self):
        self.wt.close();self.support.close();self.db.close()


class ModelStore(Mapping):
    def __init__(self,root,data=None):
        self.root=Path(root);self.data=data
        self.wt=mmap_npz(self.root/'model_arrays.npz','wt')
        self.support=mmap_npz(self.root/'model_arrays.npz','support')
        with duckdb.connect(str(self.root/'models.duckdb'),read_only=True) as db:
            self.rows={str(r[0]):r[1:] for r in db.execute('SELECT * FROM models ORDER BY event_index').fetchall()}

    def __iter__(self):return iter(self.rows)
    def __len__(self):return len(self.rows)
    def __getitem__(self,e):
        raw,wo,wr,wc,so,sr,sc=self.rows[str(e)];m=json.loads(raw)
        if wr:
            m['WT_prediction']=self.wt[wo:wo+wr*wc].reshape(wr,wc)
            m['support']['reference_inputs']=self.support[so:so+sr*sc].reshape(sr,sc)
        if self.data is not None:
            i=int(e);m['observed_denominators']=np.asarray(self.data['C'][:,self.data['offsets'][i]:self.data['offsets'][i+1]].sum(1)).ravel()
        return m


def save_structures(root,structures,n_cells):
    meta={};matrices=[];offset=0
    for gene,s in structures.items():
        m=dict(s);counts=sparse.csr_matrix(m.pop('counts'))
        m['column_start']=offset;m['column_end']=offset+counts.shape[1]
        offset+=counts.shape[1];matrices.append(counts);meta[gene]=m
    matrix=sparse.hstack(matrices,format='csr') if matrices else sparse.csr_matrix((n_cells,0))
    sparse.save_npz(Path(root)/'structure_counts.npz',matrix)
    dump(Path(root)/'structures_meta.json',meta)


class StructureStore(Mapping):
    def __init__(self,root):
        root=Path(root);self.meta=json.loads((root/'structures_meta.json').read_text())
        self.counts=sparse.load_npz(root/'structure_counts.npz').tocsr()
    def __len__(self):return len(self.meta)
    def __iter__(self):return iter(self.meta)
    def __getitem__(self,gene):
        m=dict(self.meta[gene]);m['counts']=self.counts[:,m['column_start']:m['column_end']].toarray()
        return m
