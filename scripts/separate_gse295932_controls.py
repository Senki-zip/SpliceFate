"""Separate explicit WT samples and barcode-matched non-targeting controls.

Source files remain untouched. No unlabeled cell is assumed WT. The matched
ONT/PacBio WT libraries are kept separate and are not summed as independent cells.
"""
from pathlib import Path
import collections
import csv
import gzip
import hashlib
import json
import time

ROOT=Path(__file__).resolve().parents[1]/'test/GSE295932'


def sha(path):
    with path.open('rb') as stream:return hashlib.file_digest(stream,'sha256').hexdigest()


def write_cells(path,barcodes,sample,group,observed=None):
    with path.open('w') as out:
        writer=csv.writer(out,delimiter='\t')
        writer.writerow(['cell','barcode','sample','condition','observed_count_records'])
        for b in sorted(barcodes):
            writer.writerow([sample+':'+b,b,sample,group,observed[b] if observed is not None else ''])


def main():
    validation=json.loads((ROOT/'validation.json').read_text())
    if not validation['passed']:raise ValueError('Source download validation failed')
    out=ROOT/'WT';out.mkdir(exist_ok=True)
    status={'status':'running','started_at':time.strftime('%Y-%m-%dT%H:%M:%S%z')}
    (out/'status.json').write_text(json.dumps(status,indent=2))
    reports=[]
    for sample,table in [('GSM8961552','GSM8961552_FT_CRISPR_barcode.tsv.gz'),
                         ('GSM8961554','GSM8961554_FT_stim_CRISPR_barcode.tsv.gz')]:
        with gzip.open(ROOT/'files'/table,'rt') as stream:
            assignment=list(csv.DictReader(stream,delimiter='\t'))
        labels={r['barcode']:r['target'] for r in assignment}
        if len(labels)!=len(assignment):raise ValueError('Duplicate barcode labels')
        wanted={b for b,target in labels.items() if target=='Nontarget'}
        samples=[('GSM8961552','GSM8961552_FT_sc_iso_count.txt.gz'),
                 ('GSM8961553','GSM8961553_FT_target_sc_iso_count.txt.gz')] if sample=='GSM8961552' else [
                 ('GSM8961554','GSM8961554_FT_stim_sc_iso_count.txt.gz'),
                 ('GSM8961555','GSM8961555_FT_stim_target_sc_iso_count.txt.gz')]
        for accession,name in samples:
            folder=out/'Nontarget_controls'/accession;folder.mkdir(parents=True,exist_ok=True)
            source=ROOT/'files'/name;dest=folder/(name.removesuffix('.txt.gz')+'.Nontarget.txt.gz')
            temporary=dest.with_name(dest.name+'.part');observed=collections.Counter();all_cells=set()
            rows=selected=0;selected_count=0.
            with gzip.open(source,'rt') as stream,gzip.open(temporary,'wt') as target:
                header=next(stream);fields=header.strip().replace('"','').split()
                if fields[0]!='cell':raise ValueError('Expected first field cell')
                target.write(header);count_index=fields.index('count')
                for line in stream:
                    if not line.strip():continue
                    rows+=1;tokens=line.split();barcode=tokens[0].strip('"');all_cells.add(barcode)
                    if barcode in wanted:
                        target.write(line);observed[barcode]+=1;selected+=1
                        selected_count+=float(tokens[count_index].strip('"'))
            temporary.replace(dest)
            # Reopen fully to check CRC, selection identity and record counts.
            verified=collections.Counter()
            with gzip.open(dest,'rt') as stream:
                next(stream)
                for line in stream:
                    barcode=line.split()[0].strip('"')
                    if barcode not in wanted:raise ValueError('KO/unassigned barcode in control export')
                    verified[barcode]+=1
            if observed!=verified:raise ValueError('Filtered counts changed during export')
            write_cells(folder/'cells.tsv',observed,accession,'Nontarget',observed)
            write_cells(folder/'assigned_but_not_observed.tsv',wanted-set(observed),accession,'Nontarget')
            report={'sample':accession,'group':'Nontarget','barcode_label_sample':sample,
                    'stimulated':sample=='GSM8961554','source':str(source.relative_to(ROOT)),
                    'output':str(dest.relative_to(ROOT)),'source_rows':rows,'source_observed_cells':len(all_cells),
                    'assigned_Nontarget_cells':len(wanted),'observed_Nontarget_cells':len(observed),
                    'selected_rows':selected,'sum_source_count_field':selected_count,
                    'assigned_not_observed':len(wanted-set(observed)),
                    'source_sha256':sha(source),'output_sha256':sha(dest),'gzip_CRC_passed':True}
            reports.append(report);print(json.dumps(report),flush=True)
    # Explicit wild-type ONT, kept in original format via a relative symlink.
    name='GSM8961556_Jurkat_UnconjGSP1_sc_iso_count.txt.gz';source=ROOT/'files'/name
    folder=out/'explicit_WT/GSM8961556_ONT';folder.mkdir(parents=True,exist_ok=True)
    observed=collections.Counter();rows=0
    with gzip.open(source,'rt') as stream:
        reader=csv.DictReader(stream,delimiter='\t')
        for row in reader:observed[row['cell']]+=1;rows+=1
    link=folder/name
    if not link.exists():link.symlink_to(Path('../../../files')/name)
    if link.resolve()!=source.resolve():raise ValueError('WT link resolves to unexpected source')
    write_cells(folder/'cells.tsv',observed,'GSM8961556','WT',observed)
    reports.append({'sample':'GSM8961556','group':'explicit_WT','platform':'ONT',
                    'observed_cells':len(observed),'source_rows':rows,'source':str(source.relative_to(ROOT)),
                    'output':str(link.relative_to(ROOT)),'source_sha256':sha(source)})
    # PacBio matrix: use exactly its barcode axis and report empty columns.
    from scipy.io import mmread
    import numpy as np
    name='GSM8961557_Jurkat_UnconjGSP1.hifi.isoseq_barcodes.tsv.gz'
    with gzip.open(ROOT/'files'/name,'rt') as stream:barcodes=[line.strip() for line in stream if line.strip()]
    if len(set(barcodes))!=len(barcodes):raise ValueError('Duplicate PacBio barcodes')
    matrixfile=ROOT/'files/GSM8961557_Jurkat_UnconjGSP1.hifi.isoseq_matrix.mtx.gz'
    with gzip.open(matrixfile,'rb') as stream:matrix=mmread(stream).tocsc()
    if matrix.shape[1]!=len(barcodes):raise ValueError('PacBio barcode axis mismatch')
    valid=np.asarray(matrix.sum(0)).ravel()>0
    folder=out/'explicit_WT/GSM8961557_PacBio';folder.mkdir(parents=True,exist_ok=True)
    for p in sorted((ROOT/'files').glob('GSM8961557_*')):
        link=folder/p.name
        if not link.exists():link.symlink_to(Path('../../../files')/p.name)
        if link.resolve()!=p.resolve():raise ValueError('PacBio WT link mismatch')
    write_cells(folder/'cells.tsv',barcodes,'GSM8961557','WT')
    ont=set(observed);pac_normalized={b.removesuffix('-1') for b in barcodes}
    reports.append({'sample':'GSM8961557','group':'explicit_WT','platform':'PacBio','matrix_cells':len(barcodes),
                    'nonzero_matrix_cells':int(valid.sum()),'empty_matrix_cells':int((~valid).sum()),
                    'matrix_shape':list(matrix.shape),'nonzero_entries':matrix.nnz,
                    'same_barcode_overlap_with_explicit_WT_ONT':len(ont&pac_normalized),
                    'source_matrix_sha256':sha(matrixfile)})
    summary={'status':'complete','finished_at':time.strftime('%Y-%m-%dT%H:%M:%S%z'),
             'definition':'Explicit author-labeled wild-type libraries separate from Nontarget CRISPR controls',
             'human_hg38':True,'source_files_modified':False,'samples':reports,
             'do_not_sum_matched_modalities':'ONT/PacBio and full-transcriptome/target-enriched libraries can cover the same cells; counts are per library, not independent totals',
             'no_unassigned_as_WT':True,'no_control_KO_mixture':True}
    (out/'summary.json').write_text(json.dumps(summary,indent=2))
    (out/'status.json').write_text(json.dumps({'status':'complete','finished_at':summary['finished_at']},indent=2))
    print(json.dumps(summary,indent=2),flush=True)


if __name__=='__main__':main()
