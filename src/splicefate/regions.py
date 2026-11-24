"""Region incidence from continuous long-read alignment blocks."""
from bisect import bisect_right
import collections
import gzip
import re


def merge_intervals(intervals):
    result = []
    for left, right in sorted(intervals):
        if result and left <= result[-1][1]:
            result[-1] = (result[-1][0], max(right, result[-1][1]))
        else:
            result.append((left, right))
    return result


class Regions:
    def __init__(self, intervals):
        self.intervals = merge_intervals(intervals)
        self.ends = [b for _, b in self.intervals]

    def overlap(self, left, right):
        total = 0
        i = bisect_right(self.ends, left)
        while i < len(self.intervals) and self.intervals[i][0] < right:
            a, b = self.intervals[i]
            total += max(0, min(b, right)-max(a, left))
            i += 1
        return total

    def subtract(self, intervals):
        result = []
        for left, right in merge_intervals(intervals):
            i = bisect_right(self.ends, left)
            cursor = left
            while i < len(self.intervals) and self.intervals[i][0] < right:
                a, b = self.intervals[i]
                if a > cursor:
                    result.append((cursor, min(a, right)))
                cursor = max(cursor, b)
                i += 1
            if cursor < right:
                result.append((cursor, right))
        return result


def reference_regions(gtf):
    """Annotated intron unions minus every annotated exon on either strand."""
    transcripts = collections.defaultdict(list)
    genes = {}
    global_exons = collections.defaultdict(list)
    opener = gzip.open if str(gtf).endswith('.gz') else open
    with opener(gtf, 'rt') as stream:
        for line in stream:
            if line.startswith('#'):
                continue
            row = line.rstrip().split('\t')
            if row[2] != 'exon':
                continue
            attributes = dict(re.findall(r'(\w+) "([^"]*)"', row[8]))
            gene = attributes['gene_id'].split('.')[0]
            transcript = attributes['transcript_id']
            interval = (int(row[3])-1, int(row[4]))
            identity = (row[0], row[6])
            entry = genes.setdefault(gene, {'chromosome':row[0], 'strand':row[6],
                                           'exons':[], 'annotated_introns':set()})
            if (entry['chromosome'], entry['strand']) != identity:
                raise ValueError('Gene spans multiple reference identities: '+gene)
            entry['exons'].append(interval)
            transcripts[gene, transcript].append(interval)
            global_exons[row[0]].append(interval)
    masks = {chrom:Regions(intervals) for chrom, intervals in global_exons.items()}
    for (gene, _), intervals in transcripts.items():
        ordered = merge_intervals(intervals)
        genes[gene]['annotated_introns'].update((a[1], b[0]) for a, b in zip(ordered, ordered[1:]))
    for entry in genes.values():
        entry['exons'] = Regions(entry['exons'])
        entry['annotated_introns'] = sorted(entry['annotated_introns'])
        entry['introns'] = Regions(masks[entry['chromosome']].subtract(entry['annotated_introns']))
    return genes


def classify_blocks(blocks, reference, min_exon_bp=10, min_intron_bp=20, retention_flank=10):
    """Count a cluster once per gene/region type; gaps never supply I coverage."""
    exon_bp = sum(reference['exons'].overlap(a, b) for a, b in blocks)
    intron_bp = sum(reference['introns'].overlap(a, b) for a, b in blocks)
    exonic = exon_bp >= min_exon_bp
    intronic = intron_bp >= min_intron_bp
    retained = intronic and any(
        any(a <= left-retention_flank and b >= right+retention_flank for a, b in blocks)
        for left, right in reference['annotated_introns'])
    # A continuous block spans an annotated exon/intron boundary. Only that
    # block's genuine intronic overlap counts toward the boundary support.
    boundary = intronic and any(
        reference['exons'].overlap(a, b) >= min_exon_bp
        and reference['introns'].overlap(a, b) >= min_intron_bp for a, b in blocks)
    return exonic, intronic, boundary, retained, exon_bp, intron_bp
