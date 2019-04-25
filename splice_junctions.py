#!/usr/bin/env python3
import sys
from operator import itemgetter

from utils import reverse_complement, parse_rest, parse_fasta, merge_intervals, collapse_N

"""
Usage:
python3 splice_junctions.py annotation.gtf reads.fasta
"""

def parse_gff(path, file_format):
    # Who the hell came up with this stupid file format?
    genes = {}
    rest = []
    with open(path, 'r') as fh:
        for line in fh:
            # Skip comments
            if line.startswith('#'):
                continue
            data = line.split('\t')

            fields = {
                'scaffold': data[0],
                'feature': data[2],
                'start': data[3],
                'end': data[4],
                'strand': data[6],
            }
            if fields['feature'] in ['exon', 'UTR', 'stop_codon', 'start_codon']:
                fields['rest'] = parse_rest(data[8], file_format)
                gene = fields['rest']['gene_id']
                if gene not in genes:
                    genes[gene] = {
                            'name': gene,
                            'exons': [fields],
                            'scaffold': fields['scaffold'],
                            'strand': fields['strand']
                    }
                else:
                    genes[gene]['exons'].append(fields)
    return genes

def find_splice_junctions(gff_path, fasta_path, file_format, od='.', k=31):
    genes = parse_gff(gff_path, file_format)
    print('gff parsed')
    scaffolds = parse_fasta(fasta_path)
    print('scaffolds parsed')
    starts = []
    ends = []

    wrong_scaffolds = 0
    for gene, data in genes.items():
        starts = []
        ends = []
        intervals = []
        if data['scaffold'] in scaffolds:
            sequence = scaffolds[data['scaffold']]
        else:
            print(f'{data["scaffold"]} not in FASTA file {fasta_path}')
            wrong_scaffolds += 1
            continue
        for exon in data['exons']:
            intervals.append((int(exon['start']) - k, int(exon['end']) + k - 2))

        intervals = merge_intervals(intervals)
        sjs = []
        for iv in intervals:
            if iv[1] - iv[0] < 3 * k - 3:
                # If the length of the exon is < k-1
                sjs.append(iv)
            else:
                sjs.append((iv[0], iv[0] + 2 * k - 2))
                sjs.append((iv[1] - (2 * k) + 2, iv[1]))

        splice_junctions = [
            collapse_N(sequence[iv[0]:iv[1]].upper()) for iv in sjs
        ]

        if data['strand'] == '-':
            splice_junctions = map(reverse_complement, splice_junctions)
        with open(f'{od}/splice_junctions.fasta', 'a') as fh:
            for i, sj in enumerate(splice_junctions):
                fh.write(f'>{gene}:{i}\n')
                fh.write(f'{sj}\n')
    print(f'{wrong_scaffolds} scaffolds not found')

if __name__ == '__main__':
    gff_path = sys.argv[1]
    file_format = gff_path[gff_path.rfind('.'):]
    if file_format not in ['.gff', '.gtf', '.gff3']:
        print(f'File format {file_format} not recognized')
    else:
        find_splice_junctions(sys.argv[1], sys.argv[2], file_format)
