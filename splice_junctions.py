#!/usr/bin/env python3
import sys
from operator import itemgetter

"""
Usage:
python3 splice_junctions.py annotation.gtf reads.fasta
"""

def reverse_complement(string):
    comp = {
        'A': 'T',
        'T': 'A',
        'G': 'C',
        'C': 'G',
        'N': 'N',
    }
    return ''.join(list(map(lambda c: comp[c], string))[::-1])

def parse_rest(rest, file_format='.gff3'):
    # Parses the 'Attributes' field of the gff
    if file_format == '.gff3':
        return dict(map(lambda x: x.rstrip('\n').split('='), rest.split(';')))
    elif file_format == '.gtf':
        # Sorry, this is a disgusting file format
        return dict(map(lambda x: x.rstrip('\n')\
                                   .lstrip(' ')\
                                   .replace('"', '')\
                                   .split(' ', 1),
                        rest.rstrip(';').split(';')[:-1]))

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

def parse_fasta(path):
    scaffolds = {}
    sequences = []
    current = ''
    with open(path, 'r') as fh:
        for line in fh:
            line = line.rstrip('\n')
            if line.startswith('>'):
                if len(current) > 0:
                    scaffolds[current] = ''.join(sequences)
                sequences = []
                current = line.split()[0][1:]
            else:
                sequences.append(line)

        scaffolds[current] = ''.join(sequences)
    return scaffolds

def merge_intervals(intervals):
    sorted_intervals = sorted(intervals, key=itemgetter(0))
    if not sorted_intervals:  # no intervals to merge
        return
    low, high = sorted_intervals[0]
    for iv in sorted_intervals[1:]:
        if iv[0] <= high:
            high = max(high, iv[1])
        else:
            yield low, high
            low, high = iv
    yield low, high

def collapse_N(string, k=31):
    # Collapses arbitrary length regions of N-bases into regions of length k
    curr = False
    mod = []
    counter = 0
    for c in string:
        if c == 'N' and curr and counter >= k:
            continue
        elif c == 'N':
            curr = True
            counter += 1
        else:
            curr = False
            counter = 0
        mod.append(c)
    return ''.join(mod)


def find_splice_junctions(gff_path, fasta_path, file_format, k=31):
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
        with open('splice_junctions.fasta', 'a') as fh:
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
