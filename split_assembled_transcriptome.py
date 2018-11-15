#!/usr/bin/env python3
import sys
from operator import itemgetter

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

    # Deal in transcripts in GTF.
    # TODO: Deal in transcripts in GFF3
    trs = {}
    tr2g = []
    rest = []
    with open(path, 'r') as fh:
        for line in fh:
            data = line.split('\t')

            fields = {
                'scaffold': data[0],
                'feature': data[2],
                'start': data[3],
                'end': data[4],
                'strand': data[6],
            }
            if file_format == '.gff3':
                fields['rest'] = parse_rest(data[8].split(' ')[0], file_format)
                fields['rest']['Name'] = fields['rest']['Name'].rstrip('\n')

                if fields['feature'] == 'gene':
                    genes[fields['rest']['Name']] = fields
                    genes[fields['rest']['Name']]['exons'] = []
                else:
                    rest.append(fields)
            elif file_format == '.gtf':
                if fields['feature'] in ['exon', 'UTR', 'stop_codon']:
                    # fields['rest'] = parse_rest(data[8], file_format)
                    # gene = fields['rest']['gene_id']
                    # if gene not in genes:
                        # genes[gene] = {'name': gene, 'exons': [fields]}
                        # genes[gene]['scaffold'] = fields['scaffold']
                        # genes[gene]['strand'] = fields['strand']
                    # else:
                        # genes[gene]['exons'].append(fields)
                    fields['rest'] = parse_rest(data[8], file_format)
                    tr = fields['rest']['transcript_id']
                    if tr not in trs:
                        trs[tr] = {
                                'name': tr,
                                'exons': [fields],
                                'scaffold': fields['scaffold'],
                                'strant': fileds['strand']
                        }
                        tr2g.append(f'{tr}\t{fields["rest"]["gene_id"]}')
                    else:
                        trs[tr]['exons'].append(fields)


    for line in rest:
        if line['feature'] == 'exon':
            genes[line['rest']['Target']]['exons'].append(line)

    with open('exons_per_gene.csv', 'w') as fh:
        fh.write('\n'.join(map(lambda v: f'{v[0]}\t{len(v[1]["exons"])}', genes.items())))
        # for key, value in genes.items():
            # print(f'{key}\t{len(value["exons"])}', file=fh)

    with open('tr2g', 'w') as fh:
        fh.write('\n'.join(tr2g))

    if file_format == '.gff3':
        return genes
    elif file_format == '.gtf':
        return trs

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
                current = line[1:]
            else:
                sequences.append(line)

        scaffolds[current] = ''.join(sequences)
    return scaffolds

def merge_intervals(intervals):
    # TODO:
    # Confirm with Lior or Jase whether this is the correct thing to do
    sorted_intervals = sorted(intervals, key=itemgetter(0))

    if not sorted_intervals:  # no intervals to merge
        return

    # low and high represent the bounds of the current run of merges
    low, high = sorted_intervals[0]

    for iv in sorted_intervals[1:]:
        if iv[0] <= high:  # new interval overlaps current run
            high = max(high, iv[1])  # merge with the current run
        else:  # current run is over
            yield low, high  # yield accumulated interval
            low, high = iv  # start new run

    yield low, high  # end the final run

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


def split_assembled_genome(gff_path, fasta_path, file_format):
    genes = parse_gff(gff_path, file_format)
    print('gff parsed')
    scaffolds = parse_fasta(fasta_path)
    print('scaffolds parsed')
    starts = []
    ends = []
    for gene, data in genes.items():
        starts = []
        ends = []
        intervals = []
        sequence = scaffolds[data['scaffold']]
        for exon in data['exons']:
            intervals.append((int(exon['start']) - 1, int(exon['end']) -1))
            starts.append(int(exon['start']) - 1)
            ends.append(int(exon['end']) - 1)

        intervals = list(merge_intervals(intervals))

        proc = ''.join(map(lambda iv: sequence[iv[0]:iv[1]], intervals))
        print('=============================')
        print(max(ends), min(ends))
        print(max(ends) - min(starts))
        print(intervals)

        unproc = sequence[min(starts):max(ends)]
        proc = collapse_N(proc).upper()
        unproc = collapse_N(unproc).upper()
        if data['strand'] == '-':
            proc= reverse_complement(proc)
            unproc= reverse_complement(unproc)
        with open('processed.fasta', 'a') as fh:
            fh.write(f'>{gene}\n')
            fh.write('\n'.join([proc[i:i+80] for i in range(0, len(proc), 80)]))
            fh.write('\n')
        with open('unprocessed.fasta', 'a') as fh:
            fh.write(f'>{gene}\n')
            fh.write('\n'.join([unproc[i:i+80] for i in range(0, len(unproc), 80)]))
            fh.write('\n')

if __name__ == '__main__':
    gff_path = sys.argv[1]
    file_format = gff_path[gff_path.rfind('.'):]
    if file_format not in ['.gff', '.gtf', '.gff3']:
        print(f'File format {file_format} not recognized')
    else:
        split_assembled_genome(sys.argv[1], sys.argv[2], file_format)
