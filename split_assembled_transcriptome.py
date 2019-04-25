#!/usr/bin/env python3
import sys
from operator import itemgetter


from utils import reverse_complement, parse_rest, parse_fasta, merge_intervals, collapse_N

"""
Usage:
python3 split_assembled_transcriptome.py annotation.gtf3 reads.fasta
"""


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
                                'strand': fields['strand']
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

def split_assembled_genome(gff_path, fasta_path, file_format, od='.'):
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
            intervals.append((int(exon['start']) - 1, int(exon['end']) -1))
            starts.append(int(exon['start']) - 1)
            ends.append(int(exon['end']) - 1)

        intervals = list(merge_intervals(intervals))

        proc = ''.join(map(lambda iv: sequence[iv[0]:iv[1]], intervals))

        unproc = sequence[min(starts):max(ends)]
        proc = collapse_N(proc).upper()
        unproc = collapse_N(unproc).upper()
        if data['strand'] == '-':
            proc= reverse_complement(proc)
            unproc= reverse_complement(unproc)
        with open(f'{od}/processed.fasta', 'a') as fh:
            fh.write(f'>{gene}\n')
            fh.write('\n'.join([proc[i:i+80] for i in range(0, len(proc), 80)]))
            fh.write('\n')
        with open(f'{od}/unprocessed.fasta', 'a') as fh:
            fh.write(f'>{gene}\n')
            fh.write('\n'.join([unproc[i:i+80] for i in range(0, len(unproc), 80)]))
            fh.write('\n')
    print(f'{wrong_scaffolds} scaffolds not found')

if __name__ == '__main__':
    gff_path = sys.argv[1]
    file_format = gff_path[gff_path.rfind('.'):]
    if file_format not in ['.gff', '.gtf', '.gff3']:
        print(f'File format {file_format} not recognized')
    else:
        split_assembled_genome(sys.argv[1], sys.argv[2], file_format)
