#!/usr/bin/env python3
"""
Author:
    Kristján Eldjárn Hjörleifsson
    keldjarn@caltech.edu

Usage:
python3 split_assembled_transcriptome.py annotation.gtf3 reads.fasta

Parses a gff/gtf file into the following collections
==============================================================================
    Non-ambiguous:
        a) Processed
            a.i)    Exon-exon splice junctions
        b) Unprocessed
            b.i)    Non-retained exon-intron splice junctions
            b.ii)   Intron intersection
    Ambiguous:
        c) Exon intersection
        d) Retained introns (i.e. pairwise exon set difference)
        e) Retained exon-intron splice junctions
==============================================================================

Returns the corresponding sequences from the corresponding fasta scaffolds,
correcting for strandedness.
"""
import sys
from operator import itemgetter
from itertools import chain

from utils import (reverse_complement, parse_rest, parse_fasta,
                   merge_intervals, collapse_N)

def process_gene(gene, scaffolds):
    exons = {
        tr: [(e['start'], e['end']) for e in exons['exons']]\
             for tr, exons in gene['trs'].items()
    }
    #===============================#
    #   b.ii) Intron intersection   #
    #===============================#
    # ivs = list(chain(*exons.values()))
    exon_union = merge_intervals(list(chain(*exons.values())))
    itr_intersect = [(i[0][1], i[1][0]) for i in zip(exon_union[:-1],
                                                     exon_union[1:])]

    #==========================#
    #   c) Exon intersection   #
    #==========================#
    itr_union =



def parse_gff(path, scaffolds, od, file_format):
    # Who the hell came up with this stupid file format?
    gene = {}
    tr2g = []
    with open(path, 'r') as fh:
        prev_gene_id = None
        for line in fh:
            # Skip comments
            if line.startswith('#'):
                continue
            data = line.split('\t')

            fields = {
                'scaffold': data[0],
                'feature': data[2],
                'start': int(data[3]) - 1, # Scaffold is 1-indexed
                'end': int(data[4]) - 1,
                'strand': data[6],
                'rest': parse_rest(data[8], file_format)
            }
            gene_id = fields['rest']['gene_id']

            # Process current gene and start collecting exons for
            # next gene
            if gene_id != prev_gene_id:
                if prev_gene_id is not None:
                    process_gene(gene, scaffolds)
                gene = {
                        'name': gene_id,
                        'trs': {},
                        'scaffold': fields['scaffold'],
                        'strand': fields['strand']
                }
                prev_gene_id = gene_id

            if file_format == '.gtf':
                if fields['feature'] in ['exon', 'UTR', 'stop_codon']:
                    tr = fields['rest']['transcript_id']
                    itr = fields['rest']['gene_id']
                    if tr not in gene['trs']:
                        gene['trs'][tr] = {
                                'name': tr,
                                'exons': [fields]
                        }
                        tr2g.append(f'{tr}\t{gene_id}')
                    else:
                        gene['trs'][tr]['exons'].append(fields)

        # Process last gene in file
        process_gene(gene, scaffolds)

    with open(f'{od}/tr2g', 'w') as fh:
        fh.write('\n'.join(tr2g))

def split_assembled_genome(gff_path, fasta_path, file_format, od='.'):
    scaffolds = parse_fasta(fasta_path)
    print('scaffolds parsed')
    trs, itrs = parse_gff(gff_path, scaffolds, od, file_format)
    print('gff parsed')

    """
    wrong_scaffolds = 0
    for tr, data in trs.items():
        if data['scaffold'] in scaffolds:
            sequence = scaffolds[data['scaffold']]
        else:
            print(f'{data["scaffold"]} not in FASTA file {fasta_path}')
            wrong_scaffolds += 1
            continue
        intervals = [(int(e['start']) - 1, int(e['end']) -1)) for e in data['exons']]
        inv_intervals = [(i[0][1], i[1][0]) for i in zip(intervals[:-1], intervals[1:])]

        intervals = list(merge_intervals(intervals))
        inv_intervals = list(merge_intervals(inv_intervals))

        proc = ''.join(map(lambda iv: sequence[iv[0]:iv[1]], intervals))
        itrs = ''.join(map(lambda iv: sequence[iv[0]:iv[1]], intervals))

        proc = collapse_N(proc.upper())
        itrs = collapse_N(itrs.upper())

        if data['strand'] == '-':
            proc = reverse_complement(proc)
            itrs = reverse_complement(itrs)
        with open(f'{od}/processed.fasta', 'a') as fh:
            fh.write(f'>{tr}\n')
            fh.write('\n'.join([proc[i:i+80] for i in range(0, len(proc), 80)]))
            fh.write('\n')
        with open(f'{od}/intron_union.fasta', 'a') as fh:
            fh.write(f'>{itr}\n')
            fh.write('\n'.join([itrs[i:i+80] for i in range(0, len(itrs), 80)]))
            fh.write('\n')
    print(f'{wrong_scaffolds} scaffolds not found')
    """

if __name__ == '__main__':
    gff_path = sys.argv[1]
    file_format = gff_path[gff_path.rfind('.'):]
    if file_format not in ['.gtf']:
        print(f'File format {file_format} not recognized')
    else:
        split_assembled_genome(gff_path, sys.argv[2], file_format, sys.argv[3])
