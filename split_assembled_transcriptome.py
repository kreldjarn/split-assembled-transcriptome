#!/usr/bin/env python3

"""
Author:
    Kristján Eldjárn Hjörleifsson
    keldjarn@caltech.edu

Usage:
python3 split_assembled_transcriptome.py annotation.gtf3 reads.fasta output_directory L

Parses a gff/gtf file into the following collections
===============================================================================
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
===============================================================================

Returns the corresponding sequences from the corresponding fasta scaffolds,
correcting for strandedness.
"""

import sys
from operator import itemgetter
from itertools import chain

from utils import (reverse_complement, parse_rest, parse_fasta, collapse_N,
                   merge_intervals)

def append_to_fasta(tr, data, path):
    with open(path, 'a') as fh:
        fh.write(f'>{tr}\n')
        fh.write('\n'.join([data[i:i+80] for i in range(0, len(data), 80)]))
        fh.write('\n')

def collapse_data(intervals, sequence):
    data = ''.join(map(lambda iv: sequence[iv[0]:iv[1]], intervals))
    data = collapse_N(data.upper())
    if gene['strand'] == '-':
        data = reverse_complement(data)
    return data

def process_gene(gene, scaffolds, od, L):
    sequence = scaffolds[gene['scaffold']]
    transcripts = {
        tr: [(e['start'], e['end']) for e in exons['exons']]\
             for tr, exons in gene['trs'].items()
    }
    ivs = list(chain(*transcripts.values()))
    start, end = min(ivs, key=lambda i: i[0]), max(ivs, key=lambda i: i[1])

    exon_union = merge_intervals(ivs)
    itrs = [(i[0][1], i[1][0]) for i in zip(exons[:-1], exons[1:])\
                               for _, exons in gene['trs'].items()]
    for tr, exons in transcripts.items():
        if exons[0]['start'] > start:
            itrs.append((start, exons[0]['start']))
        if exons[-1]['end'] < end:
            itrs.append((exons[-1]['end'], end))
    itrs = sorted(itrs, key=itemgetter(0))
    itr_union = merge_intervals(itrs)

    #===================================#
    #   a) Exon-exon splice junctions   #
    #===================================#
    for tr, exons in transcripts.items():
        ee_sjs = [(ee[0][1] - L + 1, ee[1][0] + L - 1)\
                  for ee in zip(exons[:-1], exons[1:])]
        append_to_fasta(tr, collapse_data(ee_sjs, sequence), f'{od}/a.fasta')

    #====================================================#
    #   b.i) Non-retained exon-intron splice junctions   #
    #====================================================#
    sjs = []
    for tr, exons in transcripts.items():
        for e in exons:
            sjs.append((e['start'] - L + 1, e['start'] + L - 1))
            sjs.append((e['end'] - L + 1, e['end'] + L - 1))

    retained_sjs = []
    nonretained_sjs = []
    for sj in sjs:
        retained = False
        for iv in itr_union + exon_union:
            if iv[0] <= sj[0] and iv[1] >= sj[1]:
                retained = True
                break
        if retained:
            retained_sjs.append(sj)
        else:
            nonretained_sjs.append(sj)

    append_to_fasta(gene['name'], collapse_data(nonretained_sjs, sequence),
                    f'{od}/bi.fasta')

    #===============================#
    #   b.ii) Intron intersection   #
    #===============================#
    itr_intersect = [(i[0][1], i[1][0]) for i in zip(exon_union[:-1],
                                                     exon_union[1:])]
    append_to_fasta(gene['name'], collapse_data(itr_intersect, sequence),
                    f'{od}/bii.fasta')

    #==========================#
    #   c) Exon intersection   #
    #==========================#
    exon_intersect = [(i[0][1], i[1][0]) for i in zip(itr_union[:-1],
                                                      itr_union[1:])]
    # If the intron union does not stretch to gene boundaries, we add leading
    # and/or trailing end of gene to exon intersection
    if start < min(itr_union, key=lambda i: i[0]):
        exon_intersect = [(start, itr_union[0][0])] + exon_intersect
    if end > max(itr_union, key=lambda i: i[1]):
        exon_intersect = exon_intersect + [(itr_union[-1][1], end)]
    append_to_fasta(gene['name'], collapse_data(exon_intersect, sequence),
                    f'{od}/c.fasta')

    #=========================#
    #   d) Retained introns   #
    #=========================#

    #==============================================#
    #   e) Retained exon-intron splice junctions   #
    #==============================================#
    append_to_fasta(gene['name'], collapse_data(retained_sjs, sequence),
                    f'{od}/e.fasta')

        # If current isoform doesn't start at gene boundaries, we add those as
        # retained introns
        # if start < min(exons, key=lambda i: i[0]):
            # itrs = [(start, exons[0][0])] + itrs
        # if end > max(exons, key=lambda i: i[1]):
            # itrs = itrs + [(exons[-1][1], end)]




def parse_gff(path, scaffolds, file_format, od, L):
    # Who the hell came up with this hecking file format?
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
                'start': int(data[3]) - 1, # Scaffold is 0-indexed
                'end': int(data[4]) - 1,
                'strand': data[6],
                'rest': parse_rest(data[8], file_format)
            }
            gene_id = fields['rest']['gene_id']

            # Process current gene and start collecting exons for next gene
            if gene_id != prev_gene_id:
                if prev_gene_id is not None:
                    process_gene(gene, scaffolds, L)
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
        process_gene(gene, scaffolds, od, L)

    with open(f'{od}/tr2g', 'w') as fh:
        fh.write('\n'.join(tr2g))

def split_assembled_genome(gff_path, fasta_path, file_format, od='.', L=150):
    scaffolds = parse_fasta(fasta_path)
    print('scaffolds parsed')
    trs, itrs = parse_gff(gff_path, scaffolds, file_format, od, L)
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
