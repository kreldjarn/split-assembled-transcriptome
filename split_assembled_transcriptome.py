#!/usr/bin/env python3

"""
Author:
    Kristján Eldjárn Hjörleifsson
    keldjarn@caltech.edu

Usage:
./split_assembled_transcriptome.py annotation.gtf3 reads.fasta output_directory L

Parses a gff/gtf file into the following collections
===============================================================================
    Unambiguous:
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
                   merge_intervals, overlap, len_overlap)

def append_to_fasta(tr, data, path):
    with open(path, 'a') as fh:
        fh.write(f'>{tr}\n')
        fh.write('\n'.join([data[i:i+80] for i in range(0, len(data), 80)]))
        fh.write('\n')

def collapse_data(intervals, sequence, strand):
    data = ''.join(map(lambda iv: sequence[iv[0]:iv[1]], intervals))
    data = collapse_N(data.upper())
    if strand == '-':
        data = reverse_complement(data)
    return data

def process_gene(gene, scaffolds, od, L):
    try:
        sequence = scaffolds[gene['scaffold']]
    except KeyError as e:
        print(f'Scaffold {gene["scaffold"]} not found.')
        return
    transcripts = {
        tr: [(e['start'], e['end']) for e in exons['exons']]\
             for tr, exons in gene['trs'].items()
    }
    ivs = list(chain(*transcripts.values()))
    start = min(ivs, key=lambda i: i[0])[0]
    end = max(ivs, key=lambda i: i[1])[1]

    exon_union = list(merge_intervals(ivs))
    itrs = [(e1[1], e2[0]) for _, exons in transcripts.items()\
                           for e1, e2 in zip(exons[:-1], exons[1:])]

    # Look for retained introns at beginning/end of gene
    for tr, exons in transcripts.items():
        if exons[0][0] > start:
            itrs.append((start, exons[0][0]))
        if exons[-1][1] < end:
            itrs.append((exons[-1][1], end))
    itrs = sorted(itrs, key=itemgetter(0))
    itr_union = list(merge_intervals(itrs))

    #===================================#
    #   a) Exon-exon splice junctions   #
    #===================================#
    for tr, exons in transcripts.items():
        ee_sjs = [(ee[0][1] - L + 1, ee[1][0] + L - 1)\
                  for ee in zip(exons[:-1], exons[1:])]
        append_to_fasta(tr, collapse_data(ee_sjs, sequence, gene['strand']),
                        f'{od}/a_exon_exon_splice_junctions.fasta')

    #====================================================#
    #   b.i) Non-retained exon-intron splice junctions   #
    #====================================================#
    sjs = []
    retained_sjs = []
    nonretained_sjs = []
    for e in ivs:
        ie = (e[0] - L + 1, e[0] + L - 1)
        ei = (e[1] - L + 1, e[1] + L - 1)
        ie_retained, ei_retained = False, False
        for iv in itr_union:
            if len_overlap(iv, (e[0], ie[1])) > 0:
                ie_retained = True
            if len_overlap(iv, (ei[0], e[1])) > 0:
                ei_retained = True
        for iv in exon_union:
            if len_overlap(iv, (ie[0], e[0])) > 0:
                ie_retained = True
            if len_overlap(iv, (e[1], ei[1])) > 0:
                ei_retained = True
        if ie_retained:
            retained_sjs.append(ie)
        else:
            nonretained_sjs.append(ie)
        if ei_retained:
            retained_sjs.append(ei)
        else:
            nonretained_sjs.append(ei)

    append_to_fasta(gene['name'], collapse_data(nonretained_sjs,
                                                sequence,
                                                gene['strand']),
                    f'{od}/bi_nonretained_exon_intron_splice_junctions.fasta')

    #===============================#
    #   b.ii) Intron intersection   #
    #===============================#
    itr_intersect = [(i[0][1], i[1][0]) for i in zip(exon_union[:-1],
                                                     exon_union[1:])]
    append_to_fasta(gene['name'], collapse_data(itr_intersect,
                                                sequence,
                                                gene['strand']),
                    f'{od}/bii_intron_intersection.fasta')

    #==========================#
    #   c) Exon intersection   #
    #==========================#
    exon_intersect = [(i[0][1], i[1][0]) for i in zip(itr_union[:-1],
                                                      itr_union[1:])]
    # If the intron union does not stretch to gene boundaries, we add leading
    # and/or trailing end of gene to exon intersection
    try:
        if start < min(itr_union, key=lambda i: i[0])[0]:
            exon_intersect = [(start, itr_union[0][0])] + exon_intersect
        if end > max(itr_union, key=lambda i: i[1])[1]:
            exon_intersect = exon_intersect + [(itr_union[-1][1], end)]
    except ValueError as e:
        # Single exon isoform
        pass
    append_to_fasta(gene['name'], collapse_data(exon_intersect,
                                                sequence,
                                                gene['strand']),
                    f'{od}/c_exon_intersection.fasta')

    #=========================#
    #   d) Retained introns   #
    #=========================#
    retained_introns = []
    for itr in itr_union:
        for ex in exon_union:
            ret = overlap(itr, ex)
            if ret != (0, 0):
                retained_introns.append(ret)
    retained_introns = sorted(retained_introns, key=itemgetter(0))
    retained_introns = list(merge_intervals(retained_introns))
    append_to_fasta(gene['name'], collapse_data(retained_introns,
                                                sequence,
                                                gene['strand']),
                    f'{od}/d_retained_introns.fasta')

    #==============================================#
    #   e) Retained exon-intron splice junctions   #
    #==============================================#
    append_to_fasta(gene['name'], collapse_data(retained_sjs,
                                                sequence,
                                                gene['strand']),
                    f'{od}/e_retained_exon_intron_splice_junctions.fasta')

def parse_gtf(path, scaffolds, file_format, od, L):
    # Who the heck came up with this hecking file format?
    genes = {}
    tr2g = []
    # Apparently, gtf files are not necessarily ordered by gene, so we cannot
    # do this in a single pass-through
    with open(path, 'r') as fh:
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

            if gene_id not in genes:
                genes[gene_id] = {
                        'name': gene_id,
                        'trs': {},
                        'scaffold': fields['scaffold'],
                        'strand': fields['strand']
                }

            if file_format == '.gtf':
                if fields['feature'] in ['exon', 'UTR', 'stop_codon']:
                    tr = fields['rest']['transcript_id']
                    if tr not in genes[gene_id]['trs']:
                        genes[gene_id]['trs'][tr] = {
                                'name': tr,
                                'exons': [fields]
                        }
                        tr2g.append(f'{tr}\t{gene_id}')
                    else:
                        genes[gene_id]['trs'][tr]['exons'].append(fields)

    for _, gene in genes.items():
        process_gene(gene, scaffolds, od, L)

    with open(f'{od}/tr2g', 'w') as fh:
        fh.write('\n'.join(tr2g))

def split_assembled_genome(gtf_path, fasta_path, file_format, od='.', L=150):
    scaffolds = parse_fasta(fasta_path)
    parse_gtf(gtf_path, scaffolds, file_format, od, L)

if __name__ == '__main__':
    gtf_path = sys.argv[1]
    file_format = gtf_path[gtf_path.rfind('.'):]
    if file_format not in ['.gtf']:
        print(f'File format {file_format} not recognized')
    else:
        split_assembled_genome(gtf_path, sys.argv[2], file_format, sys.argv[3])
