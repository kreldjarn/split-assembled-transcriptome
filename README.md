# Split Assembled Transcriptome

Usage:
```
./split_assembled_transcriptome.py annotation.gff scaffolds.fasta output_directory
```

Parses a gff/gtf file into the following categories
- Unambiguous:
    * a) Processed
        + a.i)    Exon-exon splice junctions
    * b) Unprocessed
        + b.i)    Non-retained exon-intron splice junctions
        + b.ii)   Intron intersection
- Ambiguous:
    + c) Exon intersection
    + d) Retained introns (i.e. pairwise exon set difference)
    + e) Retained exon-intron splice junctions

Outputs a separate fasta file for each category

# A Major Caveat

Please note that the phrase 'The `.GTF/.GFF3` standard' is an oxymoron. There is no standard. This means that in order to use the script with an annotation file output by a program it hasn't seen before, it is very likely that it needs to be modified.
