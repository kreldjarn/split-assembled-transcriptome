# Split Assembled Transcriptome

Given an annotation file on the `.GTF` or `.GFF3` formats and a .`fasta` file containing the corresponding scaffolds, outputs two fasta files containing the processed and the unprocessed versions of the transcriptome, respectively.

# A Major Caveat

Please note that the phrase 'The `.GTF/.GFF3` standard' is an oxymoron. There is no standard. This means that in order to use the script with an annotation file output by a program it hasn't seen before, it is very likely that it needs to be modified.
