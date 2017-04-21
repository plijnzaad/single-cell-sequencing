# single-cell-sequencing

A script for quality control of a scRNA-Seq run on one or more 384 wells
plates using the CEL-Seq2 protocol and the scripts found in the
https://bitbucket.org/princessmaximacenter/scseq repository.

This is a fork of https://github.com/MauroJM/single-cell-sequencing by
Mauro Muraro (van Oudenaarden group, Hubrecht Institute, Utrecht).

The script now substantionally deviates from the original. Among the
changes are some refactoring (plate.plot now plots exactly one plate),
graphs showing the number of unmapped reads (for that, the `*.cout?.csv`
files produced by
https://github.com/plijnzaad/scseq/blob/master/process_sam_cel384v2.pl
have been extended to include the number of unmapped reads per well),
and graphs aimed at judging whether additional sequencing is likely to
be meaningful (for this, `process_sam_cel384v2.pl` has been extended to
output the number of genes and umis 'seen' during traversal of the
`.bam` file)

The current fork will remain frozen here, but is currently developed
further on our bitbucket
acccount. [Contact me](mailto:p.lijnzaad@gmailcom) if you want access.

