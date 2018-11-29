library(testthat)

### base

test_file('test_3UTR_alignments.R')
test_file('test_CDS_alignments.R')
test_file('test_bedgraph_averaging.R')

setwd('./gtf_to_bed')

test_file('test_gtf_to_bed.R')

setwd('../targetscan')

test_file('test_convert_fasta_to_ts.R')
test_file('test_get_miRNA_families.R')
test_file('test_get_mir_for_context.R')
test_file('test_prep_alt_mir.R')
