library(testthat)

### base

test_file('test_3UTR_alignments.R')
test_file('test_CDS_alignments.R')

setwd('./gtf_to_bed')

test_file('test_gtf_to_bed.R')

setwd('../targetscan')

test_file('test_convert_fasta_to_ts.R')
test_file('test_prep_alt_mir.R')

setwd('./main')

test_file('test_contextpp_output.R')
