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
test_file('test_prep_alt_mir.R')
#test_file('targetscan/test_AIR_retrieval.R')
test_file('test_get_mir_for_context.R')

setwd('./main')

test_file('test_fix_ts_output.R')
test_file('test_contextpp_output.R')

setwd('../../reannotate_3utrs')

test_file('test_bed_reannotation.R')
#test_file('test_creation_of_longest_bed.R')

setwd('../cumulative_plots')

test_file('test_get_utr_lengths.R')
test_file('test_get_utr_sites.R')
