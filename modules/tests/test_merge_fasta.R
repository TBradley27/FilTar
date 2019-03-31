library(tidyverse)
library(plyr)

source('../get_utr_and_cds/MergeFasta_function.R')

output = MergeFasta('unmerged_fasta.fa')

print(output)

unmerged = read_tsv('unmerged_fasta.fa', col_names=c("fasta_unicolumn"))
print(unmerged)

context('test that fasta merging works')

test_that("sequences care properly concatenated", {
  expect_equal(output[2], paste('AGAGAGAGGAGGACCTGCTGTCACA','AGAGCCTCTCGCGCTCTCT',sep=''))
  expect_equal(output[4], paste('GCGCGCGCGCGTATATATGCCTGC','CGCGCGCTATATATATGCGTTATA',sep=''))
  expect_equal(output[6], paste('CGCGCTATATAGAGCTTCCGCTAGATTATATAGGC','AAATATATGGATATTAGAGA',sep=''))
})

test_that("headers are preserved", {
  expect_equal(output[1], '>seq1')
  expect_equal(output[3], '>seq2')
  expect_equal(output[5], '>seq3')
})

