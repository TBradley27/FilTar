library(testthat)
library(plyr)
library(tidyverse)
library(filtar)

output = get_utr_lengths('mock.bed')
output$tx_id = as.character(output$tx_id)

expected_output = read_tsv('expected_utr_lengths.tsv', col_names=c('tx_id','utr_length'))
print(expected_output)

context('table matching')
test_that('actual output and expected output is equal',{
	expect_equal(output,expected_output)
})
