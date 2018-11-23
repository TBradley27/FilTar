library(tidyverse)
library(plyr)
library(testthat)

x = read_tsv('../../../results/mature_mirna_seed.tsv', col_names=FALSE)

context('dimensions: table')

test_that('output has two columans', {
	expect_equal( 2, dim(x)[2] )
})

#test_that('output has two columns', {
#        expect_equal( 3, dim(x)[1] )
#})

context('content: header column')

test_that('Seed sequences do not contain invalid characters', {
        expect_match( x$X1, "^[a-z][a-z][a-z]-(miR|let|lin)-?[0-9]{1,5}[a-z]?(-[0-9]-)?((-3|5)p)?")
})

context('content: seed column')

test_that('Seed sequences do not contain invalid characters', {
        expect_match( x$X2, '^[ACGU]{7}$' )
})

