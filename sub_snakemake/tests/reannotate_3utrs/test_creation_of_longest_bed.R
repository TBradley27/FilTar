library(plyr)
library(tidyverse)
library(testthat)

source('/gpfs/afm/moxon/thomas2/APAtrap/sub_snakemake/reannotate_3utrs/merge_multiple_bed_files_function.R')

output = merge_multiple_bed_files('mock1.bed','mock_bed')

#print(output)

context('single exon forward strand - variable end codons')

test_that("tx2 has the correct start codon annotation", {
                expect_equal(350, output$start[output$id == 'tx2'])
})

test_that("tx2 has the correct end codon annotation", {
                expect_equal(400, output$end[output$id == 'tx2'])
})

context('single exon forward strand - variable start codons')

test_that("tx9 has the correct start codon annotation", {
                expect_equal(690, output$start[output$id == 'tx9'])
})

test_that("tx9 has the correct end codon annotation", {
                expect_equal(750, output$end[output$id == 'tx9'])
})

context('single exon forward strand - variable start and end codons')

test_that("tx6 has the correct start codon annotation", {
                expect_equal(540, output$start[output$id == 'tx6'])
})

test_that("tx6 has the correct end codon annotation", {
                expect_equal(600, output$end[output$id == 'tx6'])
})

context('single exon negative strand - variable end codons')

test_that("tx3 has the correct start codon annotation", {
                expect_equal(366, output$start[output$id == 'tx3'])
})

test_that("tx3 has the correct end codon annotation", {
                expect_equal(450, output$end[output$id == 'tx3'])
})

context('single exon negative strand - variable start codons')

test_that("tx4 has the correct start codon annotation", {
                expect_equal(450, output$start[output$id == 'tx4'])
})

test_that("tx4 has the correct end codon annotation", {
                expect_equal(500, output$end[output$id == 'tx4'])
})

context('single exon negative strand - variable start and end codons')

test_that("tx5 has the correct start codon annotation", {
                expect_equal(490, output$start[output$id == 'tx5'])
})

test_that("tx5 has the correct end codon annotation", {
                expect_equal(599, output$end[output$id == 'tx5'])
})

context('multi exon - positive strand - different mock2')

test_that("tx1 has the correct start codon annotation", {
                expect_equal(250, output$start[output$id == 'tx1'][1])
})

test_that("tx1 has the correct start codon annotation", {
                expect_equal(300, output$start[output$id == 'tx1'][2])
})

test_that("tx1 has the correct initial end codon annotation", {
                expect_equal(300, output$end[output$id == 'tx1'][1])
})

test_that("tx1 has the correct terminal end codon annotation", {
                expect_equal(350, output$end[output$id == 'tx1'][2])
})

context('multi exon - negative strand - different mock2')

test_that("tx8 has the correct start codon annotation", {
                expect_equal(650, output$start[output$id == 'tx8'][1])
})

test_that("tx8 has the correct start codon annotation", {
                expect_equal(600, output$start[output$id == 'tx8'][2])
})

test_that("tx8 has the correct initial end codon annotation", {
                expect_equal(700, output$end[output$id == 'tx8'][1])
})

test_that("tx8 has the correct terminal end codon annotation", {
                expect_equal(630, output$end[output$id == 'tx8'][2])
})


