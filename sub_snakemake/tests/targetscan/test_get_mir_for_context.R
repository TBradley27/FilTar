library(tidyverse)
source('../../target_prediction/targetscan/get_mirna_context_function.R')

output = get_mirna_context('mock3.tsv','mock2.tsv','mmu')
print(output)

context('pattern matching')

test_that ('The identifier column only contains 1-4 repetitions of 0-9',{
	expect_match( output$identifier, '^[0-9]{1,4}$' )
})

test_that('The first column is monotonic', {
        expect_true(all(output$identifier == cummax(output$identifier)))
})

test_that('The tax id column is correct', {
        expect_match(as.character(output$tax_id.x), '10090')
})

test_that('Seed sequences do not contain invalid characters', {
        expect_match( output$miRNA, "^mmu-(miR|let|lin)-?[0-9]{1,5}[a-z]?(-[0-9]-)?((-3|5)p)?")
})

test_that('mature seqeunce has the right length and characters', {
        expect_match( output$miRNA_sequence, '^[ACTGU]{18,24}$' )
})

context('main logic')

test_that('A single seed can lead to multiple mature sequences', {
	expect_true(output$identifier %>% table %>% max > 1)
})

input_seeds = read_tsv('mock3.tsv', col_names=TRUE) %>% filter(tax_id=='10090')

test_that('every reference seed sequence has a corresponding mature sequence', {
	expect_match(input_seeds$seq %in% substring(output$miRNA_sequence, 2, 8) %>% as.character, as.character(TRUE))
})

test_that('identifiers map to the correct sequence name', {
	expect_equal(as.character(output$identifier[output$miRNA == 'mmu-miR-4-5p']),'1')
})

test_that('identifiers map to the correct sequence', {
        expect_equal(as.character(output$identifier[output$miRNA_sequence == 'CAAAAAAAGCGCGCGCGCGT']),'1')
})

test_that('sequence maps to the correct sequence name', {
        expect_equal(as.character(output$miRNA[output$miRNA_sequence == 'CAAAAAAAGCGCGCGCGCGT']),'mmu-miR-4-5p')
})

test_that('output does not contain non-reference families',{
	expect_match('GATTGCG' %in% output$seq %>% as.character, 'FALSE')
})

context('complete match')

expected_output = tibble(identifier=c('1','1','2'),tax_id.x=c(as.integer(10090),as.integer(10090),as.integer(10090)),miRNA=c('mmu-miR-4-5p','mmu-miR-4b-5p','mmu-miR-1-3p'),
			miRNA_sequence=c('CAAAAAAAGCGCGCGCGCGT','CAAAAAAAGCGCGCGTTCGT','ACCTCTCGCGAGAGAGGCGC'))

test_that('expected output matches actual output', {
        expect_equal(expected_output, output)
})
