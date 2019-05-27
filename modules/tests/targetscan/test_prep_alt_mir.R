library(testthat)
library(tidyverse)

run_tests = function(input, species) {
	output = read_tsv(input, col_names=FALSE, col_types='cccc')
	#print(output)	 
	
	test_that('First and third columns are equal', {
		expect_equal( output$X1, output$X3)
	})

	test_that('tax_id contain 4-6 repetitions of allowed characters', {
		expect_match( output$X2, '^[0-9]{4,6}$' )
	})

	test_that('mature sequence has the right length and characters', {
		expect_match( output$X4, '^[ACTGU]{18,24}$' )
	})

	
	if (species == 'hsa') {

		context('shape')


		test_that('Output table has the correct number of rows', {
		expect_equal( dim(output)[1], 4)
		})

		context('pattern matching')

		test_that('Seed sequences do not contain invalid characters', {
		expect_match( output$X1, "^hsa-(miR|let|lin)-?[0-9]{1,5}[a-z]?(-[0-9]-)?((-3|5)p)?")
		})

		context('header sequence correspondence')

		test_that('First header matches sequence', {
			expect_equal( 'AAGCGCGCCCTCTCGCGAGA', output$X4[output$X1 == 'hsa-miR-1a-3p'])
		})

	} else if (species == 'mmu') {

		context('shape')

		test_that('Output table has the correct number of rows', {
		expect_equal( dim(output)[1], 3)
		})

		context('pattern matching')

		test_that('Seed sequences do not contain invalid characters', {
		expect_match( output$X1, "^mmu-(miR|let|lin)-?[0-9]{1,5}[a-z]?(-[0-9]-)?((-3|5)p)?")
		})

		context('header sequence correspondence')

		test_that('First header matches sequence', {
			expect_equal( 'CAAAAAAAGCGCGCGCGCGT', output$X4[output$X1 == 'mmu-miR-4-5p'])
		})

	}

}

run_tests('hsa_mock2.tsv', 'hsa')
run_tests('mmu_mock2.tsv', 'mmu' )
