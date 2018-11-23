library(tidyverse)

output = read_tsv('mock2.tsv', col_names=FALSE, col_types='cccc')
print(output)

context('shape')

test_that('Output table has the correct number of rows', {
        expect_equal( dim(output)[1], 2)
})

context('pattern matching')
 
test_that('Seed sequences do not contain invalid characters', {
        expect_match( output$X1, "^mmu-(miR|let|lin)-?[0-9]{1,5}[a-z]?(-[0-9]-)?((-3|5)p)?")
})

test_that('First and third columns are equal', {
        expect_equal( output$X1, output$X3)
})

test_that('tax_id contain 4-6 repetitions of allowed characters', {
        expect_match( output$X2, '^[0-9]{4,6}$' )
})

test_that('mature seqeunce has the right length and characters', {
        expect_match( output$X4, '^[ACTGU]{18,24}$' )
})

context('header sequence correspondence')

test_that('First header matches sequence', {
        expect_equal( 'CAAAAAAAGCGCGCGCGCGT', output$X4[output$X1 == 'mmu-miR-4-5p'])
})
