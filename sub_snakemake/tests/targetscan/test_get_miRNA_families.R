library(plyr)
library(tidyverse)
source('../../../sub_snakemake/target_prediction/targetscan/get_mirna_family_function.R')

#input = read_tsv('../../../results/mature_mirna_seed.tsv', col_names=c('identifier','seq'))

#print(input)

output = get_mirna_family('mock.tsv', 'mmu')

print(output)

write.table(output, 'mock3.tsv', col.names=TRUE,row.names=FALSE, quote=FALSE, sep="\t")

context('pattern matching')

test_that('id column contains 1-4 numbers', {
        expect_match(as.character(output$identifier), '^[0-9]{1,4}$')
})

test_that('Seed sequences contain 7 repetitions of allowed characters', {
        expect_match( as.character(output$seq), '^[ACGTU]{7}$' )
})

test_that('tax_id contain 4-6 repetitions of allowed characters', {
        expect_match( output$tax_id, '^[0-9]{4,6}$' )
})

context('main logic')

test_that('The first column is monotonic', {
        expect_true(all(output$identifier == cummax(output$identifier)))
})

test_that('The combined seq and tax_id columns contain no duplicates', {
        expect_match(output[,c('seq','tax_id')] %>% duplicated %>% as.character, 'FALSE')
})

tmp =  output[,c('seq','tax_id')] %>% filter(seq=='CCTCTCG')
print(tmp)
tmp2 = c('9606','10090') %in% tmp$tax_id

test_that('The output factorises orthologous seeds as belonging to the same family', {
	expect_match(tmp2 %>% as.character, '^TRUE$')
})

test_that('The output does not contain seed sequences without a reference orthologue', {
        expect_equal(dim(output %>% filter(seq=='AAAAAAA'))[1], 0)
})
