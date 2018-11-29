library(plyr)
library(tidyverse)

source('../../../validation/sub_snakemake/cumulative_plots/alt_sites_function.R')

output = get_alt_sites('mock_sites.tsv','mock_utr_lengths.tsv')
print(output)

test_that ('A site which is stops one short of the final nucleotide of the 3UTR is permitted', {
	expect_match('ENST00000327669.4' %in% output$a_Gene_ID %>% as.character,'TRUE') 
})

test_that ('A site which is stops st the final nucleotide of the 3UTR is not permitted', {
        expect_match('ENST000001.1' %in% output$a_Gene_ID %>% as.character,'FALSE')
})
