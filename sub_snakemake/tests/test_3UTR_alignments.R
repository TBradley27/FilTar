library(tidyverse)

all_tests = function (species, taxonomic_id) {

	context(str_interp('${species}'))	

	alignment = read_tsv( str_interp('../../results/msa/no_reannotation/${species}_chrY_3UTR_msa.tsv'), col_names=c('tx_id','tax_id','seq'))
	canon_tax_ids = read_tsv('../../config/species_basic2.tsv', col_names=c('tax_id'))
	bed_file = read_tsv( str_interp('../../results/bed/no_reannotation/${species}_3UTR.chrY.bed'), col_names=FALSE)

	tax_ids = alignment$tax_id %>% unique()

	#print(canon_tax_ids)

	test_that("The tabular data has three columns", {
	  expect_equal(dim(alignment)[2], 3) # test num columns
	})

	#test_that("The tabular data has 6362 rows", {
	#  expect_equal(dim(alignment)[1], 6362) # test num rows
	#})

	test_that("The tax_id column only contains recognised taxonomic IDs from Agarwal et al. 84-way alignment", {
	  #expect_true(grepl('^[0-9]*$', tax_ids) %>% all()) # test that expression only contains numerical characters
	  #expect_true(nchar(tax_ids) %>% min() == 4) # test that there are no taxonomic IDs which are shorter than 4 digits
	  #expect_true(nchar(tax_ids) %>% max() == 7) # test that there are no taxonomic IDs which are longer than 7 digits
	  expect_true(tax_ids %in% canon_tax_ids$tax_id %>% all())
	})

#	test_that("The tax_id uses the same 84 species used in Targetscan7 (Agarwal et al. 2015)", {
#          #expect_true(grepl('^[0-9]*$', tax_ids) %>% all()) # test that expression only contains numerical characters
#          #expect_true(nchar(tax_ids) %>% min() == 4) # test that there are no taxonomic IDs which are shorter than 4 digits
#          #expect_true(nchar(tax_ids) %>% max() == 7) # test that there are no taxonomic IDs which are longer than 7 digits
#          expect_true(canon_tax_ids$tax_id %in% tax_ids %>% all())
#        })

	test_that("The tx_id column only contains recognised transcript IDs", {
	  expect_true(alignment$tx_id %in% bed_file$X5 %>% all())
	})

	#print(alignment$tx_id[!alignment$tx_id %in% bed_file$X5])

	test_that("There is an output for every single transcript", {
	  expect_true(alignment$tx_id %>% unique() %>% length() == bed_file$X5 %>% unique() %>% length())
	})

	tmp =  !alignment %>% filter(tax_id == taxonomic_id) %>% select(tx_id,tax_id) %>% duplicated()
	#print (tmp)

	test_that("There is only one entry per transcript per species for the reference record", {
	  expect_true(tmp %>% all())
	})

	seq_lengths_df = alignment %>% filter(tax_id == taxonomic_id) %>% select(tx_id,seq) %>% mutate(seq=gsub('-','',seq)) %>% mutate(seq_length = nchar(seq))

	#seq_lengths_df[order(seq_lengths_df$seq_length),] %>% print()

	bed_seq_lengths = bed_file %>% mutate(exon_length = X3 - X2) %>% group_by(X5) %>% summarise(tx_seq_total = sum(exon_length)) %>% unique()
	#print('foo')
	#bed_seq_lengths[order(bed_seq_lengths$tx_seq_total),] %>% print()

	y = merge(seq_lengths_df, bed_seq_lengths, by.x='tx_id', by.y='X5') %>% as.tibble()

	test_that("reference sequences are all of the correct length", {
		expect_equal(y$seq_length, y$tx_seq_total)	
	})

	biomart_data = read_tsv(str_interp('${species}_mart_export2.tsv'), col_names=c('tx_id','seq'))
	alignment_human = alignment %>% filter(tax_id == taxonomic_id) %>% mutate(tx_id=gsub('\\..*','',tx_id)) %>% mutate(seq=gsub('-','',seq)) %>% mutate(seq=toupper(seq))

	merged = merge(biomart_data, alignment_human, by='tx_id') %>% as.tibble()

	#print(merged)
        #print(merged$tx_id[1])
	#print(merged$seq.x[1])
	#print(merged$seq.y[1])
        #print(merged$tx_id[70])
	#print(merged$seq.x[70])
        #print(merged$seq.y[70])

	test_that("alignment reference sequences are the same as those directly downloaded from biomart", {
		expect_equal(merged$seq.x, merged$seq.y, ignore.case=TRUE)	
})
}

all_tests('hsa','9606')
all_tests('mmu','10090')
