library(tidyverse)

output_3UTR = read_tsv('hsa_mock_3UTR.bed', col_names=c('chrom','start','end','strand','id'))
output_CDS = read_tsv('hsa_mock_CDS.bed', col_names=c('chrom','start','end','strand','id'))

output_3UTR_mmu = read_tsv('mmu_mock_3UTR.bed', col_names=c('chrom','start','end','strand','id'))
output_CDS_mmu = read_tsv('mmu_mock_CDS.bed', col_names=c('chrom','start','end','strand','id'))

gtf = read_tsv('hsa_mock.gtf', col_names=FALSE, comment='#')
gtf_mmu = read_tsv('mmu_mock.gtf', col_names=FALSE, comment='chr#')

test_main_logic = function (tx_id, transcript_feature, bed_table, gtf, species) { # test all 3UTR bed records are accurate with respect to the GTF

	context('dimensions')
	test_that('bed has five columns', {
		expect_equal(dim(bed_table)[2], 5)
	})

	context('pattern matching')

	test_that('first column only contains first chromosome', {
		expect_match(bed_table$chrom %>% as.character, '1')
	})

	test_that('start column only contains numbers - no leading zeroes', {
		expect_match(bed_table$start %>% as.character, '[1-9][0-9]+')
	})

	test_that('end column only contains numbers - no leading zeroes', {
		expect_match(bed_table$end %>% as.character, '[1-9][0-9]+')
	})

	test_that('strand column is either one or minus one', {
		expect_match(bed_table$strand %>% as.character, '(1|-1)')
	})

	if (species == 'hsa') {

	test_that('id column has the correct pattern', {
		expect_match(bed_table$id %>% as.character, '^ENST[0-9]+\\.[1-9][0-9]?$')
	})

	} else {
			
	test_that('id column has the correct pattern', {
                expect_match(bed_table$id %>% as.character, '^ENSMUST[0-9]+\\.[1-9][0-9]?$')
        })

	}

	context('main logic')
	
	gtf = filter(gtf, X3 == transcript_feature)
	gtf = separate(
		data=gtf,
		col=X9, 
		sep=';',
		into=c(
			"gene_id", "gene_version", "transcript_id", "transcript_version",
			"gene_name", "gene_source", "gene_biotype", "transcript_name",
			"transcript_source", "transcript_biotype", "tag", "transcript_support"
		)
	)
	gtf$transcript_id = gsub('[^A-Z0-9]','',gtf$transcript_id)
	gtf$transcript_version = gsub('[^0-9]','',gtf$transcript_version)
	gtf$transcript_id_full = paste(gtf$transcript_id, gtf$transcript_version, sep='.')
	gtf = filter(gtf, transcript_id_full == tx_id)
	gtf$X7 = gsub('\\+','1', gtf$X7)
	gtf$X7 = gsub('\\-','-1', gtf$X7)		

	test_that ('The strand is correct', {
		expect_equal(gtf$X7, bed_table$strand[bed_table$id == tx_id] %>% as.character)
	})

	test_that ('The chromosomes are correct', {
                expect_equal(gsub('chr','',gtf$X1) %>% as.character, bed_table$chrom[bed_table$id == tx_id] %>% as.character)
        })

	test_that ('The start codons are correct', {
                expect_equal( (gtf$X4 - 1) %>% as.character, bed_table$start[bed_table$id == tx_id] %>% as.character) #NB: The shift to the correct co-ordinate system
        })

	test_that ('The end codons are correct', {
                expect_equal( gtf$X5 %>% as.character, bed_table$end[bed_table$id == tx_id] %>% as.character) #NB: The shift to the correct co-ordinate system
        })
}

tx_ids_3UTR = output_3UTR$id %>% unique()
tx_ids_CDS =  output_CDS$id %>% unique()

tx_ids_3UTR_mmu = output_3UTR_mmu$id %>% unique()
tx_ids_CDS_mmu =  output_CDS_mmu$id %>% unique()

map(tx_ids_3UTR, test_main_logic, 'three_prime_utr', output_3UTR, gtf, 'hsa')
map(tx_ids_CDS, test_main_logic, 'CDS', output_CDS, gtf, 'hsa')
map(tx_ids_3UTR_mmu, test_main_logic, 'three_prime_utr', output_3UTR_mmu, gtf_mmu, 'mmu')
map(tx_ids_CDS_mmu, test_main_logic, 'CDS', output_CDS_mmu, gtf_mmu, 'mmu')
