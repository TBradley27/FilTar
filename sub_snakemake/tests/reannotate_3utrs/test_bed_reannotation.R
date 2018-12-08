library(tidyverse)
library(plyr)
library(testthat)

source('/gpfs/afm/moxon/thomas2/APAtrap/sub_snakemake/reannotate_3utrs/extend_bed_function2.R')

x = get_full_bed('chr10_hsa.bed6','hsa_HeLa_chr10.utr.bed','hsa_all_transcripts.txt')

#dim(x) %>% print()

context('single exon forward strand - updated 3UTR annotation')

y = x %>% filter(id=='ENST00000224237.9')

test_that("ENST00000224237.9 has the correct start codon annotation", {
                expect_equal(17237230, y$start[1], )
})

test_that("ENST00000224237.9 has the correct end codon annotation", {
                expect_equal(17237588, y$stop[1], )
})

context('multi exon positive strand - conflicting 3UTR annotations - keep canonical, discard APAtrap annotation')

y = x %>% filter(id=='ENST00000464969.6')

test_that("ENST00000464969.6 has the correct start codon annotation", {
                expect_equal(102157494, y$start[1], )
})

test_that("ENST00000464969.6 has the correct end codon annotation", {
                expect_equal(102162599, y$stop[10], )
})

context('single exon negative strand - updated 3UTR annotation')

y = x %>% filter(id=='ENST00000298510.3')

test_that("ENST00000298510.3 has the correct start codon annotation", {
                expect_equal(119167702, y$start[1], )
})

test_that("ENST00000298510.3 has the correct end codon annotation", {
                expect_equal(119168533, y$stop[1], )
})

context('multi exon negative strand - conflicting 3UTR annotations - keep canonical, discard APAtrap annotation')

y = x %>% filter(id=='ENST00000463743.5')

test_that("ENST00000463743.5 has the correct start codon annotation", {
                expect_equal(93323077, y$start[1], )
})

test_that("ENST00000463743.5 has the correct end codon annotation", {
                expect_equal(93307001, y$stop[7], )
})

context('single exon positive strand - tx previously without a 3UTR annotation')

y = x %>% filter(id=='ENST00000363306.1')

test_that("ENST00000363306.1 has the correct start codon annotation", {
                expect_equal(86889569, y$start[1], )
})

test_that("ENST00000363306.1 has the correct end codon annotation", {
                expect_equal(86889682, y$stop[1], )
})

context('single exon negative strand - tx previously without a 3UTR annotation')

y = x %>% filter(id=='ENST00000378952.7')

test_that("ENST00000378952.7 has the correct start codon annotation", {
                expect_equal(12167673, y$start[1], )
})

test_that("ENST00000378952.7 has the correct end codon annotation", {
                expect_equal(12167811, y$stop[1], )
})

context('single exon positive strand - unchanged 3UTR annotation')

y = x %>% filter(id=='ENST00000381604.8')

test_that("ENST00000378952.8 has the correct start codon annotation", {
                expect_equal(252470, y$start[1], )
})

test_that("ENST00000378952.8 has the correct end codon annotation", {
                expect_equal(254626, y$stop[1], )
})

context('multi exon positive strand - unchanged 3UTR annotation')

y = x %>% filter(id=='ENST00000474119.5')

test_that("ENST00000474119.5 has the correct start codon annotation", {
                expect_equal(4847234, y$start[1], )
})

test_that("ENST00000474119.5 has the correct end codon annotation", {
                expect_equal(4848062, y$stop[2], )
})

context('single exon negative strand - unchanged 3UTR annotation')

y = x %>% filter(id=='ENST00000564130.2')

test_that("ENST00000564130.2 has the correct start codon annotation", {
                expect_equal(46891, y$start[1], )
})

test_that("ENST00000564130.2 has the correct end codon annotation", {
                expect_equal(47056, y$stop[1], )
})

context('multi exon negative strand - unchanged 3UTR annotation')

y = x %>% filter(id=='ENST00000567466.1')

test_that("ENST00000567466.1 has the correct start codon annotation", {
                expect_equal(48424, y$start[1], )
})

test_that("ENST00000567466.1 has the correct end codon annotation", {
                expect_equal(48114, y$stop[2], )
})

context('bed file has the correct dimensions')

test_that('bed file has the correct number of columns', {
	expect_equal(dim(x)[2], 5)
})

normal_bed = read_tsv('chr10_hsa.bed6', col_names=c('chromosome','start','stop','strand','id'))
extended_utrs = read_tsv('hsa_HeLa_chr10.utr.bed', col_names=c('chromosome','start','stop','id','dummy','strand'))

normal_bed = separate(normal_bed, id, into=c('id','version'))
extended_utrs = separate(extended_utrs, id, into=c('id','dummy2','chrom_dup','strand_dup'))

normal_bed$version = NULL

extended_utrs = extended_utrs[,c('chromosome','start','stop','strand','id')]
extended_utrs$chromosome = stringr::str_replace_all(extended_utrs$chromosome, 'chr','')
normal_bed$chromosome = stringr::str_replace_all(normal_bed$chromosome, 'chr','')

old_records = extended_utrs[extended_utrs$id %in% normal_bed$id,]

test_that('bed file has the correct number of rows', {
        expect_equal(dim(x)[1], dim(normal_bed)[1] + dim(extended_utrs)[1] - dim(old_records)[1])
})

#context('bed file has the correct ordering of transcripts within a chromosome')

#z = x

#tx_ids = z$id %>% unique()
#tmp = map(tx_ids, function (x) {filter(z, id == x )[1,]   })
#tmp = ldply(tmp, data.frame) %>% as.tibble() 
#print(tmp, n=Inf)

#print(all(tmp$start[1:10] == cummax(tmp$start[1:10]))) 

# Discard test as strictly speaking, I don't think transcript start IDs should be monotonic






