#!/bin/bash Rscript

library(readr)
library(dplyr)

utr_bed = read_tsv(
	snakemake@input[[1]],
	col_names=FALSE,
	col_types = 'ciiic'
)

utr_bed = utr_bed %>% dplyr::select(X5,X2,X3)

utr_bed$utr_length = utr_bed$X3 - utr_bed$X2

# apply function over a factor - in this case the transcript ID i.e. X5
utr_bed$total_utr_length = ave(utr_bed$utr_length, as.factor(utr_bed$X5), FUN = sum)

utr_bed$tissue = snakemake@wildcards[["tissue"]]

utr_bed = utr_bed %>% dplyr::select(X5,tissue,total_utr_length)

utr_bed = unique(utr_bed)

write.table(utr_bed, file=snakemake@output[[1]], sep="\t", quote=FALSE, col.names=FALSE, row.names=FALSE)
