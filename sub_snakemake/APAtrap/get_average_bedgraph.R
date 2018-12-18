#!/bin/env Rscript

united_bedgraph = readr::read_tsv(
	file=snakemake@input[[1]],
	col_names=FALSE,
	col_types=readr::cols(.default = 'd', X1 = 'c', X2 = 'i', X3 = 'i')
)

united_bedgraph = filtar::AvgBedgraph(united_bedgraph)

united_bedgraph = united_bedgraph[,c('X1','X2','X3','avg')]

write.table(
	x=united_bedgraph,
	file=snakemake@output[[1]],
	row.names=FALSE,
	col.names=FALSE,
	quote=FALSE,
	sep="\t"
)
