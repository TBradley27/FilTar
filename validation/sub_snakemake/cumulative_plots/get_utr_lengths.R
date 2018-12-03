library(plyr)
library(tidyverse)
source('sub_snakemake/cumulative_plots/get_utr_lengths_function.R')

utr_lengths = get_utr_lengths(snakemake@input[[1]]) 

write.table(
	x=utr_lengths,
	file=snakemake@output[[1]],
	sep="\t",
	quote=FALSE,
	row.names=FALSE,
	col.names=TRUE
)	
