library(plyr)
library(tidyverse)
library(filtar)

utr_lengths = get_utr_lengths(snakemake@input[[1]]) 

write.table(
	x=utr_lengths,
	file=snakemake@output[[1]],
	sep="\t",
	quote=FALSE,
	row.names=FALSE,
	col.names=TRUE
)
