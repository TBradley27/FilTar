library(plyr)
library(tidyverse)
source('sub_snakemake/cumulative_plots/alt_sites_function.R')

alt_targets = get_alt_sites(snakemake@input$targets,snakemake@input$utr_lens)

write.table(
	x=alt_targets,
	file=snakemake@output[[1]],
	row.names=FALSE,
	col.names=TRUE,
	quote=FALSE,
	sep="\t"
)	
