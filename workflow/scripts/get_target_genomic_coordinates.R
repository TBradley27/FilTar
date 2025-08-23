#!/bin/env R

library(tidyverse)

# load in target predictions

target_predictions = readr::read_tsv(snakemake@input[['targets']])

utr_genomic_positions = readr::read_tsv(snakemake@input[['bed_file']], col_names=c('chromosome','start','stop','strand','transcript'), 
	col_types='ciicc')

utr_genomic_positions = utr_genomic_positions[!duplicated(utr_genomic_positions$transcript),]

print(target_predictions)
print(utr_genomic_positions)

target_predictions = dplyr::inner_join(target_predictions, utr_genomic_positions, by=c(`Gene ID`='transcript'))
target_predictions$genomic_start = target_predictions$start + target_predictions$`UTR start`
target_predictions$genomic_end = target_predictions$start + target_predictions$`UTR end`

print(target_predictions %>% select(`Gene ID`,`Mirbase ID`,`UTR start`,`UTR end`,start,stop,genomic_start,genomic_end), width=Inf)

write.table(x=target_predictions, file=snakemake@output[[1]], sep='\t', quote=FALSE, col.names=TRUE, row.names=FALSE)
