library(tidyverse)
library(tximport)

### read in the 3UTR length data ###

utr_lengths = read_tsv(snakemake@input[['canonical']], col_names=TRUE) 
pc_transcripts = read_tsv(snakemake@input[['pc_transcripts']], col_names='tx_id')

utr_lengths = utr_lengths[utr_lengths$tx_id %in% pc_transcripts$tx_id,]

### read in the expression data ###

real = strsplit(snakemake@input$quant, split='/')
real_accessions = c()

for (i in 1:length(real)) {
  real_accessions = c(real_accessions, real[[i]][3])
}

files = c(snakemake@input$quant)
names(files) = real_accessions

txi <- tximport(files, type="kallisto", txOut=TRUE)

print('foo')

TPMs = as.data.frame(txi$abundance)
TPMs$average = rowMeans(TPMs)
TPMs$names = rownames(txi$abundance)

#### expression filtering stats ######

medium_expression = dplyr::filter(TPMs, average >= 0.1)
utr_lengths_filtered = utr_lengths[utr_lengths$tx_id %in% medium_expression$names,]

## tx_stat ##

print('START expression filtering stats')

length(utr_lengths$tx_id) %>% print()
length(utr_lengths_filtered$tx_id) %>% print()

transcripts_lost = length(utr_lengths$tx_id) - length(utr_lengths_filtered$tx_id)

proportion = transcripts_lost / length(utr_lengths$tx_id)

percentage = proportion * 100

### nucleotide stats ####

sum(utr_lengths$utr_length) %>% print()
sum(utr_lengths_filtered$utr_length) %>% print()

bases_lost = sum(utr_lengths$utr_length) - sum(utr_lengths_filtered$utr_length)

proportion_bases = bases_lost / sum(utr_lengths$utr_length)

percentage_bases = proportion_bases * 100

data = tibble(
        bases_lost = bases_lost / 1000000,
        percentage_bases_lost = percentage_bases,
	tx_lost=transcripts_lost,
	percentage_tx_lost=percentage,
	)

write.table(
	x=data,
	file=snakemake@output[[1]],
	row.names=FALSE,
	quote=FALSE,
	sep='\t'
)




