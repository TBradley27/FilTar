library(tidyverse)
library(tximport)

## read in the data

original = read_tsv(snakemake@input[['canonical']])
new = read_tsv(snakemake@input[['new']])
pc_transcripts = read_tsv(snakemake@input[['pc_transcripts']], col_names='tx_id')

merged = merge(original,new, by = 'tx_id') %>% as_tibble()
merged = merged[merged$tx_id %in% pc_transcripts$tx_id,]

merged$diff = merged$utr_length.y - merged$utr_length.x # get change in 3UTR length

### read in the transcript quant data
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

high_expression = dplyr::filter(TPMs, average >= 5.0)

merged_filtered = merged[merged$tx_id %in% high_expression$names,] # remove lowly expressed transcripts

print(merged)
print(merged_filtered)

###

num_new_bases = sum(merged$diff[merged$diff>0])
num_old_bases = sum(merged_filtered$diff[merged_filtered$diff<0]) %>% abs()

# reannotation integer vector

merged$int_diff = merged$diff
merged_filtered$int_diff = merged_filtered$diff

merged = merged %>% mutate(int_diff = replace(int_diff, int_diff > 0, 1))
merged_filtered = merged_filtered %>% mutate(int_diff = replace(int_diff, int_diff < 0, -1))

num_tx_elongated = sum(merged$int_diff[merged$int_diff>0])
num_tx_truncated = sum(merged_filtered$int_diff[merged_filtered$int_diff<0]) %>% abs()

percent_bases_plus = ( num_new_bases / sum(merged$utr_length.x) ) * 100
percent_bases_minus = ( num_old_bases / sum(merged$utr_length.x) ) * 100 %>% abs()

percent_tx_elongated = ( num_tx_elongated / length(merged$tx_id ) ) * 100
percent_tx_truncated = ( num_tx_truncated / length(merged$tx_id ) ) * 100


data = tibble(
	sample=snakemake@wildcards$tissue, 
	num_new_bases = num_new_bases, 
	num_old_bases = num_old_bases, 
	percent_bases_gained = signif(percent_bases_plus,3),
	percent_bases_lost = signif(percent_bases_minus,3), 
	num_tx_elongated = num_tx_elongated, 
	num_tx_truncated = num_tx_truncated,
	percent_tx_elongated = signif(percent_tx_elongated,3),
	percent_tx_truncated = signif(percent_tx_truncated,3)
	)

write.table(
	x=data,
	file=snakemake@output[[1]],
	row.names=FALSE,
	quote=FALSE,
	sep='\t'
)




