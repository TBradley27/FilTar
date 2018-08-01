#!/bin/env Rscript

library(purrr)
library(readr)
library(plyr)
library(tibble)

args = commandArgs(trailingOnly = TRUE)

data = read_tsv(args[1])
#data$Name = NULL
#transposed_data = t(data)
num_transcripts = dim(data)[1]

read_sample_mappings = read_tsv('data/read_sample_mapping.tsv')

print('bar')

read_sample_mappings = read_sample_mappings[read_sample_mappings$Run %in% colnames(data),]

print ('car')

get_sample_averages = function (sample) {

	read_sample_mappings = dplyr::filter(read_sample_mappings, secondary_sample_accession == sample)
	data = dplyr::select(data, read_sample_mappings$Run)	
	transposed_data = t(data)

	avg_values = map(1:num_transcripts, function (x) mean(transposed_data[,x]) )
	avg_values = plyr::ldply(avg_values, data.frame) %>% round(digits=2) %>% as.tibble()

	return (avg_values)
} 

get_all_sample_averages = function (vec_of_samples) { #NB: there is some recursion here which I am not taking advantage of.
        tmp = as.tibble(matrix(ncol=0, nrow=num_transcripts))
	for (sample in vec_of_samples) {
		tmp = cbind(tmp,get_sample_averages(sample))
	}	

	transposed_data = t(tmp)

	avg_values = map(1:num_transcripts, function (x) mean(transposed_data[,x]) )
        avg_values = plyr::ldply(avg_values, data.frame) %>% round(digits=2) %>% as.tibble()

	return(avg_values)
}

sample_vec = levels(as.factor(read_sample_mappings$secondary_sample_accession))

print(sample_vec)

avg_tpm_column = get_all_sample_averages(sample_vec)

print ('foo')

data_with_avg = cbind(data,avg_tpm_column)
colnames(data_with_avg) = c(colnames(data),'avg')

write.table(data_with_avg, file=args[2], sep="\t", quote=FALSE, row.names=FALSE)
