#!/bin/env Rscript

library(plyr)
library(tidyverse)
library(filtar)

print(snakemake@input)

print('foo')

full_set_sorted = get_full_bed (snakemake@input$normal_bed, snakemake@input$extended_bed, snakemake@input$all_transcripts)

write.table(full_set_sorted, file=snakemake@output[[1]], quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t")
