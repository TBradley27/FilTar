#!/bin/env Rscript

library(plyr)
library(tidyverse)

source('../sub_snakemake/reannotate_3utrs/extend_bed_function2.R')

print(snakemake@input)

full_set_sorted = get_full_bed (snakemake@input$normal_bed, snakemake@input$extended_bed)

write.table(full_set_sorted, file=snakemake@output[[1]], quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t")
