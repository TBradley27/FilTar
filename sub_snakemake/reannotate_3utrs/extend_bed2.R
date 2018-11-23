#!/bin/env Rscript

library(tidyverse)
library(plyr)

source('extended_bed_function2.R')

full_set_sorted = get_full_bed (snakemake@input['normal_bed'], snakemake@input['extended_bed'])

write.table(full_set_sorted, file=snakemake@output[[1]], quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t")
