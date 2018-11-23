#!/bin/env Rscript

library(plyr)
library(tidyverse)
source('sub_snakemake/target_prediction/targetscan/get_mirna_family_function.R')

print(snakemake@wildcards)

mirna_seeds = get_mirna_family(snakemake@input[[1]], snakemake@wildcards$species)

write.table(mirna_seeds, snakemake@output[[1]], sep="\t", quote=FALSE, col.names=FALSE, row.names=FALSE)
