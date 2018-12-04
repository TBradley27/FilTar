library(tidyverse)
source('../sub_snakemake/target_prediction/targetscan/get_mirna_context_function.R')

test = get_mirna_context(snakemake@input$mirna_seed, snakemake@input$mature_mirnas, snakemake@wildcards$species)

write_tsv(test, snakemake@output[[1]], col_names=FALSE)
