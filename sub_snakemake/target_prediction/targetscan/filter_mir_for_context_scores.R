
filtar::filter_mature_mirs(snakemake@input[[1]], paste(snakemake@wildcards[["species"]],snakemake@wildcards[["miRNA"]],sep='-'), snakemake@output[[1]])
