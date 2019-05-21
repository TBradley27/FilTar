mirna_seeds = filtar::get_mirna_family(snakemake@input[[1]], snakemake@wildcards$species)

write.table(mirna_seeds, snakemake@output[[1]], sep="\t", quote=FALSE, col.names=FALSE, row.names=FALSE)
