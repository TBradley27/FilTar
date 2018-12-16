AIR_file = filtar::get_AIR_file(snakemake@input[[1]], snakemake@input[[2]])

write.table(AIR_file, snakemake@output[[1]], sep="\t", col.names=FALSE, row.names=FALSE, quote=FALSE)
