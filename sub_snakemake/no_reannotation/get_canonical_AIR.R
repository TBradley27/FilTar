AIRs = filtar::get_canonical_AIRs(snakemake@input[[1]])

write.table(AIRs, snakemake@output[[1]], quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t")
