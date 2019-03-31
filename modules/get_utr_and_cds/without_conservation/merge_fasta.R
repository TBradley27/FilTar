#!/usr/bin/env Rscript

output = filtar::MergeFasta(snakemake@input[[1]])
write.table(output, snakemake@output[[1]], quote=FALSE, sep="\t", col.names=FALSE, row.names=FALSE)
