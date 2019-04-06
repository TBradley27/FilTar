#!/usr/bin/env Rscript

if ( file.info(snakemake@input[[1]])$size == 0 ) {
        file.copy(from=snakemake@input[[1]], to=snakemake@output[[1]])
} else {
	output = filtar::MergeFasta(snakemake@input[[1]])
	write.table(output, snakemake@output[[1]], quote=FALSE, sep="\t", col.names=FALSE, row.names=FALSE)
}
