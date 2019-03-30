#!/bin/env Rscript

if ( file.info(snakemake@input$extended_bed)$size == 0 ) {
	file.copy(from=snakemake@input$normal_bed, to=snakemake@output[[1]])	
} else {
	full_set_sorted = filtar::get_full_bed(snakemake@input$normal_bed, snakemake@input$extended_bed, snakemake@input$all_transcripts, snakemake@input$tx_quant)
        write.table(full_set_sorted, file=snakemake@output[[1]], quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t")

}

write.table(full_set_sorted, file=snakemake@output[[1]], quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t")
