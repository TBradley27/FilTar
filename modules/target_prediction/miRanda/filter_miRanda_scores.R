if ( file.info(snakemake@input[['miRanda_scores']])$size == 0  ) {
	file.create(snakemake@output[[1]])	
} else {

	filtered_miRanda_scores = filtar::filter_miRanda_scores(snakemake@input[['miRanda_scores']], snakemake@input[['expression_values']],snakemake@params[['tpm_expression_threshold']])
	write.table(filtered_miRanda_scores, snakemake@output[[1]], row.names=FALSE, col.names=TRUE, sep="\t", quote=FALSE)
}
