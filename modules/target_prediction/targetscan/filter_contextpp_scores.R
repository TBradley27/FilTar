filtered_contextpp_scores = 
	filtar::filter_contextpp_scores(
					snakemake@input[['contextpp_scores']], 
					snakemake@input[['expression_values']],
					snakemake@params[['tpm_expression_threshold']]
	)

write.table(filtered_contextpp_scores, snakemake@output[[1]], row.names=FALSE, col.names=TRUE, sep="\t", quote=FALSE)
