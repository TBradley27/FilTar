input = readr::read_tsv(snakemake@input[[1]], col_names=FALSE, col_types='ciiic')

if (length(snakemake@config[['transcripts']]) == 0) {

	file.copy(from=snakemake@input[[1]],to=snakemake@output[[1]])

} else {

	for (transcript in snakemake@config[['transcripts']]) {
		if (!transcript %in% input$X5) {
			write(stringr::str_interp("The transcript identifer '${transcript}' is not a valid identifier for the selected species for ensembl release ${snakemake@config[['ensembl_release']]}"), stderr())
		}
	}

	filtered_input = input[input$X5 %in% snakemake@config[['transcripts']],]

	if (dim(filtered_input)[1] == 0) {
		write("No valid transcript identifiers in input. Halting execution.", stderr())
		quit(save="no", status=1)
	}

	write.table(filtered_input,snakemake@output[[1]], col.names=FALSE, row.names=FALSE, quote=FALSE, sep="\t")
}


