input = readr::read_tsv(snakemake@input[[1]], col_names=FALSE, col_types='ciiciciiiicc')

if (length(snakemake@config[['transcripts']]) == 0) {
        file.copy(from=snakemake@input[[1]],to=snakemake@output[[1]])
} else {
	transcripts = gsub('\\..*','', snakemake@config[['transcripts']])

	filtered_input = input[input$X4 %in% transcripts,]

	write.table(filtered_input,snakemake@output[[1]], col.names=FALSE, row.names=FALSE, quote=FALSE, sep="\t")
}
