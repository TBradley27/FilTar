input = readr::read_tsv(snakemake@input[[1]], col_names=FALSE)

print(snakemake@config[['transcripts']])

print(input)

if (length(snakemake@config[['transcripts']]) == 0) {

	file.copy(from=snakemake@input[[1]],to=snakemake@output[[1]])

} else {

	filtered_input = input[input$X5 %in% snakemake@config[['transcripts']],]
	print(filtered_input)
	write.table(filtered_input,snakemake@output[[1]], col.names=FALSE, row.names=FALSE, quote=FALSE, sep="\t")
}


