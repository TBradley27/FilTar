data = readr::read_tsv(snakemake@input[[1]], col_names=FALSE)
data = tidyr::separate(data, X5, into=c('X5.1','X5.2'), remove=TRUE)
data = tidyr::separate(data, X6, into=c('X6.1','X6.2'), remove=TRUE)

colnames(data) = c('miRNA_ID','transcript_ID','score','energy(kCal/Mol)','miRNA_start','miRNA_end','3UTR_start','3UTR_end','alignment_length','percent_matches', 'percent_matches_and_wobbles')

write.table(
	x=data,
	file=snakemake@output[[1]],
	quote=FALSE,
	row.names=FALSE,
	col.names=TRUE,
	sep="\t"
)
