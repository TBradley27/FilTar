num_lines = length(readLines(snakemake@input[[1]]))

if ( num_lines  == 1 ) { # test if the file contains any records 

	file.copy(from=snakemake@input[[1]], to=snakemake@output[[1]])

} else {

	ts_sites = filtar::fix_ts_output(snakemake@input[[1]])

	write.table(ts_sites,
	file=snakemake@output[[1]],
	sep="\t",
	row.names=FALSE,
	quote=FALSE)

}
