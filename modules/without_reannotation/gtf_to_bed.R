
data = readr::read_tsv(
	snakemake@input[[1]],
	comment='##',
	col_types='ccciicccc', 
	col_names=c('chromosome','source','feature','start','stop','unknown1','strand','unknown2','info')
	)

data = dplyr::filter(data, feature=='three_prime_UTR')
data = tidyr::separate(
	data=data,
	col='info',
	into=c('ID','Parent','gene_id','transcript_id','gene_type','gene_name','transcript_type','transcript_name','exon_number','exon_id','protein_id'),
	sep=';',
	remove=TRUE
	)
data = dplyr::select(data, c('chromosome','start','stop','strand','transcript_id'))
data$strand = gsub('\\+',1,data$strand)
data$strand = gsub('-',-1,data$strand)
data$transcript_id = gsub('transcript_id=','',data$transcript_id)

print(data)

write.table(
	x=data,
	file=snakemake@output[[1]],
	sep="\t",
	col.names=FALSE,
	row.names=FALSE,
	quote=FALSE
	)
