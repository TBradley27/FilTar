get_canonical_AIRs = function(utr_lengths) {
	AIRs = read_tsv(utr_lengths, col_names=TRUE)
	AIRs$start = 1
	AIRs$AIR = 100
	AIRs = AIRs[,c('tx_id','start','utr_length','AIR')]

	AIRs$tx_id = gsub('\\.[0-9][0-9]?','',AIRs$tx_id) # remove version number - we may change this later

	return(AIRs)
}
