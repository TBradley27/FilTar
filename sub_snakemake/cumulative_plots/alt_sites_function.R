get_alt_sites = function(targets, utr_lens) {

	targets = read_tsv(
		file = targets,
		col_names = TRUE,
		col_types = 'ccciiiiiccccci',
		trim_ws = FALSE
		)

	alt_utr_lens =  read_tsv(
		file =utr_lens,
		col_names = TRUE,
		col_types = 'ci',
		n_max=Inf
		)

	filter_targets = function (id, utr_length) {
		sub_targets = filter(targets, a_Gene_ID == id)
		sub_targets = filter(sub_targets, UTR_end < utr_length)

		return (sub_targets)
	}

	alt_targets = map2(alt_utr_lens$tx_id, alt_utr_lens$utr_length, filter_targets)
	alt_targets = ldply(alt_targets, data.frame) %>% as.tibble()

	return(alt_targets)
}	
