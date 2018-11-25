fix_ts_output = function (ts_sites_output) {

	ts_sites = read_tsv(ts_sites_output)

	sort_species = function (string) {
	  if (grepl(' ', string)) {
		 sorted_string = string %>% 
		 strsplit(split=' ', fixed=TRUE) %>% 
		 unlist %>% 
		 as.numeric %>% 
		 unique %>%
		 sort(method='quick') %>% 
		 paste(collapse=" ")
		 return(sorted_string)
	  } else {
		return (string)
	 }
	}

	# reference in order to sort site-types into correct ordering
	y = c('7mer-1a','7mer-m8','8mer-1a','6mer')

	remove_redundant_site_types = function (record) {
	#  if (grepl('+',record)) {
		new_record = record %>% strsplit('+', fixed=TRUE) %>% unlist %>% unique 
		new_record =  new_record[order(match(new_record,y))] %>%
		paste(collapse="+")
		return (new_record)
	#  }
	#  else {
	# 	return (record)
	#  }
	}

	dummy = vector(length=dim(ts_sites)[1])
	dummy2 = vector(length=dim(ts_sites)[1])
	dummy3 = vector(length=dim(ts_sites)[1])

	for (i in 1:dim(ts_sites)[1] ) {
		print(  (i/dim(ts_sites)[1]) * 100 ) # progress bar
		dummy[i] = sort_species(ts_sites$Species_in_this_group[i])
		dummy2[i] = remove_redundant_site_types(ts_sites$Group_type[i])

		if ( grepl('\\+', dummy2[i]) ) {
			dummy3[i] = 
			ts_sites$Species_in_this_group_with_this_site_type[i] %>% 
			strsplit(' ', fixed=TRUE) %>% 
			unlist %>% 
			as.numeric %>% 
			unique %>%
			sort() %>% 
			paste(collapse=" ")
		} else {
			dummy3[i]=''
		}

	}

	ts_sites$Species_in_this_group = dummy
	ts_sites$Group_type = dummy2
	ts_sites$Species_in_this_group_with_this_site_type = dummy3

	return(ts_sites)
}

