library(plyr)
library(tidyverse)

merge_multiple_bed_files = function(normal_bed_file, specific_bed_file_dir) {

	# read in canonical annotation
	tx_ids = read_tsv(
		normal_bed_file, col_types='ciiic',
		col_names=c('chrom','start','end','strand','id')
		)

	# split into positive and negative strand subsets
	tx_pos = tx_ids %>% dplyr::filter(strand == 1) %>% dplyr::select(id) %>% unique()
	tx_pos_starts = tx_ids %>% dplyr::filter(strand == 1) %>% dplyr::select(id) %>% unique()
	tx_neg = tx_ids %>% dplyr::filter(strand == -1) %>% dplyr::select(id) %>% unique()
	tx_neg_starts = tx_ids %>% dplyr::filter(strand == -1) %>% dplyr::select(id) %>% unique()	

	print(tx_pos)

	# Get a list of available bed files
	bed_files = dir(path = specific_bed_file_dir, pattern="*.bed$")
	bed_files = c(bed_files, normal_bed_file)
	print(bed_files)
	#bed_files = bed_files[!bed_files == snakemake@input[[1]] ] # remove canonical annotation

	# create space in objects for more annotations
	for (i in bed_files) {tx_pos[[i]] = NA; tx_neg[[i]] = NA; tx_pos_starts[[i]] = NA ; tx_neg_starts[[i]] = NA}

	# remove version numbers from the canonical annotations
	tx_pos$id = gsub('\\..*','',tx_pos$id)
	tx_neg$id = gsub('\\..*','',tx_neg$id)
	tx_pos_starts$id = gsub('\\..*','',tx_pos_starts$id)
	tx_neg_starts$id = gsub('\\..*','',tx_neg_starts$id)

	print(tx_pos)
	print(tx_neg)

	# iterate over all available bed files
	for (i in bed_files) {

	  print(i)
	  
	  file_path = paste(specific_bed_file_dir,'/',str_interp("${i}"),sep="")

	  if (i == normal_bed_file) { 
	  	tmp_bed = tx_ids
	  	tmp_bed$id = gsub('\\..*','',tmp_bed$id)
	  } else {
	  	tmp_bed = read_tsv( file_path, col_types='ciiic', col_names=c('chrom','start','end','strand','id'))
	  }
	  
	  # reformat tx identifier variable
	  tmp_bed$id = gsub('\\|.*','', tmp_bed$id)

	  # iterate over every row in the bed file
	  # populate tables with additional end co-ordinates  
	  for (j in 1:dim(tmp_bed)[1]) { 
	    if (tmp_bed$strand[j] == '1') {
	      tx_pos[[i]][tx_pos$id == tmp_bed$id[j] ] = as.numeric(tmp_bed$end[j])
	      tx_pos_starts[[i]][tx_pos_starts$id == tmp_bed$id[j] ] = as.numeric(tmp_bed$start[j])
	      } 
	    else {
	      tx_neg[[i]][tx_neg$id == tmp_bed$id[j] ] = as.numeric(tmp_bed$start[j])
	      tx_neg_starts[[i]][tx_neg_starts$id == tmp_bed$id[j] ] = as.numeric(tmp_bed$end[j])
	      #print(tx_neg[,c(1,2)])
	      }
	    }
	}

	print('foo')

	print(tx_pos)
        print(tx_neg)

	print ('tx_pos_start')
	print (tx_pos_starts)

	# make sure to exclude the first column - get the max or min of each row
	tx_pos_starts$start <- apply(tx_pos_starts[,-c(1)], 1, function(x) min(x, na.rm = TRUE))
	tx_pos$end <- apply(tx_pos[,-c(1)], 1, function(x) max(x, na.rm = TRUE))
	tx_neg_starts$start <- apply(tx_neg_starts[,-c(1)], 1, function(x) max(x, na.rm = TRUE))
	tx_neg$end <- apply(tx_neg[,-c(1)], 1, function(x) min(x, na.rm = TRUE))

	tx_ids$id = gsub('\\..*','',tx_ids$id)

	print('bar')

	reannotate = function (tx_id) {

		x = filter(tx_ids, id == tx_id)

		if (dim(x)[1] > 1) {
			print('noo')
			return(x)
		} else {
			if (x$strand[1] == 1) {
				x$start = tx_pos_starts$start[ tx_pos_starts$id == tx_id ]
				x$end = tx_pos$end[ tx_pos$id == tx_id ]
			}

			else if (x$strand[1] == -1) {
				x$start = tx_neg$end[ tx_neg$id == tx_id ]
				x$end = tx_neg_starts$start[ tx_neg_starts$id == tx_id ]
			}
			print('goo')

			return(x)
		}
	}

	tx_ids_unique = tx_ids$id %>% unique

	new_annotations = map(tx_ids_unique, reannotate)

	new_annotations = ldply(new_annotations, data.frame) %>% as.tibble()

	return(new_annotations)
}
