get_AIR_file = function(APA_file,transcript_lengths_file) {
	APAtrap_output = read.table(APA_file, sep="\t",
				    header=TRUE) %>% as.tibble()

	#APAtrap_output = APAtrap_output[1:1268,] # remove incomplete rows - I have no idea what this line is for

	reposition_last_APA_site = function (APA_df) {
	  
	 APA_df$strand = str_sub(APA_df$Gene, -1, -1) # retrieve the last character in the column which is the transcript strand
	 APA_df$last_APA = ''                         # create space for the last APA position
	 APA_df$start_pos = ''                        # create space for what I assume is the transcript start position
	  
	 # Account for strandedness
	 for (i in 1:dim(APA_df)[1]) {                # iterate over rows
	   if (APA_df$strand[i] == '+') {
	      
	     APA_df$last_APA[i] = sub('.*-', '', APA_df$Loci[i]) # remove everything before the hyphen
	     APA_df$start_pos[i] = sub('[0-9A-Z]*:', '', APA_df$Loci[i]) %>% 
		sub(pattern='-.*',replacement='') # capture everything after the colon and before the hyphen
	      
	    } else if (APA_df$strand[i] == '-')
	    {
	      
	      APA_df$last_APA[i] = sub('[0-9A-Z]*:', '', APA_df$Loci[i]) %>% 
		sub(pattern='-.*',replacement='')
	      APA_df$start_pos[i] = sub('.*-', '', APA_df$Loci[i])
	      
	    }
	  }
	  APA_df$Predicted_APA = str_glue_data(APA_df,"{Predicted_APA},{last_APA}") # collect all APA sites together
	  
	  # clean up
	  APA_df$last_APA = NULL
	  return (APA_df)
	}

	APAtrap_output = reposition_last_APA_site(APAtrap_output)

#	print(APAtrap_output)
	#convert_to_percentages = function (exp_string, total) {
	#  expression_values = str_split(exp_string, ',') %>% unlist %>% as.numeric
	  #percentage_values = ( expression_values / total ) * 100
	#  return (percentage_values)
	#}

	# I think this function is to remove APA sites with a quoted abundance of 0.00 which in effect would mean that they are redundant - probably initially identified because they have the correct motif

	remove_unused_APAsites = function (APA_df) {
	  converted_exp_column = map(APA_df$Group_1_1_Separate_Exp, function (x)
	    str_split(x, ',') %>% unlist %>% as.numeric) # convert to numerical vector
	  
	  converted_loci_column = map(APA_df$Predicted_APA, function (x)
	    str_split(x, ',') %>% unlist) # convert to numerical vector
	  
	  for (i in 1:length(converted_exp_column)) { # loops over each row of the data frame
	    index = converted_exp_column[[i]] != 0 # A Boolean vector
	    converted_exp_column[[i]] = converted_exp_column[[i]][index]
	    converted_loci_column[[i]] = converted_loci_column[[i]][index]
	  }
	  
	  APA_df$a = converted_exp_column
	  APA_df$b = converted_loci_column
	  
	  return (APA_df)
	  
	}

	APAtrap_output = remove_unused_APAsites(APAtrap_output)
#	print(APAtrap_output)

	get_relative_abundances = function (exp_string, total, tx_id) {
	  
	  # function applied recursively
	  get_cumulative_depreciation = function (percent_vec) {
	    cumulative_depreciation = c() #initialise vector 
	    cumulative_depreciation = c(cumulative_depreciation, 100) # initialise vector 
	    for (i in 1:length(percent_vec) - 1) { 
	      percent_fall = percent_vec[i]
	      cumulative_depreciation = c(
		cumulative_depreciation, 
		tail(cumulative_depreciation, n=1) - percent_fall
	      )
	    }
	    return (cumulative_depreciation)
	  }
	    
	  #expression_values = str_split(exp_string, ',') %>% unlist %>% as.numeric
	  percentage_values = ( exp_string / total ) * 100

	  return (get_cumulative_depreciation(percentage_values) )
	  
	}

	# get relative abundances

	y = purrr::map2(APAtrap_output$a,
		    APAtrap_output$Group_1_1_Total_Exp,
		    get_relative_abundances)

	names(y) = APAtrap_output$Gene
	y = unlist(y) %>% as.data.frame
	y$id = rownames(y)
	rownames(y) = NULL
	y$id = gsub('\\|.*','',y$id)
	y = y[,c(2,1)] # reorder columns

#	print(y) # gene names and expression values without utr loci

	new = APAtrap_output
	new$Gene = gsub('\\|.*','',APAtrap_output$Gene) 
	new = new[,c('Gene','b','Loci','strand','start_pos')]

	y = merge(y, new, by.x='id', by.y='Gene')

	y$Loci = NULL
	#y$strand = NULL

	#print(y) # end position missing

	tx_ids = y$id %>% as.factor %>% levels

	get_rel_APA_position = function(tx_id) {
	  records = subset(y, y$id == tx_id)
	  
	  absolute_start = as.numeric(records$start_pos[1])
	  
	  for (i in 1:dim(records)[1] ) {
	    
	    if (i == 1) {
		if (records$strand[i] == '+') {
	      		records$rel_start_pos[i] = 1 
	      		records$rel_end_pos[i] = as.numeric(records$b[[1]][i]) - absolute_start + 1 # initial relative end equals first APA site - absolute start
		} else if (records$strand[i] == '-') {
			records$rel_start_pos[i] = 1
			records$rel_end_pos[i] = absolute_start - as.numeric(records$b[[1]][i]) + 1
		}
	    }   else {
			if (records$strand[i] == '+') {
	      			records$rel_start_pos[i] = records$rel_end_pos[i-1] + 1
	      			records$rel_end_pos[i] = as.numeric(records$b[[1]][i]) - absolute_start + 1
		} else if (records$strand[i] == '-') {
				records$rel_start_pos[i] = records$rel_end_pos[i-1] + 1
                                records$rel_end_pos[i] = absolute_start - as.numeric(records$b[[1]][i]) + 1
			}
	    }
	  }
	  #print(absolute_start) 
	  #print(records) 
	  records$b = NULL
	  records$start = NULL
	  records = records[,c(1,5,6,2)]
	  colnames(records) = c('id','rel_start_pos','rel_end_pos','AIR')
	  
	  return (records)
	}

	APA_records = map(tx_ids, get_rel_APA_position) %>% plyr::ldply(data.frame) %>% as.tibble

	all_transcripts = read_tsv(transcript_lengths_file, col_names=TRUE)
	print(APA_records)
	non_updated_tx = all_transcripts[!all_transcripts$tx_id %in% APA_records$id,]
	print(non_updated_tx)
	non_updated_tx$start = 1
	non_updated_tx$AIR = 100
	non_updated_tx = non_updated_tx[,c('tx_id','start','utr_length','AIR')]
	colnames(non_updated_tx) = c('id','rel_start_pos','rel_end_pos','AIR')
	print(non_updated_tx)

	APA_records = rbind(APA_records,non_updated_tx)
	
	return(APA_records)
}
