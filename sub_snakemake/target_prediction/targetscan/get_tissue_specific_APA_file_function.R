library(plyr)
library(tidyverse)

APAtrap_output = read.table(snakemake@input[[1]], sep="\t",
                            header=TRUE) %>% as.tibble()

APAtrap_output$Group_1_2_Separate_Exp = NULL # for brevity
APAtrap_output$Group_1_2_Total_Exp = NULL

APAtrap_output = APAtrap_output[1:1268,] # remove incomplete rows

reposition_last_APA_site = function (APA_df) {
  
  APA_df$strand = str_sub(APA_df$Gene, -1, -1)
  APA_df$last_APA = ''
  APA_df$start_pos = ''
  
  # Account for strandedness
  for (i in 1:dim(APA_df)[1]) {
    if (APA_df$strand[i] == '+') {
      
      APA_df$last_APA[i] = sub('.*-', '', APA_df$Loci[i])
      APA_df$start_pos[i] = sub('[0-9A-Z]*:', '', APA_df$Loci[i]) %>% 
        sub(pattern='-.*',replacement='')
      
    } else if (APA_df$strand[i] == '-')
    {
      
      APA_df$last_APA[i] = sub('[0-9A-Z]*:', '', APA_df$Loci[i]) %>% 
        sub(pattern='-.*',replacement='')
      APA_df$start_pos[i] = sub('.*-', '', APA_df$Loci[i])
      
    }
  }
  APA_df$Predicted_APA = str_glue_data(APA_df,"{Predicted_APA},{last_APA}")
  
  # clean up
  APA_df$last_APA = NULL
  return (APA_df)
}

APAtrap_output = reposition_last_APA_site(APAtrap_output)

#convert_to_percentages = function (exp_string, total) {
#  expression_values = str_split(exp_string, ',') %>% unlist %>% as.numeric
  #percentage_values = ( expression_values / total ) * 100
#  return (percentage_values)
#}

remove_unused_APAsites = function (APA_df) {
  converted_exp_column = map(APA_df$Group_1_1_Separate_Exp, function (x)
    str_split(x, ',') %>% unlist %>% as.numeric)
  
  converted_loci_column = map(APA_df$Predicted_APA, function (x)
    str_split(x, ',') %>% unlist)
  
  for (i in 1:length(converted_exp_column)) {
    index = converted_exp_column[[i]] != 0
    converted_exp_column[[i]] = converted_exp_column[[i]][index]
    converted_loci_column[[i]] = converted_loci_column[[i]][index]
  }
  
  APA_df$a = converted_exp_column
  APA_df$b = converted_loci_column
  
  return (APA_df)
  
}

APAtrap_output = remove_unused_APAsites(APAtrap_output)

get_relative_abundances = function (exp_string, total, tx_id) {
  
  get_cumulative_depreciation = function (percent_vec) {
    cumulative_depreciation = c()
    cumulative_depreciation = c(cumulative_depreciation, 100)
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

new = APAtrap_output
new$Gene = gsub('\\|.*','',APAtrap_output$Gene) 
new = new[,c('Gene','b','Loci','strand','start_pos')]

y = merge(y, new, by.x='id', by.y='Gene')

y$Loci = NULL
#y$strand = NULL

tx_ids = y$id %>% as.factor %>% levels

get_rel_APA_position = function(tx_id) {
  records = subset(y, y$id == tx_id)
  
  absolute_start = as.numeric(records$start_pos[1])
  
  for (i in 1:dim(records)[1] ) {
    
    if (i == 1) {
      records$rel_start_pos[i] = 1 
      records$rel_end_pos[i] = as.numeric(records$b[[1]][i]) - absolute_start + 1
    } else {
      records$rel_start_pos[i] = records$rel_end_pos[i-1] + 1
      records$rel_end_pos[i] = as.numeric(records$b[[1]][i]) - absolute_start + 1
    }
  }
  
  records$b = NULL
  records$start = NULL
  records = records[,c(1,3,4,2)]
  colnames(records) = c('id','rel_start_pos','rel_end_pos','AIR')
  
  return (records)
}

APA_records = map(tx_ids, get_rel_APA_position) %>% plyr::ldply(data.frame) %>% as.tibble

write.table(APA_records, snakemake@output[[1]], sep="\t", col.names=FALSE, row.names=FALSE, quote=FALSE)
