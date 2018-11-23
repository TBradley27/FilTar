library(readr)
library(plyr)
library(dplyr)
library(stringr)
library(purrr)

ts_sites = read_tsv(snakemake@input[[1]])

sort_species = function (string) {
  sorted_string = string %>% 
    strsplit(split=' ', fixed=TRUE) %>% 
    unlist %>% 
    as.numeric %>% 
    unique %>%
    sort(method='quick') %>% 
    paste(collapse=" ")
  return(sorted_string)
}

#ts_sites$Species_in_this_group = map(
#  ts_sites$Species_in_this_group, 
#  sort_species
#  ) %>% 
#  plyr::ldply(data.frame) %>% unlist

# reference in order to sort site-types into correct ordering
y = c('7mer-1a','7mer-m8','8mer-1a','6mer')

remove_redundant_site_types = function (record) {
  new_record = record %>% strsplit('+', fixed=TRUE) %>% unlist %>% unique 
  new_record =  new_record[order(match(new_record,y))] %>%
    paste(collapse="+")
  return (new_record)
}


dummy = vector(length=dim(ts_sites)[1])
dummy2 = vector(length=dim(ts_sites)[1])
dummy3 = vector(length=dim(ts_sites)[1])

for (i in 1:dim(ts_sites)[1] ) {
        print(  (i/dim(ts_sites)[1]) * 100 )

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


#ts_sites$Group_type = map(ts_sites$Group_type, remove_redundant_site_types) %>%
#  plyr::ldply(data.frame) %>% unlist %>% as.character

#for (i in 1:dim(ts_sites)[1]) {
#  if (grepl('\\+', ts_sites$Group_type[i]) ) {
#    ts_sites$Species_in_this_group_with_this_site_type[i] = 
#      ts_sites$Species_in_this_group_with_this_site_type[i] %>% 
#      str_split(' ') %>% 
#      unlist %>% 
#      as.numeric %>% 
#      unique %>%
#      sort %>% 
#      paste(collapse=" ")
#    
#  } else { 
#    
#    ts_sites$Species_in_this_group_with_this_site_type[i] = ''
#    print(ts_sites)
#    }
#}


write.table(ts_sites,
file=snakemake@output[[1]],
sep="\t",
row.names=FALSE,
quote=FALSE)
