library(readr)
library(plyr)
library(dplyr)
library(stringr)
library(purrr)

f = function() {

ts_sites = read_tsv("results/hsa_chrY_msa.sites.testing.tsv")

sort_species = function (string) {
  sorted_string = string %>% 
    str_split(' ') %>% 
    unlist %>% 
    as.numeric %>% 
    unique %>%
    sort %>% 
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
  new_record = record %>% str_split('\\+') %>% unlist %>% unique 
  new_record =  new_record[order(match(new_record,y))] %>%
    paste(collapse="+")
  return (new_record)
}


for (i in 1:dim(ts_sites)[1] ) {
        print(  (i/dim(ts_sites)[1]) * 100 )
        ts_sites$Species_in_this_group[i] = sort_species(ts_sites$Species_in_this_group[i])
	ts_sites$Group_type[i] = remove_redundant_site_types(ts_sites$Group_type[i])
}

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
    
#  } else { 
    
#    ts_sites$Species_in_this_group_with_this_site_type[i] = ''
#    print(ts_sites)
#    }
#}


write.table(ts_sites,
file="results/hsa_chrY_msa.sites2.tsv",
sep="\t",
row.names=FALSE,
quote=FALSE)

}

f()




