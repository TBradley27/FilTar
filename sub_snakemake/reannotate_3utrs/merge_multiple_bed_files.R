library(plyr)
library(tidyverse)

#setwd("~/Documents/tmp/APAtrap/merge_multiple_bed_files")

# read in canonical annotation
tx_ids = read_tsv(snakemake@input[[1]], col_types='ciiic',
	col_names=c('chrom','start','end','strand','id'))

# split into positive and negative strand subsets
tx_pos = tx_ids %>% dplyr::filter(strand == 1) %>% dplyr::select(id) %>% unique()
tx_neg = tx_ids %>% dplyr::filter(strand == -1) %>% dplyr::select(id) %>% unique()

print(tx_pos)

# Get a list of available bed files
bed_files = dir(path = snakemake@params[[1]], pattern="*full.bed$")
bed_files = c(bed_files, snakemake@input[[1]])
print(bed_files)
#bed_files = bed_files[!bed_files == snakemake@input[[1]] ] # remove canonical annotation

# create space in objects for more annotations
for (i in bed_files) {tx_pos[[i]] = NA; tx_neg[[i]] = NA}

# remove version numbers from the canonical annotations
tx_pos$id = gsub('\\..*','',tx_pos$id)
tx_neg$id = gsub('\\..*','',tx_neg$id)

print(tx_pos)
print(tx_neg)

# iterate over all available bed files
for (i in bed_files) {

  print(i)
  
  file_path = paste(snakemake@params[[1]],'/',str_interp("${i}"),sep="")

  if (i == snakemake@input[[1]]) { 
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
      } 
    else {
      tx_neg[[i]][tx_neg$id == tmp_bed$id[j] ] = as.numeric(tmp_bed$start[j])
      #print(tx_neg[,c(1,2)])
      }
    }
}


# make sure to exclude the first column - get the max or min of each row
tx_pos$end <- apply(tx_pos[,-c(1)], 1, function(x) max(x, na.rm = TRUE))
tx_neg$start <- apply(tx_neg[,-c(1)], 1, function(x) min(x, na.rm = TRUE))

tx_ids$id = gsub('\\..*','',tx_ids$id)

reannotate = function(tx_id) {

x = filter(tx_ids, id == tx_id)

if (x$strand[1] == 1) {
	x$end[length(x$strand)] = tx_pos$end[ tx_pos$id == tx_id ] # amend final exon record
}

else if (x$strand[1] == -1) {
	x$start[length(x$strand)] = tx_neg$start[ tx_neg$id == tx_id ] # amend final exon record
}

print(x)

return(x)

}

tx_ids_unique = tx_ids$id %>% unique


print(tx_pos[,c(5,6)])
print(tx_neg[,c(5,6)])

print(tx_pos$end)
print(tx_neg$start)

new_annotations = map(tx_ids_unique, reannotate)

new_annotations = ldply(new_annotations, data.frame) %>% as.tibble()

#for (i in 1:length(tx_pos$id)) {
#  tx_pos$start[i] = tx_ids$start[ tx_ids$id == tx_pos$id[i] ]
#  tx_pos$chrom[i] = tx_ids$chrom[ tx_ids$id == tx_pos$id[i] ]
#}

#for (i in 1:length(tx_neg$id)) {
#  tx_neg$end[i] = tx_ids$end[ tx_ids$id == tx_neg$id[i] ]
#  tx_neg$chrom[i] = tx_ids$chrom[ tx_ids$id == tx_neg$id[i] ]
#}

###

#tx_pos$strand = '1'
#tx_neg$strand = '-1'

#tx_pos = tx_pos[,c('chrom','start','end','id','strand')]
#tx_neg = tx_neg[,c('chrom','start','end','id','strand')]

#tx_all = rbind(tx_pos,tx_neg)
#tx_all = tx_all[,c('chrom','start','end','strand','id')]

write.table(new_annotations, file=snakemake@output[[1]], sep="\t", quote=FALSE, 
            row.names=FALSE, col.names = FALSE)
