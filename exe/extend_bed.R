#!/bin/env Rscript

library(tidyverse)

args = commandArgs(trailingOnly=TRUE)

# read in the data
normal_bed = read_tsv(args[1], col_names=c('chromosome','start','stop','strand','id'), col_types=list('c','i','i','c','c'))
extended_utrs = read_tsv(args[2], col_names=c('chromosome','start','stop','id','dummy','strand'))

# tidy the data
normal_bed = separate(normal_bed, id, into=c('id','version'))
extended_utrs = separate(extended_utrs, id, into=c('id','dummy2','chrom_dup','strand_dup'))

normal_bed$version = NULL

extended_utrs = extended_utrs[,c('chromosome','start','stop','strand','id')]
extended_utrs$chromosome = stringr::str_replace_all(extended_utrs$chromosome, 'chr','')

# identify the transcripts which do and do not have APAtrap records
not_extended_utrs = normal_bed[!(normal_bed$id %in% extended_utrs$id),]
extended_utrs_normal = normal_bed[(normal_bed$id %in% extended_utrs$id),]


# get the terminal exons of all transcripts with APAtrap records

get_terminal_exons_only = function (string, df, id_column) {

  x = subset(df, df$id == string)

  if (x$strand[1] == '1') {
    y = subset(x, x$stop ==max(x$stop)) }
  else if (x$strand[1] == '-1') {
    y = subset(x, x$start == min(x$start))
  }
  return (y)
}

transcript_IDs = extended_utrs_normal$id %>% unique()

terminal_exons = map(transcript_IDs, get_terminal_exons_only, extended_utrs_normal)
terminal_exons = ldply(terminal_exons, data.frame) %>% as.tibble() # merge into one data frame

# Determine if APAtrap results in an extension of specified transcripts
# Extend the original annotation were applicable

compare_gencode_with_APAtrap = function (string) {
  x = dplyr::filter(terminal_exons, id == string)
  y = dplyr::filter(extended_utrs, id == string)
  
  if (x$strand == "1") {
    if (y$stop > x$stop) {x$stop = y$stop}
    else if (y$stop <= x$stop) {}
  } else if (x$strand == "-1") {
    if (y$start < x$start) {x$start = y$start}
    else if (y$start >= x$start) {}
  }
  
  return (x)
}

term_exons_reannotated = map(transcript_IDs, compare_gencode_with_APAtrap) 
term_exons_reannotated = ldply(term_exons_reannotated, data.frame) %>% as.tibble() # merge into one data frame

nonterminal_exons = dplyr::anti_join(normal_bed, terminal_exons) # get non-terminal exons
reannotated_utrs = rbind(nonterminal_exons, term_exons_reannotated)

print('testing')

# Determine if APAtrap results in a shortening of specified transcripts

shorten_UTRs = function (string) {
   # subet data for a particular transcript of interest
   renannoted_utrs = dplyr::filter(reannotated_utrs, id == string)
   extended_utrs = dplyr::filter(extended_utrs, id == string)

   # if APAtrap has not produced new annotations for this record - then return the bed records as they exists
   if (nrow(extended_utrs) == 0) {
	return (renannoted_utrs) } else {

   # if APAtrap records for this transcript exist - then update the bed records
   if(renannoted_utrs$strand[1] == "1") {
   	renannoted_utrs = subset(renannoted_utrs, renannoted_utrs$start < extended_utrs$stop)  # remove redundant trailing exons
        renannoted_utrs$stop[length(renannoted_utrs$stop)] = extended_utrs$stop # truncate the terminal exon
        }
   else if (renannoted_utrs$strand[1] == "-1") {
	renannoted_utrs = subset(renannoted_utrs, renannoted_utrs$stop > extended_utrs$start) 
        renannoted_utrs$start[length(renannoted_utrs$start)] = extended_utrs$start
        }
        print ("extended_utrs")
        print (extended_utrs)
        print ("normal_utrs")
	print (renannoted_utrs)
	return(renannoted_utrs)
	}
}
   
   # loop over all exons for the GENCODE annotation - find which exon the extended UTR lies
   #for (i in 1:length(normal_bed$id)) {
   #  if(normal_bed$strand == "1") {
   #    if (extended_utrs$stop[1] > normal_bed$start[i] & extended_utrs$stop[1] < normal_bed$stop[i]) {
   #       print (i)
   #       print (extended_utrs$stop[1])
   #       print (normal_bed$start[i])
   #       print (normal_bed$stop[i])
   #} else {}
   #}
#}
#   return(normal_bed)
#}

print('testing2')

extended_utrs_normal = reannotated_utrs[(reannotated_utrs$id %in% extended_utrs$id),]
transcript_IDs = extended_utrs_normal$id %>% unique()

shortened_utrs = map(transcript_IDs, shorten_UTRs)
shortened_utrs = ldply(shortened_utrs, data.frame) %>% as.tibble()

nonshortened_utrs = dplyr::anti_join(reannotated_utrs, extended_utrs_normal) 
reannotated_utrs = rbind(nonshortened_utrs, shortened_utrs)

# reformat bed file

reannotated_utrs = reannotated_utrs[order(reannotated_utrs$id, decreasing=FALSE),]
reannotated_utrs = reannotated_utrs[order(reannotated_utrs$start, decreasing=FALSE),]

desired_order = c('1','2','3','4','5','6','7','X','8','9','11','10','12','13','14','15','16','17','18','20','19','Y','22','21','KI270728.1',
'KI270727.1',
'GL000009.2',
'GL000194.1',
'GL000205.2',
'GL000195.1',
'KI270734.1',
'GL000213.1',
'GL000218.1',
'KI270731.1',
'KI270721.1',
'KI270711.1',
'KI270713.1')

reannotated_utrs$chromosome <- factor( as.character(reannotated_utrs$chromosome), levels=desired_order )

reannotated_utrs = reannotated_utrs[order(reannotated_utrs$chromosome, decreasing=FALSE),]

reannotated_utrs$chromosome = as.character(reannotated_utrs$chromosome)

tx_IDs = reannotated_utrs$id %>% unique()

reorder_bed_files = function (tx_ID) {
	transcript_specific_bed_records = subset(reannotated_utrs, reannotated_utrs$id == tx_ID)
        #print(transcript_specific_bed_records)
        if (transcript_specific_bed_records$strand[1] == '1') {
           transcript_specific_bed_records = transcript_specific_bed_records[order(transcript_specific_bed_records$start, decreasing=FALSE),]
	} 
        else if (transcript_specific_bed_records$strand[1] == '-1' ) {
           transcript_specific_bed_records = transcript_specific_bed_records[order(transcript_specific_bed_records$start, decreasing=TRUE),]
	}
        #print(transcript_specific_bed_records)
	return (transcript_specific_bed_records)
	
}

reannotated_utrs_sorted = map(tx_IDs, reorder_bed_files)
reannotated_utrs_sorted = ldply(reannotated_utrs_sorted, data.frame) %>% as.tibble()

#print(reannotated_utrs_sorted)

# write
write.table(reannotated_utrs_sorted, file=args[3], quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t")
