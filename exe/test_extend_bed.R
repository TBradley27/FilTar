#!/bin/env Rscript

library(readr)
library(stringr)
library(tidyr)
library(plyr)
library(dplyr)
library(tibble)
library(purrr)
library(crayon)

get_terminal_exons_only = function (string, df, id_column) {
  x = subset(df, df$id == string)

  if (x$strand[1] == '1') {
    y = subset(x, x$stop ==max(x$stop)) }
  else if (x$strand[1] == '-1') {
    y = subset(x, x$start == min(x$start))
  }
  return (y)
}

extend_UTRs = function (string) {
  x = dplyr::filter(terminal_exons, id == string)
  y = dplyr::filter(example_APAtrap_annotations, id == string)
  
  if (x$strand == "1") {
    if (y$stop > x$stop) {x$stop = y$stop}
    else if (y$stop <= x$stop) {}
  } else if (x$strand == "-1") {
    if (y$start < x$start) {x$start = y$start}
    else if (y$start >= x$start) {}
  }
  
  return (x)
}

shorten_UTRs = function (string) {
   # subet data for a particular transcript of interest
   renannoted_utrs = dplyr::filter(example_GENCODE_exons, id == string)
   extended_utrs = dplyr::filter(example_APAtrap_annotations, id == string)

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
	return(renannoted_utrs)
	}
}


# example exon set

#test condition 1: contains both single and multiple exons per transcript
#test condition 2: contains both forward and reverse strand transcripts
#test condition 3: All permutations of conditions 1 and 2

example_GENCODE_exons = tibble(           
	id=c('A','A','B','B','C','C','D','E','E','F','F'),
        start=c(5,10,20,15,25,50,40,200,300,550,530),
        stop=c(10,15,25,20,30,60,50,250,350,600,545),
	strand=c('1','1','-1','-1','1','1','-1','1','1','-1','-1')
)

transcript_IDs = example_GENCODE_exons$id %>% unique()

terminal_exons = map(transcript_IDs, get_terminal_exons_only, example_GENCODE_exons)
terminal_exons = ldply(terminal_exons, data.frame) %>% as.tibble()

cat(cyan('example_GENCODE_exons\n'))

print(example_GENCODE_exons)

cat(cyan('example_GENCODE_terminal_exons\n'))

print(terminal_exons)

# example APAtrap annotations

# test condition 1: contain both extensions and truncations
# test condition 2: contain both forward and reverse strand exons
# test condition 3: contain unaltered annotations
# test conditions n: all combined perumtations of the previous conditions

example_APAtrap_annotations = tibble(
        id=c('A','B','C','D','E','F'),
        start=c(5,5,25,43,200,536),
        stop=c(250,25,27,50,323,545),
        strand=c('1','-1','1','-1','1','-1')
)

term_exons_reannotated = map(transcript_IDs, extend_UTRs) 
term_exons_reannotated = ldply(term_exons_reannotated, data.frame) %>% as.tibble() # merge into one data frame

cat(cyan('example_APAtrap_annotations\n'))
print(example_APAtrap_annotations)
cat(cyan('example_APAtrap_terminal_exon_extensions\n'))
print(term_exons_reannotated)

shortened_utrs = map(transcript_IDs, shorten_UTRs)
shortened_utrs = ldply(shortened_utrs, data.frame) %>% as.tibble()

cat(cyan('UTR_shortenings'))
print(shortened_utrs)




