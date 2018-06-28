#!/bin/env Rscript

library(readr)
library(stringr)
library(tidyr)
library(plyr)
library(dplyr)
library(tibble)
library(purrr)

args = commandArgs(trailingOnly=TRUE)

print(args)

print ('foo')
normal_bed = read_tsv(args[1], col_names=c('chromosome','start','stop','strand','id'), col_types=list('c','i','i','c','c'))
print ('bar')

print ('foo')
extended_utrs = read_tsv(args[2], col_names=c('chromosome','start','stop','id','dummy','strand'))
print ('bar')

normal_bed = separate(normal_bed, id, into=c('id','version'))
extended_utrs = separate(extended_utrs, id, into=c('id','dummy2','chrom_dup','strand_dup'))

#extended_utrs$dummy = NULL
#extended_utrs$dummy2 = NULL
#extended_utrs$chrom_dup = NULL
#extended_utrs$strand_dup = NULL
normal_bed$version = NULL

extended_utrs = extended_utrs[,c('chromosome','start','stop','strand','id')]
extended_utrs$chromosome = stringr::str_replace_all(extended_utrs$chromosome, 'chr','')

not_extended_utrs = normal_bed[!(normal_bed$id %in% extended_utrs$id),]

extended_utrs_normal = normal_bed[(normal_bed$id %in% extended_utrs$id),]

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
terminal_exons = ldply(terminal_exons, data.frame) %>% as.tibble()

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
term_exons_reannotated = ldply(term_exons_reannotated, data.frame) %>% as.tibble()

nonterminal_exons = dplyr::anti_join(normal_bed, terminal_exons)

reannotated_utrs = rbind(nonterminal_exons, term_exons_reannotated)

write.table(reannotated_utrs, file=args[3], quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t")
