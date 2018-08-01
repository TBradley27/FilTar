#!/usr/bin/env Rscript

library(tibble)
library(plyr)
library(readr)
library(dplyr, warn.conflicts = FALSE)
library(tidyr)
library(stringr)

MergeFasta = function(file) {
  
  # A function which concatenates sequences (depending on order of appearance in fasta file) which are 
  # attributable to the same fasta header
  
  # Input:
  #   file: An unmered fasta file
  
  # Output:
  #   fasta: A merged fasta file
  
  fasta = read_tsv(file, col_names=c("fasta_unicolumn"))
  
  header_records = grepl(">", fasta$fasta_unicolumn) # Boolean vector
  
  fasta_header = fasta[header_records,]
  fasta_seq = fasta[!header_records,]
  
  fasta = dplyr::bind_cols(fasta_header, fasta_seq) # Convert to a n x 2 data frame
  colnames(fasta) = c("header","seq")
  
  fasta = plyr::ddply(fasta, .(header), function(dat) { # The main functional code
    paste(dat$seq, collapse='')
  }  )
  
  colnames(fasta) = c("header","seq")

  fasta = as.tibble(fasta)

  fasta = c(t(fasta))   # Transpose data frame, so that data can be written in the correct format
  
  return(fasta)
  
}

args = commandArgs(trailingOnly=TRUE)
fasta = MergeFasta(args[1])
write.table(fasta, args[2], quote=FALSE, col.names=FALSE, row.names = FALSE)
