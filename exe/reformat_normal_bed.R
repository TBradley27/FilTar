#!/bin/env Rscript

library(readr)
library(stringr)
library(tidyr)
library(plyr)
library(dplyr)
library(tibble)
library(purrr)

args = commandArgs(trailingOnly=TRUE)

bed_file = read_tsv(args[1], col_names=c('chromosome','start','stop','strand','id'), col_types=list('c','i','i','c','c'))

bed_file = separate(bed_file, id, into=c('id','version'))

bed_file$version = NULL

bed_file = bed_file[order(bed_file$id, decreasing=FALSE),]

bed_file = bed_file[order(bed_file$start, decreasing=FALSE),]

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

bed_file$chromosome <- factor( as.character(bed_file$chromosome), levels=desired_order )

bed_file = bed_file[order(bed_file$chromosome, decreasing=FALSE),]

bed_file$chromosome = as.character(bed_file$chromosome)

tx_IDs = bed_file$id %>% unique()

reorder_bed_files = function (tx_ID) {
        transcript_specific_bed_records = subset(bed_file, bed_file$id == tx_ID)
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

bed_file_sorted = map(tx_IDs, reorder_bed_files)
bed_file_sorted = ldply(bed_file_sorted, data.frame) %>% as.tibble()

#print(bed_file_sorted)

write.table(bed_file_sorted, file=args[2], quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t")
