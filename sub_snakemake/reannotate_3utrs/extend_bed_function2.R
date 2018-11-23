#!/bin/env Rscript

get_full_bed = function (normal_bed_file, extended_bed_file) {

# read in the data
normal_bed = read_tsv(normal_bed_file, col_names=c('chromosome','start','stop','strand','id'), col_types=list('c','i','i','c','c'))
extended_utrs = read_tsv(extended_bed_file, col_names=c('chromosome','start','stop','id','dummy','strand'))

# tidy the data
normal_bed = separate(normal_bed, id, into=c('id','version'))
extended_utrs = separate(extended_utrs, id, into=c('id','dummy2','chrom_dup','strand_dup'))

normal_bed$version = NULL

extended_utrs = extended_utrs[,c('chromosome','start','stop','strand','id')]
extended_utrs$chromosome = stringr::str_replace_all(extended_utrs$chromosome, 'chr','')
normal_bed$chromosome = stringr::str_replace_all(normal_bed$chromosome, 'chr','')

##

reannotate = function(tx_id) {

canonical_tx_record = filter(normal_bed, id==tx_id)
apatrap_tx_record = filter(extended_utrs, id==tx_id)

if (dim(apatrap_tx_record)[1] == 0) { # no record
	return (canonical_tx_record)
} else {

	if (dim(canonical_tx_record)[1] == 1) { # single exon

		new_tx_record = canonical_tx_record
		new_tx_record$start = apatrap_tx_record$start
		new_tx_record$stop = apatrap_tx_record$stop

		return(new_tx_record)

	} else {   # I think it is best to do this as APAtrap does not seem to acknowledge multi-exon 3UTRs 

		return(canonical_tx_record)

	#	if (apatrap_tx_record$strand[1] == '+') {

	#		new_tx_record = canonical_tx_record
        #        	new_tx_record$start[1] = apatrap_tx_record$start
        #        	new_tx_record$stop[dim(canonical_tx_record)[1]] = apatrap_tx_record$stop

	#		return (new_tx_record)

#}	#	else if (apatrap_tx_record$strand[1] == '-') {

	#		new_tx_record = canonical_tx_record
        #        	new_tx_record$stop[1] = apatrap_tx_record$stop[1]
        #        	new_tx_record$start[dim(canonical_tx_record)[1]] = apatrap_tx_record$start

	#		return (new_tx_record)

		}
	}
}

print(normal_bed)
print(extended_utrs)

tx_ids = normal_bed[normal_bed$id %in% extended_utrs$id,]
tx_ids = tx_ids$id %>% unique()

old_records_changed = map(tx_ids, reannotate)

old_records_changed = ldply(old_records_changed, data.frame) %>% as.tibble()

print('old records changed')
old_records_changed %>% select(id) %>% table() %>% print()

new_records = extended_utrs[!extended_utrs$id %in% normal_bed$id,]
#new_records$chromosome = paste('chr',new_records$chromosome,sep='')
new_records$strand = gsub('\\+','1',new_records$strand)
new_records$strand = gsub('-','-1',new_records$strand)

print('new records')
print(new_records)

#new_records %>% select(id) %>% table() %>% print()

old_records_unchanged = normal_bed[!normal_bed$id %in% extended_utrs$id,]

print('old records unchanged')
print(old_records_unchanged)

full_set = rbind(old_records_changed, new_records, old_records_unchanged)

# reformat bed file

full_set = full_set[order(full_set$id, decreasing=FALSE),]
#full_set = full_set[order(full_set$start, decreasing=FALSE),]

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

full_set$chromosome <- factor( as.character(full_set$chromosome), levels=desired_order )

full_set = full_set[order(full_set$chromosome, decreasing=FALSE),] # ordering of chromosomes within the BED file
full_set$chromosome = as.character(full_set$chromosome)

full_set = full_set[order(full_set$start, decreasing=FALSE),]

print('just after tx start position ordering')

tx_IDs = full_set$id %>% unique()

reorder_bed_files = function (tx_ID) { # ordering of exons within a transcript
        transcript_specific_bed_records = subset(full_set, full_set$id == tx_ID)
        #print(transcript_specific_bed_records)
        if (transcript_specific_bed_records$strand[1] == '1') {
           txs = transcript_specific_bed_records[order(transcript_specific_bed_records$start, decreasing=FALSE),]
        }
        else if (transcript_specific_bed_records$strand[1] == '-1' ) {
           txs = transcript_specific_bed_records[order(transcript_specific_bed_records$start, decreasing=TRUE),]
        }
        #print(transcript_specific_bed_records)
        return (txs)
}

full_set_sorted = map(tx_IDs, reorder_bed_files)
print(full_set_sorted[1:30])
full_set_sorted = ldply(full_set_sorted, data.frame) %>% as.tibble()
print(full_set_sorted[1:30,], n=Inf)
return(full_set_sorted)

}
