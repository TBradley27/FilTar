library(magrittr)

getwd()
pattern = stringr::str_interp("${snakemake@wildcards[['species']]}_[A-Za-z]+_3UTR.chr${snakemake@wildcards[['chrom']]}.bed$")
print(pattern)
file_list = list.files(path='results/bed/', pattern=pattern, full.names=TRUE)

files=list()

for (i in 1:length(file_list)) {
	files[[i]] = readr::read_tsv(file_list[i], col_names=c('chromosome','start','stop','strand','transcript_ID'), col_types='ciiic')
}

reference = files[[1]]

# identify 3UTRs with multiple exons
tmp = table(reference$transcript_ID) %>% as.data.frame()
colnames(tmp) = c('transcript_ID','frequency')
multi_exons = subset(tmp, frequency > 1 )

reference = reference[!reference$transcript_ID %in% multi_exons$transcript_ID,]

single_pos_exons = subset(reference, reference$strand == 1)
single_neg_exons = subset(reference, reference$strand == -1)

for (i in 1:length(file_list)) {
	tmp_pos = files[[i]][files[[i]]$transcript_ID %in% single_pos_exons$transcript_ID,]
	tmp_pos$stop = as.numeric(tmp_pos$stop)
	single_pos_exons = cbind(single_pos_exons, tmp_pos$stop)
	tmp_neg = files[[i]][files[[i]]$transcript_ID %in% single_neg_exons$transcript_ID,]
	tmp_neg$start = as.numeric(tmp_neg$start)
        single_neg_exons = cbind(single_neg_exons, tmp_neg$start)
}

single_pos_exons$chromosome = NULL; single_neg_exons$chromosome = NULL;
single_pos_exons$start = NULL; single_neg_exons$start = NULL;
single_pos_exons$stop = NULL; single_neg_exons$stop = NULL;
single_pos_exons$strand = NULL; single_neg_exons$strand = NULL;

single_pos_exons$max = apply(single_pos_exons[,-1],1,max) # -1 in this context means ignore the first column i.e. the transcript_id column
single_neg_exons$min = apply(single_neg_exons[,-1],1,min)

# swap in new bed values into the main reference table
for (i in 1:length(single_pos_exons$transcript_ID)) {
	df1_index = grep(single_pos_exons$transcript_ID[i],reference$transcript_ID)
	df2_index = grep(single_pos_exons$transcript_ID[i],single_pos_exons$transcript_ID)
	reference$stop[df1_index] = single_pos_exons$max[df2_index]
}

for (i in 1:length(single_neg_exons$transcript_ID)) {
        df1_index = grep(single_neg_exons$transcript_ID[i],reference$transcript_ID)
        df2_index = grep(single_neg_exons$transcript_ID[i],single_neg_exons$transcript_ID)
        reference$stop[df1_index] = single_neg_exons$min[df2_index]
}


write.table(
	x=reference,
	file=snakemake@output[[1]],
	quote=FALSE,
	row.names=FALSE,
	col.names=FALSE
)


