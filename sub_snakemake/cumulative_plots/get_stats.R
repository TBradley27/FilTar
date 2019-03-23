library(tidyverse)
library(tximport)

###

real = strsplit(snakemake@input$quant, split='/')
real_accessions = c()

for (i in 1:length(real)) {
  real_accessions = c(real_accessions, real[[i]][3])
}

files = c(snakemake@input$quant)
names(files) = real_accessions

txi <- tximport(files, type="kallisto", txOut=TRUE)

print('foo')

TPMs = as.data.frame(txi$abundance)
TPMs$average = rowMeans(TPMs)
TPMs$names = rownames(txi$abundance)

high_expression = dplyr::filter(TPMs, average >= 5.0)
medium_expression = dplyr::filter(TPMs, average >= 0.1)

print(high_expression)

#high_expression = dplyr::filter(kallisto, tpm >= 5.0)

###

cerebellum_targets = read_tsv(snakemake@input[['tissue']])
cerebellum_targets_names = cerebellum_targets$a_Gene_ID %>% unique()

pc_transcripts = read_tsv(
  snakemake@input[['pc_transcripts']],
  col_names=FALSE
  )

cerebellum_targets_names = 
  cerebellum_targets_names[cerebellum_targets_names %in% pc_transcripts$X1]
cerebellum_targets = cerebellum_targets[cerebellum_targets$a_Gene_ID %in% pc_transcripts$X1,]


canonical_targets = read_tsv(snakemake@input[['control']])
canonical_targets_names = canonical_targets$a_Gene_ID %>% unique()

canonical_targets_names = canonical_targets_names[canonical_targets_names %in% pc_transcripts$X1]
canonical_targets = canonical_targets[canonical_targets$a_Gene_ID %in% pc_transcripts$X1,]

# gained target transcripts

a=cerebellum_targets_names[!cerebellum_targets_names %in% canonical_targets_names] %>% length() #%>% print()

# lost target transcripts

old_targets=canonical_targets_names[!canonical_targets_names %in% cerebellum_targets_names]

b=old_targets[old_targets %in% high_expression$names] %>% length() #%>% print()

# expression filtered targets

filtered_targets = canonical_targets[!canonical_targets$a_Gene_ID %in% medium_expression$names,]
e = dim(filtered_targets)[1]

###

cerebellum_targets$Group_num = NULL
canonical_targets$Group_num = NULL

# gained target sites

c=setdiff(cerebellum_targets,canonical_targets) #%>% print()
c=dim(c)[1]

# lost target sites

old_target_sites=setdiff(canonical_targets,cerebellum_targets)

old_target_sites = old_target_sites[old_target_sites$a_Gene_ID %in% high_expression$names,]

d=dim(old_target_sites)[1] #%>% print()

output = tibble(
	new_transcripts=a,
        old_transcripts=b,
        new_target_sites=c,
        old_target_sites=d,
        filtered_target_sites=e
       )

print(output)

write.table(
	output,
	snakemake@output[[1]],
	quote=FALSE,
	sep="\t",
	row.names=FALSE
)


