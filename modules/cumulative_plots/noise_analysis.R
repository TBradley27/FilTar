#!/usr/bin/env Rscript

library(plyr)
library(tidyverse)
library(DESeq2)
library(tximport)
library(BiocParallel)

register(MulticoreParam(snakemake@threads[[1]]))

## read in transcript data

pc_transcripts = read_tsv(
	file = snakemake@input$pc_transcripts,
	col_names= c('tx_id'),
	col_types= 'c'
	)

## read in and preprocess target data

canon_targets = read_tsv( # read in targets derived from canonical UTR annotations
        file = snakemake@input$canon_targets,
        col_names = TRUE,
        col_types = 'ccciiiiiccccci',
        trim_ws = FALSE
        )

## filter targets for the correct site types and the correct miRNA


species_three_letters = snakemake@wildcards$species
species_tax_id = snakemake@config$tax_ids[[species_three_letters]]

#miRNA_name_with_prefix = paste(species_three_letters,snakemake@wildcards$miRNA,sep='-')
#miRNA_family = dplyr::filter(miRNA_table, mature_miRNA_name == miRNA_name_with_prefix)
#miRNA_family = miRNA_family$family_code[1]

canon_targets = filter(canon_targets, Site_type %in% snakemake@params$nontarget_site_types)
#canon_targets = filter(canon_targets, miRNA_family_ID == miRNA_family)
canon_targets = filter(canon_targets, species_ID == species_tax_id)

# retrieve run accessions

mock = strsplit(snakemake@input$quant_mock, split='/')
real = strsplit(snakemake@input$quant_real, split='/')
mock_accessions = c()
real_accessions = c()

for (i in 1:length(mock)) {
  mock_accessions = c(mock_accessions, mock[[i]][3])
}

for (i in 1:length(real)) {
  real_accessions = c(real_accessions, real[[i]][3])
}

# create the experimental design matrix
samples = data.frame(
  run=c(mock_accessions, real_accessions),
#  treatment = factor(c('negative_control','negative_control','negative_control','miRNA','miRNA'),
  treatment = factor(rep(c("negative_control","miRNA"),each=length(mock)),
                       ordered=FALSE)
)
  
rownames(samples) = samples$run

# kallisto or salmon data to be read in
files = c(snakemake@input$quant_mock, snakemake@input$quant_real)
names(files) = samples$run

### DESeq2

txi <- tximport(files, type="kallisto", txOut=TRUE)
  
ddsTxi <- DESeqDataSetFromTximport(txi, # big decision whether to use txi or txi$counts here - see the txi vignette
                                     colData = samples,
                                     design = ~ treatment)

#ddsTxi <- DESeqDataSetFromMatrix(round(txi[['counts']]),
#                                     colData = samples,
#                                     design = ~ treatment)
  
ddsTxi$treatment <- factor(ddsTxi$treatment, 
                          levels=c("negative_control","miRNA") )

ddsTxi = estimateSizeFactors(ddsTxi)

### filter by expression level

keep <- rowMeans(counts(ddsTxi, normalized=TRUE)[,1:length(mock)]) >= 0 # filter by normalised count levels
ddsTxi <- ddsTxi[keep,]
dds <- DESeq(ddsTxi, parallel=TRUE)

### shrink by log2 fold change

x = "treatment_miRNA_vs_negative_control"

resLFC <- lfcShrink(dds, coef=x, 
                      type="normal", parallel=TRUE)
  
results = cbind(resLFC@rownames, as_tibble(resLFC@listData))
  
exp_data = dplyr::filter(results, is.na(log2FoldChange) == FALSE)

exp_data = exp_data[exp_data$`resLFC@rownames` %in% pc_transcripts$tx_id,]

### expression filter by TPM

#TPMs = as.data.frame(txi$abundance)
#TPMs$average = rowMeans(TPMs)
#TPMs$names = rownames(txi$abundance)
#print(TPMs[1:10,])
#TPMs = dplyr::filter(TPMs, average >= 0)
#exp_data = filter(results, baseMean >= snakemake@params$exp_threshold)

#print(dim(results))
#exp_data = results[results$`resLFC@rownames` %in% TPMs$names,]
#print(dim(exp_data))

#print(exp_data[1:100,])

### Split the lfc data into distinct sets

non_targets_exp = exp_data$log2FoldChange[!exp_data$`resLFC@rownames` %in% canon_targets$a_Gene_ID]
non_targets_exp2 = non_targets_exp - median(non_targets_exp)

canon_targets = filter(canon_targets, Site_type %in% snakemake@params$target_site_types)
canon_targets_exp = exp_data$log2FoldChange[exp_data$`resLFC@rownames` %in% canon_targets$a_Gene_ID]

canon_targets_exp = canon_targets_exp - median(non_targets_exp)

### build ggplot df

nontargets = tibble(fc=non_targets_exp2)
nontargets$legend = stringr::str_interp("No seed site (n=${length(non_targets_exp)})")

canon_targets = tibble(fc=canon_targets_exp)
canon_targets$legend = stringr::str_interp("Seed site (n=${length(canon_targets_exp)})")

all_zero_exp = rep(0, 10000)
all_zero = tibble(fc=all_zero_exp)
all_zero$legend = stringr::str_interp("no noise (n=${length(all_zero_exp)})")

signal_targets = sqrt( sum ( (canon_targets_exp**2) ) / length(canon_targets_exp) )
signal_nontargets = sqrt( sum ( (non_targets_exp2**2) ) / length(non_targets_exp2) )
signal_noise_ratio = signal_targets / signal_nontargets 

print(signal_targets)
print(signal_nontargets)
print(signal_noise_ratio)

results = tibble(signal_targets=signal_targets, signal_nontargets=signal_nontargets, signal_noise_ratio=signal_noise_ratio)

write.table(
x=results,
file = snakemake@output[[1]],
col.names=TRUE, row.names=TRUE, sep='\t', quote=FALSE) 


#old_targets$legend = stringr::str_interp("old targets (n=${length(old_targets_exp)})")

#ggplot_df = rbind(nontargets,canon_targets, all_zero)
#
#ggplot_df$legend = factor(ggplot_df$legend, levels = c(
##	stringr::str_interp("Added seed site (n=${length(new_targets_exp)})"),
#        stringr::str_interp("Seed site (n=${length(canon_targets_exp)})"),
#	stringr::str_interp("No seed site (n=${length(non_targets_exp)})")					
#	)	
#)
#
##wilcox.test(new_targets_exp, non_targets_exp, alternative='less')
#
##ks.test(new_targets_exp, non_targets_exp, alternative='less')
#
#### begin plotting
#
#ggplot_object = ggplot(
#  ggplot_df, aes(x=fc,color=legend)
#    ) +
#  geom_step(aes(y=..y..), stat="ecdf") +
#  theme_classic() +
#  labs(
#        title=
#        bquote(
#                .(str_interp("${snakemake@wildcards$miRNA}")) ~ 'transfection' ~ .(str_interp("(${snakemake@wildcards$cell_line})"))
#        ),
#        y="Cumulative Proportion",
#	tag=expression(bold("")),
#        x=expression('log'[2]*'(mRNA Fold Change)')
#        subtitle=as.expression(bquote(~ p %~~% .(format (p_value$p.value, nsmall=3, digits=3) ) ) )
#	)  +
#  theme(legend.title=element_blank(), legend.position=c(0.75,0.20)) +
#  scale_color_manual(
#		values=c("dodgerblue","black","darkorange"),
#                breaks=c(    # change legend order
#                        stringr::str_interp("Seed site (n=${length(canon_targets_exp)})"),
#                        stringr::str_interp("No seed site (n=${length(non_targets_exp)})"),
#                        )    
#	            ) +
#  coord_cartesian(xlim = c(-snakemake@params$x_lim,snakemake@params$x_lim))
#
#
#### save
#
#saveRDS(ggplot_object, file = 'foo.rds')
#
##ggsave(snakemake@output[[1]])
