#!/usr/bin/env Rscript

library(plyr)
library(tidyverse)
library(DESeq2)
library(tximport)
library(BiocParallel)

register(MulticoreParam(10))

## read in and preprocess target data

cl_targets = read_tsv(	# read in targets derived from cell-line or tissue-specific annotations
	file = snakemake@input$cl_targets,
	col_names = TRUE,
	col_types = 'ccciiiiiccccci',
	trim_ws = FALSE
	)

canon_targets = read_tsv( # read in targets derived from canonical UTR annotations
        file = snakemake@input$canon_targets,
        col_names = TRUE,
        col_types = 'ccciiiiiccccci',
        trim_ws = FALSE
        )

## filter targets for the correct site types and the correct miRNA

miRNA_table = readr::read_tsv(
        file = snakemake@input$miRNA_dict,
        col_names = c('family_code','tax_id','mature_miRNA_name','mature_miRNA_sequence'),
        col_types = 'cccc',
        )


species_three_letters = snakemake@wildcards$species
species_tax_id = snakemake@config$tax_ids[[species_three_letters]]

miRNA_name_with_prefix = paste(species_three_letters,snakemake@wildcards$miRNA,sep='-')
miRNA_family = dplyr::filter(miRNA_table, mature_miRNA_name == miRNA_name_with_prefix)
miRNA_family = miRNA_family$family_code[1]

canon_targets = filter(canon_targets, Site_type %in% snakemake@params$nontarget_site_types)
canon_targets = filter(canon_targets, miRNA_family_ID == miRNA_family)
canon_targets = filter(canon_targets, species_ID == species_tax_id)

cl_targets = filter(cl_targets, Site_type %in% snakemake@params$nontarget_site_types) 
cl_targets = filter(cl_targets, miRNA_family_ID == miRNA_family)
cl_targets = filter(cl_targets, species_ID == species_tax_id)

# retrieve run accessions

mock = strsplit(snakemake@input$quant_mock, split='/')
real = strsplit(snakemake@input$quant_real, split='/')
mock_accessions = c()
real_accessions = c()

for (i in 1:length(mock)) {
  mock_accessions = c(mock_accessions, mock[[i]][3])
  real_accessions = c(real_accessions, real[[i]][3])
}

# create the experimental design matrix
samples = data.frame(
  run=c(mock_accessions, real_accessions),
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

### PCA plot ###

png('pca.png')
vsd = vst(dds, blind=FALSE)
plotPCA(vsd, intgroup=c("treatment"))
dev.off()

### shrink by log2 fold change

x = "treatment_miRNA_vs_negative_control"

resLFC <- lfcShrink(dds, coef=x, 
                      type="normal", parallel=TRUE)
  
results = cbind(resLFC@rownames, as_tibble(resLFC@listData))
  
exp_data = filter(results, is.na(log2FoldChange) == FALSE)

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

non_targets_exp = exp_data$log2FoldChange[!exp_data$`resLFC@rownames` %in% cl_targets$a_Gene_ID]
non_targets_exp = non_targets_exp - median(non_targets_exp)

cl_targets = filter(cl_targets, Site_type %in%  snakemake@params$target_site_types)
cl_targets_exp = exp_data$log2FoldChange[exp_data$`resLFC@rownames` %in% cl_targets$a_Gene_ID]

cl_targets_exp = cl_targets_exp - median(non_targets_exp)

canon_targets = filter(canon_targets, Site_type %in% snakemake@params$target_site_types)
canon_targets_exp = exp_data$log2FoldChange[exp_data$`resLFC@rownames` %in% canon_targets$a_Gene_ID]

canon_targets_exp = canon_targets_exp - median(non_targets_exp)

new_targets_names = cl_targets[!cl_targets$a_Gene_ID %in% canon_targets$a_Gene_ID,]
new_targets_exp = exp_data$log2FoldChange[exp_data$`resLFC@rownames` %in% new_targets_names$a_Gene_ID]

new_targets_exp = new_targets_exp - median(non_targets_exp)

old_targets_names = canon_targets[!canon_targets$a_Gene_ID %in% cl_targets$a_Gene_ID,]
old_targets_exp = exp_data$log2FoldChange[exp_data$`resLFC@rownames` %in% old_targets_names$a_Gene_ID]

old_targets_exp = old_targets_exp - median(non_targets_exp)

### build ggplot df

nontargets = tibble(fc=non_targets_exp)
nontargets$legend = stringr::str_interp("No seed binding (n=${length(non_targets_exp)})")

cl_targets = tibble(fc=cl_targets_exp)
cl_targets$legend = stringr::str_interp("${snakemake@wildcards$cell_line} targets (n=${length(cl_targets_exp)})")

canon_targets = tibble(fc=canon_targets_exp)
canon_targets$legend = stringr::str_interp("canonical targets (n=${length(canon_targets_exp)})")

new_targets = tibble(fc=new_targets_exp)
new_targets$legend = stringr::str_interp("new targets (n=${length(new_targets_exp)})")

old_targets = tibble(fc=old_targets_exp)
old_targets$legend = stringr::str_interp("old targets (n=${length(old_targets_exp)})")

ggplot_df = rbind(nontargets,canon_targets, new_targets, old_targets)

p_value = ks.test(new_targets_exp, non_targets_exp, alternative='less')

### begin plotting

ggplot(
  ggplot_df, aes(x=fc,color=legend)
    ) +
  geom_step(aes(y=..y..), stat="ecdf") +
  theme_classic() +
  labs(
        title=
        bquote(
                .(str_interp("${snakemake@wildcards$miRNA}")) ~ .(substitute(italic(vs.))) ~ 'mock transfection' ~ .(str_interp("(${snakemake@wildcards$cell_line})"))
        ),
        y="Cumulative Fraction",
        x=expression('log'[2]*'(mRNA Fold Change)'), 
        subtitle=as.expression(bquote(~ p %~~% .(format (p_value, nsmall=3, digits=3) ) ) ) )  +
  theme(legend.title=element_blank()) +
  coord_cartesian(xlim = c(-snakemake@params$x_lim,snakemake@params$x_lim))

ggsave(snakemake@output[[1]])
