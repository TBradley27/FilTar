#!/usr/bin/env Rscript

library(plyr)
library(tidyverse)
library(DESeq2)
library(tximport)

# read in the data

cl_targets = read_tsv(
	file = snakemake@input$cl_targets,
	col_names = TRUE,
	col_types = 'ccciiiiiccccci',
	trim_ws = FALSE
	)

print(cl_targets)

canon_targets = read_tsv(
        file = snakemake@input$canon_targets,
        col_names = TRUE,
        col_types = 'ccciiiiiccccci',
        trim_ws = FALSE
        )

mock1 = strsplit(snakemake@input$quant_mock[1], split='/')[[1]][3]
mock2 = strsplit(snakemake@input$quant_mock[2], split='/')[[1]][3]
real1 = strsplit(snakemake@input$quant_real[1], split='/')[[1]][3]
real2 = strsplit(snakemake@input$quant_real[2], split='/')[[1]][3]

samples = data.frame(
  run=c(mock1,mock2,real1,real2),
  treatment = factor(rep(c('miRNA',"negative_control"),each=2),
                       ordered=FALSE)
)
  
rownames(samples) = samples$run
  
files = c(
  paste(snakemake@input$quant_mock[1]),
  paste(snakemake@input$quant_mock[2]),
  paste(snakemake@input$quant_real[1]),
  paste(snakemake@input$quant_real[2])
)
names(files) = samples$run

# filter targets for 8mers and the correct miRNA

cl_targets = filter(cl_targets, Site_type %in% snakemake@params$nontarget_site_types) 
cl_targets = filter(cl_targets, miRNA_family_ID == snakemake@wildcards$miRNA)

canon_targets = filter(canon_targets, Site_type %in% snakemake@params$nontarget_site_types)
canon_targets = filter(canon_targets, miRNA_family_ID == snakemake@wildcards$miRNA)

### DESeq2

txi <- tximport(files, type="kallisto", txOut=TRUE)
  
ddsTxi <- DESeqDataSetFromTximport(txi,
                                     colData = samples,
                                     design = ~ treatment)
  
ddsTxi$treatment <- factor(ddsTxi$treatment, 
                          levels=c("negative_control","miRNA") )
  
dds <- DESeq(ddsTxi)

print(resultsNames(dds))
x = "treatment_miRNA_vs_negative_control"
print(x)

resLFC <- lfcShrink(dds, coef=x, 
                      type="normal")
  
print(resLFC)
  
results = cbind(resLFC@rownames, as.tibble(resLFC@listData))
  
results = filter(results, is.na(log2FoldChange) == FALSE)
results$`resLFC@rownames` = gsub('\\..*','',results$`resLFC@rownames`)

# remove lowly expressed transcripts

exp_data = filter(results, baseMean >= snakemake@params$exp_threshold)

# subset the expression data

non_targets_exp = exp_data$log2FoldChange[!exp_data$`resLFC@rownames` %in% canon_targets$a_Gene_ID]

cl_targets = filter(cl_targets, Site_type %in%  snakemake@params$target_site_types)  
cl_targets_exp = exp_data$log2FoldChange[exp_data$`resLFC@rownames` %in% cl_targets$a_Gene_ID]

canon_targets = filter(canon_targets, Site_type %in% snakemake@params$target_site_types)
canon_targets_exp = exp_data$log2FoldChange[exp_data$`resLFC@rownames` %in% canon_targets$a_Gene_ID]

new_targets_names = cl_targets[!cl_targets$a_Gene_ID %in% canon_targets$a_Gene_ID,]
new_targets_exp = exp_data$log2FoldChange[exp_data$`resLFC@rownames` %in% new_targets_names$a_Gene_ID]

old_targets_names = canon_targets[!canon_targets$a_Gene_ID %in% cl_targets$a_Gene_ID,]
old_targets_exp = exp_data$log2FoldChange[exp_data$`resLFC@rownames` %in% old_targets_names$a_Gene_ID]

#print(length(exp_data$Name))

# build ggplot df

nontargets = tibble(fc=non_targets_exp)
nontargets$legend = stringr::str_interp("No seed binding (n=${length(non_targets_exp)})")

cl_targets = tibble(fc=cl_targets_exp)
cl_targets$legend = stringr::str_interp("${snakemake@wildcards$cell_line} targets (n=${length(cl_targets_exp)})")

canon_targets = tibble(fc=canon_targets_exp)
canon_targets$legend = stringr::str_interp("ensembl targets (n=${length(canon_targets_exp)})")

new_targets = tibble(fc=new_targets_exp)
new_targets$legend = stringr::str_interp("new targets (n=${length(new_targets_exp)})")

old_targets = tibble(fc=old_targets_exp)
old_targets$legend = stringr::str_interp("old targets (n=${length(old_targets_exp)})")

ggplot_df = rbind(nontargets,cl_targets, new_targets, old_targets)

print(ggplot_df)

# begin plotting

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
        x=expression('log'[2]*'(mRNA Fold Change)')
	) +
#        subtitle=as.expression(bquote(~ p %~~% .(format (p_value, nsmall=3, digits=3) ) ) )
  theme(legend.title=element_blank()) +
  xlim(-snakemake@params$x_lim,snakemake@params$x_lim)

ggsave(snakemake@output[[1]])
