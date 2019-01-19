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
#mock3 = strsplit(snakemake@input$quant_mock[3], split='/')[[1]][3]
#mock4 = strsplit(snakemake@input$quant_mock[4], split='/')[[1]][3]
real1 = strsplit(snakemake@input$quant_real[1], split='/')[[1]][3]
real2 = strsplit(snakemake@input$quant_real[2], split='/')[[1]][3]
#real3 = strsplit(snakemake@input$quant_real[3], split='/')[[1]][3]
#real4 = strsplit(snakemake@input$quant_real[4], split='/')[[1]][3]

samples = data.frame(
  run=c(mock1,mock2,real1,real2),
  treatment = factor(rep(c('negative_control',"miRNA"),each=2),
                       ordered=FALSE)
)

rownames(samples) = samples$run

print(samples)
  
files = c(
  paste(snakemake@input$quant_mock[1]),
  paste(snakemake@input$quant_mock[2]),
  paste(snakemake@input$quant_real[1]),
  paste(snakemake@input$quant_real[2])
)
names(files) = samples$run

# filter targets for 8mers and the correct miRNA

canon_targets = filter(canon_targets, Site_type %in% snakemake@params$nontarget_site_types)
canon_targets = filter(canon_targets, miRNA_family_ID == snakemake@wildcards$miRNA)

canon_targets$a_Gene_ID = gsub('\\..*','', canon_targets$a_Gene_ID)

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

print(as.data.frame(resLFC))
  
results = cbind(resLFC@rownames, as.tibble(resLFC@listData))
  
results = filter(results, is.na(log2FoldChange) == FALSE)
results$`resLFC@rownames` = gsub('\\..*','',results$`resLFC@rownames`)


# remove lowly expressed transcripts

exp_data = results

# subset the expression data

non_targets_exp = exp_data$log2FoldChange[!exp_data$`resLFC@rownames` %in% canon_targets$a_Gene_ID]

non_targets_exp = non_targets_exp - median(non_targets_exp)

canon_targets = filter(canon_targets, Site_type %in% snakemake@params$target_site_types)
canon_targets_exp = exp_data$log2FoldChange[exp_data$`resLFC@rownames` %in% canon_targets$a_Gene_ID]

canon_targets_exp = canon_targets_exp - median(non_targets_exp)

filt_results = filter(exp_data, baseMean > snakemake@params$exp_threshold)
filt_targets_exp = filt_results$log2FoldChange[filt_results$`resLFC@rownames` %in% canon_targets$a_Gene_ID]

filt_targets_exp = filt_targets_exp - median(non_targets_exp)

#print(length(exp_data$Name))

# build ggplot df

nontargets = tibble(fc=non_targets_exp)
nontargets$legend = stringr::str_interp("No seed binding (n=${length(non_targets_exp)})")

canon_targets = tibble(fc=canon_targets_exp)
canon_targets$legend = stringr::str_interp("seed targets (n=${length(canon_targets_exp)})")

filt_targets = tibble(fc=filt_targets_exp)
filt_targets$legend = stringr::str_interp("filtered seed targets (n=${length(filt_targets_exp)})")

ggplot_df = rbind(nontargets,canon_targets, filt_targets)

p_value = ks.test(filt_targets_exp,canon_targets_exp, alternative='less')

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
        x=expression('log'[2]*'(mRNA Fold Change)'), 
        subtitle=as.expression(bquote(~ p %~~% .(format (p_value, nsmall=3, digits=3) ) ) )) +
  theme(legend.title=element_blank()) +
  coord_cartesian(xlim = c(-snakemake@params$x_lim,snakemake@params$x_lim))

ggsave(snakemake@output[[1]])
