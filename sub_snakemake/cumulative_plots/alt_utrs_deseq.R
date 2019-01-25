#!/usr/bin/env Rscript

library(plyr)
library(tidyverse)
library(DESeq2)
library(tximport)

# read in the data

cl_targets = read_tsv(	# read in targets derived from cell-line or tissue-specific annotations
	file = snakemake@input$cl_targets,
	col_names = TRUE,
	col_types = 'ccciiiiiccccci',
	trim_ws = FALSE
	)

print(cl_targets)

canon_targets = read_tsv( # read in targets derived from canonical UTR annotations
        file = snakemake@input$canon_targets,
        col_names = TRUE,
        col_types = 'ccciiiiiccccci',
        trim_ws = FALSE
        )

# retrieve run accessions

mock1 = strsplit(snakemake@input$quant_mock[1], split='/')[[1]][3] 
mock2 = strsplit(snakemake@input$quant_mock[2], split='/')[[1]][3]
#mock3 = strsplit(snakemake@input$quant_mock[3], split='/')[[1]][3]
#mock4 = strsplit(snakemake@input$quant_mock[4], split='/')[[1]][3]
real1 = strsplit(snakemake@input$quant_real[1], split='/')[[1]][3]
real2 = strsplit(snakemake@input$quant_real[2], split='/')[[1]][3]
#real3 = strsplit(snakemake@input$quant_real[3], split='/')[[1]][3]
#real4 = strsplit(snakemake@input$quant_real[4], split='/')[[1]][3]


# create the experimental design matrix
samples = data.frame(
  run=c(mock1,mock2,real1,real2),
  treatment = factor(rep(c("negative_control","miRNA"),each=2),
                       ordered=FALSE)
)

print(samples)
  
rownames(samples) = samples$run

# kallisto or salmon data to be read in
files = c(
  paste(snakemake@input$quant_mock[1]),
  paste(snakemake@input$quant_mock[2]),
#  paste(snakemake@input$quant_mock[3]),
#  paste(snakemake@input$quant_mock[4]),
  paste(snakemake@input$quant_real[1]),
  paste(snakemake@input$quant_real[2])
#  paste(snakemake@input$quant_real[3]),
#  paste(snakemake@input$quant_real[4])
)
names(files) = samples$run

# filter targets for 8mers and the correct miRNA

cl_targets = filter(cl_targets, Site_type %in% snakemake@params$nontarget_site_types) 
cl_targets = filter(cl_targets, miRNA_family_ID == 113)
cl_targets = filter(cl_targets, species_ID == 9606)

canon_targets = filter(canon_targets, Site_type %in% snakemake@params$nontarget_site_types)
canon_targets = filter(canon_targets, miRNA_family_ID == 113)
canon_targets = filter(canon_targets, species_ID == 9606)


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

# remove lowly expressed transcripts

exp_data = filter(results, baseMean >= snakemake@params$exp_threshold)

print(exp_data[1:100,])

# subset the expression data

non_targets_exp = exp_data$log2FoldChange[!exp_data$`resLFC@rownames` %in% cl_targets$a_Gene_ID]

#non_targets_exp = non_targets_exp - median(non_targets_exp)

print(cl_targets)
cl_targets = filter(cl_targets, Site_type %in%  snakemake@params$target_site_types)
print(cl_targets) 
cl_targets_exp = exp_data$log2FoldChange[exp_data$`resLFC@rownames` %in% cl_targets$a_Gene_ID]

#cl_targets_exp = cl_targets_exp - median(non_targets_exp)

canon_targets = filter(canon_targets, Site_type %in% snakemake@params$target_site_types)
canon_targets_exp = exp_data$log2FoldChange[exp_data$`resLFC@rownames` %in% canon_targets$a_Gene_ID]

#canon_targets_exp = canon_targets_exp - median(non_targets_exp)

new_targets_names = cl_targets[!cl_targets$a_Gene_ID %in% canon_targets$a_Gene_ID,]
new_targets_exp = exp_data$log2FoldChange[exp_data$`resLFC@rownames` %in% new_targets_names$a_Gene_ID]

#new_targets_exp = new_targets_exp - median(non_targets_exp)

old_targets_names = canon_targets[!canon_targets$a_Gene_ID %in% cl_targets$a_Gene_ID,]
old_targets_exp = exp_data$log2FoldChange[exp_data$`resLFC@rownames` %in% old_targets_names$a_Gene_ID]

#old_targets_exp = old_targets_exp - median(non_targets_exp)

#print(length(exp_data$Name))

# build ggplot df

nontargets = tibble(fc=non_targets_exp)
nontargets$legend = stringr::str_interp("No seed binding (n=${length(non_targets_exp)})")

print(cl_targets_exp)

cl_targets = tibble(fc=cl_targets_exp)
cl_targets$legend = stringr::str_interp("${snakemake@wildcards$cell_line} targets (n=${length(cl_targets_exp)})")

print(cl_targets)

canon_targets = tibble(fc=canon_targets_exp)
canon_targets$legend = stringr::str_interp("canonical targets (n=${length(canon_targets_exp)})")

new_targets = tibble(fc=new_targets_exp)
new_targets$legend = stringr::str_interp("new targets (n=${length(new_targets_exp)})")

old_targets = tibble(fc=old_targets_exp)
old_targets$legend = stringr::str_interp("old targets (n=${length(old_targets_exp)})")

ggplot_df = rbind(nontargets,canon_targets, new_targets, old_targets)

print(ggplot_df)

p_value = ks.test(new_targets_exp, non_targets_exp, alternative='less')

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
        subtitle=as.expression(bquote(~ p %~~% .(format (p_value, nsmall=3, digits=3) ) ) ) )  +
  theme(legend.title=element_blank()) +
  coord_cartesian(xlim = c(-snakemake@params$x_lim,snakemake@params$x_lim))

ggsave(snakemake@output[[1]])
