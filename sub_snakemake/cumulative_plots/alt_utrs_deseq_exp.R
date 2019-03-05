#!/usr/bin/env Rscript

library(plyr)
library(tidyverse)
library(DESeq2)
library(tximport)
library(BiocParallel)

register(MulticoreParam(20))

# read in the data

canon_targets = read_tsv(
        file = snakemake@input$canon_targets,
        col_names = TRUE,
        col_types = 'ccciiiiiccccci',
        trim_ws = FALSE
        )

# filter targets for 8mers and the correct miRNA

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
print(canon_targets)
canon_targets = filter(canon_targets, miRNA_family_ID == miRNA_family)
print(canon_targets)
canon_targets = filter(canon_targets, species_ID == species_tax_id)

## kallisto data 

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


samples = data.frame(
  run=c(mock_accessions, real_accessions),
  #treatment = factor(c('negative_control','negative_control','negative_control','negative_control','miRNA','miRNA')
  treatment = rep(c("negative_control","miRNA"),each=length(mock),
                       ordered=FALSE)
)

rownames(samples) = samples$run
files = c(snakemake@input$quant_mock, snakemake@input$quant_real)  
names(files) = samples$run

#canon_targets$a_Gene_ID = gsub('\\..*','', canon_targets$a_Gene_ID)

### DESeq2

txi <- tximport(files, type="kallisto", txOut=TRUE)
  
ddsTxi <- DESeqDataSetFromTximport(txi,
                                     colData = samples,
                                     design = ~ treatment)
  
ddsTxi$treatment <- factor(ddsTxi$treatment, 
                          levels=c("negative_control","miRNA") )
  
dds <- DESeq(ddsTxi, parallel=TRUE)

print(resultsNames(dds))
x = "treatment_miRNA_vs_negative_control"

resLFC <- lfcShrink(dds, coef=x, 
                      type="normal", parallel=TRUE)

results = cbind(resLFC@rownames, as_tibble(resLFC@listData))
 
#results$log2FoldChange[is.na(results$log2FoldChange)] <- 0
 
results = filter(results, is.na(log2FoldChange) == FALSE)
#results$`resLFC@rownames` = gsub('\\..*','',results$`resLFC@rownames`)


# remove lowly expressed transcripts

exp_data = results

# subset the expression data

non_targets_exp = exp_data$log2FoldChange[!exp_data$`resLFC@rownames` %in% canon_targets$a_Gene_ID]
non_targets_exp = non_targets_exp - median(non_targets_exp)

canon_targets = filter(canon_targets, Site_type %in% snakemake@params$target_site_types)
canon_targets_exp = exp_data$log2FoldChange[exp_data$`resLFC@rownames` %in% canon_targets$a_Gene_ID]
canon_targets_exp = canon_targets_exp - median(non_targets_exp)

print('median baseMean')
median_baseMean = median(exp_data$baseMean)
print(median_baseMean)

print(exp_data[1:100,])

filt_results = filter(exp_data, baseMean > median_baseMean)

filt_targets_exp = filt_results$log2FoldChange[filt_results$`resLFC@rownames` %in% canon_targets$a_Gene_ID]
filt_targets_exp = filt_targets_exp - median(non_targets_exp)

filt_non_targets_exp = filt_results$log2FoldChange[!filt_results$`resLFC@rownames` %in% canon_targets$a_Gene_ID]
filt_non_targets_exp = filt_non_targets_exp - median(non_targets_exp)

#print(length(exp_data$Name))

# build ggplot df

nontargets = tibble(fc=non_targets_exp)
nontargets$legend = stringr::str_interp("No seed binding (n=${length(non_targets_exp)})")

canon_targets = tibble(fc=canon_targets_exp)
canon_targets$legend = stringr::str_interp("seed targets (n=${length(canon_targets_exp)})")

filt_targets = tibble(fc=filt_targets_exp)
filt_targets$legend = stringr::str_interp("seed targets (filtered) (n=${length(filt_targets_exp)})")

filt_non_targets = tibble(fc=filt_non_targets_exp)
filt_non_targets$legend = stringr::str_interp("No seed binding (filtered) (n=${length(filt_non_targets_exp)})")

ggplot_df = rbind(nontargets,canon_targets, filt_targets)

p_value = ks.test(filt_targets_exp,canon_targets_exp, alternative='greater')

print(ggplot_df)

# begin plotting

ggplot_object = ggplot(
  ggplot_df, aes(x=fc,color=legend)
    ) +
  geom_step(aes(y=..y..), stat="ecdf") +
  theme_classic() +
  labs(
        title=
        bquote(
                .(str_interp("${snakemake@wildcards$miRNA}")) ~ 'transfection' ~ .(str_interp("(${snakemake@wildcards$cell_line})"))
        ),
        y=NULL,
	tag=expression(bold("A")),
        x=NULL, 
        subtitle=as.expression(bquote(~ p %~~% .(format (p_value$p.value, nsmall=3, digits=3) ) ) )
	) +
  theme(legend.title=element_blank(), legend.position=c(0.8,0.2)) +
  scale_color_manual(values=c("black", "sienna2", "red")) +
  coord_cartesian(xlim = c(-snakemake@params$x_lim,snakemake@params$x_lim))

## save ggplot object

saveRDS(ggplot_object, file = snakemake@output[[1]])

#ggsave(snakemake@output[[1]])
