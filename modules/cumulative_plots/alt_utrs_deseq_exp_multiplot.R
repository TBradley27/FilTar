#!/usr/bin/env Rscript

library(plyr)
library(tidyverse)
library(DESeq2)
library(tximport)
library(BiocParallel)

register(MulticoreParam(snakemake@threads[[1]]))

# read in protein-coding transcripts

pc_transcripts = read_tsv(
        file = snakemake@input$pc_transcripts,
        col_names= c('tx_id'),
        col_types= 'c'
        )

# read in the data

canon_targets = read_tsv(
        file = snakemake@input$canon_targets,
        col_names = TRUE,
        col_types = 'ccciiiiiccccci',
        trim_ws = FALSE
        )

# filter targets for 8mers and the correct miRNA

#miRNA_table = readr::read_tsv(
#        file = snakemake@input$miRNA_dict,
#        col_names = c('family_code','tax_id','mature_miRNA_name','mature_miRNA_sequence'),
#        col_types = 'cccc',
#        )


species_three_letters = snakemake@wildcards$species
species_tax_id = snakemake@config$tax_ids[[species_three_letters]]

#miRNA_name_with_prefix = paste(species_three_letters,snakemake@wildcards$miRNA,sep='-')
#miRNA_family = dplyr::filter(miRNA_table, mature_miRNA_name == miRNA_name_with_prefix)
#miRNA_family = miRNA_family$family_code[1]

canon_targets = filter(canon_targets, Site_type %in% snakemake@params$nontarget_site_types)
print(canon_targets)
#canon_targets = filter(canon_targets, miRNA_family_ID == miRNA_family)
#print(canon_targets)
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
#  treatment = factor(c('negative_control','negative_control','negative_control','miRNA','miRNA')
  treatment = rep(c("negative_control","miRNA"),each=length(mock),
                       ordered=FALSE
	)
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
exp_data = exp_data[exp_data$`resLFC@rownames` %in% pc_transcripts$tx_id,]

# subset the expression data

non_targets_exp = exp_data$log2FoldChange[!exp_data$`resLFC@rownames` %in% canon_targets$a_Gene_ID]
non_targets_exp2 = non_targets_exp - median(non_targets_exp)

canon_targets = filter(canon_targets, Site_type %in% snakemake@params$target_site_types)
canon_targets_exp = exp_data$log2FoldChange[exp_data$`resLFC@rownames` %in% canon_targets$a_Gene_ID]
canon_targets_exp = canon_targets_exp - median(non_targets_exp)

#print('median baseMean')
#median_baseMean = median(exp_data$baseMean)
#print(median_baseMean)

#print(exp_data[1:100,])

TPMs = as.data.frame(txi$abundance)
TPMs$average = rowMeans(TPMs[,1:length(mock)])
TPMs$names = rownames(txi$abundance)

##########

for (item in 1:length(snakemake@params$exp_threshold)){
	TPMs = dplyr::filter(TPMs, average >= snakemake@params$exp_threshold[item])
	filt_results = exp_data[exp_data$`resLFC@rownames` %in% TPMs$names,]
	filt_targets_exp = filt_results$log2FoldChange[filt_results$`resLFC@rownames` %in% canon_targets$a_Gene_ID]
	filt_targets_exp = filt_targets_exp - median(non_targets_exp)
	filt_targets_tmp = paste("filt_targets",item,sep="_")
	assign(filt_targets_tmp, tibble(fc=filt_targets_exp, legend=stringr::str_interp("Seed site (${snakemake@params$exp_threshold[item]} TPM) (n=${length(filt_targets_exp)})")))
	print(filt_targets_tmp)
}

print(filt_targets_1)
print(filt_targets_2)
print(filt_targets_3)
print(filt_targets_4)
print(filt_targets_5)

##########

filt_non_targets_exp = filt_results$log2FoldChange[!filt_results$`resLFC@rownames` %in% canon_targets$a_Gene_ID]
filt_non_targets_exp = filt_non_targets_exp - median(non_targets_exp)

#x = filt_targets_exp - filt_nontargets_exp
#y = canon_targets_exp - non_targets_exp

#print(length(exp_data$Name))

# build ggplot df

nontargets = tibble(fc=non_targets_exp2)
nontargets$legend = stringr::str_interp("No seed site (n=${length(non_targets_exp)})")

canon_targets = tibble(fc=canon_targets_exp)
canon_targets$legend = stringr::str_interp("Seed site (0.0 TPM) (n=${length(canon_targets_exp)})")

#filt_targets = tibble(fc=filt_targets_exp)
#filt_targets$legend = stringr::str_interp("Seed site (filtered) (n=${length(filt_targets_exp)})")

filt_non_targets = tibble(fc=filt_non_targets_exp)
filt_non_targets$legend = stringr::str_interp("No seed binding (filtered) (n=${length(filt_non_targets_exp)})")

#filtered = tibble(fc=x)
#filtered$legend = stringr::str_interp("No seed binding (n=${length(x)})")

#not_filtered = tibble(fc=y)
#not_filtered$legend = stringr::str_interp("No seed binding (n=${length(y)})")

ggplot_df = rbind(nontargets,canon_targets, filt_targets_1,filt_targets_2,filt_targets_3,filt_targets_4,filt_targets_5)

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
        y="Cumulative Proportion",
	tag=expression(bold("A")),
        x=expression('log'[2]*'(mRNA Fold Change)') 
        #subtitle=as.expression(bquote(~ p %~~% .(format (p_value$p.value, nsmall=3, digits=3) ) ) )
	) +
  theme(legend.title=element_blank(), legend.position=c(0.77,0.27)) +
  scale_color_manual(
	values=c("black","darkorange","limegreen","lightgoldenrod","darkorchid1","aquamarine","hotpink1"),
        breaks=c(    # change legend order
			stringr::str_interp("Seed site (10 TPM) (n=2224)"),
			stringr::str_interp("Seed site (5 TPM) (n=3730)"),
			stringr::str_interp("Seed site (1 TPM) (n=7253)"),
			stringr::str_interp("Seed site (0.5 TPM) (n=8579)"),
                        stringr::str_interp("Seed site (0.1 TPM) (n=11192)"),
                        stringr::str_interp("Seed site (0.0 TPM) (n=${length(canon_targets_exp)})"),
                        stringr::str_interp("No seed site (n=${length(non_targets_exp)})")
                        )
	) +
  coord_cartesian(xlim = c(-snakemake@params$x_lim,snakemake@params$x_lim))

## save ggplot object

saveRDS(ggplot_object, file = snakemake@output[[1]])

#ggsave(snakemake@output[[1]])
