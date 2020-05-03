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


gg = ggplot(results, aes(x=log2FoldChange, y=-log10(padj))) + geom_point() +
 geom_point(data = filter(
	results,
	log2FoldChange > 2.0 & padj < 0.05 | log2FoldChange < -2.0 & padj < 0.05
	), colour='red') +
 geom_hline(yintercept=-log10(0.05)) + geom_vline(xintercept=c(-2,2)) +
  theme_classic() +
  labs(
        title=
        bquote(
               bold(.(str_interp("${snakemake@wildcards$miRNA}")) ~ 'transfection' ~ .(str_interp("(${snakemake@wildcards$cell_line})"))
        ))
	)
#plot(results$log2FoldChange, -log10(results$padj),
#     main="Volcano plot", xlab="Effect size: log2(fold-change)", ylab="-log10(adjusted p-value)",
#     pch=20, cex=0.6)

#abline(v=c(-2,2), col="brown")
#abline(h=-log10(0.05), col="brown")

#with(subset(results, padj<0.05 & abs(log2FoldChange)>2), points(log2FoldChange, -log10(padj), pch=20, col="red", cex=0.5))

saveRDS(object=gg, file=snakemake@output[[1]])
