library(gridExtra)
library(grid)
library(ggplot2)

A549 = readRDS(file = snakemake@input[['A549']])
HeLa = readRDS(file = snakemake@input[['HeLa']])
NMuMG = readRDS(file= snakemake@input[['NMuMG']])
ESCs = readRDS(file= snakemake@input[['ESCs']])

g = arrangeGrob(A549, HeLa, NMuMG, ESCs,
             ncol=2, nrow=2,
             bottom=textGrob(
               expression(bold('log'[2]*'(mRNA Fold Change)'))
               ), 
             left=textGrob(
               'Cumulative Proportion', 
               rot=90,
               gp=gpar(fontface="bold")
               )
             )
ggsave(snakemake@output[[1]],plot = g, scale=1.25)
