library(gridExtra)
library(grid)
library(ggplot2)

HeLa_9 = readRDS(file = snakemake@input[['HeLa_9']])
HeLa_11 = readRDS(file = snakemake@input[['HeLa_11']])
HeLa_1 = readRDS(file =snakemake@input[['HeLa_1']])
U343 = readRDS(file=snakemake@input[['U343']])

g = arrangeGrob(U343, HeLa_9, HeLa_11, HeLa_1,
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
