library(gridExtra)
library(grid)
library(ggplot2)

## open device

a=readRDS(snakemake@input[['HeLa_21']])
b=readRDS(snakemake@input[['HeLa_22']])
c=readRDS(snakemake@input[['HeLa_23']])
d=readRDS(snakemake@input[['HeLa_24']])
e=readRDS(snakemake@input[['HeLa_1']])
f=readRDS(snakemake@input[['HeLa_2']])
g=readRDS(snakemake@input[['HeLa_3']])
h=readRDS(snakemake@input[['HeLa_4']])
i=readRDS(snakemake@input[['HeLa_17']])
j=readRDS(snakemake@input[['HeLa_18']])
k=readRDS(snakemake@input[['HeLa_19']])
l=readRDS(snakemake@input[['HeLa_20']])

g = arrangeGrob(a,b,c,d,e,f,g,h,i,j,k,l,
             ncol=3, nrow=4,
             bottom=textGrob(
               expression(bold(''))
               ),
             left=textGrob(
               '',
               rot=90,
               gp=gpar(fontface="bold")
               )
             )
ggsave(snakemake@output[[1]],plot = g, scale=2.00, dev='png', dpi=400)

dev.off()
