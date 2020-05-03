library(gridExtra)
library(grid)
library(ggplot2)

## open device

a=readRDS(snakemake@input[['A549']])
b=readRDS(snakemake@input[['HeLa_6']])
c=readRDS(snakemake@input[['HeLa_7']])
d=readRDS(snakemake@input[['HeLa_8']])
e=readRDS(snakemake@input[['HeLa_9']])
f=readRDS(snakemake@input[['HeLa_10']])
g=readRDS(snakemake@input[['HeLa_11']])
h=readRDS(snakemake@input[['HeLa_12']])
i=readRDS(snakemake@input[['HeLa_13']])
j=readRDS(snakemake@input[['HeLa_14']])
k=readRDS(snakemake@input[['HeLa_15']])
l=readRDS(snakemake@input[['HeLa_16']])
#m=readRDS(snakemake@input[['HeLa_17']])
#n=readRDS(snakemake@input[['HeLa_18']])
#o=readRDS(snakemake@input[['HeLa_19']])
#p=readRDS(snakemake@input[['HeLa_20']])

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
