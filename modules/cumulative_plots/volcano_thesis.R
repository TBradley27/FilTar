library(gridExtra)
library(grid)
library(ggplot2)

## open device

a=readRDS(snakemake@input[['U251']])
b=readRDS(snakemake@input[['U343']])
c=readRDS(snakemake@input[['Du145']])
d=readRDS(snakemake@input[['HBE14o']])
e=readRDS(snakemake@input[['U20S_1']])
f=readRDS(snakemake@input[['U20S_2']])
g=readRDS(snakemake@input[['NMuMG_1']])
h=readRDS(snakemake@input[['NMuMG_2']])
i=readRDS(snakemake@input[['NMuMG_2']])
j=readRDS(snakemake@input[['CD4_1']])
k=readRDS(snakemake@input[['CD4_2']])
l=readRDS(snakemake@input[['CD4_3']])

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
