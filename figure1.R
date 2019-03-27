library(gridExtra)
library(grid)
library(ggplot2)

A549 = readRDS(file = "results/plots/hsa_PRJNA304643_miR-1343-3p_A549_exp.rds")
HeLa = readRDS(file = "results/plots/hsa_PRJNA512378_miR-16-5p_HeLa_exp.rds")
#U343 = readRDS(file='results/plots/hsa_PRJNA231155_miR-137-3p_U343_exp.rds')
#CD4 = readRDS(file='results/plots/mmu_PRJNA309441_miR-23a-3p_CD4_exp.rds')
NMuMG = readRDS(file='results/plots/mmu_PRJNA340017_miR-1199-5p_NMuMG_exp.rds')
ESCs = readRDS(file='results/plots/mmu_PRJNA270999_miR-294-3p_ESCs_exp.rds')

#png('results/plots/figure_1_alt.png', width=700, height=700)

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
ggsave('results/plots/figure_1_alt.png',plot = g, scale=1.25)
