library(gridExtra)
library(grid)
library(ggplot2)

HeLa = readRDS(file = "results/plots/hsa_PRJNA512378_miR-16-5p_HeLa_alt_utr_old.rds")
A549 = readRDS(file = "results/plots/hsa_PRJNA304643_miR-1343-3p_A549_alt_utr_old.rds")
#HBE14 = readRDS(file = 'results/plots/hsa_PRJNA304643_miR-1343-3p_16HBE14o_alt_utr_old.rds')
#U343 = readRDS(file='results/plots/hsa_PRJNA231155_miR-137-3p_U343_alt_utr_old.rds')
#Du145 = readRDS(file='results/plots/hsa_PRJNA292016_miR-141-3p_Du145_alt_utr_old.rds')
#CD4 = readRDS(file='results/plots/mmu_PRJNA309441_miR-23a-3p_CD4_alt_utr.rds')
NMuMG = readRDS(file='results/plots/mmu_PRJNA340017_miR-1199-5p_NMuMG_alt_utr_old.rds')
ESCs = readRDS(file='results/plots/mmu_PRJNA270999_miR-294-3p_ESCs_alt_utr_old.rds')

g = arrangeGrob(grid.arrange(A549, HeLa, NMuMG, ESCs,
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
)

ggsave('results/plots/figure_3_alt.png',plot = g, scale=1.50)
