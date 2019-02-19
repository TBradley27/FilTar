library(gridExtra)
library(grid)
library(ggplot2)

HeLa = readRDS(file = "results/plots/hsa_PRJNA512378_let-7c-5p_HeLa_exp.rds")
HBE4o = readRDS(file = "results/plots/hsa_PRJNA304643_miR-1343-3p_16HBE14o_exp.rds")
A549 = readRDS(file = "results/plots/hsa_PRJNA304643_miR-1343-3p_A549_exp.rds")
U343 = readRDS(file='results/plots/hsa_PRJNA231155_miR-137-3p_U343_exp.rds')
U251 = readRDS(file='results/plots/hsa_PRJNA231155_miR-137-3p_U251_exp.rds')
Huh7 = readRDS(file='results/plots/hsa_PRJNA229375_miR-124-3p_Huh7_exp.rds')
ESCs = readRDS(file='results/plots/mmu_PRJNA270999_miR-294-3p_ESCs_exp.rds')

png('foo.png', width=900, height=900)

grid.arrange(HeLa, A549, U343, U251, Huh7, ESCs,
             ncol=2, nrow=3,
             bottom=textGrob(
               expression(bold('log'[2]*'(mRNA Fold Change)'))
               ), 
             left=textGrob(
               'Cumulative Proportion', 
               rot=90,
               gp=gpar(fontface="bold")
               )
             )
dev.off()
