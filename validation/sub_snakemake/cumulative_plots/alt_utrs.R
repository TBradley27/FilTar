#!/usr/bin/env Rscript

library(plyr)
library(tidyverse)

# read in the data

cl_targets = read_tsv(
	file = snakemake@input$cl_targets,
	col_names = TRUE,
	col_types = 'ccciiiiiccccci',
	trim_ws = FALSE
	)

print(cl_targets)

canon_targets = read_tsv(
        file = snakemake@input$canon_targets,
        col_names = TRUE,
        col_types = 'ccciiiiiccccci',
        trim_ws = FALSE
        )

mock = list()
for (i in 1:length(snakemake@input$quant_mock)) {

print(snakemake@input$quant_mock[i])

mock[[i]] = read_tsv(
        file = snakemake@input$quant_mock[i],
        col_names = TRUE,
        )
}

real = list()
for (i in 1:length(snakemake@input$quant_real)) {

print(snakemake@input$quant_real[i])

real[[i]] = read_tsv(
        file = snakemake@input$quant_real[i],
        col_names = TRUE,
        )
}

# filter targets for 8mers and the correct miRNA

cl_targets = filter(cl_targets, Site_type %in% snakemake@params$site_types) 
cl_targets = filter(cl_targets, miRNA_family_ID == snakemake@wildcards$miRNA)

canon_targets = filter(canon_targets, Site_type %in% snakemake@params$site_types)
canon_targets = filter(canon_targets, miRNA_family_ID == snakemake@wildcards$miRNA)

# get average values between replicates

#colnames(mock1) = colnames(mock1) %>% toupper

print(mock)
print(real)

make_uppercase = function (x) {
	colnames(x) = colnames(x) %>% toupper
	return (x)
}

mock = map(mock, make_uppercase) # for dual salmon/kallisto compatibility
real = map(real, make_uppercase)

print(mock)
print(real)

# get mean TPM values

mean_mock_TPM = vector(length=length(mock[[1]]$TPM))
mean_real_TPM = vector(length=length(real[[1]]$TPM))

for (i in 1:length(mock)) {
mean_mock_TPM = mean_mock_TPM + mock[[i]]$TPM
}
mean_mock_TPM = mean_mock_TPM / length(mock)

for (i in 1:length(real)) {
mean_real_TPM = mean_real_TPM + real[[i]]$TPM
}
mean_real_TPM = mean_real_TPM / length(real)

exp_data = tibble(
	Name=mock[[1]][[1]], # workaround: get the first column of either the kallisto or salmon output file
	mock_TPM = mean_mock_TPM,
	real_TPM = mean_real_TPM
)

print(exp_data)

# remove lowly expressed transcripts

exp_data = filter(exp_data, mock_TPM >= snakemake@params$exp_threshold)
exp_data = filter(exp_data, real_TPM >= snakemake@params$exp_threshold)

# add pseudocount

exp_data$mock_TPM = exp_data$mock_TPM + snakemake@params$pseudocount
exp_data$real_TPM = exp_data$real_TPM + snakemake@params$pseudocount

# get fold changes

exp_data$fold_changes = exp_data$real_TPM / exp_data$mock_TPM
exp_data$log2_fc = exp_data$fold_changes %>% log2

# subset the expression data

print(cl_targets)

exp_data$Name = gsub('\\..*','',exp_data$Name)

non_targets_exp = exp_data$log2_fc[!exp_data$Name %in% canon_targets$a_Gene_ID]

cl_targets = filter(cl_targets, Site_type %in% c('8mer-1a') )  
cl_targets_exp = exp_data$log2_fc[exp_data$Name %in% cl_targets$a_Gene_ID]

canon_targets = filter(canon_targets, Site_type %in% c('8mer-1a') )
canon_targets_exp = exp_data$log2_fc[exp_data$Name %in% canon_targets$a_Gene_ID]

new_targets_names = cl_targets[!cl_targets$a_Gene_ID %in% canon_targets$a_Gene_ID,]
new_targets_exp = exp_data$log2_fc[exp_data$Name %in% new_targets_names$a_Gene_ID]

old_targets_names = canon_targets[!canon_targets$a_Gene_ID %in% cl_targets$a_Gene_ID,]
old_targets_exp = exp_data$log2_fc[exp_data$Name %in% old_targets_names$a_Gene_ID]

#print(length(exp_data$Name))

# build ggplot df

nontargets = tibble(fc=non_targets_exp)
nontargets$legend = stringr::str_interp("No seed binding (n=${length(non_targets_exp)})")

cl_targets = tibble(fc=cl_targets_exp)
cl_targets$legend = stringr::str_interp("${snakemake@wildcards$cell_line} targets (n=${length(cl_targets_exp)})")

canon_targets = tibble(fc=canon_targets_exp)
canon_targets$legend = stringr::str_interp("ensembl targets (n=${length(canon_targets_exp)})")

new_targets = tibble(fc=new_targets_exp)
new_targets$legend = stringr::str_interp("new targets (n=${length(new_targets_exp)})")

old_targets = tibble(fc=old_targets_exp)
old_targets$legend = stringr::str_interp("old targets (n=${length(old_targets_exp)})")

ggplot_df = rbind(nontargets,canon_targets, cl_targets, new_targets, old_targets)

print(ggplot_df)

# begin plotting

ggplot(
  ggplot_df, aes(x=fc,color=legend)
    ) +
  geom_step(aes(y=..y..), stat="ecdf") +
  theme_classic() +
  labs(
        title=
        bquote(
                .(str_interp("${snakemake@wildcards$miRNA}")) ~ .(substitute(italic(vs.))) ~ 'mock transfection' ~ .(str_interp("(${snakemake@wildcards$cell_line})"))
        ),
        y="Cumulative Fraction",
        x=expression('log'[2]*'(mRNA Fold Change)')
	) +
#        subtitle=as.expression(bquote(~ p %~~% .(format (p_value, nsmall=3, digits=3) ) ) )
  theme(legend.title=element_blank()) +
  xlim(-snakemake@params$x_lim,snakemake@params$x_lim)

ggsave(snakemake@output[[1]])
