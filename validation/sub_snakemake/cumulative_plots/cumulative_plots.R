#!/usr/bin/env Rscript

library(readr)
library(dplyr)
library(purrr)
library(stringr)
library(ggplot2)

# read in the data

targets = read_tsv(
	file = snakemake@input$targets,
	col_names = FALSE,
	trim_ws = FALSE # this is important to ensure alignments can be viewed accurately
	)

mock1 = read_tsv(
        file = snakemake@input$quant_mock[1],
        col_names = TRUE,
        )

mock2 = read_tsv(
        file = snakemake@input$quant_mock[2],
        col_names = TRUE,
        )

real1 = read_tsv(
        file = snakemake@input$quant_real[1],
        col_names = TRUE,
        )

real2 = read_tsv(
        file = snakemake@input$quant_real[2],
        col_names = TRUE,
        )

# filter targets for 8mers and the correct miRNA

targets = filter(targets, X4 %in% c('8mer-1a','7mer-m8')) 

if (snakemake@wildcards$miRNA == 'miR-124') {
	targets = filter(targets, X3 == 'miR-124-3p')
} else if (snakemake@wildcards$miRNA == 'miR-155') {
	targets = filter(targets, X3 == 'miR-155-5p')
}

# get average values between replicates

exp_data = tibble(
	Name=mock1$Name,
	mock_TPM = ( mock1$TPM + mock2$TPM ) / 2,
	real_TPM = ( real1$TPM + real2$TPM ) / 2
)

# remove lowly expressed transcripts

exp_data = filter(exp_data, mock_TPM > 0)

# add pseudocount

exp_data$mock_TPM = exp_data$mock_TPM + 1
exp_data$real_TPM = exp_data$real_TPM + 1

# get fold changes

exp_data$fold_changes = exp_data$real_TPM / exp_data$mock_TPM
exp_data$log2_fc = exp_data$fold_changes %>% log2

## get 13-16 nt window

# convert spaces to dashes in the alignment just so that it is easier to view

targets$supp_window = gsub(' ','-', targets$X34)

targets$supp_window = targets$supp_window %>% strsplit('')
targets$supp_window = map(targets$supp_window, function(x) { rev(x) %>% paste(collapse='') %>% substr(13,16) } )

supp_targets = subset(targets, targets$supp_window == "||||")


# subset target data

targets_fc = exp_data$log2_fc[exp_data$Name %in% targets$X1]
supp_targets_fc = exp_data$log2_fc[exp_data$Name %in% supp_targets$X1]
nontargets_fc = exp_data$log2_fc[!exp_data$Name %in% targets$X1]

p_value = ks.test(supp_targets_fc, targets_fc, alternative='greater')$p.value

# plot

#png(filename=snakemake@output[[1]])

#plot(ecdf(nontargets_fc), cex=0.1)
#plot(ecdf(targets_fc), add=T, col='red', cex=0.1)
#plot(ecdf(supp_targets_fc), add=T, col='green', cex=0.1)

nontargets = tibble(fc=nontargets_fc)
nontargets$legend = stringr::str_interp("No seed binding (n=${length(nontargets_fc)})")

targets = tibble(fc=targets_fc)
targets$legend = stringr::str_interp("Seed binding (n=${length(targets_fc)})")

supp_targets = tibble(fc=supp_targets_fc)
supp_targets$legend = stringr::str_interp("Seed with supp. binding (n=${length(supp_targets_fc)})")

ggplot_df = rbind(nontargets,targets,supp_targets)

ggplot(
  ggplot_df, aes(x=fc,color=legend)
    ) + 
  geom_step(aes(y=..y..), stat="ecdf") +
  theme_classic() + 
  labs(
	title=
	bquote(
		.(str_interp("${snakemake@wildcards$miRNA}")) ~ .(substitute(italic(vs.))) ~ 'mock transfection' ~ .(str_interp("(${snakemake@wildcards$cell_line})"))
	)
	,
	y="Cumulative Fraction",
	x=expression('log'[2]*'(mRNA Fold Change)'),
	subtitle=as.expression(bquote(~ p %~~% .(format (p_value, nsmall=3, digits=3) ) ) )
) + 
  theme(legend.title=element_blank()) +
  xlim(-3,3)

ggsave(snakemake@output[[1]])

