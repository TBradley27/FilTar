library(sleuth)
library(tidyverse)

### read in and filter target data 

cl_targets = read_tsv(  # read in targets derived from cell-line or tissue-specific annotations
        file = snakemake@input$cl_targets,
        col_names = TRUE,
        col_types = 'ccciiiiiccccci',
        trim_ws = FALSE
        )

print(cl_targets)

canon_targets = read_tsv( # read in targets derived from canonical UTR annotations
        file = snakemake@input$canon_targets,
        col_names = TRUE,
        col_types = 'ccciiiiiccccci',
        trim_ws = FALSE
        )

miRNA_table = readr::read_tsv(
        file = snakemake@input$miRNA_dict,
        col_names = c('family_code','tax_id','mature_miRNA_name','mature_miRNA_sequence'),
        col_types = 'cccc',
        )


species_three_letters = snakemake@wildcards$species
species_tax_id = snakemake@config$tax_ids[[species_three_letters]]

miRNA_name_with_prefix = paste(species_three_letters,snakemake@wildcards$miRNA,sep='-')
miRNA_family = dplyr::filter(miRNA_table, mature_miRNA_name == miRNA_name_with_prefix)
miRNA_family = miRNA_family$family_code[1]

canon_targets = filter(canon_targets, Site_type %in% snakemake@params$nontarget_site_types)
canon_targets = filter(canon_targets, miRNA_family_ID == miRNA_family)
canon_targets = filter(canon_targets, species_ID == species_tax_id)

cl_targets = filter(cl_targets, Site_type %in% snakemake@params$nontarget_site_types)
cl_targets = filter(cl_targets, miRNA_family_ID == miRNA_family)
cl_targets = filter(cl_targets, species_ID == species_tax_id)

## prepare data

kal_dirs = c(snakemake@input[['quant_mock']], snakemake@input[['quant_real']]) %>% unlist()

kal_dirs = gsub(pattern='/abundance.tsv',replacement='',x=kal_dirs)

print(kal_dirs)

samples = c(
	snakemake@config[[snakemake@wildcards$project]][[snakemake@wildcards$cell_line]]$mock,
	snakemake@config[[snakemake@wildcards$project]][[snakemake@wildcards$cell_line]][[snakemake@wildcards$miRNA]]
	)
 
s2c = data.frame(
	sample=samples,
	condition=rep(c("negative_control","miRNA"),each=2),
	path=kal_dirs,
	stringsAsFactors=FALSE
)

print(s2c)

## sleuth

my_filter <- function(row, min_reads = 1, min_prop = 0.10) { mean(row >= min_reads) >= min_prop }

so = sleuth_prep(s2c, extra_bootstrap_summary = TRUE, filter_fun = my_filter) # load kallisto data into the object
so = sleuth_fit(so, ~condition, 'full') # estimate parameters for sleuth 'full' model

#so = sleuth_fit(so, ~1, 'reduced') # estimate parameters for 'reduced' model in which equal abundances are assumed
#so = sleuth_lrt(so, 'reduced', 'full') # differential expression analysis

so = sleuth_wt(so, which_beta='conditionnegative_control', which_model='full')

saveRDS(so, file = "my_data.rds")

models(so) %>% print
tests(so) %>% print

sleuth_table <- sleuth_results(so, 'conditionnegative_control', 'wt', show_all = FALSE)
sleuth_table$b = as.numeric(sleuth_table$b)
sleuth_table$b = -sleuth_table$b


print(sleuth_table[1:100,])

print(sleuth_table$b[1:100])

foo = extract_model(so, which_model='full')
print(foo)

# divide sleuth output into relevant distributions

non_targets_exp = sleuth_table$b[!sleuth_table$target_id %in% cl_targets$a_Gene_ID]
#non_targets_exp = non_targets_exp[!is.na(non_targets_exp)]

non_targets_exp = non_targets_exp - median(non_targets_exp)

cl_targets = filter(cl_targets, Site_type %in%  snakemake@params$target_site_types)
cl_targets_exp = sleuth_table$b[sleuth_table$target_id %in% cl_targets$a_Gene_ID]
#cl_targets_exp = cl_targets_exp[!is.na(cl_targets_exp)]

cl_targets_exp = cl_targets_exp - median(non_targets_exp)

#print(cl_targets_exp)

canon_targets = filter(canon_targets, Site_type %in% snakemake@params$target_site_types)
canon_targets_exp = sleuth_table$b[sleuth_table$target_id %in% canon_targets$a_Gene_ID]
#canon_targets_exp =canon_targets_exp[!is.na(canon_targets_exp)]

canon_targets_exp = canon_targets_exp - median(non_targets_exp)

#print(canon_targets_exp)

new_targets_names = cl_targets[!cl_targets$a_Gene_ID %in% canon_targets$a_Gene_ID,]
new_targets_exp = sleuth_table$b[sleuth_table$target_id %in% new_targets_names$a_Gene_ID]
#new_targets_exp = new_targets_exp[!is.na(new_targets_exp)]

new_targets_exp = new_targets_exp - median(non_targets_exp)

#print(new_targets_exp)

old_targets_names = canon_targets[!canon_targets$a_Gene_ID %in% cl_targets$a_Gene_ID,]
old_targets_exp = sleuth_table$b[sleuth_table$target_id %in% old_targets_names$a_Gene_ID]
#old_targets_exp = old_targets_exp[!is.na(old_targets_exp)]

old_targets_exp = old_targets_exp - median(non_targets_exp)


#print(old_targets_exp)

# build ggplot df

nontargets = tibble(fc=non_targets_exp)
nontargets$legend = stringr::str_interp("No seed binding (n=${length(non_targets_exp)})")

print(cl_targets_exp)

cl_targets = tibble(fc=cl_targets_exp)
cl_targets$legend = stringr::str_interp("${snakemake@wildcards$cell_line} targets (n=${length(cl_targets_exp)})")

print(cl_targets)

canon_targets = tibble(fc=canon_targets_exp)
canon_targets$legend = stringr::str_interp("canonical targets (n=${length(canon_targets_exp)})")

new_targets = tibble(fc=new_targets_exp)
new_targets$legend = stringr::str_interp("new targets (n=${length(new_targets_exp)})")

old_targets = tibble(fc=old_targets_exp)
old_targets$legend = stringr::str_interp("old targets (n=${length(old_targets_exp)})")

ggplot_df = rbind(nontargets,canon_targets, new_targets, old_targets)

print(ggplot_df)

p_value = ks.test(new_targets_exp, non_targets_exp, alternative='less')

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
        x=expression('log'[2]*'(mRNA Fold Change)'),
        subtitle=as.expression(bquote(~ p %~~% .(format (p_value, nsmall=3, digits=3) ) ) ) )  +
  theme(legend.title=element_blank()) +
  coord_cartesian(xlim = c(-snakemake@params$x_lim,snakemake@params$x_lim))

ggsave(snakemake@output[[1]])







