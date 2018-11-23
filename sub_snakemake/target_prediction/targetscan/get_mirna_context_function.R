get_mirna_context = function (mirna_seed_file, mature_mirnas_file, species) {

TaxID = list(hsa="9606", mmu="10090")

specific_tax_id = TaxID[[species]]

mirna_seeds = read_tsv(mirna_seed_file, col_names=c("identifier", "seq", "tax_id"))
print(mirna_seeds)

mirnas = read_tsv(mature_mirnas_file, col_names = c('miRNA_family','tax_id','miRNA','miRNA_sequence'))
print(mirnas)
mirnas$seed = stringr::str_sub(mirnas$miRNA_sequence, 2, 8)

test = merge(mirnas, mirna_seeds, by.x = 'seed', by.y = 'seq')
test = subset(test, test$tax_id.y == specific_tax_id)

test = test[!duplicated(test), ]
test = test[,c('identifier','tax_id.x', 'miRNA', 'miRNA_sequence')]

return(test)
}
