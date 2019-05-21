test = filtar::get_mirna_context(snakemake@input$mirna_seed, snakemake@input$mature_mirnas, snakemake@wildcards$species)

readr::write_tsv(test, snakemake@output[[1]], col_names=FALSE)
