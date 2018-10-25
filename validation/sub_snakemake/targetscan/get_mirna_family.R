#!/bin/env Rscript

library(readr)
library(dplyr)
library(stringr)
library(purrr)

mirna_seeds = read_table2(snakemake@input[[1]], col_names=c("identifier", "seq"))

# group miRNA families together
mirna_seeds = mirna_seeds[order(mirna_seeds$seq),]

# return the three letter identifier for each species
mirna_seeds$species = stringr::str_sub(mirna_seeds$identifier, 1,3)

TaxID = list(hsa="9606", ptr="9598", ggo="9595",   pab="9601",   nle="9581",    rma="9544",   mfa="9541",  
  pan="9557",   csa="60711",   cja="9483",   sbo="27679",   oga="30611",   tch="246437",   str="43179",   jja="51337",  
  moc="79684",   cgr="10029",   mau="10036",   mmu="10090",   rno="10116",   hgl="10181",   cpo="10141",   cla="34839",  
  ode="10160",   ocu="9986",   opr="9978",   ssc="9823",   vpa="30538",   cfe="419612",   ttr="9739",   oor="9733",  
  pod="59538",   bta="9913",   oar="9940",   chi="9925",   eca="9796",   csi="9807",   fca="9685",   cfa="9615",  
  mfu="9669",   ame="9646",   oro="9708",   lwe="9713",   pal="9402",   pva="132908",   mda="225400",   mlu="59463", 
  efu="29078",   eeu="9365",   sar="42254",   ccr="143302",   laf="9785",   eed="28737",   tma="127582",   cas="185453", 
  ete="9371",   oaf="1230840",   dno="9361",   mdo="13616",   sha="9305",   meu="9315",   oan="9258",   fch="345164", 
  fpe="9854",   fal="59894",   zal="44394",   gfo="48883",   tgu="59729",   phu="181119",   mun="13146",   avi="241585", 
  ama="176014",   cli="8932",  apl="8839",   gga="9031",  mga="9103",   ami="8496",   cmy="8469",  cpi="8478", 
  psi="13735",   asp="55534",   aca="28377",   xtr="8364",   lch="7897")

mirna_seeds = mirna_seeds[mirna_seeds$species %in% names(TaxID),]

map_ids = function(string) {
  tax_id = TaxID[[string]]
  
  return (tax_id)
}

mirna_seeds$tax_id = map(mirna_seeds$species, map_ids)
mirna_seeds$species = NULL

mirna_seeds$tax_id = as.character(mirna_seeds$tax_id)
mirna_seeds$seq = as.factor(mirna_seeds$seq)
mirna_seeds$identifier = as.integer(mirna_seeds$seq)

write.table(mirna_seeds, snakemake@output[[1]], sep="\t", quote=FALSE, col.names=FALSE, row.names=FALSE)







