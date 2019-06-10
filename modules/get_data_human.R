library(tidyverse)

human_data = list (
  # ENCODE project - 161 
  PRJNA30709 = readr::read_tsv(url("https://www.ebi.ac.uk/gxa/experiments-content/E-GEOD-26284/resources/experiment-design")),
  # Skeletal Muscle - 26
  PRJNA252429 = readr::read_tsv(url("https://www.ebi.ac.uk/gxa/experiments-content/E-GEOD-58387/resources/experiment-design")),
  # The Human Protein Atlas - 200
  PRJEB4337 = readr::read_tsv(url("https://www.ebi.ac.uk/gxa/experiments-content/E-MTAB-2836/resources/experiment-design")),
  # Multiple Organs - 19
  PRJNA143627 = readr::read_tsv(url("https://www.ebi.ac.uk/gxa/experiments-content/E-MTAB-3716/resources/experiment-design")),
  # ENCODE project - multiple organs - strand-specific - 25 - strange identifiers
  #PRJNA30709 = readr::read_tsv(url("https://www.ebi.ac.uk/gxa/experiments-content/E-MTAB-4344/resources/experiment-design")),
  # Bronchus epithelium
  PRJEB13935 = readr::read_tsv(url("https://www.ebi.ac.uk/gxa/experiments-content/E-MTAB-4729/resources/experiment-design")),
  #Illumina Body Map - tissues and organs - 16
  PRJEB2445 = readr::read_tsv(url("https://www.ebi.ac.uk/gxa/experiments-content/E-MTAB-513/resources/experiment-design"))

  # Haemopoetic lineages - blueprint project - 96 - Managed access
  #PRJNA143627 = readr::read_tsv(url("https://www.ebi.ac.uk/gxa/experiments-content/E-MTAB-3819/resources/experiment-design"))
  #GTEx - Raw data not available
  #HipSci Project - cannot find data
)

print('foo')

human_data$PRJEB6971 = human_data$PRJEB4337[172:200,]
human_data$PRJEB4337 = human_data$PRJEB4337[1:171,]

RemoveColumn = function(df, columnName) {
  if ( columnName %in% colnames(df)) {
    df = df[, !names(df) %in% columnName]
    return(df)
  }
  else {
    return(df)
  }
}

RenameColumn = function(df, columnName, newcolumnName) {
  if ( columnName %in% colnames(df)) {
    names(df)[names(df) == columnName] <- newcolumnName
    return(df)
  }
  else {
    return(df)
  }
}

# This is obviously bad code
human_data = purrr::map(human_data, RemoveColumn, "Sample Characteristic Ontology Term[age]")
human_data = purrr::map(human_data, RemoveColumn, "Sample Characteristic Ontology Term[cell type]")
human_data = purrr::map(human_data, RemoveColumn, "Sample Characteristic Ontology Term[organism]")
human_data = purrr::map(human_data, RemoveColumn, "Sample Characteristic Ontology Term[organism part]")
human_data = purrr::map(human_data, RemoveColumn, "Sample Characteristic Ontology Term[sex]")
human_data = purrr::map(human_data, RemoveColumn, "Sample Characteristic Ontology Term[strain]")
human_data = purrr::map(human_data, RemoveColumn, "Factor Value[cell type]")
human_data = purrr::map(human_data, RemoveColumn, "Factor Value Ontology Term[cell type]")
human_data = purrr::map(human_data, RemoveColumn, "Analysed")
human_data = purrr::map(human_data, RemoveColumn, "Factor Value[organism part]")
human_data = purrr::map(human_data, RemoveColumn, "Factor Value Ontology Term[organism part]")
human_data = purrr::map(human_data, RemoveColumn, "Sample Characteristic Ontology Term[developmental stage]")
human_data = purrr::map(human_data, RemoveColumn, "Sample Characteristic Ontology Term[genotype]")
human_data = purrr::map(human_data, RemoveColumn, "Factor Value[strain]")
human_data = purrr::map(human_data, RemoveColumn, "Factor Value Ontology Term[strain]")
human_data = purrr::map(human_data, RemoveColumn, "Factor Value[developmental stage]")
human_data = purrr::map(human_data, RemoveColumn, "Factor Value Ontology Term[developmental stage]")
human_data = purrr::map(human_data, RemoveColumn, "Sample Characteristic[biosource provider]")
human_data = purrr::map(human_data, RemoveColumn, "Sample Characteristic Ontology Term[biosource provider]")
human_data = purrr::map(human_data, RemoveColumn, "Sample Characteristic Ontology Term[strain or line]")
human_data = purrr::map(human_data, RemoveColumn, "Sample Characteristic Ontology Term[individual]")
human_data = purrr::map(human_data, RemoveColumn, "Sample Characteristic Ontology Term[individual]")
human_data = purrr::map(human_data, RemoveColumn, "Sample Characteristic Ontology Term[cell line]")
human_data = purrr::map(human_data, RemoveColumn, "Sample Characteristic Ontology Term[cellular component]")
#human_data = purrr::map(human_data, RemoveColumn, "Sample Characteristic Ontology Term[disease state]")
#human_data = purrr::map(human_data, RemoveColumn, "Sample Characteristic Ontology Term[karyotype]")
human_data = purrr::map(human_data, RemoveColumn, "Sample Characteristic Ontology Term[biopsy site]")
human_data = purrr::map(human_data, RemoveColumn, "Sample Characteristic Ontology Term[clinical information]")
#human_data = purrr::map(human_data, RemoveColumn, "Sample Characteristic Ontology Term[disease]")
human_data = purrr::map(human_data, RemoveColumn, "Sample Characteristic Ontology Term[ethnic group]")
human_data = purrr::map(human_data, RemoveColumn, "Sample Characteristic[sex]")
human_data = purrr::map(human_data, RemoveColumn, "Sample Characteristic[biopsy site]")
human_data = purrr::map(human_data, RemoveColumn, "Sample Characteristic[ethnic group]")
human_data = purrr::map(human_data, RemoveColumn, "Sample Characteristic[clinical information]")
human_data = purrr::map(human_data, RemoveColumn, "Factor Value Ontology Term[cell line]")
human_data = purrr::map(human_data, RemoveColumn, "Factor Value Ontology Term[RNA]")
human_data = purrr::map(human_data, RemoveColumn, "Factor Value[sex]")
human_data = purrr::map(human_data, RemoveColumn, "Factor Value Ontology Term[sex]")
human_data = purrr::map(human_data, RemoveColumn, "Factor Value Ontology Term[cellular component]")
human_data = purrr::map(human_data, RemoveColumn, "Sample Characteristic[individual]")
human_data = purrr::map(human_data, RenameColumn, "Sample Characteristic[developmental stage]","Sample Characteristic[age]" )
human_data = purrr::map(human_data, RenameColumn, "Sample Characteristic[strain or line]","Sample Characteristic[strain]" )


# Remove PRJNA229389 because study says samples were taken at postnatal day 7
# PRJNA267539 is very iffy. Something to do with a bone marrow transplant. 
#   Mice less < 8 weeks old
# PRJNA283850 - does not say - I may need to read the paper more carefully
# PRJEB6329 - As far as I can tell, age has not been specified
#  PRJNA177791 = animals of breeding age

# add mising age data
# http://science.sciencemag.org/content/sci/suppl/2012/12/19/338.6114.1593.DC1/Merkin.SM.pdf
#human_data$PRJNA177791$`Sample Characteristic[age]` = 'breeding age'

#human_data[['PRJNA187158']][['Project']] = 'PRJNA187158'

add_project_column = function(df,string) {
  df$Project = string
  return (df)
}

# incorporate project names into the data frames
human_data = map2(human_data, names(human_data), add_project_column )

merged_data = dplyr::bind_rows(human_data)

# Delete aberrant records
merged_data = merged_data[
  merged_data$`Sample Characteristic[cellular component]` %in% 
    c("whole cell",NA),]

colnames(merged_data) %>% print()

merged_data = merged_data[
  merged_data$`Sample Characteristic[disease state]` %in% 
    c("normal",NA),]

#merged_data = merged_data[
#  merged_data$`Sample Characteristic[disease]` %in% 
#    c("normal",NA),]

#merged_data = merged_data[
#  merged_data$`Sample Characteristic[karyotype]` %in% 
#    c("normal",NA),]

#merged_data = merged_data[
#  merged_data$`Factor Value[RNA]` %in% 
#    c("long polyA RNA",NA),]

#merged_data = merged_data[
#  !merged_data$`Sample Characteristic[age]` %in% 
#    c(NA),]

print(merged_data)

# information hard-coded from browsing ENA for relevant information
#single_end_projects = c('PRJNA187158', 'PRJNA187161', 'PRJNA193409', 
#                        'PRJNA143627', 'PRJEB13590')
paired_end_projects = c('PRJEB4337', 'PRJEB6971')

mixed_single_and_paired = c('PRJEB2445')

merged_data$Run_Type = map(
  merged_data$Project, function(x) if (x %in% paired_end_projects) 
  {'paired'} 
  else 
  {'mixed_single_and_paired'}
) %>% as.character()

# warning: This section is hard-coded from ENA web data
merged_data[172:187,"Run_Type"] = 'paired' # FROM ERR030872 - ERR030887
merged_data[188:203,"Run_Type"] = 'single'

cleaned_data = merged_data[,c(
  'Project',
  'Run',
  'Sample Characteristic[organism part]',
  'Run_Type'
)]

print('bar')

write.table(cleaned_data, file = "human_experiments.tsv", 
            quote = FALSE, sep='\t', col.names=TRUE, row.names=FALSE)

#downloaded = readr:::read_tsv("mouse_downloaded_datasets.txt", col_names=FALSE)
#downloaded$X1 = stringr::str_replace(downloaded$X1, '.fastq.gz', '')
#downloaded$X1 = stringr::str_replace(downloaded$X1, '_1', '')
#downloaded$X1 = stringr::str_replace(downloaded$X1, '_2', '')
#downloaded = downloaded[!duplicated(downloaded$X1),]

#missing = merged_data[!merged_data$Run %in% downloaded$X1 ,]

#library(yaml)

#as.yaml(list(foo=1:10, bar=c("test1", "test2"))) %>% cat()

#df <- data.frame(name = c("bob","joe"),
#                 target = c(c("yellow","blue"), "grey"),
#                 code1 = c("fly", "walk"),
#                 code2 = c("jump", "run"))

#list(
#  samples=split(replace(df, "name", NULL), 
#                df$name)
#     )

#as.yaml(list(samples=split(replace(df, "name", NULL), df$name))) %>% cat()

#cleaned_data$Project = NULL

#cleaned_data$Run_Type = unlist(cleaned_data$Run_Type)

#list(
#  tissues=split(
#    cleaned_data$Run_Type, 
#    cleaned_data$`Sample Characteristic[organism part]`)
#) %>% as.yaml() %>% cat
