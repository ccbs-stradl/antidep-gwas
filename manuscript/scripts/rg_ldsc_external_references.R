# create .cols sidecar meta data file for
# rg_ldsc_external_references.csv
source(here::here("manuscript/scripts/supplementary_tables_excell_create_cols_meta_FUN.R"))

file_name <- here::here("manuscript/tables/rg_ldsc_external_references.csv")

external_datasets <- read.csv(file_name)

colname_descriptions <- c(
  "Abbreviation" = "Abbreviation of the phenotype",
  "Phenotype" = "Phenotype name",
  "Reference" = "Publication reference",
  "Publication" = "Link to publication",
  "Source" = "Cohort source",
  "Download" = "Link to download the dataset"
  )

create_cols_meta(
  file_name = file_name,
  table_variable_name = external_datasets,
  colname_descriptions = colname_descriptions
)
