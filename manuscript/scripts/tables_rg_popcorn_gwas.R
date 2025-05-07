# collate tables of Popcorn rg for input GWAS

library(dplyr)
library(tidyr)
library(readr)
library(stringr)


# popcorn log directory
log_dir = here::here("results/models/popcorn/gwas")

# metaset to analyse
metaset <- read_csv(here::here("inputs/metasets/antidep-2501.csv"))

# popcorn output files
popcorn_files <- list.files(log_dir,
  pattern = "*.correlations.txt",
  full.names = TRUE)

names(popcorn_files) <- basename(popcorn_files)

popcorn_tables <- lapply(popcorn_files, read_tsv)

# parse dataset information from filename
popcorns <- bind_rows(popcorn_tables, .id = "filename") |>
    rename(estimate = ...1) |>
    mutate(datasets = str_remove(filename, pattern = ".correlations.txt")) |>
    separate(datasets, into = c("p1", "p2"), sep = "--") |>
    separate(p1, into = c("p1_cohort", "p1_pheno", "p1_cluster", "p1_version"), sep = "-") |>
    separate(p2, into = c("p2_cohort", "p2_pheno", "p2_cluster", "p2_version"), sep = "-") |>
    filter(estimate == "pgi") |>
    select(-filename, -estimate) |>
    select(starts_with("p1"), starts_with("p2"), pgi = `Val (obs)`, everything())

# keep combinations that are in the metaset, cross-cluster only
cohorts_versions <- metaset |>
  mutate(cohort_version = str_c(cohort, version, sep = "-")) |>
  pull(cohort_version)

popcorns_keep <- popcorns |>
  filter(str_c(p1_cohort, p1_version, sep = "-") %in% cohorts_versions,
         str_c(p2_cohort, p2_version, sep = "-") %in% cohorts_versions) |>
         filter(p1_cluster != p2_cluster) |>
         filter(p1_pheno != "N06AX", p2_pheno != "N06AX")

file_name <- here::here("manuscript/tables/rg_popcorn_gwas.csv")
write_csv(popcorns_keep, file_name)

popcorns_keep <- fread(file_name)

# create .cols sidecar meta data file
colname_descriptions <- c("p1_cohort" = "Name of first cohort",
                          "p1_pheno" = "Anti-depressant phenotype of first cohort",
                          "p1_cluster" = "Ancestry cluster of first cohort",
                          "p1_version" = "Version of first cohort",
                          "p2_cohort" = "Name of second cohort",
                          "p2_pheno" = "Anti-depressant phenotype of second cohort",
                          "p2_cluster" = "Ancestry cluster of second cohort",
                          "p2_version" = "Version of second cohort",
                          "pgi" = "Genetic impact correlation",
                          "SE" = "Standard error of genetic impact correlation",
                          "Z" = "Z statistic of genetic impact correlation",
                          "P (Z)" = "p-value of Z statistic of genetic impact correlation"
)


source(here::here("manuscript/scripts/supplementary_tables_excell_create_cols_meta_FUN.R"))

create_cols_meta(
  file_name = file_name,
  table_variable_name = popcorns_keep,
  colname_descriptions = colname_descriptions
)
