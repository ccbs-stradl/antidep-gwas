# collate tables of Popcorn rg for meta-analyses

library(dplyr)
library(tidyr)
library(readr)
library(stringr)


# popcorn log directory
log_dir = here::here("results/models/popcorn/meta")

popcorn_files <- list.files(log_dir,
                            pattern = "antidep-2501.+antidep-2501",
                            full.names = TRUE)

names(popcorn_files) <- basename(popcorn_files)

popcorn_tables <- lapply(popcorn_files, read_tsv)

# parse dataset information from filename
popcorns <- bind_rows(popcorn_tables, .id = "filename") |>
    rename(estimate = ...1) |>
    mutate(datasets = str_remove(filename, pattern = ".correlations.txt")) |>
    separate(datasets, into = c("p1", "p2"), sep = "--") |>
    separate(p1, into = c("p1_meta", "p1_version", "p1_method", "p1_pheno", "p1_cluster"), sep = "-") |>
    separate(p2, into = c("p2_meta", "p2_version", "p2_method", "p2_pheno", "p2_cluster"), sep = "-") |>
    filter(estimate == "pgi") |>
    select(-filename, -estimate) |>
    select(starts_with("p1"), starts_with("p2"), pgi = `Val (obs)`, everything())

# keep cross-cluster comparisons
popcorns_keep <- popcorns |>
  filter(p1_cluster != p2_cluster)

file_name <- here::here("manuscript/tables/rg_popcorn_meta.csv")
write_csv(popcorns_keep, file_name)

# create .cols sidecar meta data file
popcorns_keep <- fread(file_name)

# create .cols sidecar meta data file
colname_descriptions <- c("p1_meta" = "Name of first meta-analysed phenotype",
                          "p1_version" = "Version of first meta-analysed phenotype", 
                          "p1_method" = "Method of meta-analysis for first phenotype", 
                          "p1_pheno" = "Anti-depressant phenotype of first meta-analysed phenotype", 
                          "p1_cluster" = "Ancestry cluster of first meta-analysed phenotype", 
                          "p2_meta" = "Name of second meta-analysed phenotype", 
                          "p2_version" = "Version of second meta-analysed phenotype", 
                          "p2_method" = "Method of meta-analysis for second phenotype", 
                          "p2_pheno" = "Anti-depressant phenotype of second meta-analysed phenotype", 
                          "p2_cluster" = "Ancestry cluster of second meta-analysed phenotype", 
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
