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

write_csv(popcorns_keep, here::here("manuscript/tables/rg_popcorn_meta.csv"))
