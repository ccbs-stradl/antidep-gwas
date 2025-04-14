# collate tables of LDSC rg for meta

library(dplyr)
library(tidyr)
library(readr)
library(stringr)


# between meta ldsc log directory
log_dir = here::here("results/models/rg/meta")

# between ldsc log files
ldsc_log_files <- list.files(log_dir,
  pattern = "antidep-2501.+antidep-2501.+log",
  full.names = TRUE)

# external ldsc log directory
ext_log_dir <- here::here("results/models/rg/meta/external")

ext_ldsc_log_files <- list.files(ext_log_dir,
  pattern = "antidep-2501.+log",
  full.names = TRUE)

# parse results table from ldsc file
parse_ldsc_table <- function(ldsc_file) {
  ldsc_log <- readr::read_lines(ldsc_file)
  # index location of table
  start_idx <- str_which(ldsc_log, pattern = "Summary of Genetic Correlation Results") + 1
  end_idx <- str_which(ldsc_log, pattern = "Analysis finished at") - 1

  # extract lines and concatenate back into a string, so we can parse it as a table
  ldsc_table_lines <- ldsc_log[start_idx:end_idx]
  ldsc_table_str <- str_c(ldsc_table_lines, sep = "\n")
  ldsc_table <- read_table(ldsc_table_str)
  return(ldsc_table)
}

ldsc_tables <- bind_rows(lapply(ldsc_log_files, parse_ldsc_table))
ext_ldsc_tables <- bind_rows(lapply(ext_ldsc_log_files, parse_ldsc_table))

# parse meta information from filename
ldsc_datasets <- ldsc_tables |>
  separate(p1, into = c("p1_meta", "p1_version", "p1_method", "p1_pheno", "p1_cluster", "p1_ext"),
           sep = "[-.]+", extra = "merge") |>
  separate(p2, into = c("p2_meta", "p2_version", "p2_method", "p2_pheno", "p2_cluster", "p2_ext"),
           sep = "[-.]+", extra = "merge") |>
  # remove column with file extension
  select(-ends_with('_ext'))

ext_ldsc_datasets <- ext_ldsc_tables |>
  separate(p1, into = c("p1_meta", "p1_version", "p1_method", "p1_pheno", "p1_cluster", "p1_ext"),
           sep = "[-.]+", extra = "merge") |>
  separate(p2, into = c("p2_pheno", "p2_ext"),
           sep = "[-.]+", extra = "merge") |>
  # remove column with file extension
  select(-ends_with('_ext'))

write_csv(ldsc_datasets, here::here("manuscript/tables/rg_ldsc_meta.csv"))
write_csv(ext_ldsc_datasets, here::here("manuscript/tables/rg_ldsc_meta_external.csv"))
