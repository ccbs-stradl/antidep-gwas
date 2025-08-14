# collate tables of LDSC rg for input GWAS

library(dplyr)
library(tidyr)
library(readr)
library(stringr)
library(glue)


# ldsc log directory
log_dir = here::here("results/models/rg/gwas")

# metaset to analyse
metaset <- read_csv(here::here("inputs/metasets/antidep-2501.csv"))

# ldsc log files
ldsc_log_files <- list.files(log_dir,
  pattern = "*.log",
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

# parse meta information from filename
ldsc_datasets <- ldsc_tables |>
  separate(p1, into = c("p1_cohort", "p1_pheno", "p1_cluster", "p1_version", "ext"),
           sep = "[-.]+", extra = "merge") |>
  separate(p2, into = c("p2_cohort", "p2_pheno", "p2_cluster", "p2_version", "ext"),
           sep = "[-.]+", extra = "merge")

# keep combinations that are in the metaset
cohorts_versions <- metaset |>
  mutate(cohort_version = str_c(cohort, version, sep = "-")) |>
  pull(cohort_version)

ldsc_datasets_keep <- ldsc_datasets |>
  filter(str_c(p1_cohort, p1_version, sep = "-") %in% cohorts_versions,
         str_c(p2_cohort, p2_version, sep = "-") %in% cohorts_versions) |>
  filter(p1_pheno != "N06AX", p2_pheno != "N06AX") |>
  select(-ext) |>
  select(starts_with("p1"), starts_with("p2"), everything()) |>
  mutate(ci = str_c("(", round(rg + qnorm(0.025) * se, 3), ", ", round(rg + qnorm(0.975) * se, 3), ")")) |>
  relocate(ci, .after = se)

file_name <- here::here("manuscript/tables/rg_ldsc_gwas.csv")
write_csv(ldsc_datasets_keep, file_name)

# create .cols sidecar meta data file
colname_descriptions <- c("p1_cohort" = "Name of first cohort",
                          "p1_pheno" = "Anti-depressant phenotype of first cohort",
                          "p1_cluster" = "Ancestry cluster of first cohort",
                          "p1_version" = "Version of first cohort",
                          "p2_cohort" = "Name of second cohort",
                          "p2_pheno" = "Anti-depressant phenotype of second cohort",
                          "p2_cluster" = "Ancestry cluster of second cohort",
                          "p2_version" = "Version of second cohort",
                          "rg" = "Genetic correlation",
                          "se" = "Standard error of genetic correlation",
                          "ci" = "95% confidence interval of genetic correlation",
                          "z" = "Z-score of genetic correlation",
                          "p" = "p-value of genetic correlation", 
                          "h2_obs" = "Observed scale heritability for second cohort",
                          "h2_obs_se" = "Standard error of observed scale heritability for second cohort",
                          "h2_int" = "Single-trait Linkage Disequilibrium Score regression intercept for second cohort",
                          "h2_int_se" = "Standard error of single-trait Linkage Disequilibrium Score regression intercept for second cohort",
                          "gcov_int" = "Cross-trait Linkage Disequilibrium Score regression intercept",
                          "gcov_int_se" = "Standard error of cross-trait Linkage Disequilibrium Score regression intercept"
)


source(here::here("manuscript/scripts/supplementary_tables_excell_create_cols_meta_FUN.R"))

create_cols_meta(
  file_name = file_name,
  table_variable_name = ldsc_datasets_keep,
  colname_descriptions = colname_descriptions
)

