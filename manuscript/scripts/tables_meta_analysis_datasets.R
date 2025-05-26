# Compile inputs for supplementary table of GWAS datasets and metasets
# Similar inputs as `meta_analysis_cohort_summary_table` but in
# long format with additional version/build information and
# column descriptions

library(readr)
library(dplyr)
library(tidyr)
library(stringr)
library(tools)

# directory with meta-analysis outputs
meta_dir <- here::here("results/meta/antidep-2501")

# get list of datasets in each meta analysis
datasets_files <- list.files(meta_dir, pattern = ".+\\.csv", full.names = TRUE)
# name list with basename w/o extension
names(datasets_files) <- file_path_sans_ext(basename(datasets_files))

# read in all datasets
datasets_list <- lapply(datasets_files, read_csv)
datasets <- bind_rows(datasets_list, .id = "metaanalysis")

# get input GWAS as list of distinct datasets
datasets_gwas <- datasets |>
  select(-metaanalysis) |>
  distinct() |>
  relocate(version, build, dataset, .after = last_col()) |>
  arrange(cohort, pheno, cluster)

# input dataset table column descriptions
datasets_gwas_descriptions <- c(
  cohort = "Cohort name",
  pheno = "ATC subgroup",
  cluster = "Genetic similarity (ancestry) cluster",
  cases = "Number of cases",
  controls = "Number of controls",
  neff = "Effective sample size",
  version = "Sumstats version identifier",
  build = "Sumstats original genome build",
  dataset = "Sumstats filename"
)

gwas_descriptions_table <- tibble(
  column = names(datasets_gwas_descriptions),
  description = datasets_gwas_descriptions
)

if (any(gwas_descriptions_table$column != colnames(datasets_gwas))) {
  stop(glue("Column names in datasets_gwas are not all described"))
}

write_csv(datasets_gwas, here::here("manuscript/tables/datasets_gwas.csv"))
write_csv(gwas_descriptions_table, here::here("manuscript/tables/datasets_gwas.csv.cols"))

# datasets included in each fixed and mr-mega meta-analyses
datasets_meta <- datasets |>
  rename(cohort_cluster = cluster) |>
  select(-pheno, -version, -build, -dataset) |>
  separate(
    metaanalysis,
    into = c("meta", "date", "method", "pheno", "meta_cluster"),
    remove = FALSE
  ) |>
  unite("version", meta, date, sep = "-") |>
  mutate(metaset = str_glue("{metaanalysis}.vcf.gz")) |>
  select(
    method, pheno, meta_cluster,
    cohort, cohort_cluster,
    cases, controls, neff,
    version, metaset
  ) |>
  arrange(method, pheno, meta_cluster, cohort, cohort_cluster)

datasets_meta_descriptions <- c(
  "method" = "Meta analysis method",
  "pheno" = "ATC subgroup",
  "meta_cluster" = "Meta-analysis genetic similirity (ancestry) cluster",
  "cohort" = "Cohort name",
  "cohort_cluster" = "Cohort genetic similirity (ancestry) cluster",
  "cases" = "Number of cases",
  "controls" = "Numer of controls",
  "neff" = "Effective sample size",
  "version" = "Meta-analysis version identifier",
  "metaset" = "Meta-analysis filename"
)

meta_descriptions_table <- tibble(
  column = names(datasets_meta_descriptions),
  description = datasets_meta_descriptions
)

if (any(meta_descriptions_table$column != colnames(datasets_meta))) {
  stop(str_glue("Column names in datasets_meta are not all described"))
}

write_csv(datasets_meta, here::here("manuscript/tables/datasets_meta_gwas.csv"))
write_csv(meta_descriptions_table, here::here("manuscript/tables/datasets_meta_gwas.csv.cols"))

# subtotal and total case/control/neff counts
datasets_meta_fixed <- datasets_meta |>
  filter(method == "fixed")

datasets_meta_mrmega <- datasets_meta |>
  filter(method == "mrmega")

datasets_meta_fixed_totals <- datasets_meta_fixed |>
  group_by(method, pheno, cohort_cluster, version, metaset) |>
  summarize(across(c(cases, controls, neff), sum)) |>
  ungroup() |>
  rename(cluster = cohort_cluster)

datasets_meta_mrmega_subtotals <- datasets_meta_mrmega |>
  group_by(method, pheno, cohort_cluster, version, metaset) |>
  summarize(across(c(cases, controls, neff), sum)) |>
  ungroup() |>
  rename(cluster = cohort_cluster)

datasets_meta_mrmega_totals <- datasets_meta_mrmega_subtotals |>
  group_by(method, pheno, version, metaset) |>
  summarize(across(c(cases, controls, neff), sum)) |>
  mutate(cluster = "DIV")

datasets_meta_totals <- bind_rows(
  datasets_meta_fixed_totals,
  datasets_meta_mrmega_subtotals,
  datasets_meta_mrmega_totals
) |>
  select(
    method, pheno, cluster,
    cases, controls, neff,
    version, metaset
  )

totals_description <- c(
  "method" = "Meta-analysis method",
  "pheno" = "ATC subgroup",
  "cluster" = "Genetic similirity (ancestry) cluster",
  "cases" = "Number of cases",
  "controls" = "Number of controls",
  "neff" = "Effective sample size",
  "version" = "Meta-analysis version",
  "metaset" = "Meta-analysis filename"
)

totals_descriptions_table <- tibble(
  column = names(totals_description),
  description = totals_description
)

if (any(totals_descriptions_table$column != colnames(datasets_meta_totals))) {
  stop(str_glue("Column names in datasets_meta_totals are not all described"))
}

write_csv(datasets_meta_totals, here::here("manuscript/tables/datasets_meta_totals.csv"))
write_csv(totals_descriptions_table, here::here("manuscript/tables/datasets_meta_totals.csv.cols"))
