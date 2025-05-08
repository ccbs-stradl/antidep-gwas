# Create supplementary tables in excel spreadsheets for the following:
# Clumping and finemapping
# Gene mapping mBAT, cross method
# GWAS catalog
# LDSC/popcorn: between gwas, between meta, with external traits, reference table
# SMR table (already collated)
# Full drug targetor
#
# --------------------------------------------
# Load libraries
library(data.table)
library(glue)
library(here)
library(openxlsx)
library(dplyr)
library(tidyr)
library(stringr)
library(tools)

# Load in functions (move below functions to this script when finished)
source(here("manuscript/scripts/supplementary_tables_excel_functions.R"))

###############################################
#### Use main() to run create_table functions #
###############################################
main <- function() {
  # Set table_index to 1
  table_index <- update_table_index(0)

  # Create the supplementary for datasets and meta-analysis samples
  create_table(
    paths = rep("manuscript/tables", 3),
    regex = c(
      "datasets_gwas",
      "datasets_meta_gwas",
      "datasets_meta_totals"
    ),
    sheet_names = c(
      "GWAS datasets",
      "Meta datasets",
      "Meta totals"
    ),
    excel_file_name = here::here(
      glue("manuscript/tables/S{table_index}_datasets.xlsx")
    ),
    table_index,
    legend_title = "GWAS datasets, meta-analysis inputs, and total sample sizes",
    legend_text_prefix = "Datasets and sample sizes are listed for ",
    legend_text_sections = c(
      "per cohort and ancestry input GWAS",
      "datasets included in each meta-analysis",
      "subtotal and total sample sizes for each meta-analyis"
          ),
    cell_title_width = 30,
    cell_title_height = 50
  )

  # Update table index
    table_index <- update_table_index(table_index)

  # Create the supplementary table for the clumps and fine mapping results
  create_table(paths = rep("manuscript/tables", 4),
               regex = c("clumps_fixed_antidep-2501.clumps",
                         "clumps_mrmega_antidep-2501.clumps",
                         "susiex_significant_summary",
                         "susiex_significant_cs"),
               sheet_names = c("clumps fixed",
                               "clumps MR-MEGA",
                               "SuSiEx summary",
                               "SuSiEx credible sets"),
               excel_file_name = here::here(glue("manuscript/tables/S{table_index}_clumps_finemap.xlsx")),
               table_index,
               legend_title = "Clumping and fine mapping results for the meta-analysis of the antidepressant GWAS.",
               legend_text_prefix = "Results are divided into ",
               legend_text_sections = c("fixed clumping results across all ancestries and antidepressant phenotypes",
                                        "MR-MEGA clumping results",
                                         "significant SuSiEx summary statistics",
                                         "significant SuSiEx credible sets"),
               cell_title_width = 30,
               cell_title_height = 50)

  # Update table index
    table_index <- update_table_index(table_index)

  # Create the supplementary table for the gene mapping results
  create_table(paths = "manuscript/tables",
               regex = c("mBAT-combo.csv"),
               sheet_names = c("mBAT-combo"),
               excel_file_name = here::here(glue("manuscript/tables/S{table_index}_gene_mapping.xlsx")),
               table_index,
               legend_title = "Positional mapping results for EUR, AFR and SAS fixed meta-analyses of the antidepressant GWAS (N06A, N06AA and N06AB).",
               legend_text_prefix = "",
               legend_text_sections = c("Results shown for Bonferroni corrected mBAT-combo p-value < 0.05"),
               cell_title_width = 39,
               cell_title_height = 49)
  
  # Update table index
    table_index <- update_table_index(table_index)

  # Create the second supplementary table for the gene mapping results
  create_table(
    paths = "manuscript/tables",
    regex = c("mBAT-combo.csv"),
    sheet_names = c("mBAT-combo"),
    excel_file_name = here::here(glue("manuscript/tables/S{table_index}_gene_mapping.xlsx")),
    table_index,
    legend_title = "Positional mapping results for EUR, AFR and SAS fixed meta-analyses of the antidepressant GWAS (N06A, N06AA and N06AB).", # nolint
    legend_text_prefix = "",
    legend_text_sections = c("Results shown for Bonferroni corrected mBAT-combo p-value < 0.05"),
    cell_title_width = 39,
    cell_title_height = 49
  )

  # Update table index
    table_index <- update_table_index(table_index)

  # Create the second supplementary table for the GWAS catalog
  create_table(
    paths = rep("manuscript/tables", 11),
    regex = c(
      "gwascat_fixed_table_N06A-AFR.csv",
      "gwascat_fixed_table_N06A-EAS.csv",
      "gwascat_fixed_table_N06A-EUR.csv",
      "gwascat_fixed_table_N06A-SAS.csv",
      "gwascat_fixed_table_N06AA-AFR.csv",
      "gwascat_fixed_table_N06AA-EUR.csv",
      "gwascat_fixed_table_N06AB-AFR.csv",
      "gwascat_fixed_table_N06AB-EUR.csv",
      "gwascat_fixed_table_N06AB-MID.csv",
      "gwascat_fixed_table_N06AB-SAS.csv",
      "gwascat_mrmega_antidep-2501.csv"
    ),
    sheet_names = c(
      "N06A-AFR",
      "N06A-EAS",
      "N06A-EUR",
      "N06A-SAS",
      "N06AA-AFR",
      "N06AA-EUR",
      "N06AB-AFR",
      "N06AB-EUR",
      "N06AB-MID",
      "N06AB-SAS",
      "MRMEGA"
    ),
    excel_file_name = here::here(glue("manuscript/tables/S{table_index}_gwas_cat.xlsx")),
    table_index,
    legend_title = "NHGRI-EBI GWAS catalogue lookup for the antidepressant meta-analysis GWAS",
    legend_text_prefix = "Results split by ancestries and phenotypes for fixed and MRMEGA meta-analysis GWAS.",
    legend_text_sections = c(
      "N06A-AFR",
      "N06A-EAS",
      "N06A-EUR",
      "N06A-SAS",
      "N06AA-AFR",
      "N06AA-EUR",
      "N06AB-AFR",
      "N06AB-EUR",
      "N06AB-MID",
      "N06AB-SAS",
      "MRMEGA"
    ),
    cell_title_width = 39,
    cell_title_height = 49
  )

  # Update table index
    table_index <- update_table_index(table_index)

  # Create the fourth supplementary table for the LDSC/popcorn results
  create_table(
    paths = rep("manuscript/tables", 6),
    regex = c(
      "rg_ldsc_gwas.csv",
      "rg_ldsc_meta.csv",
      "rg_ldsc_meta_external",
      "rg_ldsc_external_references",
      "rg_popcorn_gwas",
      "rg_popcorn_meta"
    ),
    sheet_names = c(
      "LDSC gwas",
      "LDSC meta",
      "LDSC meta external",
      "GWAS references",
      "Popcorn gwas",
      "Popcorn meta"
    ),
    excel_file_name = here::here(glue("manuscript/tables/S{table_index}_ldsc_popcorn.xlsx")),
    table_index,
    legend_title = "Cross-ancestry genetic correlations (LDSC and popcorn) between cohorts, meta-analyses and external traits", # nolint
    legend_text_prefix = "Results are divided into ",
    legend_text_sections = c(
      "LDSC of input cohort GWASs",
      "LDSC of meta-analysed GWASs",
      "LDSC of meta-analysed GWASs with external traits",
      "References for external traits used in the LDSC analysis",
      "Popcorn of input cohort GWASs",
      "Popcorn of meta-analysed GWASs"
    ),
    cell_title_width = 47,
    cell_title_height = 32
  )
}

# ---------------------------------------------
# Run the main function
main()

# ---------------------------------------------
# Run lintr and styler
# lint(here("manuscript/scripts/supplementary_tables_excell.R"))
# style_file(here("manuscript/scripts/supplementary_tables_excell.R"))
