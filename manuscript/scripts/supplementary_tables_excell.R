# Create supplementary tables in excel spreadsheets for the following:
# Datasets and meta-analysis sample sizes
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

  # Create the supplementary table for the LDSC/popcorn results
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
  
  # Update table index
  table_index <- update_table_index(table_index)
  
  create_table(
    paths = "manuscript/tables",
    regex = "manuscript/tables/across_methods_and_mdd_gwas_antidep_subset",
    sheet_names = "MDD GWAS",
    excel_file_name = here::here(glue("manuscript/tables/S{table_index}_mdd_comparison.xlsx")),
    table_index,
    legend_title = "Comparison of results to major depression GWAS.",
    legend_text_prefix = "",
    legend_text_sections = c(
      "Genes identified in antidepressant GWAS meta-analysis and MDD GWAS (Adams et al. 2025). Rows are TRUE if the gene was identified with the method described in the column header."
    ),
    cell_title_width = 39,
    cell_title_height = 49
  )
  
}

# ---------------------------------------------
# Run the main function
main()

# ---------------------------------------------
# Run lintr and styler
# lint(here("manuscript/scripts/supplementary_tables_excell.R"))
# style_file(here("manuscript/scripts/supplementary_tables_excell.R"))
