# Create supplementary tables in excel spreadsheets for the following:
# Clumping / finemapping
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
main <- function(){
  # Set table_index to 1
  table_index <- update_table_index(0)
  
  # Create the first supplementary table for the clumps and fine mapping results
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
  
  # Set table index to 2
  table_index <- update_table_index(table_index)
  
  # Create the second supplementary table for the gene mapping results
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
  
  # Set table index to 3
  table_index <- update_table_index(table_index)
  
}

# ---------------------------------------------
# Run the main function
main()

