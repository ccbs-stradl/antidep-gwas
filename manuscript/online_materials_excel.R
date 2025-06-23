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


# ---------------------------------------
# Overide functions to change "Table S1" to "Online Material 01" etc.

# function to add legend title
add_legend_title <- function(table_index, legend_title, wb) {
  # Update legend_title
  legend_title <- as.character(glue("Online Material 0{table_index}. {legend_title}"))
  
  # Write the table legend title
  writeData(wb, sheet = "README", legend_title, startRow = 1, startCol = 1)
}

create_table <- function(paths, regex,
                         sheet_names,
                         excel_file_name,
                         table_index,
                         legend_title,
                         legend_text_prefix, legend_text_sections,
                         cell_title_width, cell_title_height,
                         bold_cols = NULL,
                         bold_condition = NULL,
                         bold_threshold = NULL) {
  # Check inputs for paths, regex and sheet names are all the same length
  if (length(paths) != length(regex) | length(paths) != length(sheet_names)) {
    stop("paths, regex and sheet_names must all be the same length")
  }
  
  # Get supplementary table character index, eg. A, B, C etc.
  letter_index <- LETTERS[1:length(sheet_names)]
  
  # Create sheet names with suffix for table index and A, B, C etc.
  sheet_names <- paste0(glue("{letter_index} "), sheet_names)
  
  # Excel sheet names must be less than 31 characters
  # If they are, stop with error, stating which sheet name is more than 31 chars
  if (any(nchar(sheet_names) > 31)) {
    stop(paste0(
      "Excel sheet names must be less than 31 characters. ",
      "The following sheet names are too long: ",
      paste(sheet_names[nchar(sheet_names) > 31], collapse = ", ")
    ))
  }
  
  # Read in full paths for where results are stored
  # get .csv/.tsv and sidecar .cols file (with column name descriptions)
  files <- mapply(get_main_file_names, paths, regex) %>%
    unlist(use.names = FALSE)
  
  # Read these tables and store as a list of data frames
  # where each dataframe has an assocaited sidecar .cols file
  # the names of this list are the file names
  results <- tables_list(files, sheet_names)
  
  # Make excel spreadsheet
  make_excel(
    results, excel_file_name, table_index, letter_index, legend_title,
    legend_text_prefix, legend_text_sections,
    cell_title_width, cell_title_height, 
    bold_cols, bold_condition, bold_threshold
  )
}
# --------------------------------------

main <- function() {
  # Set table_index to 1
  table_index <- update_table_index(0)


# Create the supplementary table for the clumps and fine mapping results
create_table(
  paths = rep("manuscript/tables", 4),
  regex = c(
    "clumps_fixed_antidep-2501.clumps",
    "clumps_mrmega_antidep-2501.clumps",
    "susiex_significant_summary",
    "susiex_significant_cs"
  ),
  sheet_names = c(
    "clumps fixed",
    "clumps MR-MEGA",
    "SuSiEx summary",
    "SuSiEx credible sets"
  ),
  excel_file_name = here::here(glue("manuscript/tables/OM{table_index}_clumps_finemap.xlsx")),
  table_index,
  legend_title = "Clumping and fine mapping results for the meta-analysis of the antidepressant GWAS.",
  legend_text_prefix = "Results are divided into ",
  legend_text_sections = c(
    "fixed clumping results across all ancestries and antidepressant phenotypes",
    "MR-MEGA clumping results",
    "significant SuSiEx summary statistics",
    "significant SuSiEx credible sets"
  ),
  cell_title_width = 30,
  cell_title_height = 50
)

# Update table index
table_index <- update_table_index(table_index)

  # Create the supplementary table for the gene mapping results
  create_table(
    paths = "manuscript/tables",
    regex = c("mBAT-combo.csv"),
    sheet_names = c("mBAT-combo"),
    excel_file_name = here::here(glue("manuscript/tables/OM{table_index}_gene_mapping.xlsx")),
    table_index,
    legend_title = "Positional mapping results for EUR, AFR and SAS fixed meta-analyses of the antidepressant GWAS (N06A, N06AA and N06AB).",
    legend_text_prefix = "",
    legend_text_sections = c("Results shown for Bonferroni corrected mBAT-combo p-value < 0.05"),
    cell_title_width = 39,
    cell_title_height = 49
  )

  # Update table index
  table_index <- update_table_index(table_index)

  # Create the supplementary table for the GWAS catalog
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
    excel_file_name = here::here(glue("manuscript/tables/OM{table_index}_gwas_cat.xlsx")),
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

  # Create the supplementary table for the SMR results
  create_table(
    paths = rep("manuscript/tables/", 6),
    regex = c(
      "blood_trait_eSMR",
      "blood_trait_mSMR",
      "blood_trait_pSMR",
      "brainmeta_trait_eSMR",
      "brainmeta_trait_mSMR",
      "brainmeta_trait_sSMR"
    ),
    sheet_names = c(
      "blood eSMR",
      "blood mSMR",
      "blood pSMR",
      "brain eSMR",
      "brain mSMR",
      "brain sSMR"
    ),
    excel_file_name = here::here(glue("manuscript/tables/OM{table_index}_smr.xlsx")),
    table_index,
    legend_title = "SMR analysis in blood and brain across multi-omic data types in antidepressant GWAS meta-analysis (N06A) (EUR ancestry).",
    legend_text_prefix = "",
    legend_text_sections = c(
      "blood eSMR",
      "blood mSMR",
      "blood pSMR",
      "brain eSMR",
      "brain mSMR",
      "brain sSMR"
    ),
    cell_title_width = 39,
    cell_title_height = 49,
    bold_cols = c("p_SMR_Bonferroni", "p_HEIDI"),
    bold_condition = c("<", ">"),
    bold_threshold = c(0.05, 0.05)
  )

  # Update table index
  table_index <- update_table_index(table_index)

  # Create the supplementary table for the full drug targetor results
  create_table(
    paths = rep("manuscript/tables", 2),
    regex = c(
      "antidep-2501-fixed-N06A-EUR.pathway.output_drugsAllp.drugclass_with4.selp.csv",
      "antidep-2501-fixed-N06A-EUR.pathway.output_drugsAllp.pathway.output_drugsAllp.csv"
    ),
    sheet_names = c(
      "Gene targets",
      "Gene sets"
    ),
    excel_file_name = here::here(glue("manuscript/tables/OM{table_index}_drug_targetor.xlsx")),
    table_index,
    legend_title = "Drug targetor results for fixed meta-analysis of N06A in EUR ancestry",
    legend_text_prefix = "FDR Q-value (BH) < 0.05 are highlighted in bold. Results are split by ",
    legend_text_sections = c(
      "Enrichments of gene targets",
      "Enrichments of gene sets"
    ),
    cell_title_width = 39,
    cell_title_height = 49,
    bold_cols = c("q_valueBH"),
    bold_condition = c("<"),
    bold_threshold = c(0.05)
  )
}

# Run the main function
main()