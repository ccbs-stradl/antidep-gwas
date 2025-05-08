# Read in drug targetor results
# Create a .cols file for each drug targetor results file
# Save the drug targetor results as csv files ready for the excel pipeline
# --------------------------------------
# Load required libraries
library(readxl)

# Load function that creates .cols file
# The function create_cols_meta(file_name, table_variable_name, colname_descriptions)
source(here::here("manuscript/scripts/supplementary_tables_excell_create_cols_meta_FUN.R"))

# --------------------------------------
# Define colname descriptions for .col meta sidecar
colname_descriptions_drugsAllp <- c(
  "COMP_P" = "P-value",
  "CODE" = "Drug class",
  "NAME" = "Drug name",
  "N" = "Number of genes",
  "AUC" = "AUC",
  "N_FOUND" = "Number of genes found",
  "q_valueBH" = "FDR q-value (BH)",
  "q_valueBY" = "FDR q-value (BY)",
  "p_valueBF" = "Bonferroni p-value"
)
colname_descriptions_drugsAllpEUR <- c(
  "COMP_P" = "P-value",
  "NAME" = "Drug name",
  "NGENES" = "Number of genes",
  "q_valueBH" = "FDR q-value (BH)",
  "q_valueBY" = "FDR q-value (BY)",
  "p_valueBF" = "Bonferroni p-value"
)


############ FUNS #####################
# Read in drug targetor results
read_results <- function(xlsx_path) {
  # Read the xlsx file
  results <- read_excel(xlsx_path, sheet = 1)

  # Return the combined results
  results
}

# Save a .cols file for each drug targetor results file
make_cols <- function(results, colname_descriptions, file_name) {
  create_cols_meta(
    file_name = file_name,
    table_variable_name = results,
    colname_descriptions = colname_descriptions
  )
}

# Save the drug targetor results as csv files ready for the excel pipeline
save_results <- function(file_name, results) {
  # Save the results as a csv file
  write.csv(
    results,
    file = here::here(file_name),
    row.names = FALSE,
    quote = FALSE
  )
}

# --------------------------------------
main <- function(xlsx_paths, colname_descriptions, file_names) {
  # Read in drug targetor results
  results_list <- lapply(xlsx_paths, read_results)

  # Save a .cols file for each drug targetor results file
  Map(make_cols, results_list, colname_descriptions, file_names)

  # Save the drug targetor results as csv files
  Map(save_results, file_names, results_list)
}

########################################

main(
  list(
    "manuscript/tables/antidep-2501-fixed-N06A-EUR.pathway.output_drugsAllp.drugclass_with4.selp.xls",
    "manuscript/tables/antidep-2501-fixed-N06A-EUR.pathway.output_drugsAllpEUR.pathway.output_drugsAllp.xls"
  ),
  list(
    colname_descriptions_drugsAllp,
    colname_descriptions_drugsAllpEUR
  ),
  list(
    "manuscript/tables/antidep-2501-fixed-N06A-EUR.pathway.output_drugsAllp.drugclass_with4.selp.csv",
    "manuscript/tables/antidep-2501-fixed-N06A-EUR.pathway.output_drugsAllp.pathway.output_drugsAllp.csv"
  )
)
