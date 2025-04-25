# Load libraries
library(data.table)
library(glue)
library(here)
library(openxlsx)
library(dplyr)
library(tidyr)
library(stringr)
library(tools)
library(readr)

# Load in functions (move below functions to this script when finished)
source(here("manuscript/scripts/supplementary_tables_excel_functions.R"))
# ---------------------------------------------
test_data_location <- here("manuscript/tests/excel_pipeline/test_data/csv/")
test_data_csv_names <- paste0("test_data_sheet_", 1:3, ".csv")
data_nrows <- 200
set.seed(3205722)

# create csv tables
for (i in 1:3) {
  # Create a data frame with random data
  df <- data.frame(
    CHR = sample(1:22, data_nrows, replace = TRUE),
    BP_START = sample(1:1000000, data_nrows, replace = TRUE),
    BP_END = sample(1:1000000, data_nrows, replace = TRUE),
    CS_ID = paste0("CS", sample(1:100, data_nrows, replace = TRUE)),
    SNP = paste0("rs", sample(1:10000, data_nrows, replace = TRUE)),
    BP = sample(1:1000000, data_nrows, replace = TRUE),
    CS_PIP = runif(data_nrows),
    OVRL_PIP = runif(data_nrows),
    p_val = c(runif(data_nrows/2, 0, 0.049), 
              runif(data_nrows/2, 0.05, 1) ) %>% sample()
  )
  
  # Write the data frame to a CSV file
  write.csv(df, here::here(glue('{test_data_location}{test_data_csv_names[i]}')), row.names = FALSE)
}


# create sidecar .cols tables
colname_descriptions <- c(
  "CHR" = "Chromosome",
  "BP_START" = "Start position of the finemapping region",
  "BP_END" = "End position of the finemapping region",
  "CS_ID" = "Credible Set ID",
  "SNP" = "SNP identifier",
  "BP" = "The base pair coordinate of the SNP in Gr37",
  "CS_PIP" = "Posterior inclusion probability (PIP) of the SNP in the credible set.",
  "OVRL_PIP" = "Posterior inclusion probability (PIP) of the SNP in the entire region.",
  "p_val" = "dummy p-value column for testing purposes"
)

colname_descriptions_table <- tibble(column = names(colname_descriptions), description = colname_descriptions)

lapply(test_data_csv_names, function(x){
write_tsv(colname_descriptions_table, here::here(glue('{test_data_location}{x}.cols')))
})

# ---------------------------------------------
# Create main function to run create_table for each table
# use test data generated above
main <- function(){
  # Set table_index to 1
  table_index <- update_table_index(0)
  
  # Create the first supplementary table for the clumps and fine mapping results
  create_table(paths = rep(test_data_location, 3),
               regex = c("test_data_sheet_1",
                         "test_data_sheet_2",
                         "test_data_sheet_3"),
               sheet_names = c("test sheet",
                               "test sheet signif",
                               "test sheet"),
               excel_file_name = here::here(glue("manuscript/tests/excel_pipeline/test_data/xlsx/S{table_index}_test_excel.xlsx")),
               table_index,
               legend_title = "Here is a legend title that should be numbered and in bold.",
               legend_text_prefix = "This is the start of the legend text",
               legend_text_sections = c("this is part A for the first sheet",
                                        "this is part B for the first sheet",
                                        "this is part C for the first sheet"),
               cell_title_width = 30,
               cell_title_height = 50,
               bold_significant = "p_val")
}

# ---------------------------------------------
# Run the main function
main()
