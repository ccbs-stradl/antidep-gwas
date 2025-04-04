# Format SMR results

# Create smy link to results
# system("ln -s /Volumes/GenScotDepression/data/AMBER/antidep-gwas/maps results/")

# load libraries
library(openxlsx)
library(dplyr)
library(tidyr)
library(stringr)

# Function that reads in the results
read_results <- function(rel_path){
  # Get the path to the results
  path <- here::here(rel_path)
  files <- list.files(path, full.names = TRUE, recursive = TRUE)
  
  # Read in tables and store as a list
  smr_tables <- lapply(files, function(x) {
    read.table(x, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  })
  
  # rename smr_tables with a short descriptive name: "blood_eSMR", "brain_eSMR" etc.
  names(smr_tables) <- extract_names(files)
  
  return(smr_tables)
}

# extract name from file path
extract_names <- function(files){
  # Name tables with the tissue type and omic type
  # Extract "blood" or "brain" from directory path
  tissue <- sub(".*/smr/([^/]+)/.*", "\\1", files)  # Extracts the first directory after "smr/"
  
  # Simplify "brainmeta" to "brain"
  tissue <- gsub("brainmeta", "brain", tissue)
  
  # Extract "eSMR/mSMR/pSMR" from the filename
  smr_type <- sub(".*trait_([a-zA-Z]+)\\.merged\\.tsv", "\\1", basename(files))
  
  # Combine to get "blood_eSMR", "brain_eSMR", etc.
  short_names <- paste(tissue, smr_type, sep = "_")
  
  return(short_names)
}

# Reorder rows in results
reorder_results <- function(tables_list){
  # Check if required columns exist
  sapply(tables_list, function(df) {
    if (!all(c("p_SMR", "p_HEIDI") %in% colnames(df))) {
      stop("Required columns 'p_SMR' and 'p_HEIDI' not found in the data frame.")
    }
  })
  
  # Reorder the rows in each table by the "p" column
  for (i in seq_along(tables_list)) {
    tables_list[[i]] <- tables_list[[i]] %>% head() %>%
      arrange(desc(p_HEIDI)) %>%
      arrange(p_SMR)
  }
  
  return(tables_list)
}

# Create an excell spreadsheet with a new sheet for each file
make_excell <- function(tables_list, file_name, sup_table_num){
  wb <- createWorkbook()
  
  # Add a readme for the first sheet
  add_readme(tables_list, wb, sup_table_num)

  # Add each table to a new sheet
  for (i in seq_along(tables_list)) {
    addWorksheet(wb, sheetName = names(tables_list)[i])
    writeData(wb, sheet = names(tables_list)[i], tables_list[[i]], startRow = 1, colNames = TRUE)
  }
  
  # Make rows that pass significance threshold bold
  lapply(1:length(tables_list), function(i) {
    sheet_name <- names(tables_list)[i]
    make_bold_rows(tables_list[[i]], sheet_name, wb)
  })
  
  saveWorkbook(wb, file_name, overwrite = TRUE)
}

# Make rows that pass significance threshold bold
# significance is p_SMR < 0.05 and p_HEIDI > 0.05
make_bold_rows <- function(df, sheet_name, wb){
  # Get the row indices that pass the significance threshold
  bold_rows <- which(df$p_SMR < 0.05 & df$p_HEIDI > 0.05)
  
  # Make those rows bold
  addStyle(wb, sheet = sheet_name, rows = bold_rows + 1, cols = 1:ncol(df), 
           style = createStyle(textDecoration = "bold"), gridExpand = TRUE)

}

# Add a readme for the first sheet
add_readme <- function(tables_list, wb, sup_table_num){
  # Create a readme sheet
  addWorksheet(wb, sheetName = "README")
  
  # Write the readme text
  writeData(wb, sheet = "README", "SMR results of N06A in EUR, with multi-omic data (transcription, splicing, DNA methylation) from blood (GTEx Consortium, 2020; McRae et al., 2018) and brain tissue (Qi et al., 2022). Rows highlighted in bold indicate significant results: p_SMR < 0.05 and p_HEIDI > 0.05")
  
  # Make first col width wider
  setColWidths(wb, sheet = "README", cols = 1, widths = 30)
  
  # Make first row width longer
  setRowHeights(wb, sheet = "README", rows = 1, heights = 127)
  
  # Add text wrapping style to the first cell
  wrap_style <- createStyle(wrapText = TRUE)
  addStyle(wb, sheet = "README", style = wrap_style, rows = 1, cols = 1)
  
  # List the Supplementary Table Names:
  supplementary_tables <- paste("Supplementary Table ", sup_table_num, LETTERS[1:length(tables_list)])

  writeData(wb, sheet = "README", supplementary_tables, startRow = 2, startCol = 1)
  
  # List the tissues
  writeData(wb, sheet = "README", names(tables_list), startRow = 2, startCol = 2)
  
}


main <- function(rel_path, excell_file_name, sup_table_num){
  # Read in results
  tables_list <- read_results(rel_path)
  
  # Reorder rows in results
  tables_list_ordered <- reorder_results(tables_list)
  
  # Create an excell spreadsheet with a new sheet for each file, 
  # significant rows bold, and a readme for the first sheet
  make_excell(tables_list_ordered, excell_file_name, sup_table_num)
  
}

main("results/maps/smr",
     here::here("manuscript/tables/smr.xlsx"),
     "XX") # Supplementary Table number placeholder

