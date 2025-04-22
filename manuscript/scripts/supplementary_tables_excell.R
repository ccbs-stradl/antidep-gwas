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

# ---------------------------------------------
# Function to update the counter of the supplementary table eg. S1, S2, S3
update_table_index <- function(table_index){
  # check table_index is numeric
  if(!is.numeric(table_index)){
    stop("table_index must be numeric")
  }
  table_index <- table_index + 1
  return(table_index)
}

# ---------------------------------------------
# function to get file names from a regular expression string
# function takes:
# a string for the path where csv/tsv files are located
# a string to match file name
# function returns the .csv or .tsv file names, removes .cols but checks they exist
get_main_file_names <- function(path, file_regex){
  files <- list.files(path, full.names = TRUE)
  files <- str_subset(files, file_regex)
  
  # If there are no files stop with error
  if(length(files) == 0){
    stop(paste0("No files found at: ", path, " with regex: ", regex))
  }
  
  # If files are found, check if they all end with .csv, .tsv or .cols
  # If they are not csv or tsv, stop with error
  if(!any(grepl("\\.csv$", files)) & !any(grepl("\\.tsv$", files)) & !any(grepl("\\.col$", files))){
    stop(paste0("No csv or tsv files found at: ", path, " with regex: ", regex))
  }
  
  # Check each .csv or .tsv in files has a corresponding .cols file
  # If not, stop with error
  for (file in files) {
    if (grepl("\\.csv$", file) | grepl("\\.tsv$", file)) {
      cols_file <- paste0(file, ".cols")
      if (!cols_file %in% files) {
        stop(paste0("No corresponding .cols file found for: ", file))
      }
    }
  }
  
  # Subset to only return .tsv or .csv files
  files <- str_subset(files, "\\.csv$|\\.tsv$")

  return(files)
}

# ---------------------------------------------
# Function to read in the results tables and sidecar .cols files
read_results <- function(full_path){
  # Read in results tables
  main <- fread(full_path)
  
  # Read in sidecar .cols files
  meta <- fread(paste0(full_path, ".cols"))
  
  # Combine into a list
  results <- list(main = main, meta = meta)
  
  return(results)
}

# ---------------------------------------------
# function that creates a list of results tables
# and names them by the sheet name

tables_list <- function(full_path_vector, sheet_names){
  
  # Check if full_path_vector is a character vector
  if(!is.character(full_path_vector)){
    stop("Input must be a character vector")
  }
  
  # Read in results tables and return a list
  # each item in the list is a list with 1) main results and 2) meta data
  results <- lapply(full_path_vector, read_results)
  
  # rename results with the sheet name (based on the file name
  names(results) <- sheet_names

  return(results)
}

# ---------------------------------------------
# function that creates an excel spreadshseet 
# it inserts a readme in the first sheet
# and inserts tables in list of tables for all subsequent sheets
# Create an excel spreadsheet with a new sheet for each file
make_excel <- function(results, file_name, table_index, letter_index, legend_title, 
                       legend_text_prefix, legend_text_sections,
                        cell_title_width, cell_title_height){
  wb <- createWorkbook()
  
  # Add a readme for the first sheet
  add_readme(results, wb, table_index, letter_index, legend_title, 
             legend_text_prefix, legend_text_sections,
             cell_title_width, cell_title_height)
  
  # Add each table to a new sheet
  for (i in seq_along(results)) {
    addWorksheet(wb, sheetName = names(results)[i])
    writeData(wb, sheet = names(results)[i], results[[i]]$main, startRow = 1, colNames = TRUE)
  }
  
  saveWorkbook(wb, file_name, overwrite = TRUE)
}

# ---------------------------------------------
# function that creates README for first sheet in excel spreadsheet
add_readme <- function(results, wb, table_index, letter_index, legend_title, 
                       legend_text_prefix, legend_text_sections,
                       cell_title_width, cell_title_height){
  # Create a readme sheet
  addWorksheet(wb, sheetName = "README")
  
  # Update legend_title
  legend_title <- as.character(glue("Table S{table_index}. {legend_title}"))
  
  # Write the table legend title
  writeData(wb, sheet = "README", legend_title, startRow = 1, startCol = 1)
  
  # Create legend_text, using appropriate grammar depending on number of sections
  number_of_sections <- length(legend_text_sections)
  if(number_of_sections == 1){
    legend_text <- paste0(legend_text_prefix, 
                          "(", letter_index[1], ") ",
                          legend_text_sections[1],
                          ".")
  } else if(number_of_sections == 2){
    legend_text <- paste0(legend_text_prefix, 
                          "(", letter_index[1], ") ",
                          legend_text_sections[1], 
                          " and (", letter_index[2], ") ",
                          legend_text_sections[2],
                          ".")
  } else if(number_of_sections > 2){
    legend_text <- paste0(legend_text_prefix, 
                          paste0("(", letter_index[1:(number_of_sections-1)], ") ",
                                 legend_text_sections[1:(number_of_sections)-1], collapse = ", "),
                          " and (", letter_index[number_of_sections], ") ",
                          legend_text_sections[number_of_sections],
                          ".")
  }
  
  # Write the table legend text
  writeData(wb, sheet = "README", legend_text, startRow = 2, startCol = 1)
  
  # Make first col width wider
  setColWidths(wb, sheet = "README", cols = 1, widths = cell_title_width)
  
  # Make first row width longer
  setRowHeights(wb, sheet = "README", rows = 1, heights = cell_title_height)
  
  # Add text wrapping style to the first cell (table legend)
  wrap_style <- createStyle(wrapText = TRUE)
  # use stack = TRUE to ensure multiple styles can be added, else it will overwrite
  addStyle(wb, sheet = "README", style = wrap_style, rows = 1, cols = 1, stack = TRUE) 
  
  # Add bold font style to first cell (table legend)
  bold_style <- createStyle(textDecoration = "bold")
  addStyle(wb, sheet = "README", style = bold_style, rows = 1, cols = 1, stack = TRUE)
  
  # Each item in 'results' contains a 'meta' data frame describing column names 
  # (from .cols sidecar files).
  # Extract that meta data and combine them into a single data frame:
  # add a column for which sheet the column name is from ie. results name
  col_name_descriptions <- do.call(rbind, lapply(1:length(results), function(i) {
    results[[i]]$meta %>%
      mutate(sheet_name = names(results[i])) %>%
      relocate(sheet_name)
  }))
  
  
  
  writeData(wb, sheet = "README", col_name_descriptions, startRow = 3, startCol = 1)
  
  # make row with column name and descriptions bold
  addStyle(wb, sheet = "README", style = bold_style, rows = 3, cols = 1:2, stack = TRUE)
  
}

###############################################
#### Create table #############################
###############################################
create_table <- function(paths, regex, 
                        sheet_names,
                        excel_file_name, 
                        table_index, 
                        legend_title, 
                        legend_text_prefix, legend_text_sections,
                        cell_title_width, cell_title_height){

  # Check inputs for paths, regex and sheet names are all the same length
  if(length(paths) != length(regex) | length(paths) != length(sheet_names)){
    stop("paths, regex and sheet_names must all be the same length")
  }
  
  # Get supplementary table character index, eg. A, B, C etc.
  letter_index <- LETTERS[1:length(sheet_names)]
  
  # Create sheet names with suffix for table index and A, B, C etc.
  sheet_names <- paste0(glue("Table S{paste0(table_index, letter_index)} "), sheet_names)
  
  # Excel sheet names must be less than 31 characters
  # If they are, stop with error, stating which sheet name is more than 31 chars
  if(any(nchar(sheet_names) > 31)){
    stop(paste0("Excel sheet names must be less than 31 characters. ",
                "The following sheet names are too long: ", 
                paste(sheet_names[nchar(sheet_names) > 31], collapse = ", ")))
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
  make_excel(results, excel_file_name, table_index, letter_index, legend_title, 
             legend_text_prefix, legend_text_sections,
              cell_title_width, cell_title_height)
  
}

###############################################
#### Use main() to run create_table functions #
###############################################
main <- function(){
  # Set table_index to 1
  table_index <- update_table_index(0)
  
  # Create the first supplementary table for the clumps and fine mapping results
  create_table(paths = rep("manuscript/tables", 3),
               regex = c("clumps_fixed_antidep-2501.clumps",
                         "susiex_significant_summary",
                         "susiex_significant_cs"),
               sheet_names = c("clumps fixed",
                               "SuSiEx summary",
                               "SuSiEx credible sets"),
               excel_file_name = here::here(glue("manuscript/tables/S{table_index}_clumps_finemap.xlsx")),
               table_index,
               legend_title = "Clumping and fine mapping results for the meta-analysis of the antidepressant GWAS.",
               legend_text_prefix = "Results are divided into ",
               legend_text_sections = c("fixed clumping results across all ancestries and antidepressant phenotypes",
                                         "significant SuSiEx summary statistics",
                                         "significant SuSiEx credible sets"),
               cell_title_width = 30,
               cell_title_height = 50)
  
  # Set table index to 2
  table_index <- update_table_index(table_index)
}

# ---------------------------------------------
# Run the main function
main()

