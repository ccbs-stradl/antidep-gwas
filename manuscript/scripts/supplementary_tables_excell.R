# Create supplementary tables in excell spreadsheets for the following:
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

# Load in functions
source(here("manuscript/scripts/supplementary_tables_excell_functions.R"))

# Load in descriptions for the column names that will be inserted into the README sheets
source(here("manuscript/scripts/supplementary_tables_excell_colname_descriptions.R"))

###############################################
#### Meta-analysis and fine mapping table #####
###############################################
# Update table index
update_table_index <- function(table_index){
  # check table_index is numeric
  if(!is.numeric(table_index)){
    stop("table_index must be numeric")
  }
  table_index <- table_index + 1
  return(table_index)
}

# Create excell table with 
# readme as the first sheet with table legend
# clumping tables
# fine mapping tables for significant results

# function to get file names from a regular expression string
# function takes:
# a string for the path where csv/tsv files are located
# a string to match file name
get_file_names <- function(path, file_regex){
  files <- list.files(path, full.names = TRUE)
  str_subset(files, file_regex)
}

# function to read in table of results (csv or tsv)
# function takes:
# a string for the full path of csv/tsv file
read_file <- function(path, file){
  fread(path)
}

# function that reads in tables and stores them as a list
# these tables can then be inserted into the excell sheets
# function takes:
# a vector of strings for the full path of csv/tsv files
tables_list <- function(full_path_vector){
  tables_list <- lapply(full_path_vector, read_file)
  # get basename without file extension
  sheet_name <- file_path_sans_ext(basename(full_path_vector))
  # sheet name is max 31 characters, if any are above this then truncate them
  if(any(nchar(sheet_name) > 31)){
    sheet_name <- str_trunc(sheet_name, 31)
  }
  # rename items in list by this name, which will later be the sheet name
  names(tables_list) <- sheet_name
  return(tables_list)
}

# function that gets column names from data frames in the sheets and returns
# a dataframe with two columns:
# first column are column names in sheets, second column is a description.
col_name_descriptions <- function(tables_list, colname_descriptions){
  
  # Get all unique column names from the list of data frames
  # Note this assumes the descriptions are the same for duplicate column names
  # the user needs to check this when they create the descriptions variable
  col_names <- unique(unlist(lapply(tables_list, colnames)))
  
  # Create a data frame with column names and match with descriptions
  col_name_df <- data.frame(
    `column name` = col_names,
    `description` = sapply(col_names, function(name) {
      if (name %in% names(colname_descriptions)) {
        colname_descriptions[[name]] # get the description that matches the col name
      } else {
        ""
      }
    }),
    stringsAsFactors = FALSE
  )
  
  return(col_name_df)
}
  
# function that creates README for first sheet in excell spreadsheet
add_readme <- function(tables_list, wb, table_index, legend_title, legend_text, 
                       cell_title_width, cell_title_height, colname_descriptions){
  # Create a readme sheet
  addWorksheet(wb, sheetName = "README")
  
  # Write the table legend title
  writeData(wb, sheet = "README", as.character(legend_title), startRow = 1, startCol = 1)
  
  # Write the table legend text
  writeData(wb, sheet = "README", as.character(legend_text), startRow = 2, startCol = 1)
  
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
  
  # Explain what the column names in each of the sheets mean
  # for column names that are the same across sheets collapse these into one table
  col_name_descriptions <- col_name_descriptions(tables_list, colname_descriptions)
  
  writeData(wb, sheet = "README", col_name_descriptions, startRow = 3, startCol = 1)
  
  # make row with column name and descriptions bold
  addStyle(wb, sheet = "README", style = bold_style, rows = 3, cols = 1:2, stack = TRUE)
  
  # # List the Supplementary Table Names:
  # supplementary_tables <- paste("Supplementary Table ", table_index, LETTERS[1:length(tables_list)])
  # 
  # writeData(wb, sheet = "README", supplementary_tables, startRow = 3, startCol = 1)
  # 
  # # List the names of the tables_list next to the table number
  # writeData(wb, sheet = "README", names(tables_list), startRow = 3, startCol = 2)
  
}

# function that creates an excell spreadshseet 
# it inserts a readme in the first sheet
# and inserts tables in list of tables for all subsequent sheets
# Create an excell spreadsheet with a new sheet for each file
make_excell <- function(tables_list, file_name, table_index, legend_title, legend_text,
                        cell_title_width, cell_title_height, colname_descriptions){
  wb <- createWorkbook()
  
  # Add a readme for the first sheet
  add_readme(tables_list, wb, table_index, legend_title, legend_text,
             cell_title_width, cell_title_height, colname_descriptions)
  
  # Add each table to a new sheet
  for (i in seq_along(tables_list)) {
    addWorksheet(wb, sheetName = names(tables_list)[i])
    writeData(wb, sheet = names(tables_list)[i], tables_list[[i]], startRow = 1, colNames = TRUE)
  }
  
  saveWorkbook(wb, file_name, overwrite = TRUE)
}

# Create the first supplementary table for the clumps and fine mapping results
create_meta_finemapping <- function(excell_file_name, table_index, legend_title, legend_text,
                                    cell_title_width, cell_title_height, colname_descriptions){

  # Read in full paths for where clumps and fine mapping results are stored
  clumps_path <- "results/meta/antidep-2501"
  clumps_regex <- "clumps"
  clumps <- get_file_names(clumps_path,
                           clumps_regex)
  if(length(clumps) == 0){
    stop(paste0("No files found at: ", clumps_path, " with regex: ", clumps_regex))
  }
  finemapping_path <- "manuscript/tables"
  finemapping_regex <- "susiex_significant"
  finemapping <- get_file_names(finemapping_path,
                               finemapping_regex)
  if(length(finemapping) == 0){
    stop(paste0("No files found at: ", finemapping_path, " with regex: ", finemapping_regex))
  }
  
  # Read these tables and store as a list of data frames
  # the names of this list are the file names
  results_list <- tables_list(c(clumps, finemapping))

  # Make excell spreadsheet
  make_excell(results_list, excell_file_name, table_index, legend_title, legend_text,
              cell_title_width, cell_title_height, colname_descriptions)
  
}

main <- function(){
  # Set table_index to 1
  table_index <- update_table_index(0)
  
  # Create the first supplementary table for the clumps and fine mapping results
  create_meta_finemapping(here::here(glue("manuscript/tables/S{table_index}_clumps_finemap.xlsx")),
                          table_index,
                          glue("Table S{table_index}. Clumping and fine mapping results for the meta-analysis of the antidepressant GWAS."),
                          glue("Clumping results are divided into"),
                          cell_title_width = 30,
                          cell_title_height = 50,
                          colname_descriptions = meta_finemapping_colname_descriptions)
  
  # Set table index to 2
  table_index <- update_table_index(table_index)
}

# Run the main function
main()

