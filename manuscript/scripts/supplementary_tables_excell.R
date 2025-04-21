# Create supplementary tables in excell spreadsheets for the following:
# Clumping / finemapping
# Gene mapping mBAT, cross method
# GWAS catalog
# LDSC/popcorn: between gwas, between meta, with external traits, reference table
# SMR table (already collated)
# Full drug targetor
#
# --------------------------------------------
# Create sym links to results
# meta analysis results:
system("ln -s /Volumes/GenScotDepression/data/AMBER/antidep-gwas/meta results") 
# fine mapping results:
system("ln -s /Volumes/GenScotDepression/data/AMBER/antidep-gwas/meta results") 

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

# function that creates README for first sheet in excell spreadsheet
add_readme <- function(tables_list, wb, sup_table_num){
  # Create a readme sheet
  addWorksheet(wb, sheetName = "README")
  
  # Write the readme text
  writeData(wb, sheet = "README", "Table Legend Placeholder")
  
  # Make first col width wider
  setColWidths(wb, sheet = "README", cols = 1, widths = 30)
  setColWidths(wb, sheet = "README", cols = 1, widths = cell_title_width)
  
  # Make first row width longer
  setRowHeights(wb, sheet = "README", rows = 1, heights = 127)
  setRowHeights(wb, sheet = "README", rows = 1, heights = cell_title_height)
  
  # Add text wrapping style to the first cell
  wrap_style <- createStyle(wrapText = TRUE)
  addStyle(wb, sheet = "README", style = wrap_style, rows = 1, cols = 1)
  
  # List the Supplementary Table Names:
  supplementary_tables <- paste("Supplementary Table ", sup_table_num, LETTERS[1:length(tables_list)])
  
  writeData(wb, sheet = "README", supplementary_tables, startRow = 2, startCol = 1)
  
  # List the names of the tables_list next to the table number
  writeData(wb, sheet = "README", names(tables_list), startRow = 2, startCol = 2)
  
}

# function that creates an excell spreadshseet 
# it inserts a readme in the first sheet
# and inserts tables in list of tables for all subsequent sheets
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
  
  saveWorkbook(wb, file_name, overwrite = TRUE)
}


create_meta_finemapping <- function(excell_file_name, sup_table_num){

  # Read in full paths for where clumps and fine mapping results are stored
  clumps <- get_file_names("results/meta/antidep-2501",
                 "clumps")
  finemapping <- get_file_names("manuscript/tables",
                 "susiex_significant")
  
  # Read these tables and store as a list of data frames
  # the names of this list are the file names
  results_list <- tables_list(c(clumps, finemapping))

  # Make excell spreadsheet
  make_excell(results_list, excell_file_name, sup_table_num)
  
}

main <- function(){
  create_meta_finemapping(here::here("manuscript/tables/Supplementary_Table_X_clumps_finemap.xlsx"),
                          "XX") # supplementary table number placeholder
  # Set table_index to 1
  table_index <- update_table_index(0)
  create_meta_finemapping(here::here(glue("manuscript/tables/S{table_index}_clumps_finemap.xlsx")),
                          table_index,
                          glue("Table S{table_index}. Clumping and fine mapping results for the meta-analysis of the antidepressant GWAS."),
                          glue("Clumping results are divided into"),
                          cell_title_width = 30,
                          cell_title_height = 50,
  # Set table index to 2
  table_index <- update_table_index(table_index)
}


main()
