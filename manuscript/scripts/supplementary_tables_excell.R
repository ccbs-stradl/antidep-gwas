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
library(here)
library(openxlsx)
library(dplyr)
library(tidyr)
library(stringr)

# Load in functions
source(here("manuscript/scripts/supplementary_tables_excell_functions.R"))

# Set up counter for keeping track of supplementary table number

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
list_tables <- function(full_path_vector){
  list_tables <- lapply(full_path_vector, read_file)
  names(list_tables) <- basename(full_path_vector)
  return(list_tables)
}

main <- function(){

  # Read in full paths for where clumps and fine mapping results are stored
  clumps <- get_file_names("results/meta/antidep-2501",
                 "clumps")
  finemapping <- get_file_names("manuscript/tables",
                 "susiex_significant")
  
  # Read these tables and store as a list of data frames
  # the names of this list are the file names
  results_list <- list_tables(c(clumps, finemapping))

  return(results_list)
  
}

main()
