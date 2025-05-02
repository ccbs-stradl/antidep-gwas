# Make summary tables of cohorts included in the meta analyses
# MR-MEGA table: rbind all the .csv files in
# results/meta/antidep-2501/*-mrmega* (N06A, N06AA, N06AB)
# For the fixed effect tables: create a tale from
# results/meta/antidep-2501/*-fixed-* for each ancestry and
# each phenotype (N06A, N06AA, N06AB)

# --------------------------------
# Load libraries
library(data.table) # supplies fread()
library(tools) # supplies file_path_sans_ext()
library(dplyr) # supplies operator %>%
library(tidyr) # supplies pivot_wider

# --------------------------------
# Helper function to read in files
read_files <- function(dir_name, pattern) {
  # Get full paths
  paths <- list.files(dir_name,
    pattern = paste0(pattern, ".*\\.csv$"), full.names = TRUE
  )

  # Loop through each file and read it using data.table::fread
  files <- lapply(paths, fread)

  # Use file name as list item name
  # Remove dir from file name
  file_names <- basename(paths)

  # Remove file extension
  names(files) <- file_path_sans_ext(file_names)

  files
}

# --------------------------------
# add the list item name in a new column
add_name_as_col <- function(dataframe, dataframe_name) {
  # add dataframe name as new column
  dataframe$meta_analysis <- dataframe_name

  # ensure meta-analysis name is the first column
  dataframe <- dataframe %>%
    relocate(meta_analysis)

  dataframe
}

# --------------------------------
# rbind a list of dataframes
rbind_dataframes_list <- function(dataframe_list, cols_to_rm) {
  # Add item list name as a new column
  dataframe_list_named <- lapply(names(dataframe_list), function(name) {
    # Apply function with dataframe and its name
    add_name_as_col(dataframe_list[[name]], name)
  })

  # rbind list of dataframes
  summary_table <- do.call(rbind, dataframe_list_named)

  # Remove cols we don't need (dataset & build)
  summary_table <- summary_table %>%
    select(-all_of(cols_to_rm))

  summary_table
}

# --------------------------------
# save csv
save_csv <- function(table, summary_tables_path, meta_type) {
  # Check path exists, other exit with error
  if (!dir.exists(summary_tables_path)) {
    stop(paste0("Error: ", summary_tables_path, " directory does not exist."))
  }

  write.csv(table,
    paste0(summary_tables_path,
      "/meta_analysis_cohort_summary_table_",
      meta_type, ".csv"
    ),
    quote = FALSE,
    row.names = FALSE
  )
}

# --------------------------------
make_wide <- function(summary_table, first_cols) {
  pivot_wider(summary_table,
    names_from = cohort,
    values_from = c(cases, controls, neff)
  ) %>%
    rename_with(~ gsub("(.+)_(.+)", "\\2_\\1", .), -c(first_cols)) %>% # Reorder names but exclude firstCols
    select(all_of(first_cols), order(colnames(.)[-c(1:length(first_cols))]) + length(first_cols)) # Keep firstCols first and reorder the rest
}

# --------------------------------
create_table_MRMEGA <- function(summary_table_path) {
  # Read in MR-MEGA cohort csv files
  dataframes <- read_files("results/meta/antidep-2501", "mrmega")

  # Concatenate the files together (rbind) and add column for the file name
  summary_table <- rbind_dataframes_list(dataframes,
    c("dataset", "build", "version", "pheno")
  )

  # Make the table wider, with columns for cases, controls,
  # and neff for each cohort
  # Keep ancestry as long format
  summary_table_wide <- make_wide(summary_table,
    c("meta_analysis", "cluster")
  )

  # Write csv
  save_csv(summary_table_wide, summary_table_path, "mrmega")

  summary_table_wide
}

# --------------------------------
create_table_fixed <- function(summary_table_path) {
  # Read in files
  dataframes <- read_files("results/meta/antidep-2501", "fixed")

  # Concatenate the files together (rbind) and add column for the file name
  summary_table <- rbind_dataframes_list(dataframes,
    c("dataset", "build", "pheno", "version", "cluster")
  )

  summary_table_wide <- make_wide(summary_table, "meta_analysis")

  # Write csv
  save_csv(summary_table_wide, summary_table_path, "fixed")

  summary_table_wide
}

# --------------------------------
# Create tables and save them as csv files
create_summary_tables <- function(summary_table_path) {
  create_table_MRMEGA(summary_table_path)

  create_table_fixed(summary_table_path)
}

# --------------------------------
# Call main function
# prints both summary tables in console
# saves both summary tables to "manuscript/tables"

create_summary_tables("manuscript/tables")
