# Run this script after other tables_rg_*.R scripts, if you need to add CIs
# columns to those tables. 
# Then run supplementary_tables_excell.R to update the Supplementary table file.
# Note that this script does not make the sidecar .cols files with updated CI so 
# manual edit of the README sheet in excel may be required.
# ------------------------------------------------------------------------------
# Function to read in tables and returned a named list,
# where the name is the file name without extension
# and the value is the data frame of the table
read_tables <- function(file_names) {
  tables_list <- sapply(file_names, function(file_name) {
    read_csv(here::here("manuscript/tables/", file_name))
  })
  names(tables_list) <- sapply(file_names, function(file_name) {
    tools::file_path_sans_ext(basename(file_name))
  })
  tables_list
}

# Function to add a CI column to the rg results tables
# estimate and SE columns must be named "rg" and "se" in ldsc tables
# and "pgi" and "SE" in popcorn tables
add_ci_column <- function(tables_list) {
  for (name in names(tables_list)) {
    table <- tables_list[[name]]

    if ("rg" %in% colnames(table) && "se" %in% colnames(table)) {
      table <- table %>%
        mutate(ci = str_c(
          "(", round(rg - 1.96 * se, 3),
          ", ", round(rg + 1.96 * se, 3), ")"
        )) %>%
        relocate(ci, .after = se)
    } else if ("pgi" %in% colnames(table) && "SE" %in% colnames(table)) {
      table <- table %>%
        mutate(CI = str_c(
          "(", round(pgi - 1.96 * SE, 3),
          ", ", round(pgi + 1.96 * SE, 3), ")"
        )) %>%
        relocate(CI, .after = SE)
    }

    tables_list[[name]] <- table
  }

  tables_list
}

# Resave the tables by overwriting the original files
save_tables <- function(tables_list) {
  for (name in names(tables_list)) {
    file_name <- here::here("manuscript/tables/", paste0(name, ".csv"))
    write_csv(tables_list[[name]], file_name)
  }
}


# Main function that reads in tables and add CI
# columns, then saves the updated tables
process_rg_tables <- function(file_names) {
  tables_list <- read_tables(file_names)
  tables_list <- add_ci_column(tables_list)
  save_tables(tables_list)
}


# ------------------------------------------------------------------------------
# Tables to read in:
file_names <- c(
  "rg_ldsc_gwas.csv",
  "rg_ldsc_meta.csv",
  "rg_ldsc_meta_external.csv",
  "rg_popcorn_gwas.csv",
  "rg_popcorn_meta.csv"
)

# Call main function to process the tables
process_rg_tables(file_names)


