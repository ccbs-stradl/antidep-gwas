# Functions for supplementary_tables_excell.R
# ---------------------------------------------
# Function to update the counter of the supplementary table eg. S1, S2, S3
update_table_index <- function(table_index) {
  # check table_index is numeric
  if (!is.numeric(table_index)) {
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
get_main_file_names <- function(file_path, file_regex) {
  files <- list.files(file_path, full.names = TRUE)
  files <- str_subset(files, file_regex)

  # If there are no files stop with error
  if (length(files) == 0) {
    stop(paste0("No files found at: ", file_path, " with regex: ", file_regex))
  }

  # If files are found, check if they all end with .csv, .tsv or .cols
  # If they are not csv or tsv, stop with error
  if (!any(grepl("\\.csv$", files)) & !any(grepl("\\.tsv$", files)) & !any(grepl("\\.col$", files))) {
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
read_results <- function(full_path) {
  # Sometimes fread does not read csvs correctly eg. if there is an extraneous comma at the end of a line
  # in this case read_csv() can be used instead, however it is slower
  main <- tryCatch(
    {
      # Read in results tables
      fread(full_path)
    },
    warning = function(w) {
      # Print the warning message
      message(
        "Warning: ", conditionMessage(w),
        paste0("\nUsing read_csv() instead of fread() for ", full_path)
      )
      read_csv(full_path)
    }
  )

  # Read in sidecar .cols files
  meta <- fread(paste0(full_path, ".cols"))

  # Combine into a list
  results <- list(main = main, meta = meta)

  return(results)
}

# ---------------------------------------------
# function that creates a list of results tables
# and names them by the sheet name

tables_list <- function(full_path_vector, sheet_names) {
  # Check if full_path_vector is a character vector
  if (!is.character(full_path_vector)) {
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
                       cell_title_width, cell_title_height, 
                       bold_cols, bold_condition, bold_threshold) {
  wb <- createWorkbook()

  # Add a readme for the first sheet
  add_readme(
    results, wb, table_index, letter_index, legend_title,
    legend_text_prefix, legend_text_sections,
    cell_title_width, cell_title_height
  )

  # Add each table to a new sheet
  for (i in seq_along(results)) {
    addWorksheet(wb, sheetName = names(results)[i])
    writeData(wb, sheet = names(results)[i], results[[i]]$main, startRow = 1, colNames = TRUE)
  }

  # Style each sheet with bold font based on user specified condition
  # bold_cols, bold_condition, bold_threshold = NULL by default
  bold_rows(wb, results, bold_cols, bold_condition, bold_threshold)

  # Auto-fit columns
  auto_fit_columns(wb, results)

  saveWorkbook(wb, file_name, overwrite = TRUE)
}

# ---------------------------------------------
# function that creates README for first sheet in excel spreadsheet
add_readme <- function(results, wb, table_index, letter_index, legend_title,
                       legend_text_prefix, legend_text_sections,
                       cell_title_width, cell_title_height) {
  # Create a readme sheet
  addWorksheet(wb, sheetName = "README")

  # Add legend title
  add_legend_title(table_index, legend_title, wb)

  # Add legend text
  add_legend_text(
    legend_text_sections,
    legend_text_prefix, letter_index,
    wb
  )

  # Add meta data on column names
  add_metadata(results, wb)

  # Style the readme
  style_readme(wb, cell_title_width, cell_title_height)
}

# --------------------------------------------
# function to add legend title
add_legend_title <- function(table_index, legend_title, wb) {
  # Update legend_title
  legend_title <- as.character(glue("Table S{table_index}. {legend_title}"))

  # Write the table legend title
  writeData(wb, sheet = "README", legend_title, startRow = 1, startCol = 1)
}

# --------------------------------------------
# function to add legend text
add_legend_text <- function(legend_text_sections,
                            legend_text_prefix, letter_index,
                            wb) {
  # check letter_text_sections equals length of letter_index
  if (length(legend_text_sections) != length(letter_index)) {
    stop("legend_text_sections and letter_index must be the same length")
  }

  # prepend each legend_text_section with the letter_index
  legend_text_sections <- paste0("(", letter_index, ") ", legend_text_sections)

  # concatenate the legend_text_sections with commas & "and" for the last item
  section_strings <- concatenate_comma_and(legend_text_sections)

  # concatenate the legend_text_prefix with the legend_text_sections
  legend_text <- paste0(legend_text_prefix, section_strings)

  # Write the table legend text
  writeData(wb, sheet = "README", legend_text, startRow = 2, startCol = 1)
}

# --------------------------------------------
# function that concatenates strings with commas and "and"
concatenate_comma_and <- function(strings) {
  # start with comma separators for each item
  separators <- rep(", ", times = length(strings))
  # replace last separator with "and"
  separators[length(separators)] <- " and "
  # replace first separator with empty string
  # (separators proceed strings)
  separators[1] <- ""
  # concatenate: "a, b, and c" etc
  text <- str_c(separators, strings, collapse = "")
  return(text)
}

# --------------------------------------------
# function to add column description metadata
# Each item in 'results' contains a 'meta' data frame describing column names
# (from .cols sidecar files).
# Extract that meta data and combine them into a single data frame:
# add a column for which sheet the column name is from ie. results name
add_metadata <- function(results, wb) {
  # extract the meta data from each results item
  meta_list <- lapply(results, function(x) {
    x$meta
  })

  # rename list of meta tables with names of results so it can be used by bind_rows
  names(meta_list) <- names(results)

  # combine all meta data into a single data frame and add id column
  col_name_descriptions <- bind_rows(meta_list, .id = "sheet_name")

  writeData(wb, sheet = "README", col_name_descriptions, startRow = 3, startCol = 1)
}

# --------------------------------------------
# function to set style in README file
style_readme <- function(wb, cell_title_width, cell_title_height) {
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

  # make row with column name and descriptions bold
  addStyle(wb, sheet = "README", style = bold_style, rows = 3, cols = 1:3, stack = TRUE)

  # autofit 2nd and 3rd columns to width of cell
  setColWidths(wb, sheet = "README", cols = 2:3, widths = "auto")
}

# --------------------------------------------
# function to style rows as bold if they meet user specified condition
bold_rows <- function(wb, results, bold_cols, bold_condition, bold_threshold) {
  # Use gaurd clause to check if bold_cols, bold_condition and bold_threshold are all provided
  if (is.null(bold_cols) || is.null(bold_condition) || is.null(bold_threshold)) {
    # return early from the function and don't do anything if these are all NULL
    return()
  }

  # Use guard clause to check if bold_cols, bold_condition and bold_threshold are all the same length
  if (length(bold_cols) != length(bold_condition) || length(bold_cols) != length(bold_threshold)) {
    stop("bold_cols, bold_condition, and bold_threshold must all be the same length.")
  }

  # Check if the column name exists in all sheets
  for (sheet in names(results)) {
    for (bold_col in bold_cols) {
      if (!bold_col %in% colnames(results[[sheet]]$main)) {
        stop(paste0("Column '", bold_col, "' not found in sheet '", sheet, "'"))
      }
    }
  }

  # Loop through each sheet and make significant rows bold
  for (sheet in names(results)) {
    rows_to_bold <- lapply(seq_along(bold_cols), function(i){
      col <- bold_cols[i]
      cond <- bold_condition[i]
      thresh <- bold_threshold[i]
      
      # Get the column index of the significant column
      col_index <- which(colnames(results[[sheet]]$main) == col)
  
      # Get the row index that meets the user specified condition
      op <- match.fun(cond)
      rows_matching_condition <- op(results[[sheet]]$main[[col_index]], thresh)
      row_index <- which(rows_matching_condition) + 1 # + 1 because of header
      return(row_index)
    })
    
    # Get intersection of all row indices
    row_index <- Reduce(intersect, rows_to_bold)
    
    # Make the significant rows bold
    addStyle(wb,
      sheet = sheet, style = createStyle(textDecoration = "bold"),
      rows = row_index, cols = 1:ncol(results[[sheet]]$main), stack = TRUE,
      gridExpand = TRUE
    )
  }
}

# --------------------------------------------
# auto-fit columns
auto_fit_columns <- function(wb, results) {
  # Loop through each sheet and auto-fit columns
  for (sheet in names(results)) {
    # Auto-fit all columns
    setColWidths(wb, sheet = sheet, cols = 1:ncol(results[[sheet]]$main), widths = "auto")
  }
}


# --------------------------------------------
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
  sheet_names <- paste0(glue("Table S{paste0(table_index, letter_index)} "), sheet_names)

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
