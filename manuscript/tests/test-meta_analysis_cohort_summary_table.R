# Script to test functions in "manuscript/scripts/meta_analysis_cohort_summary_table.R"
# To run tests setwd to root of repo and run in console:
# source("manuscript/scripts/test-meta_analysis_cohort_summary_table.R")
# -----------------------------
library(testthat)

# -----------------------------
# Source file with main function and helper functions we want to test
source("manuscript/scripts/meta_analysis_cohort_summary_table.R")

# -----------------------------
test_that("directory to save summary tables exist", {
  test_data <- data.frame(A = c(1,2,3))
  expect_error(save_csv(test_data, "this/path/does/not/exist", "metatype"),
               "Error: this/path/does/not/exist directory does not exist.")
})

# -----------------------------
test_that("MR-MEGA files exist in meta directory",{
  expect_no_error(read_files("meta", "mrmega"))
})

# -----------------------------
test_that("Fixed meta analysis files exist in meta directory",{
  expect_no_error(read_files("meta", "fixed"))
})

# -----------------------------
# Add tests for:
# *summary_table_fixed.csv matches values in input .csv (test for a couple of inputs)
  # this is because there is a risk when pivoting the table from long to wide and reordering the cols that something could go wrong at this step