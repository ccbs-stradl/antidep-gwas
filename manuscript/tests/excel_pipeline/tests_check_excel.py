# Run Python tests to check excel file looks as expected
# Create test data by running the same R function that creates real data
# on test input csv files
# -----------------------------------------------
# Import functions to run tests
from funs_check_excel import (check_file_exists,
                              check_readme_cell_contents_exist,
                              check_cells_are_bold,
                              check_cell_values_match_expected,
                              check_conditional_bold_cells)

# Define the file_path variable for the test data created
file_path = "test_data/xlsx/S1_test_excel.xlsx"

# Import module to source R script that generates test data
import subprocess


def setup_module(module):
    """Run this once before all tests. Generate test data."""
    subprocess.run(['Rscript', 'generate_test_data.R'])


# ----------- FILE EXISTS ---------------------
def test_sheets_exist():
    assert check_file_exists(file_path) is True


# ----------- README SHEET --------------------
# check content of cells in the readme sheet exists
# content should be in the following cells: row 1 & col 1, row 2 & col 1, row 3-10 & cols 1-3
def test_readme_cell_contents_exist():
    assert check_readme_cell_contents_exist(file_path) is True


# check cell in row 1 and col 1 ie. cell A1 in first sheet is bold (this is the legend title)
def test_legend_title_is_bold():
    assert check_cells_are_bold(file_path, 0, 'A1') is True


# check cells in row 3 and cols 1-3 are bold
def test_col_meta_headings_are_bold():
    assert check_cells_are_bold(file_path, 0, 'A3:C3') is True


# check cells in row 3 and cols 1-3 have the words: sheet_name, column and description in each cell
def test_col_meta_headings_content():
    assert check_cell_values_match_expected(file_path,
                                            0,
                                            'A3:C3',
                                            ['sheet_name', 'column', 'description']) is True


# ----------- NON README SHEETS -------------------
# check sheets (except README) contain content

# check that rows are bolded in any sheets (except README)
# when they meet the specified condition(s) in the given column(s) (when supplied by the user)
def test_conditional_bold_cells():
    assert check_conditional_bold_cells(file_path,
                                        list(range(1, 3)),
                                        ['p_val', 'CS_PIP'],
                                        ['<', '>'],
                                        [-1, 1]) is True


# check rows that should not be bold are not bold
def test_conditional_not_bold_cells():
    assert check_conditional_bold_cells(file_path,
                                        list(range(1, 3)),
                                        ['p_val'],
                                        ['>'],
                                        [0.05]) is False
# ----------- SHEETS -------------------
# check number of sheets

# check sheet names
