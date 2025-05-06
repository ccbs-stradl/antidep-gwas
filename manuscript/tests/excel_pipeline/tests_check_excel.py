
from funs_check_excel import (check_file_exists,
                              check_readme_cell_contents_exist,
                              check_cells_are_bold)

file_path = "test_data/xlsx/S1_test_excel.xlsx"

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

# ----------- NON README SHEETS -------------------
# check sheets (except README) contain content

# check that rows are bolded in any sheets (except README)
# when they meet the specified condition(s) in the given column(s) (when supplied by the user)

# ----------- SHEETS -------------------
# check number of sheets

# check sheet names