
import pytest
from funs_check_excel import (check_file_exists,
                              check_readme_cell_contents_exist,
                              check_legend_title_is_bold)

file_path = "test_data/xlsx/S1_test_excel.xlsx"

def test_sheets_exist():
    assert check_file_exists(file_path) is True

def test_readme_cell_contents_exist():
    assert check_readme_cell_contents_exist(file_path) is True

# check content of cells in the readme sheet exists
# content should be in the following cells: row 1 & col 1, row 2 & col 1, row 3-10 & cols 1-3

# check cell in row 1 and col 1 in first sheet is bold (this is the legend title)
def test_legend_title_is_bold():
    assert check_legend_title_is_bold(file_path) is True

# check cell in row 2 and col 1 in first sheet has no style

# check cells in row 3 and cols 1-3 are bold

# check cells in row 3 and cols 1-3 have the words: sheet_name, column and description in each cell