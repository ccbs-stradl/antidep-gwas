import os
import xlwings as xw
from xlwings.main import Range

os.environ["XLWINGS_LICENSE_KEY"] = "noncommercial"

# ----------- FILE EXISTS ---------------------
# check excel file with test data exists
# this file is created by running generate_test_data.R
def check_file_exists(file_path: str) -> bool:
    """Check if the file exists."""
    return os.path.isfile(file_path)

# ----------- README SHEET --------------------
# check content of cells in the readme sheet exists
# content should be in the following cells: row 1 & col 1, row 2 & col 1, row 3-10 & cols 1-3
def check_readme_cell_contents_exist(file_path: str) -> bool:
    """Check if the expected cells in the README sheet contain content"""
    wb = xw.Book(file_path)
    readme = wb.sheets[0]
    cols_to_check = ['A1:A10' , 'B3:B10', 'C3:C10']
    for col in cols_to_check:
        values = readme.range(col).value
        is_empty = any(value is None for value in values)
        if is_empty:
            return False
    return True

# ----------- BOLD STYLE --------------------
def cell_not_bold(cell) -> bool:
    """Return false if the cell is not bold."""
    cell_is_not_bold = not cell.font.bold
    return cell_is_not_bold

def get_cells(file_path: str, sheet_index: int, cell_range: str) -> Range:
    """Return cell range from given sheet."""
    wb = xw.Book(file_path)
    cells = wb.sheets[sheet_index][cell_range]
    return cells

def check_cells_are_bold(file_path: str, sheet_index: int, cell_range: str) -> bool:
    """Check if a range of cells are bold."""
    cells = get_cells(file_path, sheet_index, cell_range)
    # Check if any of the given cells are not bold
    # if any cell is not bold this will return True
    is_not_bold = any(cell_not_bold(cell) for cell in cells)
    if is_not_bold:
        return False
    return True

# check cells in row 3 and cols 1-3 have the words: sheet_name, column and description in each cell

# ----------- NON README SHEETS -------------------
# check sheets (except README) contain content

# check that rows are bolded in any sheets (except README)
# when they meet the specified condition(s) in the given column(s) (when supplied by the user)

# ----------- SHEETS -------------------
# check number of sheets

# check sheet names