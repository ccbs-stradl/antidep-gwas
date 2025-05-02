import os
import xlwings as xw

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
        value = readme.range(col).value
        if value is None or str(value).strip() == '':
            return False
        else:
            return True

# check cell in row 1 and col 1 in first sheet is bold (this is the legend title)
def check_legend_title_is_bold(file_path: str) -> bool:
    """Check if the legend title in cell A1 is bold."""
    wb = xw.Book(file_path)
    legend = wb.sheets[0]['A1']
    legend_is_bold = legend.font.bold
    return legend_is_bold

# check cell in row 2 and col 1 in first sheet has no style

# check cells in row 3 and cols 1-3 are bold

# check cells in row 3 and cols 1-3 have the words: sheet_name, column and description in each cell

# ----------- NON README SHEETS -------------------
# check sheets (except README) contain content

# check that rows are bolded in any sheets (except README)
# when they meet the specified condition(s) in the given column(s) (when supplied by the user)

# ----------- SHEETS -------------------
# check number of sheets

# check sheet names