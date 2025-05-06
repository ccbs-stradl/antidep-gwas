import os
import xlwings as xw
from xlwings.main import Range
import pandas as pd
import operator

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
    cols_to_check = ['A1:A10', 'B3:B10', 'C3:C10']
    for col in cols_to_check:
        values = readme.range(col).value
        is_empty = any(value is None for value in values)
        if is_empty:
            return False
    return True


# ----------- GET CELLS ---------------------
def get_cells(file_path: str, sheet_index: int, cell_range: str) -> Range:
    """Return cell range from given sheet."""
    wb = xw.Book(file_path)
    if cell_range != 'used_range':
        cells = wb.sheets[sheet_index][cell_range]
    else:
        cells = wb.sheets[sheet_index].used_range
    return cells


# get the cell contents for a range of cells
# check those contents match a list of strings of expected values
def get_cell_values(file_path: str, sheet_index: int, cell_range: str) -> list:
    contents = get_cells(file_path, sheet_index, cell_range)
    values = contents.value
    return values


# ----------- BOLD STYLE --------------------
def cell_not_bold(cell: Range) -> bool:
    """Return false if the cell is not bold."""
    cell_is_not_bold = not cell.font.bold
    return cell_is_not_bold


def check_cells_are_bold(file_path: str, sheet_index: int, cell_range: str) -> bool:
    """Check if a range of cells are bold."""
    cells = get_cells(file_path, sheet_index, cell_range)
    # Check if any of the given cells are not bold
    # if any cell is not bold this will return True
    is_not_bold = any(cell_not_bold(cell) for cell in cells)
    if is_not_bold:
        return False
    return True


# ----------- CELL CONTENTS ------------
def check_cell_values_match_expected(file_path: str,
                                     sheet_index: int,
                                     cell_range: str,
                                     expected_values: list) -> bool:
    values = get_cell_values(file_path, sheet_index, cell_range)
    values_match_expected = values == expected_values
    return values_match_expected


# ----------- NON README SHEETS -------------------
# check sheets (except README) contain content

# check that rows are bolded in any sheets (except README)
# when they meet the specified condition(s) in the given column(s) (when supplied by the user)

def convert_sheet_to_dataframe(file_path, sheet_index) -> pd.DataFrame:
    """Convert sheet data from Excel to a DataFrame."""
    sheet = get_cell_values(file_path, sheet_index, 'used_range')
    return pd.DataFrame(sheet[1:], columns=sheet[0])


def get_matched_rows(df, column_names, condition, threshold) -> list[int]:
    """Get rows matching the given conditions."""
    # Get a dictionary of possible operators using the operator module
    ops = {
        '<': operator.lt,
        '>': operator.gt,
        '==': operator.eq,
        '<=': operator.le,
        '>=': operator.ge,
        '!=': operator.ne
    }

    # for each item in column name, condition and threshold values
    # get the condition which is to be met for a row to be bold
    matched_rows_list = []
    for col, cond, thresh in zip(column_names, condition, threshold):
        # get the row indices which satisfy this condition
        row_indices = df.index[ops[cond](df[col], thresh)].tolist()
        # add 2 to account for starting at 0 and header row, gets actual row numbers
        excel_row_indices = [index + 2 for index in row_indices]
        matched_rows_list.append(excel_row_indices)

    # if there are multiple conditions to match,
    # ie. multiple items in column_names, condition and threshold
    # then take the intersection of these row indices
    # so that only rows that match all conditions are returned
    if len(matched_rows_list) > 1:
        matched_rows = sorted(set.intersection(*map(set, matched_rows_list)))
    else:
        matched_rows = matched_rows_list[0]

    return matched_rows

def check_conditional_bold_cells(file_path: str,
                                 sheet_indices: list,
                                 column_names: list,
                                 condition: list,
                                 threshold: list) -> bool:

    # define helper functions not needed outside this function:
    def rows_are_bold(sheet_index, matched_rows, max_col_letter) -> bool :
        """Check if rows in the matched rows list are bold."""
        row_is_bold = []
        for row in matched_rows:
            cell_range = 'A' + str(row) + ':' + max_col_letter + str(row)

            bold_row = check_cells_are_bold(file_path,
                                            sheet_index,
                                            cell_range)
            row_is_bold.append(bold_row)

        # check that all these cells are bold
        all_rows_are_bold = all(row_is_bold)

        return all_rows_are_bold

    def get_max_col_letter(df) -> str:
        """Get the letter of the final column in the excel sheet."""
        number_of_cols = len(df.columns)
        max_col_letter = chr(ord('@') + number_of_cols)
        return max_col_letter

    # create empty list to store True/False values for each sheet
    boldness_by_sheet = []

    # loop over each sheet in sheet_indices to check each sheet individually
    # for each sheet get all the values
    for sheet_index in sheet_indices:

        df = convert_sheet_to_dataframe(file_path, sheet_index)

        # for each item in column name, condition and threshold values
        # get the condition which is to be met for a row to be bold
        matched_rows = get_matched_rows(df, column_names, condition, threshold)

        # Get the letter of the final column in the excel sheet
        max_col_letter = get_max_col_letter(df)

        # iterate over each row that should be bold
        all_rows_are_bold = rows_are_bold(sheet_index,
                                          matched_rows,
                                          max_col_letter)

        # append all_rows_are_bold to the sheet value
        boldness_by_sheet.append(all_rows_are_bold)

    all_rows_in_all_sheets_are_bold = all(boldness_by_sheet)

    return all_rows_in_all_sheets_are_bold

# ----------- SHEETS -------------------
# check number of sheets

# check sheet names
