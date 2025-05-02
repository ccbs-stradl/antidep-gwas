
from funs_check_excel import check_file_exists
import pytest

def test_sheets_exist():
    file_path = "test_data/xlsx/S1_test_excel.xlsx"
    assert check_file_exists(file_path) is True