from openpyxl import Workbook
import pytest

def test_excel_exists():
    """
    Test the Excel file exists
    """
    # Load the workbook
    wb = Workbook()
    wb = load_workbook(filename = 'manuscript/tests/excel_pipeline/test_data/xlsx/S1_test_excel.xlsx')
    assert wb is not None, "Excel test data does not exist"
