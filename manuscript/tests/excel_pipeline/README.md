# Test Excel spreadsheets created by the functions in `manuscript/scripts/supplementary_tables_excell.R`.

1) create test .csv and .cols files
2) run the pipeline that generates the xlsx file in R
3) Use pytest and xlwings to read the test .xlsx files and run assert equal statements to check cell content and style.

## Running the tests

### Create a .venv and install the requirements
```bash
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

### Run the tests
This script will run the R script `generate_test_data.R` which creates test .csv and .cols files and generate the xlsx file using the R functions used to create the real xlsx files. The tests then check the xlsx test file for expected content and style.
```
pytest manuscript/tests/excel_pipeline/tests_check_excel.py
```
