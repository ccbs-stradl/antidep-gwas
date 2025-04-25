# Test to check excel spreadsheet renders as expected 

Excel spreadsheets created in `manuscript/scripts/supplementary_tables_excell.R` are tested here.

1) create test .csv and .cols files
2) run the pipeline that generates the xlsx file in R
3) Use pytest and openpyxl to read the test .xlsx files and run assert equal statements to check cell content and style.

## Running the tests

### Manually

```
Rscript manuscript/tests/generate_test_data.R
pytest manuscript/tests/e2e_test_excel.py
```

### GitHub actions
Set up GitHub actions to run tests whenever changes are made to
`manuscript/tests/excel_pipeline/` or `manuscript/scripts/supplementary_tables_excel*`
