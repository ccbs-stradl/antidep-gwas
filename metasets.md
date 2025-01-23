Antidepressant exposure GWAS meta-analysis sets
================

Sets of input summary statistics to group for meta-analysis.

``` r
library(readr)
library(dplyr)
```

``` r
datasets <- read_csv("datasets.csv")
```

    ## Rows: 57 Columns: 8
    ## ── Column specification ───────────────────────────────────────────────────────────────────────────────────────────────────────────────────
    ## Delimiter: ","
    ## chr (6): cohort, pheno, dataset, version, build, cluster
    ## dbl (2): cases, controls
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
datasets |> distinct(cohort, version)
```

    ## # A tibble: 7 × 2
    ##   cohort  version    
    ##   <chr>   <chr>      
    ## 1 AllOfUs v7_2408    
    ## 2 AllOfUs v7_2411    
    ## 3 GenScot rmUKBB_2309
    ## 4 BBJ     SK2020     
    ## 5 FinnGen R9         
    ## 6 FinnGen R12        
    ## 7 UKB     21007_2309

``` r
metaset_2408_list <- list(
  list(cohort = "AllOfUs", version = "v7_2408"),
  list(cohort = "GenScot", version = "rmUKBB_2309"),
  list(cohort = "BBJ",     version = "SK2020"),
  list(cohort = "FinnGen", version = "R9"),
  list(cohort = "UKB",     version = "21007_2309")
)
metaset_2408 <- bind_rows(lapply(metaset_2408_list, as_tibble))

metaset_2501_list <- list(
  list(cohort = "AllOfUs", version = "v7_2411"),
  list(cohort = "GenScot", version = "rmUKBB_2309"),
  list(cohort = "BBJ",     version = "SK2020"),
  list(cohort = "FinnGen", version = "R12"),
  list(cohort = "UKB",     version = "21007_2309")
)
metaset_2501 <- bind_rows(lapply(metaset_2408_list, as_tibble))

dir.create("metasets")
```

    ## Warning in dir.create("metasets"): 'metasets' already exists

``` r
write_csv(metaset_2408, "metasets/antidep-2408.csv")
write_csv(metaset_2501, "metasets/antidep-2501.csv")
```
