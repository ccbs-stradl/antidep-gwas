Between GWAS genetic correlations
================

``` r
library(readr)
library(dplyr)
library(ggplot2)
library(corrplot)
library(stringr)
```

Read in LDSC (within-ancestry) and Popcorn (cross ancestry) genetic
correlation estimates. Retain only estimates where $se <= 1/1.96 = 0.51$
(i.e., if $r_g = 1$ it could be distinguished from zero).

``` r
ldsc <- read_csv(here::here("manuscript/tables/rg_ldsc_gwas.csv"))
```

    ## Rows: 166 Columns: 18
    ## ── Column specification ───────────────────────────────────────────────────────────────────────────────────────────────────
    ## Delimiter: ","
    ## chr  (8): p1_cohort, p1_pheno, p1_cluster, p1_version, p2_cohort, p2_pheno, ...
    ## dbl (10): rg, se, z, p, h2_obs, h2_obs_se, h2_int, h2_int_se, gcov_int, gcov...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
popcorns <- read_csv(here::here("manuscript/tables/rg_popcorn_gwas.csv"))
```

    ## Rows: 168 Columns: 12
    ## ── Column specification ───────────────────────────────────────────────────────────────────────────────────────────────────
    ## Delimiter: ","
    ## chr (8): p1_cohort, p1_pheno, p1_cluster, p1_version, p2_cohort, p2_pheno, p...
    ## dbl (4): pgi, SE, Z, P (Z)
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
rgs <-
bind_rows(
    ldsc,
    rename(popcorns, rg = pgi, se = SE, z = Z, p = `P (Z)`)
) |>
filter(se <= 1/qnorm(0.975))
```

Set up dataset names for correlation matrix

``` r
p1_datasets <- rgs |>
    mutate(dataset = str_c(p1_cohort, p1_pheno, p1_cluster, sep = "-")) |>
    pull(dataset)

p2_datasets <- rgs |>
    mutate(dataset = str_c(p2_cohort, p2_pheno, p2_cluster, sep = "-")) |>
    pull(dataset)

datasets <- unique(c(p1_datasets, p2_datasets))

rg <- rgs |> pull(rg)
```

Matrix to store correlations

``` r
rg_matrix <- matrix(data = NA,
                    nrow = length(datasets),
                    ncol = length(datasets),
                    dimnames = list(datasets, datasets))

for(i in seq_along(rg)) {
    dataset1 <- p1_datasets[i]
    dataset2 <- p2_datasets[i]
    rg_matrix[dataset1, dataset2] <- rg_matrix[dataset2, dataset1] <- rg[i]
}
```

Squish genetic correlations to (-1, 1)

``` r
rg_matrix_squish <- rg_matrix
rg_matrix_squish[rg_matrix > 2] <- NA
rg_matrix_squish[rg_matrix < -2] <- NA
rg_matrix_squish[rg_matrix < 2 & rg_matrix > 1] <- 1
rg_matrix_squish[rg_matrix > -2 & rg_matrix < -1] <- -1
```

Correlation plots

``` r
corrplot(rg_matrix, is.corr = FALSE)
```

![](rg-gwas_files/figure-gfm/corrplot-1.png)<!-- -->

``` r
corrplot(rg_matrix_squish)
```

![](rg-gwas_files/figure-gfm/corrplot_squish-1.png)<!-- -->

Text to paste in manuscript results section for min/max observed
correlations.

``` r
# get unique rgs
rgs_distinct <- rgs |>
  filter(p1_cohort <= p2_cohort, p1_pheno <= p2_pheno, p1_cluster <= p2_cluster) |>
  mutate(rg = signif(rg, digits = 3))

# within ancestry, cross cohort
rg_w_anc_x_cohort <- rgs_distinct |>
  filter(p1_cluster == p2_cluster, p1_cohort != p2_cohort)
# cross ancestry, within cohort
rg_x_anc_w_cohort <- rgs_distinct |>
  filter(p1_cluster != p2_cluster, p1_cohort == p2_cohort)
# cross ancestry, cross cohort
rg_x_anc_x_cohort <- rgs_distinct |>
  filter(p1_cluster != p2_cluster, p1_cohort != p2_cohort)

correlation_text <- function(correlations, intro_text) {
  cmin <- correlations |>
    filter(rg == min(rg))
  cmax <- correlations |>
    filter(rg == max(rg))

  str_glue("{intro_text} r_g ranged from {cmin$rg} between {cmin$p1_pheno} in {cmin$p1_cohort} {cmin$p1_cluster} and {cmin$p2_pheno} in {cmin$p2_cohort} {cmin$p2_cluster} to {cmax$rg} between {cmax$p1_pheno} in {cmax$p1_cohort} {cmax$p1_cluster} and {cmax$p2_pheno} in {cmax$p2_cohort} {cmax$p2_cluster}.")
}

cat(correlation_text(rg_w_anc_x_cohort, "Within clusters and across cohorts,"), " ", correlation_text(rg_x_anc_x_cohort, "Across clusters and cohorts,"))
```

Within clusters and across cohorts, r_g ranged from 0.02 between N06AA
in AllOfUs SAS and N06AB in UKB SAS to 0.93 between N06A in GenScot EUR
and N06AB in UKB EUR. Across clusters and cohorts, r_g ranged from 0.462
between N06A in FinnGen EUR and N06AB in UKB SAS to 0.637 between N06A
in AllOfUs EUR and N06A in UKB SAS.
