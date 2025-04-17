Antidepressant exposure GWAS SNP lists for forest plots
================

Requires
[SNPlocs.Hsapiens.dbSNP150.GRCh38](https://bioconductor.org/packages/release/data/annotation/html/SNPlocs.Hsapiens.dbSNP150.GRCh38.html)
for SNP lookup

``` r
library(readr)
library(dplyr)
library(stringr)
library(SNPlocs.Hsapiens.dbSNP150.GRCh38)
```

List of SNPs from meta-analyses and SusieX

``` r
metaset <- "antidep-2501"
fixed <- read_csv(here::here(str_glue("manuscript/tables/clumps_fixed_{metaset}.clumps.csv")))
```

    ## Rows: 77 Columns: 23
    ## ── Column specification ───────────────────────────────────────────────────────────────────────────────────────────────────
    ## Delimiter: ","
    ## chr  (5): dataset, CHROM, ID, REF, ALT
    ## dbl (18): LOCUS, POS, studies, BETA, SE, CHISQ, P, Q, QP, INFO, AFCAS, AFCON...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
mrmega <- read_csv(here::here(str_glue("manuscript/tables/clumps_mrmega_{metaset}.clumps.csv")))
```

    ## Rows: 33 Columns: 33
    ## ── Column specification ───────────────────────────────────────────────────────────────────────────────────────────────────
    ## Delimiter: ","
    ## chr  (6): dataset, MarkerName, EA, NEA, Effects, CHR
    ## dbl (26): Locus, Chromosome, Position, EAF, Nsample, Ncohort, beta_0, se_0, ...
    ## lgl  (1): Comments
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
susiex <- read_csv(here::here("manuscript/tables/susiex_significant_cs.csv"))
```

    ## Rows: 4 Columns: 26
    ## ── Column specification ───────────────────────────────────────────────────────────────────────────────────────────────────
    ## Delimiter: ","
    ## chr  (7): SNP, REF_ALLELE_AFR, REF_ALLELE_SAS, REF_ALLELE_EUR, ALT_ALLELE_AF...
    ## dbl (19): CHR, BP_START, BP_END, CS_ID, BP, REF_FRQ_AFR, REF_FRQ_SAS, REF_FR...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

Get SNPs for meta-analysis (build 38)

``` r
fixed_snps <- fixed |>
  select(CHROM, POS, ID, REF, ALT)
mrmega_snps <- mrmega |>
  mutate(CHROM = if_else(Chromosome == 23, true = "chrX", false = str_c("chr", Chromosome))) |>
  select(CHROM, POS = Position, ID = MarkerName, REF = NEA, ALT = EA)
```

Lookup build 38 coordinates for SusieX results

``` r
# ref/alt alleles should all be the same, but pick AFR since it
# is available for all selected SNPs
susiex_locs <- snpsById(SNPlocs.Hsapiens.dbSNP150.GRCh38, ids = pull(susiex, SNP))
susiex_snps <- susiex_locs |>
  as.data.frame() |>
  as_tibble() |>
  left_join(susiex, by = c("RefSNP_id" = "SNP")) |>
  transmute(CHROM = str_c("chr", as.character(seqnames)),
            POS = pos, 
            ID = RefSNP_id,
            REF = REF_ALLELE_AFR, ALT = ALT_ALLELE_AFR)
```

Merge unique SNPs

``` r
snps <- bind_rows(fixed_snps, mrmega_snps, susiex_snps) |>
  distinct()

dir.create(here::here("results/forest"))
```

    ## Warning in dir.create(here::here("results/forest")):
    ## '/Users/mark/Work/antidep-gwas/results/forest' already exists

``` r
write_tsv(snps, here::here(str_glue("results/forest/{metaset}.snplist")), col_names = FALSE)
```
