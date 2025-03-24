Antidepressant exposure GWAS fixed-effects meta-analysis
================

``` r
library(dplyr)
library(readr)
library(stringr)
library(topr)
library(UpSetR)
library(plyranges)
```

## Results

### Sumstats

``` r
meta_gzs <- list.files("meta", "fixed-.+meta\\.gz", full.names = TRUE)
metas <- str_remove(meta_gzs, ".meta.gz")
names(metas) <- str_remove(basename(metas), "fixed-")
           
sumstats_paths <- lapply(metas, function(.x) str_c(.x, "meta.gz", sep = "."))

sumstats <- lapply(sumstats_paths, function(path) {
    read_tsv(path) |>
        select(CHROM = CHR, POS = BP, ID = SNP,
               REF = A2, ALT = A1, P, OR,
               AF = AFCON, NCAS, NCON, NEFF)
})
```

    ## Rows: 18446355 Columns: 21
    ## ── Column specification ──────────────────────────────────────────────────────────────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (4): CHR, SNP, A1, A2
    ## dbl (17): BP, studies, OR, SE, P, OR_R, SE_R, P_R, Q, I, INFO, AFCAS, AFCON,...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
    ## Rows: 7250885 Columns: 21
    ## ── Column specification ──────────────────────────────────────────────────────────────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (4): CHR, SNP, A1, A2
    ## dbl (17): BP, studies, OR, SE, P, OR_R, SE_R, P_R, Q, I, INFO, AFCAS, AFCON,...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
    ## Rows: 8391178 Columns: 21
    ## ── Column specification ──────────────────────────────────────────────────────────────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (4): CHR, SNP, A1, A2
    ## dbl (17): BP, studies, OR, SE, P, OR_R, SE_R, P_R, Q, I, INFO, AFCAS, AFCON,...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
    ## Rows: 10748047 Columns: 21
    ## ── Column specification ──────────────────────────────────────────────────────────────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (4): CHR, SNP, A1, A2
    ## dbl (17): BP, studies, OR, SE, P, OR_R, SE_R, P_R, Q, I, INFO, AFCAS, AFCON,...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
    ## Rows: 8842191 Columns: 21
    ## ── Column specification ──────────────────────────────────────────────────────────────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (4): CHR, SNP, A1, A2
    ## dbl (17): BP, studies, OR, SE, P, OR_R, SE_R, P_R, Q, I, INFO, AFCAS, AFCON,...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
    ## Rows: 10195596 Columns: 21
    ## ── Column specification ──────────────────────────────────────────────────────────────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (4): CHR, SNP, A1, A2
    ## dbl (17): BP, studies, OR, SE, P, OR_R, SE_R, P_R, Q, I, INFO, AFCAS, AFCON,...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
    ## Rows: 18463581 Columns: 21
    ## ── Column specification ──────────────────────────────────────────────────────────────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (4): CHR, SNP, A1, A2
    ## dbl (17): BP, studies, OR, SE, P, OR_R, SE_R, P_R, Q, I, INFO, AFCAS, AFCON,...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
    ## Rows: 10357731 Columns: 21
    ## ── Column specification ──────────────────────────────────────────────────────────────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (4): CHR, SNP, A1, A2
    ## dbl (17): BP, studies, OR, SE, P, OR_R, SE_R, P_R, Q, I, INFO, AFCAS, AFCON,...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
    ## Rows: 10206559 Columns: 21
    ## ── Column specification ──────────────────────────────────────────────────────────────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (4): CHR, SNP, A1, A2
    ## dbl (17): BP, studies, OR, SE, P, OR_R, SE_R, P_R, Q, I, INFO, AFCAS, AFCON,...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
    ## Rows: 18465873 Columns: 21
    ## ── Column specification ──────────────────────────────────────────────────────────────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (4): CHR, SNP, A1, A2
    ## dbl (17): BP, studies, OR, SE, P, OR_R, SE_R, P_R, Q, I, INFO, AFCAS, AFCON,...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
    ## Rows: 10357712 Columns: 21
    ## ── Column specification ──────────────────────────────────────────────────────────────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (4): CHR, SNP, A1, A2
    ## dbl (17): BP, studies, OR, SE, P, OR_R, SE_R, P_R, Q, I, INFO, AFCAS, AFCON,...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
    ## Rows: 8840275 Columns: 21
    ## ── Column specification ──────────────────────────────────────────────────────────────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (4): CHR, SNP, A1, A2
    ## dbl (17): BP, studies, OR, SE, P, OR_R, SE_R, P_R, Q, I, INFO, AFCAS, AFCON,...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
    ## Rows: 10199157 Columns: 21
    ## ── Column specification ──────────────────────────────────────────────────────────────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (4): CHR, SNP, A1, A2
    ## dbl (17): BP, studies, OR, SE, P, OR_R, SE_R, P_R, Q, I, INFO, AFCAS, AFCON,...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

### Loci

``` r
clumped_paths <- sumstats_paths <- lapply(metas, function(.x) str_c(.x, "clumped", sep = "."))

clumped <- lapply(clumped_paths[sapply(clumped_paths, file.exists)], read_table)
```

    ## 
    ## ── Column specification ──────────────────────────────────────────────────────────────────────────────────────────────────────────────────
    ## cols(
    ##   CHR = col_double(),
    ##   F = col_double(),
    ##   SNP = col_character(),
    ##   BP = col_double(),
    ##   P = col_double(),
    ##   TOTAL = col_double(),
    ##   NSIG = col_double(),
    ##   S05 = col_double(),
    ##   S01 = col_double(),
    ##   S001 = col_double(),
    ##   S0001 = col_double(),
    ##   SP2 = col_character()
    ## )
    ## 
    ## 
    ## ── Column specification ──────────────────────────────────────────────────────────────────────────────────────────────────────────────────
    ## cols(
    ##   CHR = col_double(),
    ##   F = col_double(),
    ##   SNP = col_character(),
    ##   BP = col_double(),
    ##   P = col_double(),
    ##   TOTAL = col_double(),
    ##   NSIG = col_double(),
    ##   S05 = col_double(),
    ##   S01 = col_double(),
    ##   S001 = col_double(),
    ##   S0001 = col_double(),
    ##   SP2 = col_character()
    ## )
    ## 
    ## 
    ## ── Column specification ──────────────────────────────────────────────────────────────────────────────────────────────────────────────────
    ## cols(
    ##   CHR = col_double(),
    ##   F = col_double(),
    ##   SNP = col_character(),
    ##   BP = col_double(),
    ##   P = col_double(),
    ##   TOTAL = col_double(),
    ##   NSIG = col_double(),
    ##   S05 = col_double(),
    ##   S01 = col_double(),
    ##   S001 = col_double(),
    ##   S0001 = col_double(),
    ##   SP2 = col_character()
    ## )
    ## 
    ## 
    ## ── Column specification ──────────────────────────────────────────────────────────────────────────────────────────────────────────────────
    ## cols(
    ##   CHR = col_double(),
    ##   F = col_double(),
    ##   SNP = col_character(),
    ##   BP = col_double(),
    ##   P = col_double(),
    ##   TOTAL = col_double(),
    ##   NSIG = col_double(),
    ##   S05 = col_double(),
    ##   S01 = col_double(),
    ##   S001 = col_double(),
    ##   S0001 = col_double(),
    ##   SP2 = col_character()
    ## )
    ## 
    ## 
    ## ── Column specification ──────────────────────────────────────────────────────────────────────────────────────────────────────────────────
    ## cols(
    ##   CHR = col_double(),
    ##   F = col_double(),
    ##   SNP = col_character(),
    ##   BP = col_double(),
    ##   P = col_double(),
    ##   TOTAL = col_double(),
    ##   NSIG = col_double(),
    ##   S05 = col_double(),
    ##   S01 = col_double(),
    ##   S001 = col_double(),
    ##   S0001 = col_double(),
    ##   SP2 = col_character()
    ## )
    ## 
    ## 
    ## ── Column specification ──────────────────────────────────────────────────────────────────────────────────────────────────────────────────
    ## cols(
    ##   CHR = col_double(),
    ##   F = col_double(),
    ##   SNP = col_character(),
    ##   BP = col_double(),
    ##   P = col_double(),
    ##   TOTAL = col_double(),
    ##   NSIG = col_double(),
    ##   S05 = col_double(),
    ##   S01 = col_double(),
    ##   S001 = col_double(),
    ##   S0001 = col_double(),
    ##   SP2 = col_character()
    ## )

``` r
ranges_paths <- sumstats_paths <- lapply(metas, function(.x) str_c(.x, "clumped.ranges", sep = "."))

clumped_ranges <- lapply(ranges_paths[sapply(ranges_paths, file.exists)], read_table)
```

    ## 
    ## ── Column specification ──────────────────────────────────────────────────────────────────────────────────────────────────────────────────
    ## cols(
    ##   CHR = col_double(),
    ##   SNP = col_character(),
    ##   P = col_double(),
    ##   N = col_double(),
    ##   POS = col_character(),
    ##   KB = col_double(),
    ##   RANGES = col_character()
    ## )

    ## 
    ## ── Column specification ──────────────────────────────────────────────────────────────────────────────────────────────────────────────────
    ## cols(
    ##   CHR = col_double(),
    ##   SNP = col_character(),
    ##   P = col_double(),
    ##   N = col_double(),
    ##   POS = col_character(),
    ##   KB = col_double(),
    ##   RANGES = col_character()
    ## )
    ## 
    ## 
    ## ── Column specification ──────────────────────────────────────────────────────────────────────────────────────────────────────────────────
    ## cols(
    ##   CHR = col_double(),
    ##   SNP = col_character(),
    ##   P = col_double(),
    ##   N = col_double(),
    ##   POS = col_character(),
    ##   KB = col_double(),
    ##   RANGES = col_character()
    ## )
    ## 
    ## 
    ## ── Column specification ──────────────────────────────────────────────────────────────────────────────────────────────────────────────────
    ## cols(
    ##   CHR = col_double(),
    ##   SNP = col_character(),
    ##   P = col_double(),
    ##   N = col_double(),
    ##   POS = col_character(),
    ##   KB = col_double(),
    ##   RANGES = col_character()
    ## )
    ## 
    ## 
    ## ── Column specification ──────────────────────────────────────────────────────────────────────────────────────────────────────────────────
    ## cols(
    ##   CHR = col_double(),
    ##   SNP = col_character(),
    ##   P = col_double(),
    ##   N = col_double(),
    ##   POS = col_character(),
    ##   KB = col_double(),
    ##   RANGES = col_character()
    ## )
    ## 
    ## 
    ## ── Column specification ──────────────────────────────────────────────────────────────────────────────────────────────────────────────────
    ## cols(
    ##   CHR = col_double(),
    ##   SNP = col_character(),
    ##   P = col_double(),
    ##   N = col_double(),
    ##   POS = col_character(),
    ##   KB = col_double(),
    ##   RANGES = col_character()
    ## )

## Manhattan plot

``` r
sumstats2 <- lapply(sumstats, function(ss) filter(ss, P <= 1e-2))
manhattan(sumstats2, legend_labels = names(sumstats), ntop = 6, sign_thresh = 5e-08, build = 38)
```

![](fixed_files/figure-gfm/manhattan-1.png)<!-- -->

## Loci

Construct genomic ranges.

``` r
clumped_ranges_grs <- lapply(clumped_ranges, function(cr) {
    cr |>
    mutate(range = str_match(POS, "chr[0-9]+:([0-9]+)\\.\\.([0-9]+)")) |>
    transmute(seqnames = str_c("chr", CHR),
              start = as.numeric(range[,2]),
              end = as.numeric(range[,3]), P, SNP) |>
    filter(P <= 5e-8) |>
    as_granges() |>
    set_genome_info(genome = "hg38")
})

loci_ranges_grs <- lapply(clumped_ranges_grs, function(gr) {
    reduce_ranges(gr, SNPs = c(SNP), Ps = c(P))
})
```

Count number of loci

``` r
sapply(loci_ranges_grs, length)
```

    ##  N06A-EAS  N06A-EUR N06AA-EUR N06AB-AFR N06AB-EUR N06AB-SAS 
    ##         1        51         3         1         2         1

## Overlaps

``` r
all_gr <- reduce_ranges(bind_ranges(loci_ranges_grs))

hits_upset <- lapply(loci_ranges_grs, function(gr) findOverlaps(all_gr, gr)@from)

upset(fromList(hits_upset), nsets = length(hits_upset), order.by='freq', text.scale=2)
```

![](fixed_files/figure-gfm/overlaps-1.png)<!-- -->

## Regions

``` r
regionplot(sumstats["N06A-EUR"], legend_labels = names(sumstats["N06A-EUR"]),
    chr = 11, xmin = 112956144, xmax = 113578846,
    build = 38, show_overview = FALSE)
```

    ## [1] "Zoomed to region:  11:112956144-113578846"

![](fixed_files/figure-gfm/chr11_11-1.png)<!-- -->
