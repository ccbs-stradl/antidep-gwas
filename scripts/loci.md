Antidepressant exposure GWAS meta-analysis loci
================

``` r
library(dplyr)
library(here)
library(readr)
library(stringr)
library(tidyr)
library(topr)
library(UpSetR)
library(plyranges)
library(ggplot2)
library(readxl)
library(easylift)
```

``` r
metaset <- "antidep-2501"
```

## Clumps

Collate all clump files together and mark which meta-analysis they are
from.

``` r
clumps_paths <- list.files(
  here::here("results", "meta", metaset),
  str_c(metaset, ".+\\.clumps\\.tsv"),
  full.names = TRUE
)
clump_prefixes <- str_remove(basename(clumps_paths), ".clumps.tsv")

pheno_cluster <- sapply(str_split(clump_prefixes, pattern = "-"), function(x) str_c(x[4], x[5], sep = "-"))
pheno_cluster <- str_replace(pheno_cluster, "DIV", "MRMEGA")

names(clumps_paths) <- pheno_cluster

clumps <- lapply(clumps_paths, read_table, col_types = cols(REF = col_character(), ALT = col_character()))
```

    ## Warning: The following named parsers don't match the column names: REF, ALT
    ## Warning: The following named parsers don't match the column names: REF, ALT
    ## Warning: The following named parsers don't match the column names: REF, ALT

## Genomics ranges

Construct genomic ranges.

``` r
clumped_ranges_grs <- lapply(clumps, function(cr) {
  cr |>
    select(
      CHROM = matches("chrom"),
      start,
      end,
      ID = matches("^(ID|MarkerName)$")) |>
    mutate(
      CHROM = if_else(
        str_detect(CHROM, "chr"),
        true = as.character(CHROM),
        false = str_c("chr", CHROM)
      )
    ) |>
    as_granges(seqnames = CHROM, start = start, end = end) |>
    set_genome_info(genome = "hg38")
})
```

Count number of loci

``` r
all_datasets <- names(clumped_ranges_grs)
all_n06a <- str_subset(names(clumped_ranges_grs), "N06A-")
all_n06aa <- str_subset(names(clumped_ranges_grs), "N06AA-")
all_n06ab <- str_subset(names(clumped_ranges_grs), "N06AB-")
fixed_datasets <- str_subset(names(clumped_ranges_grs), "MEGA", negate = TRUE)
fixed_n06a_datasets <- str_subset(names(clumped_ranges_grs), "N06A-[^M]")
fixed_n06aa_datasets <- str_subset(names(clumped_ranges_grs), "N06AA-[^M]")
fixed_n06ab_datasets <- str_subset(names(clumped_ranges_grs), "N06AB-[^M]")
mega_datasets <- str_subset(names(clumped_ranges_grs), "MEGA")
mega_n06a_datasets <- str_subset(names(clumped_ranges_grs), "N06A-MRMEGA")
mega_n06aa_datasets <- str_subset(names(clumped_ranges_grs), "N06AA-MRMEGA")
mega_n06ab_datasets <- str_subset(names(clumped_ranges_grs), "N06AB-MRMEGA")

datasets_lists <- list(
  all = all_datasets,
  all_n06a = all_n06a,
  all_n06ab = all_n06aa,
  all_n06ab = all_n06ab,
  fixed = fixed_datasets,
  fixed_n06a = fixed_n06a_datasets,
  fixed_n06aa = fixed_n06aa_datasets,
  fixed_n06ab = fixed_n06ab_datasets,
  mega = mega_datasets,
  mega_n06a = mega_n06a_datasets,
  mega_n06aa = mega_n06aa_datasets,
  mega_n06ab = mega_n06ab_datasets
)
```

``` r
sapply(datasets_lists, function(x) sum(sapply(clumped_ranges_grs[x], length)))
```

    ##         all    all_n06a   all_n06ab   all_n06ab       fixed  fixed_n06a 
    ##         111          88           9          14          77          62 
    ## fixed_n06aa fixed_n06ab        mega   mega_n06a  mega_n06aa  mega_n06ab 
    ##           7           7          34          26           2           6

## Overlaps

``` r
all_gr <- reduce_ranges(bind_ranges(clumped_ranges_grs))

hits_upset <- lapply(clumped_ranges_grs, function(gr) findOverlaps(all_gr, gr)@from)

upset(fromList(hits_upset), nsets = length(hits_upset), order.by = "freq", text.scale = 2)
```

![](loci_files/figure-gfm/overlaps-1.png)<!-- -->

Number of unique loci total

``` r
sapply(datasets_lists, function(x) length(reduce_ranges(bind_ranges(clumped_ranges_grs[x]))))
```

    ##         all    all_n06a   all_n06ab   all_n06ab       fixed  fixed_n06a 
    ##          80          64           9          11          74          62 
    ## fixed_n06aa fixed_n06ab        mega   mega_n06a  mega_n06aa  mega_n06ab 
    ##           7           7          33          26           2           6

SNPs unique to MR-MEGA analyses

``` r
filter_by_non_overlaps(
  clumped_ranges_grs[['N06AB-MRMEGA']],
  clumped_ranges_grs[['N06A-EUR']]
)
```

    ## GRanges object with 5 ranges and 1 metadata column:
    ##       seqnames              ranges strand |          ID
    ##          <Rle>           <IRanges>  <Rle> | <character>
    ##   [1]     chr2     3390659-3492791      * |  rs74766704
    ##   [2]     chr3 158112929-158799345      * |   rs1656375
    ##   [3]    chr16   83443927-83509611      * |   rs8045989
    ##   [4]    chr23   70448613-71021855      * |   rs6624508
    ##   [5]    chr23 143889504-144077530      * | rs113661605
    ##   -------
    ##   seqinfo: 5 sequences from hg38 genome; no seqlengths

``` r
filter_by_non_overlaps(
  clumped_ranges_grs[['N06AB-MRMEGA']],
  bind_ranges(
    clumped_ranges_grs[['N06A-EUR']],
    clumped_ranges_grs[['N06AB-EUR']]
  )
)
```

    ## GRanges object with 3 ranges and 1 metadata column:
    ##       seqnames              ranges strand |          ID
    ##          <Rle>           <IRanges>  <Rle> | <character>
    ##   [1]     chr2     3390659-3492791      * |  rs74766704
    ##   [2]    chr23   70448613-71021855      * |   rs6624508
    ##   [3]    chr23 143889504-144077530      * | rs113661605
    ##   -------
    ##   seqinfo: 5 sequences from hg38 genome; no seqlengths

## MDD2025

Compare to PGC MDD

``` r
mdd2025_all <- read_excel(here::here("manuscript/reference/MDD_GWAS_Online_Results_(COJO).xlsx"), sheet = 2) |>
  mutate(P = as.numeric(P))

mdd2025_eur <- read_excel(here::here("manuscript/reference/MDD_GWAS_Online_Results_(COJO).xlsx"), sheet = 3)

mdd2025_list <- list(mdd2025_all = mdd2025_all, mdd2025_eur = mdd2025_eur)

mdd2025_hg19_grs <- lapply(mdd2025_list, function(cr) {
  cr |>
    transmute(
      CHROM = if_else(CHR == 23, true = "chrX", false = str_c("chr", CHR)),
      start = Locus.Left,
      end = Locus.Right,
      SNP
    ) |>
    as_granges(seqnames = CHROM, start = start, end = end) |>
    set_genome_info(genome = "hg19")
})

chain <- here::here("reference/hg19ToHg38.over.chain.gz")

mdd2025_grs <- lapply(mdd2025_hg19_grs, function(gr) easylift(gr, "hg38", chain))

mdd2025_gr <- reduce_ranges(bind_ranges(mdd2025_grs))
```

Overlaps

``` r
length(filter_by_overlaps(all_gr, mdd2025_gr))
```

    ## [1] 50
