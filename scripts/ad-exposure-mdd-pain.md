``` r
library(readr)
library(dplyr)
library(stringr)
library(GenomicSEM)

build_all <- params$build_all
```

## Sumstats

- pain: “Back pain”, Sakaue et al
  [10.1038/s41588-021-00931-x](https://doi.org/10.1038/s41588-021-00931-x),
  [hum0197-v3-220](https://humandbs.dbcls.jp/en/hum0197-v3-220).
- depression: “Major depression”, Meng et al
  [10.1038/s41588-023-01596-4](https://www.nature.com/articles/s41588-023-01596-4),
  [mdd2023diverse](https://figshare.com/articles/dataset/mdd2023diverse/24799299);
  Adams et al
  [10.1016/j.cell.2024.12.002](https://doi.org/10.1016/j.cell.2024.12.002),
  [mdd2025](https://doi.org/10.6084/m9.figshare.27061255.v4).

### Backpain

Download

``` sh
mkdir -p reference/sumstats

curl https://humandbs.dbcls.jp/files/hum0197/hum0197.v3.BBJ.BP.v1.zip > reference/sumstats/hum0197.v3.BBJ.BP.v1.zip
curl https://humandbs.dbcls.jp/files/hum0197/hum0197.v3.EUR.BP.v1.zip > reference/sumstats/hum0197.v3.EUR.BP.v1.zip

unzip reference/sumstats/hum0197.v3.BBJ.BP.v1.zip -d reference/sumstats/
unzip reference/sumstats/hum0197.v3.EUR.BP.v1.zip -d reference/sumstats/
```

Get RSIDs

``` sh
gunzip -c reference/sumstats/hum0197.v3.BBJ.BP.v1/GWASsummary_Back_pain_Japanese_SakaueKanai2020.auto.txt.gz | awk 'OFS = "\t" {if(NR > 1) {print $2, $3}}' > reference/sumstats/back_pain_eas.tsv

gunzip -c reference/sumstats/hum0197.v3.EUR.BP.v1/GWASsummary_Back_pain_EUR_SakaueKanai2020.auto.txt.gz | awk 'OFS = "\t" {if(NR > 1) {print $2, $3}}' > reference/sumstats/back_pain_eur.tsv

cat reference/sumstats/back_pain_eas.tsv reference/sumstats/back_pain_eur.tsv | sort | uniq > reference/sumstats/back_pain_targets.tsv

# extract targets from the dbSNP VCF file and print chromosome, position, SNP ID, and alleles
bcftools view \
--targets-file reference/sumstats/back_pain_targets.tsv \
--output-type b \
reference/dbsnp.v153.b37.vcf.gz |\
bcftools query --print-header \
--format '%CHROM\t%POS\t%ID\t%REF\t%ALT{0}\n' > reference/sumstats/back_pain_chr_pos_rsid.tsv
```

Read in sumstats and SNP list

``` r
hum0197_backpain_eas <- read_tsv(here::here("reference/sumstats/hum0197.v3.BBJ.BP.v1/GWASsummary_Back_pain_Japanese_SakaueKanai2020.auto.txt.gz"))

hum0197_backpain_eur <- read_tsv(here::here("reference/sumstats/hum0197.v3.EUR.BP.v1/GWASsummary_Back_pain_EUR_SakaueKanai2020.auto.txt.gz"))

back_pain_chr_pos_rsid <- read_tsv(here::here("reference/sumstats/back_pain_chr_pos_rsid.tsv"))
```

Update IDs and approximate effective sample size

``` r
hum0197_backpain_eas_cases <- 1732
hum0197_backpain_eas_controls <- 176994
hum0197_backpain_eas_tot <- hum0197_backpain_eas_cases + hum0197_backpain_eas_controls
hum0197_backpain_eas_neff <- 4 / (1/hum0197_backpain_eas_cases + 1/hum0197_backpain_eas_controls)

hum0197_backpain_eas_sumstats <-
hum0197_backpain_eas |>
  mutate(n_count = AC_Allele2/2) |>
  mutate(n_pct = n_count / hum0197_backpain_eas_tot) |>
  mutate(N = n_pct * hum0197_backpain_eas_neff) |>
  transmute(SNP = SNPID, A1 = Allele2, A2 = Allele1, INFO = imputationInfo, FRQ = AF_Allele2, OR = exp(BETA), P = p.value, N)
```

``` r
hum0197_backpain_eur_ukb_cases <- 14893
hum0197_backpain_eur_ukb_controls <- 342765
hum0197_backpain_eur_ukb_tot <- hum0197_backpain_eur_ukb_cases + hum0197_backpain_eur_ukb_controls
hum0197_backpain_eur_ukb_neff <- 4 / (1/hum0197_backpain_eur_ukb_cases + 1/hum0197_backpain_eur_ukb_controls)

hum0197_backpain_eur_finngen_cases <- 7520   
hum0197_backpain_eur_finngen_controls <- 103091
hum0197_backpain_eur_finngen_tot <- hum0197_backpain_eur_finngen_cases + hum0197_backpain_eur_finngen_controls
hum0197_backpain_eur_finngen_neff <- 4 / (1/hum0197_backpain_eur_finngen_cases + 1/hum0197_backpain_eur_finngen_controls)

hum0197_backpain_eur_sumstats <-
hum0197_backpain_eur |>
  mutate(CHR = as.character(CHR)) |>
  left_join(back_pain_chr_pos_rsid, by = c("CHR" = "#[1]CHROM", "POS" = "[2]POS"),
            relationship = "many-to-many") |>
  filter((Allele1 == `[4]REF` & Allele2 ==  `[5]ALT`) |
         (Allele2 == `[4]REF` & Allele1 ==  `[5]ALT`)) |>
  filter(!duplicated(`[3]ID`)) |>
  mutate(cohorts = case_when(!is.na(AF_Allele2_UKB) & !is.na(AF_Allele2_FG) ~ "both",
                              is.na(AF_Allele2_UKB) & !is.na(AF_Allele2_FG) ~ "fg",
                             !is.na(AF_Allele2_UKB) & is.na(AF_Allele2_FG) ~ "ukb",
                             .default = NA_character_)) |>
  mutate(N = case_match(cohorts,
                         "both" ~ hum0197_backpain_eur_ukb_neff + hum0197_backpain_eur_finngen_neff,
                         "fg" ~ hum0197_backpain_eur_finngen_neff,
                         "ukb" ~ hum0197_backpain_eur_ukb_neff,
                         .default = NA_real_
                       ),
        FRQ = case_match(cohorts,
           "both" ~ (AF_Allele2_UKB * hum0197_backpain_eur_ukb_tot + AF_Allele2_FG * hum0197_backpain_eur_finngen_tot) / (hum0197_backpain_eur_ukb_tot + hum0197_backpain_eur_finngen_tot),
           "fg" ~ AF_Allele2_FG,
           "ukb" ~ AF_Allele2_UKB,
           .default = NA_real_
         )) |>
  transmute(SNP = coalesce(`[3]ID`, v), A1 = Allele2, A2 = Allele1, FRQ = FRQ, OR = exp(BETA), P = p.value, N = N)
```

``` r
write_tsv(hum0197_backpain_eas_sumstats, here::here("reference/sumstats/back_pain_SakaueKanai2020_eas.txt"))
write_tsv(hum0197_backpain_eur_sumstats, here::here("reference/sumstats/back_pain_SakaueKanai2020_eur.txt"))
```

### Major depression

``` sh
curl -L https://figshare.com/ndownloader/files/51487019 > reference/sumstats/pgc-mdd2025_no23andMe_eur_v3-49-24-11.tsv.gz

curl -L https://figshare.com/ndownloader/files/43621287 > reference/sumstats/mdd2023diverse_AFR_Neff.csv

curl -L https://figshare.com/ndownloader/files/43621266 > reference/sumstats/mdd2023diverse_EAS_Neff.csv

curl -L https://figshare.com/ndownloader/files/43621275 > reference/sumstats/mdd2023diverse_HIS_Neff.csv

curl -L https://figshare.com/ndownloader/files/43621242 > reference/sumstats/mdd2023diverse_SAS_Neff.csv
```

``` r
adams_md_eur <- read_tsv(here::here("reference/sumstats/pgc-mdd2025_no23andMe_eur_v3-49-24-11.tsv.gz"), comment = "##")

adams_md_eur_sumstats <- adams_md_eur |>
  transmute(SNP = ID, A1 = EA, A2 = NEA, FRQ = FCON, INFO = IMPINFO, OR = exp(BETA), P = PVAL, N = NEFF)
  
write_tsv(adams_md_eur_sumstats, here::here("reference/sumstats/md_adams_eur.txt"))
```

``` r
for(cluster in c("AFR", "EAS", "HIS", "SAS")) {
  meng_md <- read_csv(here::here(str_glue("reference/sumstats/mdd2023diverse_{cluster}_Neff.csv")))

  meng_md_sumstats <- meng_md |>
    transmute(SNP = coalesce(SNP, str_glue("{Chromosome}:{Position}:{EA}:{NEA}")),
              A1 = str_to_upper(EA), A2 = str_to_upper(NEA),
              FRQ = EAF, OR = exp(logOR), P, N = Neff)
              
  write_tsv(meng_md_sumstats, here::here(str_glue("reference/sumstats/md_meng_{cluster}.txt")))
}
```

## Munge

``` r
# hm3 SNPs
hm3 <- here::here("reference/w_hm3.snplist")

for(cluster in c("AFR", "AMR", "EAS", "EUR", "SAS")) {
  munge(files = here::here(str_glue("results/txt/meta/antidep-2501-fixed-N06A-{cluster}.txt")),
       trait.names = str_glue("AD_{cluster}"), hm3 = hm3,
       info.filter = 0.9, maf.filter = 0.01, N = NA,
       column.names = list(MAF = "FRQ"))
}
```

``` r
for(cluster in c("eas", "eur")) {
  munge(files = here::here(str_glue("reference/sumstats/back_pain_SakaueKanai2020_{cluster}.txt")),
       trait.names = str_glue("Pain_{str_to_upper(cluster)}"), hm3 = hm3,
       info.filter = 0.9, maf.filter = 0.01, N = NA,
       column.names = list(MAF = "FRQ"))
}
```

``` r
for(cluster in c("AFR", "HIS", "EAS", "SAS")) {
  munge(files = here::here(str_glue("reference/sumstats/md_meng_{cluster}.txt")),
       trait.names = str_glue("MD_{cluster}"), hm3 = hm3,
       info.filter = 0.9, maf.filter = 0.01, N = NA,
       column.names = list(MAF = "FRQ"))
}

munge(files = here::here("reference/sumstats/md_adams_eur.txt"),
  trait.names = "MD_EUR",
  hm3 = hm3, info.filter = 0.9, maf.filter = 0.01, N = NA,
  column.names = list(MAF = "FRQ"))
```

## LDSC

### LDSC reference

``` sh
curl -L https://pan-ukb-us-east-1.s3.amazonaws.com/ld_release/UKBB.ALL.ldscore.tar.gz > reference/UKBB.ALL.ldscore.tar.gz
tar xzf reference/UKBB.ALL.ldscore.tar.gz
for cluster in AFR AMR CSA EAS EUR; do
  mkdir reference/UKBB.ALL.ldscore/UKBB.${cluster}
  cp reference/UKBB.ALL.ldscore/UKBB.${cluster}.rsid.l2.ldscore.gz reference/UKBB.ALL.ldscore/UKBB.${cluster}/1.l2.ldscore.gz
  cp reference/UKBB.ALL.ldscore/UKBB.${cluster}.l2.M reference/UKBB.ALL.ldscore/UKBB.${cluster}/1.l2.M
  cp reference/UKBB.ALL.ldscore/UKBB.${cluster}.l2.M_5_50 reference/UKBB.ALL.ldscore/UKBB.${cluster}/1.l2.M_5_50
done
```

### Estimate genetic covariances

``` r
# if building, estimate all covariances and save out as R code.
# otherwise load them in
if(build_all) {
  dir.create(here::here("scripts/ad-exposure-mdd_files"))
  afr_ldsc <- ldsc(traits = c("AD_AFR.sumstats.gz", "MD_AFR.sumstats.gz"),
  trait.names = c("AD", "MD"),
  sample.prev = c(0.5, 0.5, 0.5),
  population.prev = c(0.25, 0.15),
  ld = here::here("reference/UKBB.ALL.ldscore/UKBB.AFR"),
  wld = here::here("reference/UKBB.ALL.ldscore/UKBB.AFR"),
  chr = 1
  )
  dput(afr_ldsc, here::here("scripts/ad-exposure-mdd_files/afr_ldsc.R"), control=c('all', 'digits17'))
  
  amr_ldsc <- ldsc(traits = c("AD_AMR.sumstats.gz", "MD_HIS.sumstats.gz"),
    trait.names = c("AD", "MD"),
    sample.prev = c(0.5, 0.5, 0.5),
    population.prev = c(0.25, 0.15),
    ld = here::here("reference/UKBB.ALL.ldscore/UKBB.AMR"),
    wld = here::here("reference/UKBB.ALL.ldscore/UKBB.AMR"),
    chr = 1
    )
  dput(amr_ldsc, here::here("scripts/ad-exposure-mdd_files/amr_ldsc.R"), control=c('all', 'digits17'))

  eas_ldsc <- ldsc(traits = c("AD_EAS.sumstats.gz", "MD_EAS.sumstats.gz", "Pain_EAS.sumstats.gz"),
                  trait.names = c("AD", "MD", "Pain"),
                  sample.prev = c(0.5, 0.5, 0.5),
                  population.prev = c(0.02, 0.15, 0.05),
                  ld = here::here("reference/UKBB.ALL.ldscore/UKBB.EAS"),
                  wld = here::here("reference/UKBB.ALL.ldscore/UKBB.EAS"),
                  chr = 1
                  )
  dput(eas_ldsc, here::here("scripts/ad-exposure-mdd_files/eas_ldsc.R"), control=c('all', 'digits17'))
                  
  eur_ldsc <- ldsc(traits = c("AD_EUR.sumstats.gz", "MD_EUR.sumstats.gz", "Pain_EUR.sumstats.gz"),
                    trait.names = c("AD", "MD", "Pain"),
                    sample.prev = c(0.5, 0.5, 0.5),
                    population.prev = c(0.25, 0.15, 0.05),
                    ld = here::here("reference/UKBB.ALL.ldscore/UKBB.EUR"),
                    wld = here::here("reference/UKBB.ALL.ldscore/UKBB.EUR"),
                    chr = 1
                    )
  dput(eur_ldsc, here::here("scripts/ad-exposure-mdd_files/eur_ldsc.R"), control=c('all', 'digits17'))

  sas_ldsc <- ldsc(traits = c("AD_SAS.sumstats.gz", "MD_SAS.sumstats.gz"),
                    trait.names = c("AD", "MD"),
                    sample.prev = c(0.5, 0.5, 0.5),
                    population.prev = c(0.25, 0.15),
                    ld = here::here("reference/UKBB.ALL.ldscore/UKBB.CSA"),
                    wld = here::here("reference/UKBB.ALL.ldscore/UKBB.SAS"),
                    chr = 1
                    )
  dput(sas_ldsc, here::here("scripts/ad-exposure-mdd_files/sas_ldsc.R"), control=c('all', 'digits17'))
} else {
  afr_ldsc <- dget(here::here("scripts/ad-exposure-mdd_files/afr_ldsc.R"))
  amr_ldsc <- dget(here::here("scripts/ad-exposure-mdd_files/amr_ldsc.R"))
  eas_ldsc <- dget(here::here("scripts/ad-exposure-mdd_files/eas_ldsc.R"))
  eur_ldsc <- dget(here::here("scripts/ad-exposure-mdd_files/eur_ldsc.R"))
  sas_ldsc <- dget(here::here("scripts/ad-exposure-mdd_files/sas_ldsc.R"))

}
```

## Model

Cholesky decomposition of MD and AD in EAS cohorts.

``` r
chol2_model <- "
F1 =~ NA*MD + AD
F2 =~ NA*AD

F1 ~~ 1*F1
F2 ~~ 1*F2
F1 ~~ 0*F2

MD ~~ 0*MD
MD ~~ 0*AD
AD ~~ 0*AD
"

chol2_amr_fit <- usermodel(amr_ldsc, estimation = "DWLS", model = chol2_model)
```

    ## [1] "Running primary model"
    ## [1] "Calculating CFI"
    ## [1] "Calculating Standardized Results"
    ## [1] "Calculating SRMR"
    ## elapsed 
    ##   0.162 
    ## [1] "Model fit statistics are all printed as NA as you have specified a fully saturated model (i.e., df = 0)"

``` r
chol2_amr_fit
```

    ## $modelfit
    ##    chisq df p_chisq AIC CFI SRMR
    ## df    NA  0      NA  NA  NA   NA
    ## 
    ## $results
    ##   lhs op rhs Unstand_Est         Unstand_SE STD_Genotype   STD_Genotype_SE
    ## 3  F1 =~  MD  0.20753978 0.0776749827284036   1.00000001 0.374265513893631
    ## 2  F1 =~  AD  0.00234703 0.0948429179666066   0.01559353 0.630131105953904
    ## 6  F2 =~  AD -0.15049469 0.0762924958590567   0.99987841 0.506883130204275
    ## 4  F1 ~~  F1  1.00000000                      1.00000000                  
    ## 7  F2 ~~  F2  1.00000000                      1.00000000                  
    ## 5  F1 ~~  F2  0.00000000                      0.00000000                  
    ## 9  MD ~~  MD  0.00000000                      0.00000000                  
    ## 8  MD ~~  AD  0.00000000                      0.00000000                  
    ## 1  AD ~~  AD  0.00000000                      0.00000000                  
    ##      STD_All     p_value
    ## 3 1.00000001 0.007542314
    ## 2 0.01559353 0.980257165
    ## 6 0.99987841 0.048540985
    ## 4 1.00000000          NA
    ## 7 1.00000000          NA
    ## 5 0.00000000          NA
    ## 9 0.00000000          NA
    ## 8 0.00000000          NA
    ## 1 0.00000000          NA

``` r
chol2_sas_fit <- usermodel(sas_ldsc, estimation = "DWLS", model = chol2_model)
```

    ## [1] "Running primary model"
    ## [1] "Calculating CFI"
    ## [1] "Calculating Standardized Results"
    ## [1] "Calculating SRMR"
    ## elapsed 
    ##   0.159 
    ## [1] "Model fit statistics are all printed as NA as you have specified a fully saturated model (i.e., df = 0)"
    ## [1] "The S matrix was smoothed prior to model estimation due to a non-positive definite matrix. The largest absolute difference in a cell between the smoothed and non-smoothed matrix was  0.0289914488799009 As a result of the smoothing, the largest Z-statistic change for the genetic covariances was  0.457881788622328 . We recommend setting the smooth_check argument to true if you are going to run a multivariate GWAS."

    ## Warning in usermodel(sas_ldsc, estimation = "DWLS", model = chol2_model): A
    ## difference greater than .025 was observed pre- and post-smoothing in the
    ## genetic covariance matrix. This reflects a large difference and results should
    ## be interpreted with caution!! This can often result from including low powered
    ## traits, and you might consider removing those traits from the model. If you are
    ## going to run a multivariate GWAS we strongly recommend setting the smooth_check
    ## argument to true to check smoothing for each SNP.

    ## Warning in usermodel(sas_ldsc, estimation = "DWLS", model = chol2_model): A
    ## difference greater than .025 was observed pre- and post-smoothing for
    ## Z-statistics in the genetic covariance matrix. This reflects a large difference
    ## and results should be interpreted with caution!! This can often result from
    ## including low powered traits, and you might consider removing those traits from
    ## the model. If you are going to run a multivariate GWAS we strongly recommend
    ## setting the smooth_check argument to true to check smoothing for each SNP.

``` r
chol2_sas_fit
```

    ## $modelfit
    ##    chisq df p_chisq AIC CFI SRMR
    ## df    NA  0      NA  NA  NA   NA
    ## 
    ## $results
    ##   lhs op rhs  Unstand_Est        Unstand_SE STD_Genotype   STD_Genotype_SE
    ## 3  F1 =~  MD 0.1776404949 0.101156188872328  1.000000010 0.569443300822671
    ## 2  F1 =~  AD 0.2328960178 0.298506714776198  0.999999926  1.28171659874201
    ## 6  F2 =~  AD 0.0001276785   683.23010062648  0.000499934  3216.97612520184
    ## 4  F1 ~~  F1 1.0000000000                    1.000000000                  
    ## 7  F2 ~~  F2 1.0000000000                    1.000000000                  
    ## 5  F1 ~~  F2 0.0000000000                    0.000000000                  
    ## 9  MD ~~  MD 0.0000000000                    0.000000000                  
    ## 8  MD ~~  AD 0.0000000000                    0.000000000                  
    ## 1  AD ~~  AD 0.0000000000                    0.000000000                  
    ##       STD_All    p_value
    ## 3 1.000000010 0.07907115
    ## 2 0.999999875 0.43527104
    ## 6 0.000499934 0.99999985
    ## 4 1.000000000         NA
    ## 7 1.000000000         NA
    ## 5 0.000000000         NA
    ## 9 0.000000000         NA
    ## 8 0.000000000         NA
    ## 1 0.000000000         NA

Cholesky decomposition of pain, MD, and AD in EUR cohorts.

``` r
chol_model <- "
F1 =~ NA*Pain + MD + AD
F2 =~ NA*MD + AD
F3 =~ NA*AD

F1 ~~ 1*F1
F2 ~~ 1*F2
F3 ~~ 1*F3
F1 ~~ 0*F2 + 0*F3
F2 ~~ 0*F3

Pain ~~ 0*Pain
Pain ~~ 0*MD
Pain ~~ 0*AD
MD ~~ 0*MD
MD ~~ 0*AD
AD ~~ 0*AD
"

chol_eas_fit <- usermodel(eas_ldsc, estimation = "DWLS", model = chol_model)
```

    ## [1] "Running primary model"
    ## [1] "Calculating CFI"
    ## [1] "Calculating Standardized Results"
    ## [1] "Calculating SRMR"
    ## elapsed 
    ##   0.191 
    ## [1] "Model fit statistics are all printed as NA as you have specified a fully saturated model (i.e., df = 0)"

``` r
chol_eas_fit
```

    ## $modelfit
    ##    chisq df p_chisq AIC CFI SRMR
    ## df    NA  0      NA  NA  NA   NA
    ## 
    ## $results
    ##     lhs op  rhs Unstand_Est         Unstand_SE STD_Genotype   STD_Genotype_SE
    ## 4    F1 =~ Pain  0.62045942  0.441939747916373    1.0000000 0.712278225485746
    ## 3    F1 =~   MD  0.08419025  0.107001513994632    0.4170087 0.529996759798668
    ## 2    F1 =~   AD  0.16356716  0.156382079034663    0.9803774 0.937311914996099
    ## 9    F2 =~   MD  0.18349912 0.0598534079363947    0.9089025 0.296464141181936
    ## 8    F2 =~   AD  0.02033470  0.158731878693593    0.1218807 0.951395960466297
    ## 12   F3 =~   AD  0.02584974  0.959843664636841    0.1549365  5.75303974056239
    ## 5    F1 ~~   F1  1.00000000                       1.0000000                  
    ## 10   F2 ~~   F2  1.00000000                       1.0000000                  
    ## 13   F3 ~~   F3  1.00000000                       1.0000000                  
    ## 6    F1 ~~   F2  0.00000000                       0.0000000                  
    ## 7    F1 ~~   F3  0.00000000                       0.0000000                  
    ## 11   F2 ~~   F3  0.00000000                       0.0000000                  
    ## 18 Pain ~~ Pain  0.00000000                       0.0000000                  
    ## 17 Pain ~~   MD  0.00000000                       0.0000000                  
    ## 16 Pain ~~   AD  0.00000000                       0.0000000                  
    ## 15   MD ~~   MD  0.00000000                       0.0000000                  
    ## 14   MD ~~   AD  0.00000000                       0.0000000                  
    ## 1    AD ~~   AD  0.00000000                       0.0000000                  
    ##      STD_All     p_value
    ## 4  1.0000000 0.160335013
    ## 3  0.4170087 0.431390979
    ## 2  0.9803774 0.295586115
    ## 9  0.9089025 0.002170819
    ## 8  0.1218807 0.898064130
    ## 12 0.1549365 0.978514607
    ## 5  1.0000000          NA
    ## 10 1.0000000          NA
    ## 13 1.0000000          NA
    ## 6  0.0000000          NA
    ## 7  0.0000000          NA
    ## 11 0.0000000          NA
    ## 18 0.0000000          NA
    ## 17 0.0000000          NA
    ## 16 0.0000000          NA
    ## 15 0.0000000          NA
    ## 14 0.0000000          NA
    ## 1  0.0000000          NA

``` r
chol_eur_fit <- usermodel(eur_ldsc, estimation = "DWLS", model = chol_model)
```

    ## [1] "Running primary model"
    ## [1] "Calculating CFI"
    ## [1] "Calculating Standardized Results"
    ## [1] "Calculating SRMR"
    ## elapsed 
    ##   0.197 
    ## [1] "Model fit statistics are all printed as NA as you have specified a fully saturated model (i.e., df = 0)"

``` r
chol_eur_fit
```

    ## $modelfit
    ##    chisq df p_chisq AIC CFI SRMR
    ## df    NA  0      NA  NA  NA   NA
    ## 
    ## $results
    ##     lhs op  rhs Unstand_Est          Unstand_SE STD_Genotype    STD_Genotype_SE
    ## 4    F1 =~ Pain  0.26226517  0.0115612511360334    1.0000000 0.0440822971476503
    ## 3    F1 =~   MD  0.15750059  0.0096911725225661    0.6450561  0.039690962159948
    ## 2    F1 =~   AD  0.20384757  0.0116914873718483    0.7681157 0.0440545586853024
    ## 9    F2 =~   MD  0.18657563 0.00741848769346632    0.7641352 0.0303830020791555
    ## 8    F2 =~   AD  0.14233968  0.0122295563603865    0.5363485 0.0460820504817166
    ## 12   F3 =~   AD  0.09282028  0.0152102285419718    0.3497550 0.0573134871108457
    ## 5    F1 ~~   F1  1.00000000                        1.0000000                   
    ## 10   F2 ~~   F2  1.00000000                        1.0000000                   
    ## 13   F3 ~~   F3  1.00000000                        1.0000000                   
    ## 6    F1 ~~   F2  0.00000000                        0.0000000                   
    ## 7    F1 ~~   F3  0.00000000                        0.0000000                   
    ## 11   F2 ~~   F3  0.00000000                        0.0000000                   
    ## 18 Pain ~~ Pain  0.00000000                        0.0000000                   
    ## 17 Pain ~~   MD  0.00000000                        0.0000000                   
    ## 16 Pain ~~   AD  0.00000000                        0.0000000                   
    ## 15   MD ~~   MD  0.00000000                        0.0000000                   
    ## 14   MD ~~   AD  0.00000000                        0.0000000                   
    ## 1    AD ~~   AD  0.00000000                        0.0000000                   
    ##      STD_All       p_value
    ## 4  1.0000000 6.323180e-114
    ## 3  0.6450561  2.163200e-59
    ## 2  0.7681157  4.432361e-68
    ## 9  0.7641352 1.410028e-139
    ## 8  0.5363485  2.610904e-31
    ## 12 0.3497550  1.044279e-09
    ## 5  1.0000000            NA
    ## 10 1.0000000            NA
    ## 13 1.0000000            NA
    ## 6  0.0000000            NA
    ## 7  0.0000000            NA
    ## 11 0.0000000            NA
    ## 18 0.0000000            NA
    ## 17 0.0000000            NA
    ## 16 0.0000000            NA
    ## 15 0.0000000            NA
    ## 14 0.0000000            NA
    ## 1  0.0000000            NA
