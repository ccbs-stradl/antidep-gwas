# Documentation for Reproducing GWAS Meta-analysis of Antidepressant Exposure Pipeline

This document provides step-by-step instructions for setting up the environment, downloading necessary files, and running the various Nextflow pipelines to reproduce the GWAS meta-analysis of antidepressant exposure.

## Important things to note



## Requirements

See [list of workflow process requirements](workflow_requirements.md) and [reference files](workflow_references.md).

## Workflow steps

### 0. Setup

Specify configuration file and working directory

```sh
config=config/eddie.config
workdir=/exports/eddie/scratch/${USER}/ad/work/format
```

[`eddie.config`](config/eddie.config) is a config file used for running scripts on the University of Edinburgh's High Performance Computer (called EDDIE).
Users running these Nextflow scripts outside of Eddie may need to create a new `.config` file to suit their computational environment. 


### 1. Format antidepressant exposure summary statistics

Run the Nextflow pipeline to format the GWAS summary statistics.

```sh
nextflow run workflows/format_gwas.nf -resume \
-work-dir $workdir \
-c $config
```

### 2. Convert hg38 summary statistics to VCF

Convert the formatted summary statistics to VCF format for genome build 38.

```sh
nextflow run workflows/vcf.nf -resume \
--sumstats "format/gwas/GRCh38/*.{txt,json,csv}" \
--chr 'chr#' \
--publish "vcf/gwas/GRCh38" \
--assembly "reference/Homo_sapiens_assembly38.{fasta,fasta.fai,dict}" \
--dbsnp "reference/dbsnp.v153.hg38.vcf.{gz,gz.tbi}" \
-work-dir $workdir \
-c $config
```

### 3. Convert hg19 summary statistics to VCF

Convert the formatted summary statistics to VCF format for genome build 19.

```sh
nextflow run workflows/vcf.nf -resume \
--sumstats "format/gwas/GRCh37/*.{txt,json,csv}" \
--chr '#' \
--publish "vcf/gwas/GRCh37" \
--assembly "reference/human_g1k_v37.{fasta,fasta.fai,dict}" \
--dbsnp "reference/dbsnp.v153.b37.vcf.{gz,gz.tbi}" \
-work-dir $workdir \
-c $config
```

### 4. Liftover hg19 VCFs to hg38

Liftover the VCF files from genome build 19 to genome build 38.

```sh
nextflow run workflows/liftover.nf -resume \
--sumstats "vcf/gwas/GRCh37/*.{vcf.gz,vcf.gz.tbi}" \
--source "reference/human_g1k_v37.{fasta,fasta.fai}" \
--destination "reference/Homo_sapiens_assembly38.{fasta,fasta.fai}" \
--chain reference/hg19ToHg38.over.chain.gz \
--publish liftover/gwas \
-work-dir $workdir \
-c $config
```

### 5. Move lifted over files

Move the lifted over files to the appropriate directory and rename them.

```sh
mv liftover/gwas/*.Homo_sapiens_assembly38.vcf.* vcf/gwas/GRCh38
rename Homo_sapiens_assembly38.vcf vcf vcf/gwas/GRCh38/*.Homo_sapiens_assembly38.vcf.*
cp vcf/gwas/GRCh37/*.csv vcf/gwas/GRCh38/
```

### 6. QC

Check case/control allele frequencies across all the sumstats

```sh
nextflow run workflows/af.nf -resume \
  -work-dir $workdir \
  -c $config
```

## Meta-analysis steps

### 7. Run meta-analysis

Run the meta-analysis pipeline.

```sh
nextflow run workflows/meta.nf -resume \
-work-dir $workdir \
-c $config
```

### 8. Convert meta-analysis results to VCF

Convert the meta-analysis results to VCF format.

```sh
nextflow run workflows/format_meta.nf -resume \
-work-dir $workdir \
-c $config
```

```sh
nextflow run workflows/vcf.nf -resume \
--sumstats "format/meta/GRCh38/*.{txt,json,csv}" \
--chr 'chr#' \
--publish "vcf/meta/GRCh38" \
--assembly "reference/Homo_sapiens_assembly38.{fasta,fasta.fai,dict}" \
--dbsnp "reference/dbsnp.v153.hg38.vcf.{gz,gz.tbi}" \
-work-dir $workdir \
-c $config -with-trace
```

### 9. Liftover hg38 VCFs to hg19

Liftover the meta-analysis VCF files from genome build 38 to genome build 19 to match reference builds for some downstream analyses.

```sh
nextflow run workflows/liftover.nf -resume \
--sumstats "vcf/meta/GRCh38/*.{vcf.gz,vcf.gz.tbi}" \
--source "reference/Homo_sapiens_assembly38.{fasta,fasta.fai}" \
--destination "reference/human_g1k_v37.{fasta,fasta.fai}" \
--chain reference/hg38ToHg19.over.chain.gz \
--publish liftover/meta \
-work-dir $workdir \
-c $config 
```

### 10. Move lifted over files and copy sidecar files

Move the lifted over files to the appropriate directory and rename them.

```sh
mv liftover/meta/*-fixed-*.human_g1k_v37.vcf.gz.* vcf/meta/GRCh37
rename human_g1k_v37.vcf vcf vcf/meta/GRCh38/*.human_g1k_v37.vcf.*
cp vcf/meta/GRCh37/*.csv vcf/meta/GRCh38/
```

## Downstream analyses

### 11. Run mBAT-combo on build hg19/GRCh37

Run the mBAT-combo pipeline on genome build 19.

```sh
nextflow run workflows/genes.nf -resume \
-work-dir $workdir \
-c $config \
--build 'hg19'
```

### 12. Run mBAT-combo on build hg38/GRCh38

Run the mBAT-combo pipeline on genome build 38.

```sh
nextflow run workflows/genes.nf -resume \
-work-dir $workdir \
-c $config \
--build 'hg38'
```

### 13. Run popcorn

Run the popcorn pipeline on GWAS sumstats

```sh
nextflow run workflows/workflows/popcorn.nf -resume \
--vcf "results/vcf/gwas/GRCh38/antidep-2501-fixed-*.{csv,json,vcf.gz,vcf.gz.tbi}" \
--output "gwas" \
-work-dir $workdir \
-c $config
```

Run the popcorn pipeline on fixed effects meta sumstats

```sh
nextflow run workflows/workflows/popcorn.nf -resume \
--vcf "results/vcf/meta/GRCh38/*.{csv,json,vcf.gz,vcf.gz.tbi}" \
--output "meta" \
-work-dir $workdir \
-c $config
```

### 14. Convert sumstats to text format

Convert the summary statistics to text format.

```sh
nextflow run workflows/txt.nf -resume \
-work-dir $workdir \
-c $config \
--sumstats "liftover/*.{vcf.gz,vcf.gz.tbi}"
```

### 15. Run SuSiEx on build hg19

Run the SuSiEx pipeline on genome build 19.

```sh
nextflow run workflows/fine_mapping.nf -resume \
-work-dir $workdir \
-c $config -with-dag fineMapping/fine_mapping_dag.png
```

### 16. Plot results of SuSiEx

Plot the results of SuSiEx using R.

```sh
Rscript scripts/fine_mapping_plots.R
```

### 17. Run LDSC

#### Between GWAS

Munge the sumstats
```sh
nextflow run workflows/workflows/txt.nf -resume \
 --sumstats "results/vcf/gwas/GRCh38/*.{vcf.gz,vcf.gz.tbi}" \
 --format ldsc --out gwas \
-work-dir $workdir \
-c $config
```

Estimate LDSC genetic correlations within each cluster
```sh
clusters=("AFR" "AMR" "EAS" "EUR" "SAS")
refs=("AFR" "AMR" "EAS" "EUR" "CSA")
for i in $(seq 0 4); do
  CLUSTER=${clusters[$i]}
  REF=${refs[$i]}
  nextflow run workflows/workflows/ldsc.nf -resume \
  --source "results/txt/munged/gwas/*${CLUSTER}*.sumstats.gz" \
  --target "results/txt/munged/gwas/*${CLUSTER}*.sumstats.gz" \
  --w_ld_chr "reference/UKBB.ALL.ldscore/UKBB.${REF}" \
  --out gwas \
  -work-dir $workdir \
  -c $config
done
```

#### Between meta

Munge the sumstats
```sh
nextflow run workflows/workflows/txt.nf -resume \
 --sumstats "results/vcf/meta/GRCh38/antidep-2501-fixed-*.{vcf.gz,vcf.gz.tbi}" \
 --format ldsc --out meta \
-work-dir $workdir \
-c $config
```

Estimate LDSC genetic correlations within each cluster
```sh
clusters=("AFR" "AMR" "EAS" "EUR" "SAS")
refs=("AFR" "AMR" "EAS" "EUR" "CSA")
for i in $(seq 1 5); do
  CLUSTER=${clusters[$i]}
  REF=${refs[$i]}
  nextflow run workflows/workflows/ldsc.nf -resume \
  --source "results/txt/munged/meta/*${CLUSTER}*.sumstats.gz" \
  --target "results/txt/munged/meta/*${CLUSTER}*.sumstats.gz" \
  --w_ld_chr "reference/UKBB.ALL.ldscore/UKBB.${REF}" \
  --out meta \
  -work-dir $workdir \
  -c $config
done
```

#### With external phenotypes

```sh
nextflow run workflows/workflows/ldsc.nf -resume \
--source "results/txt/munged/meta/*EUR.sumstats.gz" \
--target "reference/munged/EUR/*.sumstats.gz" \
--w_ld_chr "reference/UKBB.ALL.ldscore/UKBB.EUR" \
--out meta/external \
-work-dir $workdir \
-c $config
```

Compile tables together
```sh
Rscript manuscript/scripts/tables_rg_ldsc_meta.R
```

### 18. Prepare sumstats for drug targetor

Prepare the sumstats into a format used as input for drug targetor
```sh
nextflow run workflows/format_gwas_drugtar.nf -resume \
-work-dir $workdir \
-c $config -with-dag fineMapping/fine_mapping_dag.png
```

### 19. Clumps and forest plots

Get SNP list from clumps and finemapping
```sh
cd scripts
Rscript -e "rmarkdown::render('fixed.Rmd')"
Rscript -e "rmarkdown::render('multi.Rmd')"
Rscript -e "rmarkdown::render('forest_snplist.Rmd')"
cd ..
```

Run workflow to make plots.
```sh
nextflow run workflows/workflows/forest.nf -resume \
work-dir $workdir \
-c $config
```