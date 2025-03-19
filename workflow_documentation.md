# Documentation for Reproducing GWAS Meta-analysis of Antidepressant Exposure Pipeline

This document provides step-by-step instructions for setting up the environment, downloading necessary files, and running the various Nextflow pipelines to reproduce the GWAS meta-analysis of antidepressant exposure.

## Important things to note

[`eddie.config`](eddie.config) is a config file used for running scripts on the University of Edinburgh's High Performance Computer (called EDDIE).
Users running these Nextflow scripts outside of Eddie may need to create a new `.config` file to suit their computational environment.

## Requirements

See [list of workflow process requirements](workflow_requirements.md) and [reference files](workflow_references.md).

## Workflow steps

### 1. Format antidepressant exposure summary statistics

Run the Nextflow pipeline to format the GWAS summary statistics.

```sh
nextflow run format_gwas.nf -resume \
-work-dir /exports/eddie/scratch/${USER}/ad/work/format \
-c eddie.config
```

### 2. Convert hg38 summary statistics to VCF

Convert the formatted summary statistics to VCF format for genome build 38.

```sh
nextflow run vcf.nf -resume \
--sumstats "format/gwas/GRCh38/*.{txt,json,csv}" \
--chr 'chr#' \
--publish "vcf/gwas/GRCh38" \
--assembly "reference/Homo_sapiens_assembly38.{fasta,fasta.fai,dict}" \
--dbsnp "reference/dbsnp.v153.hg38.vcf.{gz,gz.tbi}" \
-work-dir /exports/eddie/scratch/${USER}/ad/work/vcf \
-c eddie.config
```

### 3. Convert hg19 summary statistics to VCF

Convert the formatted summary statistics to VCF format for genome build 19.

```sh
nextflow run vcf.nf -resume \
--sumstats "format/gwas/GRCh37/*.{txt,json,csv}" \
--chr '#' \
--publish "vcf/gwas/GRCh37" \
--assembly "reference/human_g1k_v37.{fasta,fasta.fai,dict}" \
--dbsnp "reference/dbsnp.v153.b37.vcf.{gz,gz.tbi}" \
-work-dir /exports/eddie/scratch/${USER}/ad/work/vcf \
-c eddie.config
```

### 4. Liftover hg19 VCFs to hg38

Liftover the VCF files from genome build 19 to genome build 38.

```sh
nextflow run liftover.nf -resume \
--sumstats "vcf/gwas/GRCh37/*.{vcf.gz,vcf.gz.tbi}" \
--source "reference/human_g1k_v37.{fasta,fasta.fai}" \
--destination "reference/Homo_sapiens_assembly38.{fasta,fasta.fai}" \
--chain reference/hg19ToHg38.over.chain.gz \
--publish liftover/gwas \
-work-dir /exports/eddie/scratch/${USER}/ad/work/liftover \
-c eddie.config
```

### 5. Move lifted over files

Move the lifted over files to the appropriate directory and rename them.

```sh
mv liftover/gwas/*.Homo_sapiens_assembly38.vcf.* vcf/gwas/GRCh38
rename Homo_sapiens_assembly38.vcf vcf vcf/gwas/GRCh38/*.Homo_sapiens_assembly38.vcf.*
cp vcf/gwas/GRCh37/*.csv vcf/gwas/GRCh38/
```

## Meta-analysis steps

### 6. Run meta-analysis

Run the meta-analysis pipeline.

```sh
nextflow run meta.nf -resume \
-work-dir /exports/eddie/scratch/${USER}/ad/work/meta \
-c eddie.config
```

### 7. Convert meta-analysis results to VCF

Convert the meta-analysis results to VCF format.

```sh
nextflow run format_meta.nf -resume \
-work-dir /exports/eddie/scratch/${USER}/ad/work/meta \
-c eddie.config
```

```sh
nextflow run vcf.nf -resume \
--sumstats "format/meta/GRCh38/*.{txt,json,csv}" \
--chr 'chr#' \
--publish "vcf/meta/GRCh38" \
--assembly "reference/Homo_sapiens_assembly38.{fasta,fasta.fai,dict}" \
--dbsnp "reference/dbsnp.v153.hg38.vcf.{gz,gz.tbi}" \
-work-dir /exports/eddie/scratch/${USER}/ad/work/vcf \
-c eddie.config -with-trace
```

### 8. Liftover hg19 VCFs to hg38

Liftover the meta-analysis VCF files from genome build 19 to genome build 38.

```sh
nextflow run liftover.nf -resume \
--sumstats "vcf/meta/GRCh38/*.{vcf.gz,vcf.gz.tbi}" \
--source "reference/Homo_sapiens_assembly38.{fasta,fasta.fai}" \
--destination "reference/human_g1k_v37.{fasta,fasta.fai}" \
--chain reference/hg38ToHg19.over.chain.gz \
--publish liftover/meta \
-work-dir /exports/eddie/scratch/${USER}/ad/work/liftover \
-c eddie.config 
```

### 9. Move lifted over files and copy sidecar files

Move the lifted over files to the appropriate directory and rename them.

```sh
mv liftover/meta/*-fixed-*.human_g1k_v37.vcf.gz.* vcf/meta/GRCh37
rename human_g1k_v37.vcf vcf vcf/meta/GRCh38/*.human_g1k_v37.vcf.*
cp vcf/meta/GRCh37/*.csv vcf/meta/GRCh38/
```

## Downstream analyses

### 10. Run mBAT-combo on build hg19/GRCh37

Run the mBAT-combo pipeline on genome build 19.

```sh
nextflow run genes.nf -resume \
-work-dir /exports/eddie/scratch/${USER}/ad/work_hg19 \
-c eddie.config \
--build 'hg19'
```

### 11. Run mBAT-combo on build hg38/GRCh38

Run the mBAT-combo pipeline on genome build 38.

```sh
nextflow run genes.nf -resume \
-work-dir /exports/eddie/scratch/${USER}/ad/work_hg38 \
-c eddie.config \
--build 'hg38'
```

### 12. Run popcorn

Run the popcorn pipeline.

```sh
nextflow run popcorn.nf -resume \
-work-dir /exports/eddie/scratch/${USER}/ad/work_hg38 \
-c eddie.config
```

### 13. Convert sumstats to text format

Convert the summary statistics to text format.

```sh
nextflow run txt.nf -resume \
-work-dir /exports/eddie/scratch/${USER}/ad/work \
-c eddie.config \
--sumstats "liftover/*.{vcf.gz,vcf.gz.tbi}"
```

### 14. Run SuSiEx on build hg19

Run the SuSiEx pipeline on genome build 19.

```sh
nextflow run fine_mapping.nf -resume \
-work-dir /exports/eddie/scratch/${USER}/ad/work \
-c eddie.config -with-dag fineMapping/fine_mapping_dag.png
```

### 15. Plot results of SuSiEx

Plot the results of SuSiEx using R.

```sh
Rscript fine_mapping_plots.R
```
