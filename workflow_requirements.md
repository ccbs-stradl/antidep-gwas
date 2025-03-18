# Workflow requirements

Workflow processes use the '[label](https://www.nextflow.io/docs/latest/reference/process.html#process-label)' directive to specify their requirements, which can then be supplied by the configuration file.

## label 'analysis'

- Used by: fine_mapping.nf, format_gwas.nf, format_meta.nf, meta.nf
- Requirements: R 4.4, Python 3.12
- R packages: dplyr, readr, stringr, [fastman](https://github.com/kaustubhad/fastman)
- Python packages: polars

## label 'gatk'

- Used by: None
- Requirments: [gatk](https://gatk.broadinstitute.org/hc/en-us)

## label 'gwas2vcf'

- Used by: vcf.nf
- Requirments: [gwas2vcf](https://github.com/MRCIEU/gwas2vcf)
- Conda: gwas2vcf.yaml

Clone repository:
```sh
mkdir vendor
git clone https://github.com/MRCIEU/gwas2vcf.git vendor/gwas2vcf
```

## label 'hail'

- Used by: swath.nf
- Requirments: [Hail](https://hail.is)

## label 'ldsc'

- Used by: ldsc.nf, txt.nf
- Requirements: [ldsc](https://github.com/bulik/ldsc) or [ldsc python3 port](https://github.com/belowlab/ldsc) (recommended)
- Conda: env/ldsc.yaml
```

## label 'mrmega'


- Used by: meta.nf
- Requirments: [MR-MEGA](https://genomics.ut.ee/en/tools)

This program can be downloaded, compiled, and placed in the `bin` directory.

```sh
mkdir workflows/bin
wget https://tools.gi.ut.ee/tools/MR-MEGA_v0.2.zip -d MR-MEGA
cd MR-MEGA
make
cd ..
mv MR-MEGA/MR-MEGA workflows/bin/
```

## label 'plink2'

- Used by: genes.nf
- Requirements: plink2

## label 'popcorn'

- Used by: popcorn.nf
- Reqruirements: [popcorn](https://github.com/brielin/Popcorn)
- Conda: popcorn.yaml

## label 'rscript'

- Used by: genes.nf vcf.nf
- Requirements: R 4.4
- Packages: cowplot, data.table, dplyr, ggplot2, tidyr, purr, stringr, [susiexR](https://github.com/AmeliaES/susiexR), [plyranges](https://tidyomics.github.io/plyranges/), [topr](https://github.com/totajuliusd/topr), UpSetR, readr, ieugwasr

## label 'tools'

- Used by: fine_mapping.nf, genes.nf, liftover.nf, meta.nf, popcorn.nf, txt.nf, vcf.nf
- Requirements: bcftools 1.20, samtools 1.9

The [bcftools score plugins](https://github.com/freeseek/score) are also required:

```sh
mkdir plugins
curl -O https://software.broadinstitute.org/software/score/score_1.20-20240505.zip
unzip -d plugins score_1.20-20240505.zip
```