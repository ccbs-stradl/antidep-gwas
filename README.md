# GWAS Meta-analysis of Antidepressant Prescribing

Cohorts 

- Biobank Japan, [ATC_N06A: Antidepressants](https://pheweb.jp/pheno/ATC_N06A)
- FINNGEN, [R9 Antidepressants](https://r9.finngen.fi/pheno/ANTIDEPRESSANTS)
- Generation Scotland, PIS

## 0. Download reference files

```sh
mkdir reference
curl "https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.{fasta,fasta.fai,dict}" -o "reference/Homo_sapiens_assembly38.#1"
curl "https://ftp.ncbi.nih.gov/snp/organisms/human_9606_b151_GRCh38p7/VCF/common_all_20180418.vcf.{gz,gz.tbi}" -o "reference/common_all_20180418.vcf.#1"
```

## 1. GWAS VCF

Convert sumstats to [GWAS VCF](https://github.com/MRCIEU/gwas-vcf-specification) format. 

```sh
nextflow run gwasvcf.nf --sumstats sumstats/*.{gz,sh} -resume
```

