# GWAS Meta-analysis of Antidepressant Prescribing

Cohorts 

- Biobank Japan, [ATC_N06A: Antidepressants](https://pheweb.jp/pheno/ATC_N06A)
- FINNGEN, [R9 Antidepressants](https://r9.finngen.fi/pheno/ANTIDEPRESSANTS)
- Generation Scotland, PIS

## 1. GWAS VCF

Convert sumstats to [GWAS VCF](https://github.com/MRCIEU/gwas-vcf-specification) format. 

```sh
nextflow run gwasvcf.nf --sumstats sumstats/*.{gz,sh} -resume
```

