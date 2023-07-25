# GWAS Meta-analysis of Antidepressant Prescribing

Cohorts 

- Biobank Japan, [ATC_N06A: Antidepressants](https://pheweb.jp/pheno/ATC_N06A)
- FINNGEN, [R9 Antidepressants](https://r9.finngen.fi/pheno/ANTIDEPRESSANTS)
- Generation Scotland, PIS

### Requirements

- [Singularity](https://docs.sylabs.io/)

## 0. Download reference files

```sh
mkdir reference
curl "https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.{fasta,fasta.fai,dict}" -o "reference/Homo_sapiens_assembly38.#1"
curl "http://fileserve.mrcieu.ac.uk/dbsnp/dbsnp.v153.hg38.vcf.{gz,gz.tbi}" -o "reference/dbsnp.v153.hg38.vcf.#1"
curl "https://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz" -o "reference/hg19ToHg38.over.chain.gz"
```

## 1. GWAS VCF

Install `gwas2vcf`(https://mrcieu.github.io/gwas2vcf/)

```sh
mkdir vendor
git clone https://github.com/MRCIEU/gwas2vcf.git vendor/gwas2vcf
```


Convert sumstats to [GWAS VCF](https://github.com/MRCIEU/gwas-vcf-specification) format. 

```sh
nextflow run gwasvcf.nf --sumstats "sumstats/*.{gz,sh}" -resume
```

