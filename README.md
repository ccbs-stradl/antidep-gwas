# GWAS Meta-analysis of Antidepressant Prescribing

Cohorts 

- Biobank Japan, [ATC_N06A: Antidepressants](https://pheweb.jp/pheno/ATC_N06A)
- FINNGEN, [R9 Antidepressants](https://r9.finngen.fi/pheno/ANTIDEPRESSANTS)
- Generation Scotland, PIS
- UK Biobank, GP prescription records
- All of Us, Drug exposure

### Requirements

- [Singularity](https://docs.sylabs.io/)

## 0. Download reference files

```sh
mkdir reference
curl "https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.{fasta,fasta.fai,dict}" -o "reference/Homo_sapiens_assembly38.#1"
curl "http://fileserve.mrcieu.ac.uk/dbsnp/dbsnp.v153.hg38.vcf.{gz,gz.tbi}" -o "reference/dbsnp.v153.hg38.vcf.#1"
curl "http://fileserve.mrcieu.ac.uk/dbsnp/dbsnp.v153.b37.vcf.{gz,gz.tbi}" -o "reference/dbsnp.v153.b378.vcf.#1"
curl "http://fileserve.mrcieu.ac.uk/ref/2.8/b37/human_g1k_v37.fasta.gz" | gunzip -c > reference/human_g1k_v37.fasta
curl "http://fileserve.mrcieu.ac.uk/ref/2.8/b37/human_g1k_v37.{fasta.fai,dict}" -o "reference/human_g1k_v37.#1"
curl "https://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz" -o "reference/hg19ToHg38.over.chain.gz"
curl "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/liftOver/hg38ToHg19.over.chain.gz" -o "reference/hg38ToHg19.over.chain.gz"
curl "https://raw.githubusercontent.com/Share-AL-work/mBAT/main/glist_ensgid_hg38_v40.txt" -o "reference/glist_ensgid_hg38_v40.txt"
curl "https://raw.githubusercontent.com/Share-AL-work/mBAT/main/glist_ensgid_hg38_v40_symbol_gene_names.txt" -o "reference/glist_ensgid_hg38_v40_symbol_gene_names.txt"
curl -L "https://www.dropbox.com/s/j72j6uciq5zuzii/all_hg38.pgen.zst?dl=1" -o reference/all_hg38.pgen.zst
plink2 --zst-decompress reference/all_hg38.pgen.zst > reference/all_hg38.pgen
rm reference/all_hg38.pgen.zst
curl -L "https://www.dropbox.com/s/ngbo2xm5ojw9koy/all_hg38_noannot.pvar.zst?dl=1" -o reference/all_hg38.pvar.zst
curl -L "https://www.dropbox.com/s/2e87z6nc4qexjjm/hg38_corrected.psam?dl=1" -o reference/all_hg38.psam
```

## 1. Format GWAS

Harmonise input sumstats formatting.

```sh
nextflow run format_gwas.nf -resume
```

## 2. GWAS VCF

Install `gwas2vcf`(https://mrcieu.github.io/gwas2vcf/)

```sh
mkdir vendor
git clone https://github.com/MRCIEU/gwas2vcf.git vendor/gwas2vcf
```

Install [bcftools score plugins](https://github.com/freeseek/score)

```sh
mkdir plugins
curl -O https://software.broadinstitute.org/software/score/score_1.20-20240505.zip
unzip -d plugins score_1.20-20240505.zip
```

Convert sumstats to [GWAS VCF](https://github.com/MRCIEU/gwas-vcf-specification) format. 

Genome build 38:
```sh
nextflow run vcf.nf \
--sumstats "format/gwas/*-GRCh38.{txt,json}" \
--assembly "reference/Homo_sapiens_assembly38.{fasta,fasta.fai,dict}" \
--dbsnp "reference/dbsnp.v153.hg38.vcf.{gz,gz.tbi}" \
-resume
```

The workflow also requires [bcftools](https://samtools.github.io/bcftools/bcftools.html), [gatk](https://gatk.broadinstitute.org), R, and [plyranges](https://sa-lee.github.io/plyranges/index.html).

### 2. Multi-ancestry meta-analysis

Meta-analyse multi-ancestry sumstats using [MR-MEGA](https://genomics.ut.ee/en/tools)

```sh
mkdir bin
wget https://tools.gi.ut.ee/tools/MR-MEGA_v0.2.zip -d MR-MEGA
cd MR-MEGA
make
cd ..
mv MR-MEGA/MR-MEGA bin/
```

Run workflow
```sh
nextflow run meta.nf -resume
```

### 3. Hail

Work with GWASVCF files in [Hail](https://hail.is). Merge VCFs and import to `MatrixTable`

```sh
nextflow run swath.nf -resume --name N06A
```
