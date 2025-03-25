# Resource and reference files

Download required reference files for the workflows.

First, create a directory for reference files

```sh
mkdir reference
```

## Genome builds

FASTA and dbSNP files.

```sh
curl "https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.{fasta,fasta.fai,dict}" -o "reference/Homo_sapiens_assembly38.#1"
curl "http://fileserve.mrcieu.ac.uk/dbsnp/dbsnp.v153.hg38.vcf.{gz,gz.tbi}" -o "reference/dbsnp.v153.hg38.vcf.#1"
curl "http://fileserve.mrcieu.ac.uk/dbsnp/dbsnp.v153.b37.vcf.{gz,gz.tbi}" -o "reference/dbsnp.v153.b37.vcf.#1"
curl "http://fileserve.mrcieu.ac.uk/ref/2.8/b37/human_g1k_v37.fasta.gz" | gunzip -c > reference/human_g1k_v37.fasta
curl "http://fileserve.mrcieu.ac.uk/ref/2.8/b37/human_g1k_v37.{fasta.fai,dict}" -o "reference/human_g1k_v37.#1"
```

## Liftover

Chain files.

```sh
curl "https://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz" -o "reference/hg19ToHg38.over.chain.gz"
curl "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/liftOver/hg38ToHg19.over.chain.gz" -o "reference/hg38ToHg19.over.chain.gz"
```

## Gene mapping

Gene lists.

```sh
curl "https://raw.githubusercontent.com/Share-AL-work/mBAT/main/glist_ensgid_hg38_v40.txt" -o "reference/glist_ensgid_hg38_v40.txt"
curl "https://raw.githubusercontent.com/Share-AL-work/mBAT/main/glist_ensgid_hg38_v40_symbol_gene_names.txt" -o "reference/glist_ensgid_hg38_v40_symbol_gene_names.txt"
```

## Reference genomes

Files from [plink2](https://www.cog-genomics.org/plink/2.0/resources#phase3_1kg).

```sh
curl -L "https://www.dropbox.com/s/j72j6uciq5zuzii/all_hg38.pgen.zst?dl=1" -o reference/all_hg38.pgen.zst
plink2 --zst-decompress reference/all_hg38.pgen.zst > reference/all_hg38.pgen
rm reference/all_hg38.pgen.zst
curl -L "https://www.dropbox.com/s/ngbo2xm5ojw9koy/all_hg38_noannot.pvar.zst?dl=1" -o reference/all_hg38.pvar.zst
curl -L "https://www.dropbox.com/s/2e87z6nc4qexjjm/hg38_corrected.psam?dl=1" -o reference/all_hg38.psam
```

## LD Scores

Files from [Pan UKBB](https://pan.ukbb.broadinstitute.org/downloads).

```sh
curl https://zenodo.org/records/7773502/files/w_hm3.snplist.gz -o reference/w_hm3.snplist.gz
gunzip reference/w_hm3.snplist.gz
curl "https://pan-ukb-us-east-1.s3.amazonaws.com/ld_release/UKBB.ALL.ldscore.tar.gz" -o reference/UKBB.ALL.ldscore.tar.gz
tar xzf reference/UKBB.ALL.ldscore.tar.gz -C reference
# make subdirectories for use with --ld-chr arguments
for cluster in AFR AMR CSA EAS EUR; do
  prefix=UKBB.${cluster}
  mkdir reference/UKBB.ALL.ldscore/${prefix}
  mv reference/UKBB.ALL.ldscore/${prefix}.l2.M reference/UKBB.ALL.ldscore/${prefix}/1.l2.M
  mv reference/UKBB.ALL.ldscore/${prefix}.l2.M_5_50 reference/UKBB.ALL.ldscore/${prefix}/1.l2.M_5_50
  mv reference/UKBB.ALL.ldscore/${prefix}.rsid.l2.ldscore.gz reference/UKBB.ALL.ldscore/${prefix}/1.l2.ldscore.gz
done
```