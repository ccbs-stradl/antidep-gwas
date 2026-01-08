/* 
  Allele frequency checks
*/

params.vcf = "results/vcf/gwas/GRCh38/*.{vcf.gz,vcf.gz.tbi,csv}"
params.metasets = "inputs/metasets/*.csv"

workflow {
  /*
      Inputs
  */
  
    // sumstats, filenames as COHORT-PHENO-CLUSTER-VERSION
    VCF_CH = Channel
      .fromFilePairs(params.vcf, size: 3)
      .map { it -> [it[0].split("-")].plus(it) }
      .map { it -> [[cohort: it[0][0], pheno: it[0][1], cluster: it[0][2], version: it[0][3]], it[1], it[2]] }
      
    
    // dataset lists for each meta-analysis
    METASET_CH = Channel
      .fromPath(params.metasets)
      .map { it -> [it.baseName, it.splitCsv(header: true)] }
      .transpose()
      .map { it -> it[1].plus([metaset: it[0]]) }

  /*
    Pre-processing
  */
    
    // select datasets that are used by a metaset
    // key to cohort and version
    VCF_CV_CH = VCF_CH
      .map { it -> [it[0].cohort, it[0].version].plus(it) }
    METASET_CV_CH = METASET_CH
      .map { it -> [it.cohort, it.version].plus(it) }
      
    // merge datasets with metasets, then gather distinct datasets
    VCF_IN_METASET_CH = METASET_CV_CH
      .combine(VCF_CV_CH, by: [0, 1])
      .map { it -> [it[3], it[4], it[5], it[2]]}
      .groupTuple(by: [0, 1, 2])

    // remove metaset information for input processing
    VCF_IN_CH = VCF_IN_METASET_CH
      .map { it -> it[0..2]}

    // track metasets
    IN_META_CH = VCF_IN_METASET_CH
      .map { it ->  it[0,3] }
        
  /*
    Extract and convert
  */

    // get allele frequency colums from VCF
    AF_CH = AF_IN(VCF_IN_CH)
    // convert tables to parquet
    AF_PARQUET_CH = PARQUET(AF_CH)

/*
  Allele frequency correlations
*/
  // process inputs into metasets
  AF_IN_CH = AF_PARQUET_CH
    .combine(IN_META_CH, by: 0)
    .transpose(by: 4)
    .map { it -> [[metaset: it[4].metaset]].plus(it) }
    .groupTuple(by: 0)
    // do not keep last element, which is the metaset tracker dictionary
    .map { it -> it[0..4] }

  CORR_CH = CORR(AF_IN_CH)

}

/* Query VCFs to extract allele frequency information
*/

process AF_IN {
  tag "${dataset}"
  label 'tools'

  cpus 1
  memory 1.GB
  time 30.m

  input:
  tuple val(details), val(dataset), path(vcf)

  output:
  tuple val(details), val(dataset), path("${dataset}.txt.gz"), path("${dataset}.csv", includeInputs: true)

  script:
  """
  # bcf query to get required columns
  bcftools norm -Ou -m- ${dataset}.vcf.gz |\
  bcftools query \
  -H -f '%CHROM %POS %ID %ALT %REF [%AFCAS] [%AFCON]\n' |\
  awk -v OFS='\t' '\$1 = \$1' |\
  gzip -c > ${dataset}.txt.gz
  """


}

/*
Convert tables to arrow/parquet
*/

process PARQUET {
  tag "${dataset}"
  label 'analysis'

  cpus 1
  memory 8.GB
  time 10.m

  input:
  tuple val(details), val(dataset), path(afs), path(csv)

  output:
  tuple val(details), val(dataset), path("${dataset}.parquet"), path(csv, includeInputs: true)

  script:
  """
  #!python3
  import os

  os.environ["POLARS_MAX_THREADS"] = "${task.cpus}"
  import polars as pl

  afs_paths = "${afs}"

  afs = pl.scan_csv(afs_paths, separator = "\\t", new_columns = ["CHROM", "POS", "ID", "ALT", "REF", "${dataset}.AFCAS", "${dataset}.AFCON"])

  afs.sink_parquet("${dataset}.parquet")
  """
}

/*
  Calculate correlations between allele frequencies
*/
process CORR {
  tag "${metaset.metaset}"
  label 'analysis'

  publishDir "manuscript/tables", mode: 'copy'

  cpus 8
  memory 31.GB
  time 30.m

  input:
  tuple val(metaset), val(details), val(datasets), path(afs), path(csvs)

  output:
  tuple val(metaset), path("${metaset.metaset}.corr.csv")

  script:
  """
  #!python3
  import os

  os.environ["POLARS_MAX_THREADS"] = "${task.cpus}"
  import polars as pl
  from functools import reduce

  afs_paths = "${afs}"
  variant_ids = ["CHROM", "POS", "ID", "ALT", "REF"]

  afs_list = [pl.scan_parquet(path) for path in afs_paths.split()]

  # join all data frame together
  afs = reduce(
    lambda left, right: left.join(right, on = ["CHROM", "POS", "ID", "ALT", "REF"], how = "inner"),
    afs_list
  )

  # select allele frequency columns
  af_cas_con = afs.select(pl.exclude(["CHROM", "POS", "ID", "ALT", "REF"]))

  af_corr = af_cas_con.collect(engine = "streaming").corr()

  af_corr.write_csv("${metaset.metaset}.corr.csv")
  """
}
