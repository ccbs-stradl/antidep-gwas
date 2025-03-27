/* 
  Format meta-analysed sumstats to be ready for GWAS VCF
  
  Inputs:
    - meta: meta/metaset-method-phenotype-cluster.{csv,gz}: Meta-analysed sumstats
 
  Outputs:
    - .txt file with the following columns
      - chr
      - pos
      - ea
      - oa
      - beta
      - se
      - pval
      - ncase
      - ncontrol
      - snp
      - eaf
      - imp_info
      - eaf_case
      - eaf_control
      - neff
    - .{json,csv} files with metaset information  
*/

import groovy.json.JsonOutput

// Input files: sumstats gz and metaset.csv
// just get main meta-analysis files that have a .gz and a .csv
params.sumstats = "meta/*.{csv,gz}"
params.build = "GRCh38"

workflow {

/* 
  Input sumstats
*/

  // meta-analysed sumstats
  SUMSTATS_CH = Channel.fromFilePairs(params.sumstats, size: 2)
  
  // split into fixed and mrmega 
  METHOD_CH = SUMSTATS_CH
    .branch {
      fixed: it[0] =~ /fixed/
      mrmega: it[0] =~ /mrmega/
    }
    
  // format fixed effects sumstats
  FORMAT_FIXED_CH = FORMAT_FIXED(METHOD_CH.fixed)
  
  // json sidecar file
  JSON_CH = FORMAT_JSON(FORMAT_FIXED_CH)
  
}

// Reformat sumstats for GWASVCF input
process FORMAT_FIXED {
  tag "${dataset}"
  
  publishDir "results/format/meta/${params.build}"

  cpus = 1
  memory = 1.GB
  time = '10m'

  input:
  tuple val(dataset), path(sumstats)

  output:
  tuple val(dataset), path("${dataset}.txt"), path("${dataset}.csv", includeInputs: true)

  shell:
  '''
  # fixed effects meta-analysis columns
  #  1	CHROM
  #  2	POS
  #  3	ID
  #  4	REF
  #  5	ALT
  #  6	studies
  #  7	BETA
  #  8	SE
  #  9	CHISQ
  # 10	P
  # 11	Q
  # 12	QP
  # 13	INFO
  # 14	AFCAS
  # 15	AFCON
  # 16	NCAS
  # 17	NCON
  # 18	NEFF
  # 19	NTOT
  
  gunzip -c !{dataset}.gz | awk 'OFS = "\t" {if(NR == 1) {print "chr", "pos", "ea", "oa", "beta", "se", "pval", "ncase", "ncontrol", "snp", "eaf", "imp_info", "eaf_case", "eaf_control", "neff"} else {print $1, $2, $5, $4, $7, $8, $10, $16, $17, $3, $15, $13, $14, $15, $18}}' > !{dataset}.txt
  
  # output columns
  # chr
  # pos
  # ea
  # oa
  # beta
  # se
  # pval
  # ncase
  # ncontrol
  # snp
  # eaf
  # imp_info
  # eaf_case
  # eaf_control
  # neff
  '''
}

// Metaset json
process FORMAT_JSON {
  tag "${dataset}"
  label 'analysis'

  publishDir "results/format/meta/${params.build}"

  cpus = 1
  memory = 1.GB
  time = '10m'

  input:
  tuple val(dataset), path(sumstats), path(csv)

  output:
  tuple val(dataset), path("${dataset}.json")

  script:
  """
  #!Rscript

  library(jsonlite)
  library(dplyr)
  library(readr)
  library(stringr)

  datasets <- read_csv("${csv}")
  summary <- datasets |>
    summarize(across(cases:neff, sum))
  
  info <- str_split("${dataset}", pattern = "-")[[1]]
  
  metaset <- list(cohort = "${dataset}",
               pheno = info[4],
               dataset = "${dataset}.gz",
               version = str_glue("{info[2]}-{info[3]}"),
               build = "${params.build}",
               cluster = info[5],
               cases = pull(summary, cases),
               controls = pull(summary, controls),
               neff = pull(summary, neff))
               
  # treat data as a singleton to avoid values getting embedded in arrays             
  write_json(unbox(as.data.frame(metaset)), "${dataset}.json")

  """
}
