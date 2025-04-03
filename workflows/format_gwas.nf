/* 
  Format cohort-level sumstats to be ready for GWAS VCF
  
  Inputs:
    - inputs/datasets.csv: meta data file with cohort, phenotype, version, build for each dataset
    - inputs/sumstats/cohort-dataset-version.gz: original GWAS sumstats
    - inputs/sumstats/cohort-version.sh: shell script to reformat sumstats to the following columns:		
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
    
*/

import groovy.json.JsonOutput

// Input files: sumstats gz, formatting scripts, and meta data csv
params.sumstats = "inputs/sumstats/*.gz"
params.scripts = "inputs/sumstats/*.sh"
params.datasets = "inputs/datasets.csv"

workflow {

/* 
  Input sumstats
*/

  // cohort-level sumstats
  SUMSTATS_CH = Channel.fromPath(params.sumstats)
    .map { it -> [it.name, it] }
  
  // formatting scripts
  SCRIPTS_CH = Channel.fromPath(params.scripts)
    .map { it -> [it.baseName, it] }

  // sumstats datasets information, keyed to filekey column
  CSV_CH = Channel
    .fromPath(params.datasets)
  DATASETS_CH = CSV_CH
    .splitCsv(header: true)

  // merge sumstats, scripts, and meta data information
  SUMSTATS_DATA_CH = DATASETS_CH
    .map { it -> [it.dataset, it] }
    .join(SUMSTATS_CH)
    .map { it -> ["${it[1].cohort}-${it[1].version}", it[1], it[2]]}
    .combine(SCRIPTS_CH, by: 0)
    .map { it -> ["${it[1].cohort}-${it[1].pheno}-${it[1].cluster}-${it[1].version}"].plus(it) }
    .map { it -> it.plus(JsonOutput.toJson(it[2])) }
    
  // run original sumstats through its reformatting script
  FORMAT_CH = FORMAT(SUMSTATS_DATA_CH)

  // datasets information table
  TABLE(SUMSTATS_DATA_CH.combine(CSV_CH))
}

// Reformat sumstats for GWASVCF input
process FORMAT {
  tag "${dataset}"
  
  publishDir "results/format/gwas/${meta.build}"

  cpus = 1
  memory = 1.GB
  time = '10m'

  input:
  tuple val(dataset), val(cohortversion), val(meta), path(sumstats), path(script), val(metajson)

  output:
  tuple val(dataset), val(meta), path("${dataset}.txt"), path("${dataset}.json")

  shell:
  """
  sh ${script} ${sumstats} ${dataset}.txt
  cat <<EOF > ${dataset}.json
  ${metajson}
  EOF
  """
}

// Dataset information table
process TABLE {
  tag "${dataset}"
  label 'analysis'

  publishDir "results/format/gwas/${meta.build}"

  cpus = 1
  memory = 1.GB
  time = '10m'

  input:
  tuple val(dataset), val(cohortversion), val(meta), path(sumstats), path(script), val(metajson), path(datasets)

  output:
  tuple val(dataset), val(meta), path("${dataset}.csv")

  script:
  """
  #!Rscript

  library(dplyr)
  library(readr)

  datasets <- read_csv("${datasets}")

  ds <- datasets |>
    filter(dataset == "${meta.dataset}") |>
    mutate(neff = 4 / (1/cases + 1/controls))

  write_csv(ds, "${dataset}.csv")
  """
}