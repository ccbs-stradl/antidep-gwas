/* 
  Format cohort-level sumstats to be ready for GWAS VCF
  
  Inputs:
    - sumstats.csv: meta data file with cohort, phenotype, version, build for each dataset
    - sumstats/cohort-dataset-version.gz: original GWAS sumstats
    - sumstats/cohort-version.sh: shell script to reformat sumstats to the following columns:		
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
    
*/

import groovy.json.JsonOutput

// Input files: sumstats gz, formatting scripts, and meta data csv
params.sumstats = "sumstats/*.gz"
params.scripts = "sumstats/*.sh"
params.meta = "sumstats.csv"

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

  // sumstats meta data, keyed to filekey column
  META_CH = Channel
    .fromPath(params.meta)
    .splitCsv(header: true)

  // merge sumstats, scripts, and meta data information
  SUMSTATS_META_CH = META_CH
    .map { it -> [it.dataset, it] }
    .join(SUMSTATS_CH)
    .map { it -> ["${it[1].cohort}-${it[1].version}", it[1], it[2]]}
    .combine(SCRIPTS_CH, by: 0)
    .map { it -> ["${it[1].cohort}-${it[1].pheno}-${it[1].cluster}-${it[1].version}-${it[1].build}"].plus(it) }
    .map { it -> it.plus(JsonOutput.toJson(it[2])) }
    //.view()
    
  // run original sumstats through its reformatting script
  FORMAT_CH = FORMAT(SUMSTATS_META_CH)
}

// Reformat sumstats for GWASVCF input
process FORMAT {
  tag "${dataset}"
  
  publishDir "format/gwas"

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