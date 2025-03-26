/* 
  Liftover sumstats
*/

params.sumstats = null
params.source = null
params.destination = null
params.chain = null
params.plugins = "plugins"
params.publish = "liftover"

workflow {
  // VCF sumstats
  SUMSTATS_CH = Channel.fromFilePairs(params.sumstats, checkIfExists: true)
  
  // source and destination fasta files
  SOURCE_CH = Channel.fromFilePairs(params.source, checkIfExists: true)
  DEST_CH = Channel.fromFilePairs(params.destination, checkIfExists: true)
  
  // liftover chain
  CHAIN_CH = Channel.fromPath(params.chain, checkIfExists: true)
  
  // bcftools plugins
  PLUGIN_CH = Channel.fromPath(params.plugins, type: "dir", checkIfExists: true)
  
  DATA_CH = SUMSTATS_CH
    .combine(SOURCE_CH)
    .combine(DEST_CH)
    .combine(CHAIN_CH)
    .combine(PLUGIN_CH)
  
  // Perform liftover  
  LIFTED_CH = LIFTOVER(DATA_CH)
  
}

process LIFTOVER {
  tag "${sumstats}"
  label 'tools'

  publishDir "${params.publish}", mode: 'copy'
  
  cpus = 1
  memory = 8.GB
  
  input:
  tuple val(sumstats), path(vcf), val(source), path(srcfasta), val(dest), path(fasta), path(chain), path(plugins)
  
  output:
  path("*.{gz,gz.tbi}")
  
  script:
  """
  export BCFTOOLS_PLUGINS="${plugins}"
  bcftools norm -m+ --output-type u ${sumstats}.vcf.gz |\
  bcftools +liftover -- \
  --src-fasta-ref ${source}.fasta \
  --fasta-ref ${dest}.fasta \
  --chain ${chain} \
  --af-tags INFO/AF,FMT/AF,FMT/AFCAS,FMT/AFCON \
  --es-tags FMT/ES |\
  bcftools norm -m- --output-type u |\
  bcftools sort --output-type z \
  --output ${sumstats}.${dest}.vcf.gz

  tabix ${sumstats}.${dest}.vcf.gz
  """
}