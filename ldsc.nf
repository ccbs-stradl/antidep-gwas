/* LD Score genetic correlations
*/

params.source = "txt/munged/*.sumstats.gz"
params.target = "reference/munged/EUR/*.sumstats.gz"
params.w_ld_chr = "reference/eur_w_ld_chr"
params.out = "meta"

workflow {
  SOURCE_CH = Channel.fromPath(params.source)
  TARGET_CH = Channel.fromPath(params.target)
  W_LD_CH = Channel.fromPath(params.w_ld_chr, type: "dir")
  
  SOURCE_TARGET_CH = SOURCE_CH
    .combine(TARGET_CH)
    .filter { it-> it[0] != it[1] }
    .combine(W_LD_CH)
    
  RG_CH = LDSC(SOURCE_TARGET_CH)
}

process LDSC {
  tag "${source.simpleName}--${target.simpleName}"
  label 'ldsc'
  
  publishDir "models/rg/${params.out}", mode: 'copy'
  
  cpu = 1
  memory = 4.GB
  
  errorStrategy 'ignore'
  
  input:
  tuple path(source), path(target), path(w_ld)
  
  output:
  tuple val(source.simpleName), val(target.simpleName), path("*.log")
  
  script:
  """
  ldsc.py \
  --rg ${source},${target} \
  --ref-ld-chr ${w_ld}/ \
  --w-ld-chr ${w_ld}/ \
  --out ${source.simpleName}--${target.simpleName}
  """
  
}