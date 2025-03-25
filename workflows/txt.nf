/* Convert VCF sumstats to plain text with SNP, A1, A2, BETA/OR, N columns
   for LDSC/GenomicSEM */
   
params.sumstats = null
params.format = "cojo"
params.merge_alleles = "reference/w_hm3.snplist"
params.out = ""

workflow {

  // get pairs of .vcf.gz/vcf.gz.tbi
  VCF_CH = Channel.fromFilePairs(params.sumstats, size: -1) { it -> it.simpleName }
  
  if(params.format == "cojo") {
    TXT_CH = COJO(VCF_CH)
  } else if(params.format == "ldsc") {
    MERGE_CH = Channel.fromPath(params.merge_alleles)
    TXT_CH = TXT(VCF_CH)
      .combine(MERGE_CH)
    MUNGE_CH = MUNGE(TXT_CH)
  } else {
    TXT_CH = TXT(VCF_CH)
  }
  
}

/* output cojo text from gwasvcf 
*/
process COJO {
  tag "${dataset}"
  label 'tools'
  
  publishDir 'results/txt/cojo/${params.out}', mode: 'copy'
  
  input:
  tuple val(dataset), path(vcf)
  
  output:
  path("*.txt.gz")
  
  shell:
  '''
  echo "SNP A1 A2 freq b se p n" > !{dataset}.txt
  bcftools query \
    -f "%ID %ALT %REF [%AFCON] [%ES] [%SE] [%LP] [%NE]\\n" \
    !{dataset}.vcf.gz | awk  '{print $1, $2, $3, $4, $5, $6, 10^(-$7), $8}' >> !{dataset}.txt
  cat !{dataset}.txt | awk -v OFS='\t' '{print $1, $2, $3, $4, $5, $6, $7, $8}' | gzip -c  > !{dataset}.txt.gz
  '''
}

/* output plain text from gwasvcf 
*/
process TXT {
  tag "${dataset}"
  label 'tools'
  
  publishDir "results/txt/${params.out}"
  
  input:
  tuple val(dataset), path(vcf)
  
  output:
  tuple val(dataset), path("*.txt")
  
  shell:
  '''
  echo "SNP\tA1\tA2\tOR\tP\tINFO\tFRQ\tN" > !{dataset}.txt
  bcftools norm -m- --output-type u !{dataset}.vcf.gz |\
  bcftools query \
  -f "%ID\\t%ALT\\t%REF\\t[%ES]\\t[%LP]\\t[%SI]\\t[%AFCON]\\t[%NE]\\n" | awk -v OFS='\t' '{if($6 == ".") $6 = 1; print $1, $2, $3,  exp($4), 10^-($5), $6, $7, $8}' >> !{dataset}.txt
  '''
}

/* munge sumstats in LDSC
*/
process MUNGE {
  tag "${dataset}"
  label 'ldsc'

  cpu = 1
  memory = 4.GB
  
  publishDir "results/txt/munged/${params.out}", mode: 'copy'
  
  input:
  tuple val(dataset), path(sumstats), path(merge)
  
  output:
  tuple val(dataset), path("${dataset}.sumstats.gz")
  
  script:
  """
  munge_sumstats.py \
    --sumstats ${sumstats} \
    --merge-alleles ${merge} \
    --out ${dataset}
  """
}