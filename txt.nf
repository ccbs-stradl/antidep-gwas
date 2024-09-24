/* Convert VCF sumstats to plain text with SNP, A1, A2, BETA/OR, N columns
   for LDSC/GenomicSEM */
   
params.sumstats = null

workflow {

  // get pairs of .vcf.gz/vcf.gz.tbi
  VCF_CH = Channel.fromFilePairs(params.sumstats, size: 2)
    
  TXT_CH = TXT(VCF_CH)
  
}

/* output plain text from gwasvcf 
*/
process TXT {
  tag "${vcf}"
  
  publishDir 'txt', mode: 'copy'
  
  input:
  tuple val(vcf), path(vcftbi)
  
  output:
  path("*.txt")
  
  shell:
  '''
  echo "SNP\tA1\tA2\tOR\tP\tINFO\tFRQ\tN" > !{vcf}.txt
  bcftools query \
  -f "%ID\\t%ALT\\t%REF\\t[%ES]\\t[%LP]\\t[%SI]\\t[%AFCON]\\t[%NE]\\n" \
  !{vcf}.gz | awk -v OFS='\t' '{print $1, $2, $3,  exp($4), 10^-($5), $6, $7, $8}' >> !{vcf}.txt
  '''
}