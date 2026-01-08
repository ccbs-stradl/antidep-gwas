/* 
  Allele frequency checks
*/

params.vcf = "results/vcf/gwas/GRCh38/*.{vcf.gz,vcf.gz.tbi,csv}"
params.metasets = "metasets/*.csv"

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

}
