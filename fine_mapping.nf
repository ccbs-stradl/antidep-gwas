nextflow.enable.dsl=2

/*
  Fine-mapping of GWAS meta using SuSiEx (hg19 only)

  Summary of processes and inputs for each process:
    * make bim, bed, fam files for each ancestry:
      - pfile 
      - IDs where for each ancestry (output from make-pgen.sh)
    * calculate effective sample size for each meta-analysis sumstats:
      - format/meta/GRCh38/antidep-2408-fixed-N06A-{ANCESTRY}.json (output from format_meta.nf)
    * process meta sumstats into correct format:
      - meta-analysed sumstats: meta/metaset-method-phenotype-cluster.csv (output from format_meta.nf)
    * clumping to identify regions for fine mapping:
      - processed sumstats (from previous process in current script)
      - pfile, per chr and ancestry, (output from make-pgen_hg19.sh)
    * load clumping results into R and determine region boundaries:
      - output from clumping process, per chr and ancestry (from previous process in current script)
    * run SuSiEx on each region:
      - output from R determining region boundaries (from previous process in current script)
      - processed sumstats (from process above)
    * explore results:
      - use susiexR package in R to format and plot output from SuSiEx

*/

params.pfile = '/exports/igmm/eddie/GenScotDepression/data/ukb/genetics/impv3_pgen/ukb_imp_v3.qc'
params.ancestry_ids = 'reference/ukb-ld-ref_ancestry.id'
params.cluster = 'EUR'
params.chr = 21

workflow {

  PFILE_CH = Channel.of(params.pfile)          
  IDS_CH = Channel.fromPath(params.ancestry_ids)
  CLUSTER_CH = Channel.of(params.cluster)      
  CHR_CH = Channel.of(params.chr)    

  REF_CH = PFILE_CH
    .combine(IDS_CH)
    .combine(CLUSTER_CH)
    .combine(CHR_CH)

  BFILE_CH = MAKE_BFILE(REF_CH)

}

process MAKE_BFILE {
  tag "$cluster:chr$chr"

  cpus = 1
  memory = 8.GB
  time = '1h'

  input:
    tuple val(pfile), path(ancestry_ids), val(cluster), val(chr)

  output:
    path "ukb_imp_v3.qc.geno02.mind02_${cluster}_${chr}.*"
 
  script:
  """
  plink2 \
    --pfile  ${pfile} \
    --keep-col-match ${ancestry_ids} ${cluster} \
    --chr ${chr} \
    --geno 0.02 \
    --maf 0.005 \
    --make-pgen 'vzs' \
    --out ukb_imp_v3.qc.geno02_${cluster}_${chr} \
    --threads 4

    plink2 \
    --pfile ukb_imp_v3.qc.geno02_${cluster}_${chr} 'vzs' \
    --mind 0.02 \
    --make-bed \
    --out ukb_imp_v3.qc.geno02.mind02_${cluster}_${chr} \
    --threads 4
  """

}




