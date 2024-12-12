nextflow.enable.dsl=2

/*
  Fine-mapping of GWAS meta using SuSiEx
*/

// genotyype build parameters
// first define which options are available
def valid_builds = ['hg19', 'hg38']

// default to hg19 if not provided
params.build = params.build ?: 'hg19'

// validate that the provided build is one of the valid options
if (!(params.build in valid_builds)) {
    error "Invalid value for --build. Please specify one of: ${valid_builds.join(', ')}"
}

// assign labels for process MA based on build:
def processMALabel = (params.build == 'hg19') ? 'tools' : (params.build == 'hg38') ? 'rscript' : null

def mapsDuplicatesDir = (params.build == 'hg19') ? 'maps_hg19_duplicates' : (params.build == 'hg38') ? 'maps_hg38_duplicates' : null

def mapsDir = (params.build == 'hg19') ? 'maps_hg19' : (params.build == 'hg38') ? 'maps_hg38' : null

// input files dependent on genotype build type:
if (params.build == 'hg19') {
  // files for hg19 build
  // input meta-analysis sumstats
  params.meta = "liftover/fixed-*.vcf.gz"

  // ld reference BED files */
  // Edit this to read in multiple different ancestries:
  params.ref = "reference/ukb_imp_v3.qc_ancestry.{pgen,psam,pvar.zst}"

} else if (params.build == 'hg38') {
  // files for hg38 build
  // input meta-analysis sumstats
  params.meta = "meta/fixed-*.meta.gz"

  // ld reference BED files */
  // params.ref = "reference/all_hg38.{pgen,psam,pvar.zst}"
  // Insert command to stop nextflow pipeline if 'hg38' is chosen
  // as separate LD BED files for each ancestry are not yet in "reference/"

}

// effective sample size QC parameter
params.neff_pct = 0.8

workflow {

  // sumstats from fixed effects meta 
  META_CH = Channel.fromPath(params.meta)
    .map { it -> [it.simpleName.split("-"), it] }
    .map { it -> [it[0][2], it[0][0], it[0][1], it[1]] }

  // QC parameters
  NEFF_CH = Channel.of(params.neff_pct)

  // format sumstats to .ma
  MA_CH = MA(META_CH, NEFF_CH)

  // run SuSiEx on all outputs of MA_CH when MA_CH is complete

}

/* Reformat to .ma for GCTA input 
  exclude SNPs with low relative NEFF
*/
process MA {
  tag "${sumstats}"
  label processMALabel

  cpus = 1
  memory = 16.GB
  time = '30m'

  input: 
  tuple val(pop), val(meta), val(pheno), path(sumstats)
  each neff_pct

  output:
  tuple val(pop), val(meta), val(pheno), path("${sumstats.simpleName}.ma")

  script:
  if (params.build == 'hg19') {
    shell:
      """
    echo -e "SNP\tA1\tA2\tfreq\tBETA\tSE\tP\tN\tCHR\tBP\tNE" > ${sumstats.simpleName}.ma
    bcftools query \
      -f "%ID\\t%ALT\\t%REF\\t[%AFCON]\\t[%ES]\\t[%SE]\\t[%LP]\\t[%SS]\\t%CHROM\\t%POS\\t[%NE]" \
      ${sumstats} | awk -v OFS='\\t' -v neff_threshold=!{params.neff_pct} '\$11 >= neff_threshold * \$11 {print \$1, \$2, \$3, \$4, \$5, \$6, 10^-(\$7), \$8, \$9, \$10, \$11}' >> ${sumstats.simpleName}.ma
      """
  } else if (params.build == 'hg38') {
    """
    #!Rscript
    if( !require('readr') ){
    install.packages('readr', repos = "https://cloud.r-project.org/")
    }
    library(readr)

    if( !require('dplyr') ){
    install.packages('dplyr', repos = "https://cloud.r-project.org/")
    }
    library(dplyr)

    if( !require('stringr') ){
    install.packages('stringr', repos = "https://cloud.r-project.org/")
    }
    library(stringr)

    sumstats <- read_tsv("${sumstats}")

    ma <- sumstats |>
      filter(NEFF >= ${neff_pct} * NEFF) |>
      transmute(SNP, A1, A2, freq = AFCON, BETA = log(OR), SE, P, N = NEFF,
          str_remove(CHR, "chr"), BP)
    
    write_tsv(ma, "${sumstats.simpleName}.ma")
    """
  }
}

/* Run SuSiEx 
  for fine mapping
*/
process FINEMAPPING {

  label susiex

  SuSiEx \
    --sst_file=EUR.sumstats.txt,AFR.sumstats.txt,SAS.sumstats.txt \
    --n_gwas=50000,50000 \
    --ref_file=EUR,AFR,SAS \
    --ld_file=EUR,AFR,SAS \
    --out_dir=./ \
    --out_name=SuSiEx.EUR.AFR.SAS.output.cs95 \
    --level=0.95 \
    --pval_thresh=1e-5 \
    --maf=0.005
    --snp_col=1,1,1 \
    --chr_col=9,9,9 \
    --bp_col=10,10,10 \
    --a1_col=2,2,2 \
    --a2_col=3,3,3 \
    --eff_col=5,5,5 \
    --se_col=6,6,6 \
    --pval_col=7,7,7 \
    --plink=./gpfs/igmmfs01/eddie/GenScotDepression/amelia/packages/plink2 \
    --mult-step=True \
    --keep-ambig=True \
    --threads=16

}





