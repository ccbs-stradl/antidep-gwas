// Format GWASs (hg19) from VCF to MA to use as input for drug targetor
/*
To run script:
nextflow run drug_targetor/format_gwas_drugtar.nf -resume \
-work-dir /exports/eddie/scratch/${USER}/ad/work \
-c eddie.config -process.executor local

*/

nextflow.enable.dsl=2

// -----  Import groovy functions ---------------
import groovy.json.JsonSlurper

// -----  Params
// MA INPUTS:
  // load in gwas meta sumstats that have been lifted over to hg19
params.meta = 'vcf/meta/GRCh37/antidep-2501-fixed-N06A-*.{csv,json,vcf.gz,vcf.gz.tbi}'


// -----  Workflow
workflow {

  def jsonSlurper = new JsonSlurper()
  META_CH = Channel.fromFilePairs(params.meta, size: 4) { it -> it.simpleName }
    .map { it -> 
        def metadata = jsonSlurper.parseText(it[1][1].text)  // Extract JSON metadata
        return [metadata.cluster, it[0], it[1], metadata]  // Explicitly store ancestry as a key
    } // last element is a json containing metadata about the meta analysis gwas, eg. pheno, ancestry, sample sizes, effective sample size

  // QC parameters
  NEFF_QC_CH = Channel.of(params.neff_pct)

  // format sumstats to .ma
  MA_CH = MA(META_CH, NEFF_QC_CH)

}

// -----  Processes
process MA {
  tag "${cluster}"
  label 'tools'

  cpus = 1
  memory = 16.GB
  time = '30m'

  publishDir "drug_targetor/input_files", mode: "copy"

  input: 
  tuple val(cluster), val(dataset), path(vcf), val(info)
  each neff_pct

  output:
  tuple val(cluster), val(dataset), val(info.neff), path("*.ma")

  script:
    """
    echo -e "SNP\tA1\tA2\tfreq\tBETA\tSE\tP\tN\tCHR\tBP\tNE" > ${dataset}.ma
    bcftools query \
      -f "%ID\\t%ALT\\t%REF\\t[%AFCON]\\t[%ES]\\t[%SE]\\t[%LP]\\t[%SS]\\t%CHROM\\t%POS\\t[%NE]" \
      ${dataset}.vcf.gz | awk -v OFS='\\t' -v neff_threshold=!{params.neff_pct} '\$11 >= neff_threshold * \$11 {print \$1, \$2, \$3, \$4, \$5, \$6, 10^-(\$7), \$8, \$9, \$10, \$11}' >> ${dataset}.ma

    Rscript -e "
      library(data.table)
      library(dplyr)

      ma_path <- '${dataset}.ma'

      # Read in sumsstats and remove rows where there are duplicate SNPs
      sumstats <- fread( ma_path ) %>%
                    distinct(SNP, .keep_all = TRUE) # keeps first distinct SNP row; .keep_all =T returns all columns

      # Overwrite sumstats
      fwrite(sumstats, ma_path, sep = '\t')

      "
    """

}



