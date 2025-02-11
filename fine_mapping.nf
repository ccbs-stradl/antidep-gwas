nextflow.enable.dsl=2

/*
// ---------------------------------------------
// ************ SCRIPT OVERVIEW ****************
// ---------------------------------------------
// Fine-mapping of GWAS meta using SuSiEx (hg19 only)

PER ANCESTRY, CHR
    * make bim, bed, fam files for each ancestry:
      - pfile 
      - IDs where for each ancestry (output from make-pgen.sh)

PER ANCESTRY
    * process meta sumstats into correct format:
      - meta-analysed sumstats: 'vcf/meta/GRCh37/antidep-2501-fixed-N06A-*.{csv,json,vcf.gz,vcf.gz.tbi}' (output from format_meta.nf)
        - these are subsetted to ancestries we have LD references for
      - json file contains meta data including the effective sample size  

PER ANCESTRY, CHR
    * clumping to identify regions for fine mapping:
      - processed sumstats (from previous process in current script)
      - bfile, per chr and ancestry, (output from make-pgen_hg19.sh)

COMBINE ALL ANCESTRIES AND ALL CHR
    * load clumping results into R and determine region boundaries:
      - output from clumping process, per chr and ancestry (from previous process in current script)

PER CHR
    * run SuSiEx on each region:
      - output from R determining region boundaries (from previous process in current script)
      - processed sumstats (from process above)
*/

// -----  Import groovy functions ---------------
import groovy.json.JsonSlurper

// ----- Get params -----------------------------
// MAKE_BFILE INPUTS:
  params.pfile = '/exports/igmm/eddie/GenScotDepression/data/ukb/genetics/impv3_pgen/ukb_imp_v3.qc.{pgen,psam,pvar}' // prefix of the pfiles to pass to --pfile in plink 
  params.ancestry_ids = 'reference/ukb-ld-ref_ancestry.id' // Output from script make-pgen.sh
  params.chr = [3, 21] // testing on chr 3 = produces results, and chr 21 = should create a null file

// MA INPUTS:
  // load in gwas meta sumstats that have been lifted over to hg19
  params.meta = 'vcf/meta/GRCh37/antidep-2501-fixed-N06A-*.{csv,json,vcf.gz,vcf.gz.tbi}'

  // effective sample size QC parameter
  params.neff_pct = 0.8 // gwas is filtered by effectve sample size, threshold is 0.8

// ---------------------------------------------
// *************** WORKFLOW ********************
// ---------------------------------------------
workflow {

/*
  ----- MAKE_BFILE process ---------------------
  make bim, bed, fam files for each ancestry/cluster
  this only needs to happen for chromosomes that will be fine mapped, 
  is there a way to get this from the sumstats?
*/
  // Make reference files in bfile format
    PFILE_CH = Channel.fromFilePairs(params.pfile, size: 3)// returns tuple: [prefix [pfile paths]]          
    IDS_CH = Channel.fromPath(params.ancestry_ids) 
    CHR_CH = Channel.fromList(params.chr)   

    // Create CLUSTER_CH to extract unique ancestries from the 3rd column of IDS_CH
    CLUSTER_CH = IDS_CH
      .splitCsv(header: false, sep: ' ') 
      .map { it[2] }                     
      .unique()

    REF_CH = PFILE_CH
     .combine(IDS_CH)
     .combine(CLUSTER_CH)
     .combine(CHR_CH)

    BFILE_CH = MAKE_BFILE(REF_CH)

/*
  ----- MA process ---------------------------- 
  process meta sumstats into correct format, 
  and collect effective sample size from json
*/

  def jsonSlurper = new JsonSlurper()
  META_CH = Channel.fromFilePairs(params.meta, size: 4) { it -> it.simpleName }
    .map { it -> 
        def metadata = jsonSlurper.parseText(it[1][1].text)  // Extract JSON metadata
        return [metadata.cluster, it[0], it[1], metadata]  // Explicitly store ancestry as a key
    } // last element is a json containing metadata about the meta analysis gwas, eg. pheno, ancestry, sample sizes, effective sample size
    .join(CLUSTER_CH) // Subset META_CH so it only contains ancestry/cluster in CLUSTER_CH 

  // QC parameters
  NEFF_QC_CH = Channel.of(params.neff_pct)

  // format sumstats to .ma
  MA_CH = MA(META_CH, NEFF_QC_CH)

/*
  ----- CLUMP process -------------------------
  clumping to identify regions for fine mapping
*/

  // get a channel of CHRs (which have significant SNPs) and their corresponding ancestry
  // to filter the clumping step to only happen for CHR and CLUSTERS with genome wide sig SNPs
    SIG_CHR_CH = MA_CH
      .map {it -> it[4]}
      .splitCsv(header: false, sep: '\t') 
      .map { it -> [ it[1], it[0] ] } // [CLUSTER, CHR] - by default numeric values are converted to strings

  // join sumstats with bfiles ready for clumping
    MA_BFILE_CH = BFILE_CH
      .combine(MA_CH, by: 0) // join all chr based on ancestries
      .map(it -> [ it[0], it[1].toString(), it[2], it[3], it[4], it[5], it[6] ]) // convert CHR to string and drop *.sig_chr
      .join(SIG_CHR_CH, by: [0,1])
      
  CLUMP_CH = CLUMP(MA_BFILE_CH)

}

// ---------------------------------------------
// ************* PROCESSES *********************
// ---------------------------------------------

process MAKE_BFILE {
  tag "chr${chr}:${cluster}"

  cpus = 1
  memory = 8.GB
  time = '1h'

  input:
    tuple val(dataset), path(pfile), path(ancestry_ids), val(cluster), val(chr)

  output:
    tuple val(cluster), val(chr), val("${dataset}.geno02.mind02_${cluster}_${chr}"), path("${dataset}.geno02.mind02_${cluster}_${chr}.*")
 
  script:
  """
  plink2 \
    --pfile  ${dataset} \
    --keep-col-match ${ancestry_ids} ${cluster} \
    --chr ${chr} \
    --geno 0.02 \
    --maf 0.005 \
    --make-pgen 'vzs' \
    --out ${dataset}.geno02_${cluster}_${chr} \
    --threads ${task.cpus} \
    --memory ${task.memory.bytes.intdiv(1000000)}

    plink2 \
    --pfile ${dataset}.geno02_${cluster}_${chr} 'vzs' \
    --mind 0.02 \
    --make-bed \
    --out ${dataset}.geno02.mind02_${cluster}_${chr} \
    --threads ${task.cpus} \
    --memory ${task.memory.bytes.intdiv(1000000)}
  """

}

process MA {
  tag "${cluster}"
  label 'tools'

  cpus = 1
  memory = 16.GB
  time = '30m'

  input: 
  tuple val(cluster), val(dataset), path(vcf), val(info)
  each neff_pct

  output:
  tuple val(cluster), val(dataset), val(info.neff), path("*.ma"), path("*.sig_chr")

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

      # Rewrite sumstats with suffix 'noZero'
      fwrite(sumstats, ma_path, sep = '\t')

      # Write a csv with a CHR column, and the number of that CHR if there is a significant SNP for it
      # if there are no significant SNPs then write an empty file

      sig_chr <- sumstats %>%
                  filter(P <= 5e-8) %>%
                  distinct(CHR) %>%
                  mutate(CLUSTER = '${cluster}')

      if ( nrow(sig_chr) > 0 ) {
        fwrite(sig_chr, paste0('${dataset}', '.sig_chr' ) , sep = '\t', col.names = FALSE)
      } else {
        # Create an empty file 
        file.create(paste0('${dataset}', '.sig_chr'))
      }

      "
    """

}

process CLUMP {
  tag "chr${chr}:${cluster}"
  
  cpus = 1
  memory = 8.GB
  time = '10m'

  input:
    tuple val(cluster), val(chr), val(bfile_prefix), path(bfiles), val(dataset), val(neff), path(ma)

  output:
    tuple val(cluster), val(chr), val(dataset), path(ma), path("${dataset}.${chr}.clumps"), path("${dataset}.${chr}.log")
 
  script:
  """
  plink2 \
    --clump ${dataset}.ma \
    --clump-id-field SNP \
    --clump-p-field P \
    --clump-p1 5e-8 \
    --clump-p2 0.05 \
    --clump-r2 0.1 \
    --clump-kb 1000 \
    --bfile ${bfile_prefix} \
    --out ${dataset}.${chr} \
    --threads ${task.cpus} \
    --memory ${task.memory.bytes.intdiv(1000000)}
  
  # Check if .clumps file exists; if not, create the file with content "NULL"
  # This is needed otherwise nextflow exits with an error because plink has not made a .clumps file (because there are no clumps in the gwas for that chromosome)
    if [ ! -f "${ma.baseName}.${chr}.clumps" ]; then
        echo "NULL" > "${ma.baseName}.${chr}.clumps"
    fi

  """
}
