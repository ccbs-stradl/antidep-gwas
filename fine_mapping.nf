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
  params.chr = [3, 6, 7, 21]

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
  ----- CLUMP process --------------------------
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

/*
  ----- CLUMP_POST process ---------------------
  load clumping results into R and determine region boundaries
*/

  // Combine all the ancestries clumping results into one channel,
  // as we need all ancestries and all chr in the R script together to 
  // find overlapping regions between ancestries.


  CLUMP_CLUSTER_CH = CLUMP_CH
    .collect(flat : false)

  CLUMP_CLUSTER_CH.view()

  CLUMP_POST_CH = CLUMP_POST(CLUMP_CLUSTER_CH)

  CLUMP_POST_CH.view()

}

// ---------------------------------------------
// ************* PROCESSES *********************
// ---------------------------------------------

process MAKE_BFILE {
  tag "chr${chr}:${cluster}"

  cpus = 4
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
  
  cpus = 2
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

process CLUMP_POST {
  label 'analysis'

  cpus = 1
  memory = 32.GB
  time = '30m'

  input:
    val(nested_clumps)

  output:
    tuple path("*.finemapRegions"), path("chr.txt")

  script:
  """
  #!Rscript
  library(stringr)
  library(data.table)
  library(dplyr)

  if (!require("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
  }

  if (!require("plyranges", quietly = TRUE)) {
    BiocManager::install("plyranges")
  }

  library(plyranges) # reduce_ranges

  nested_clumps <- "${nested_clumps}"

  # Remove the inner and outer "[[" and "]]"
  nested_clumps_trim <- gsub("[[", "", nested_clumps, fixed = TRUE) %>%
                          gsub("]]", "", ., fixed = TRUE)

  # Create a list of the collected channels, so each item in the list is a channel before they were collected here:  CLUMP_CLUSTER_CH = CLUMP_CH.collect(flat : false)
  # each item in the list corresponds to a different ancestry
  nested_clumps_list <- as.list(unlist(str_split(nested_clumps_trim, stringr::fixed("], [") )))

  # clean up white space and spit each item in the list by "," into a vector
  nested_clumps_clean <- lapply(nested_clumps_list, function(l){
                        str_split(l, ",") %>%
                        unlist() %>%
                        str_trim() %>%
                        gsub("'", "", ., fixed = TRUE)
                      })

  # start lapply here over nested_clumps_clean, get a list of GRanges for each ancestry
  # return NULL if there is no GRanges object
  granges_list <- lapply(nested_clumps_clean, function(nested_clump_clean){

    # Read in clumps results
    clumped_data <- fread(nested_clump_clean[5]) # index 5 is the .clumps file (see order in output from CLUMP process; index counting starts from 1 not 0 in R)

    if(nrow(clumped_data) == 0){
      # skip this and go to next element in lapply if there is no clumped data
      return() # returns NULL
    } else {

      loci <- clumped_data %>%
        dplyr::select(CHR= `#CHROM`, POS, SNP=ID)

      hg19 <- genome_info('hg19')
      # pull out lengths manually since seqnames uses "chrN" instead of "N"
      hg19_chr_lengths <- as_tibble(hg19) |> slice(1:23) |> pull(width)

      # add genome info for autosomes and X
      grng <- loci %>%
                arrange(CHR) %>%
                as_granges(seqnames = CHR,
                              start = POS,
                              end = POS) 

      seqlevels(grng) <- as.character(1:23)
      seqlengths(grng) <- hg19_chr_lengths
      grng <- set_genome_info(grng, 'hg19', is_circular = rep(FALSE, 23))

      # stretch each region, by 100 kb upstream and downstream, then trim back to position boundaries
      grng_stretched <- stretch(anchor_center(grng), 200000) %>% trim()

      return(grng_stretched)
    }

  })

  # For the regions that are to be fine mapped make a note of which chr they come from
  chr <- sapply(granges_list, function(l){ 
            as.data.frame(l)[["seqnames"]] %>% 
              unique() %>% 
              as.numeric() }) %>%
         unique() %>%
         list()

  # Reduce list of granges into one granges object, then reduce any overlapping regions (from different ancestries)
  # First deal with all NULLs (ie. no clumps data, and no regions to fine map for all ancestries)
  if (all(sapply(granges_list, is.null))) {
    # If the list is all NULL, create a file with NULL content
    writeLines("NULL", paste0("chr.finemapRegions"))
  } else {

  # If the first item in granges_list is a GRanges object then the following code runs:
  # do.call(c, granges_list) %>%
  #                 reduce() 
  # However, if the first item in the list is NULL then that code will result in an error
  # granges_list_null_first <- list(granges_list[[2]], granges_list[[1]], granges_list[[3]])
  #   do.call(c, granges_list_null_first) %>%
  #                   reduce()
  # Error in (function (classes, fdef, mtable)  : 
  # unable to find an inherited method for function ‘reduce’ for signature ‘"list"’

  # Therefore remove NULL items from list to avoid this error
  granges_list <- granges_list[granges_list != "NULL"]

  # Now there are no NULL items, and at least one GRanges object (because of previous "if (all(sapply(granges_list, is.null)))") 
  grng_streched <- do.call(c, granges_list) %>%
                    reduce()

  # Reduce ranges to collapse overlapping or nearby regions
  grng_reduced <- grng_streched %>%
    reduce_ranges(min.gapwidth = 5000)  %>% # What should this be set to?
    as_tibble() %>%
    mutate(WIDTH = end - start + 1) %>%
    dplyr::select(CHR = seqnames, BP_START = start, BP_END = end, WIDTH)


  # Save a table of min and max BP positions for SuSiEx, per chr
  write.table(grng_reduced, paste0("chr.finemapRegions"), row.names = F, quote = F, sep = "\t")
  }

  # Save chr as tab separated file
  fwrite(chr, 'chr.txt' , sep = '\t', col.names = FALSE)
  """
} 


