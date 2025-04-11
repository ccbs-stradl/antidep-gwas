// make forest plots from various sumstats and snplist inputs
// takes clumped and susiex summaries and creates forest plots across
// all GWAS and meta-analysis inputs


params.gwas = "results/vcf/gwas/GRCh38/*.{csv,vcf.gz,vcf.gz.tbi}"
params.meta = "results/vcf/meta/GRCh38/antidep-2501-*.{csv,json,vcf.gz,vcf.gz.tbi}"
params.metaset = "inputs/metasets/antidep-2501.csv"
params.snplist = "results/forest/antidep-2501.snplist"

workflow {
  // sumstats, filenames as COHORT-PHENO-CLUSTER-VERSION
  // parse filename to dictionary
  GWAS_CH = Channel
    .fromFilePairs(params.gwas, size: 3)
    .map { it -> [it[0].split("-")].plus(it) }
    .map { it -> [[cohort: it[0][0], pheno: it[0][1], cluster: it[0][2], version: it[0][3]], it[1], it[2]] }

  // meta-analysis filename as META-SET-METHOD-PHENO-CLUSTER
  // parse filename to dictionary
  META_CH = Channel.fromFilePairs(params.meta, size: 4)
    .map { it -> [it[0].split("-")].plus(it) }
    .map { it -> [[cohort: it[0][2], pheno: it[0][3], cluster: it[0][4], version: it[0][0] + '-' + it[0][1]], it[1], it[2]] }

  SNPLIST_CH = Channel.fromPath(params.snplist)

  // keep input GWAS that were used in one of the meta-analyses
   METASET_CH = Channel
      .fromPath(params.metaset)
   METASET_NAME_CH = METASET_CH
    .map { it -> it.simpleName }

   METASET_INFO_CH = METASET_CH
      .map { it -> [it.baseName, it.splitCsv(header: true)] }
      .transpose()
      .map { it -> it[1].plus([metaset: it[0]]) }

  // select datasets that are used by a metaset
  // key to cohort and version
  VCF_CV_CH = GWAS_CH
    .map { it -> [it[0].cohort, it[0].version].plus(it) }
  METASET_CV_CH = METASET_INFO_CH
    .map { it -> [it.cohort, it.version].plus(it) }

  // merge datasets with metasets, then gather distinct datasets
  GWAS_IN_METASET_CH = METASET_CV_CH
    .combine(VCF_CV_CH, by: [0, 1])
    .map { it -> [it[3], it[4], it[5], it[2]]}
    .groupTuple(by: [0, 1, 2])
    .map { it-> it[0, 1, 2] }
  
  // stack gwas and meta together
  GWAS_META_CH = GWAS_IN_METASET_CH
    .concat(META_CH)
    .combine(SNPLIST_CH)
  
  // extract snplist from sumstats
  EXTRACT_CH = EXTRACT(GWAS_META_CH)
    .map { it -> it[2] }
    .collect()

  // make forest plots
  FOREST_CH = FOREST(EXTRACT_CH, SNPLIST_CH, METASET_NAME_CH)
  
}

process EXTRACT {
  tag "${dataset}"
  label 'tools'


  cpus = 1
  memory = 1.GB
  time = '30m'

  input:
  tuple val(details), val(dataset), path(vcf), path(snplist)

  output:
  tuple val(details), val(dataset), path("${dataset}.tsv")

  script:
  """
  cat ${snplist} | awk -v OFS="\t" '{print \$1, \$2}' > regions.txt
  bcftools view \
    --regions-file regions.txt \
    --output-type u \
    ${dataset}.vcf.gz |\
  bcftools query \
    --format "%CHROM %POS %ID %ALT %REF [%ES] [%SE] [%LP]\\n" \
    > query.txt

  cat query.txt | awk -v OFS="\t" '{if(NR == 1) {print "CHROM", "POS", "ID", "ALT", "REF", "ES", "SE", "LP", "dataset", "cohort", "pheno", "cluster", "version"} else {print \$1, \$2, \$3, \$4, \$5, \$6, \$7, \$8, "${dataset}", "${details.cohort}", "${details.pheno}", "${details.cluster}", "${details.version}"}}' > ${dataset}.tsv
  """

}

process FOREST {
  label 'analysis'

  publishDir "results/forest/${metaset}", mode: "copy"

  cpus = 1
  memory = 4.GB
  time = '30m'

  input:
  path(sumstats)
  path(snplist)
  val(metaset)

  output:
  path("*.png")

  script:
  """
  #!Rscript

  library(dplyr)
  library(readr)
  library(stringr)
  library(ggplot2)

  # parse filenames from inputs string, read in tables, and 
  # bind together
  paths <- str_split("${sumstats}", pattern = " ")[[1]]
  sumstats_list <- lapply(paths, read_tsv)
  sumstats <- bind_rows(sumstats_list)

  # list of SNPs to be plotted
  snplist <- read_tsv("${snplist}",
                      col_names = c("CHROM", "POS", "ID",
                                    "REF", "ALT"))
  
  cohorts <- sumstats |>
    distinct(cohort) |>
    filter(cohort != 'fixed') |>
    arrange(desc(cohort)) |>
    pull(cohort)

  # match 
  # remove unused phenotype, sort by cohort then meta
  sumstats_meta <- sumstats |>
    inner_join(snplist) |>
    filter(pheno != "N06AX") |>
    mutate(Cohort = factor(cohort,
                           levels = c("fixed", cohorts),
                           labels = c("Fixed", cohorts)))

  # group rows by variant
  sumstats_variants <- sumstats_meta |>
    group_by(CHROM, POS, ID, ALT, REF)

  # plotting function on a subset of variants
  forest_plot <- function(sumstats, variant) {

    title <- variant |> mutate(title = str_glue("{CHROM}:{POS}:{REF}:{ALT} - {ID}")) |> pull(title)
    filename <- variant |> mutate(filename = str_glue("{CHROM}_{POS}_{REF}_{ALT}-{ID}.png")) |> pull(filename)

    g <- ggplot(sumstats, aes(x = Cohort, colour = cohort, y = exp(ES), ymin = exp(ES - SE), ymax = exp(ES + SE))) +
    geom_pointrange() +
    #geom_linerange(aes(y = ES, ymin = ES + qnorm(2.5e-8) * SE, ymax = ES + qnorm(1-2.5e-8) * SE), linetype = "dashed") +
    facet_grid(pheno + cluster ~ ., scales = "free", space = "free") +
    scale_y_continuous("Odds Ratio") +
    # scale_colour_manual(values = c("AllOfUs" = "#4E9CD7", "BBJ" = "#CE3F7F", "FinnGen" = "#3900CF", "GenScot" = "#0B5EB0", "UKB" = "#0C5276", "fixed" = "black")) +
    scale_colour_manual(values = c("AllOfUs" = hsv(206/360, .64, 1), "BBJ" = hsv(333/360, .69, 1), "FinnGen" = hsv(257/360, 1, 1), "GenScot" = hsv(210/360, .94, 1), "UKB" = hsv(192/360, .71, 1), "fixed" = "black")) +
    coord_flip() +
    ggtitle(title) +
    theme_minimal() +
    theme(strip.text.y = element_text(angle = 0))

    ggsave(filename, plot = g)

    return(list(variant = variant, sumstats = sumstats, forest = g))

  }

  sumstats_forest <- sumstats_variants |> group_map(.f = forest_plot)
  """
}