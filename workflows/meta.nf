nextflow.enable.dsl=2

/*
    Meta analyse GWASVCF files
*/

params.vcf = "results/vcf/gwas/GRCh38/*.{vcf.gz,vcf.gz.tbi,csv}"
params.metasets = "metasets/*.csv"
params.ref = "reference/all_hg38.{pgen,psam,pvar.zst}"

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

    REF_CH = Channel
        .fromFilePairs(params.ref, size: 3)
        
    /*
      Pre-processing
    */
      
      // select datasets that are used by a metaset
      // key to cohort and version
      VCF_CV_CH = VCF_CH
        .map { it -> [it[0].cohort, it[0].version].plus(it) }
      METASET_CV_CH = METASET_CH
        .map { it -> [it.cohort, it.version].plus(it) }
        
      // merge datasets with metasets, then gather distinct datasets
      VCF_IN_METASET_CH = METASET_CV_CH
        .combine(VCF_CV_CH, by: [0, 1])
        .map { it -> [it[3], it[4], it[5], it[2]]}
        .groupTuple(by: [0, 1, 2])

      // remove metaset information for input processing
      VCF_IN_CH = VCF_IN_METASET_CH
        .map { it -> it[0..2]}

      // track metasets
      IN_META_CH = VCF_IN_METASET_CH
        .map { it ->  it[0,3] }
        
    /*
        MR-MEGA
        https://genomics.ut.ee/en/tools
    */
    
      // process for input into MR-MEGA
      // group by phenotype and metaset
      // conduct meta-analysis if there are at least 4 cohorts
      MEGA_IN_CH = MEGA_IN(VCF_IN_CH)
        .combine(IN_META_CH, by: 0)
        .transpose(by: 4)
        .map { it -> [[metaset: it[4].metaset, pheno: it[0].pheno]].plus(it) }
        .groupTuple(by: 0)
        .filter { it -> it[1].size() >= 4 }
        .map { it -> it[0..4] }
      
	
      MEGA_CH = MEGA(MEGA_IN_CH)
     
      MEGA_POST_CH = MEGA_POST(MEGA_CH)

      MEGA_ASSOC_CH = MEGA_POST_CH
        .map { it -> it[0,1,4] }
        .transpose()

    /*
        Fixed effects (polars)
    */
    
    // process for input into fixed-effects analysis
    // group by phenotype, cluster, and metaset
    // conduct meta-analysis if there are at least 2 cohorts
    FIXED_IN_CH = FIXED_IN(VCF_IN_CH)
      .combine(IN_META_CH, by: 0)
      .transpose(by: 4)
      .map { it -> [[metaset: it[4].metaset, pheno: it[0].pheno, cluster: it[0].cluster]].plus(it) }
      .groupTuple(by: 0)
      .filter { it -> it[1].size() >= 2 }
      .map { it -> it[0..4] }
      
    FIXED_CH = FIXED(FIXED_IN_CH)
    
    FIXED_POST_CH = FIXED_POST(FIXED_CH)

    FIXED_ASSOC_CH = FIXED_POST_CH
        .map { it -> it[0,1,4] }
        .transpose()

    // /* 
    //     Post analysis
    // */
    
    ASSOC_CH = MEGA_ASSOC_CH
         .concat(FIXED_ASSOC_CH)

    MANHATTAN(ASSOC_CH)

    // count number of GWsig variants
    // further processing of results with 1 or more loci
    ASSOC_GW_CH = GW(ASSOC_CH)
      .filter { it -> it[3].toInteger() > 0 }
      .map { it -> it[0,1,2] }

    // match association and reference SNPs with CPID
    REF_CPID_CH = REF_CPID(REF_CH)

    ASSOC_REF_CH = ASSOC_GW_CH
         .combine(REF_CPID_CH)

    // clumping
    CLUMP_CH = CLUMP(ASSOC_REF_CH)
      .groupTuple(by: [0, 1])

    // reports
    POST_CH = MEGA_POST_CH
      .concat(FIXED_POST_CH)
      .combine(CLUMP_CH, by: [0, 1])
    
    POST(POST_CH)
	
}

/* Query VCF to create MR-MEGA GWA input file 
    1) MARKERNAME – snp name
    2) EA – effect allele
    3) NEA – non effect allele
    4) OR - odds ratio
    5) OR_95L - lower confidence interval of OR
    6) OR_95U - upper confidence interval of OR
    7) EAF – effect allele frequency
    8) N - sample size
    9) CHROMOSOME  - chromosome of marker
    10) POSITION - position of marker
*/
process MEGA_IN {
    tag "${dataset}"
    label 'tools'

    cpus = 1
    memory = 1.GB
    time = '30m'

    input:
    tuple val(details), val(dataset), path(vcf)

    output:
    tuple val(details), val(dataset), path("${dataset}.txt.gz"), path("${dataset}.csv", includeInputs: true)

    script:
    """
    # bcf query to get required columns
    # remove variants with MAF < 0.5%, INFO < 0.4
    # rename missing variants names to CPID
    bcftools query \
    -i 'FORMAT/AF >= 0.005 & FORMAT/AF <= 0.995 & (FORMAT/SI >= 0.4 | FORMAT/SI = ".")' \
	-f "%ID\\t%ALT\\t%REF\\t[%ES]\\t[%SE]\\t[%AF]\\t[%SS]\\t[%NC]\\t[%NE]\\t%CHROM\\t%POS\\n" \
	${dataset}.vcf.gz | awk '{if(\$1 == ".") \$1 = \$9":"\$10; print \$0}' > ${dataset}.query.txt

    # header
    echo "" | awk 'OFS = "\\t" {print "MARKERNAME", "EA", "NEA", "OR", "OR_95L", "OR_95U", "EAF", "N", "CHROMOSOME", "POSITION"}' > ${dataset}.txt

    # reformat
    # remove "chr" from chromosome names
    # turn beta effect size back into OR
    # calculate N as effective sample size
    cat ${dataset}.query.txt | awk \
    'OFS = "\\t" {gsub("chr", "", \$10); print \$1, \$2, \$3, exp(\$4), exp(\$4-1.96*\$5), exp(\$4+1.96*\$5), \$6, \$9, \$10, \$11}' >> ${dataset}.txt

    # compress
    gzip ${dataset}.txt
    """
/*
    column numbers of bcf query
     1	%ID
     2	%ALT
     3	%REF
     4	[%ES]
     5	[%SE]
     6	[%AF]
     7	[%SS]
     8	[%NC]
     9  [%NE]
    10	%CHROM
    11	%POS
*/
}

process MEGA {
    tag "${metaset.metaset}-${metaset.pheno}"

    publishDir "results/meta/${metaset.metaset}", pattern: "*.log", mode: "copy"

    cpus = 1
    memory = 48.GB
    time = '3h'

    input:
    tuple val(metaset), val(details), val(datasets), path(gwas), path(csv)

    output:
    tuple val(metaset), val(details), val("${metaset.metaset}-mrmega-${metaset.pheno}-DIV"), path("*.result"), path("*.log"), path(csv, includeInputs: true)

    script:
    """
    ls *.txt.gz > mrmega.in

    MR-MEGA --filelist mrmega.in \
    --pc ${[gwas.size() - 3, 3].min()} \
    --out ${metaset.metaset}-mrmega-${metaset.pheno}-DIV
    """
}

process MEGA_POST {
    tag "${analysis}"
    label 'analysis'

    publishDir "results/meta/${metaset.metaset}", pattern: "*.{gz,csv}", mode: 'copy'

    cpus = 1
    memory = 31.GB
    time = '30m'

    input:
    tuple val(metaset), val(details), val(analysis), path(result), path(log), path(datasets)

    output:
    tuple val(metaset), val(analysis), path("${analysis}.gz"), path("${analysis}.csv"), path("*.assoc")

    script:
    """
    #!Rscript
    library(dplyr)
    library(readr)
    library(stringr)
	
    # meta-analysis
	  mega <- read_tsv("${result}")
    # datasets info
    dataset_files <- str_split("${datasets}", pattern = " ")
    datasets <- bind_rows(lapply(dataset_files, read_csv)) |>
      arrange(cohort, cluster)

    # recalculate P-values to allow for smaller values
    mega_p <- mega |>
        filter(!is.na(beta_0)) |>
        mutate(`P-value_association` = pchisq(chisq_association, ndf_association, lower.tail=FALSE),
               `P-value_ancestry_het` = pchisq(chisq_ancestry_het, ndf_ancestry_het, lower.tail=FALSE),
               `P-value_residual_het` = pchisq(chisq_residual_het, ndf_residual_het, lower.tail=FALSE)
        ) |>
    arrange(Chromosome, Position) |>
    mutate(CHR = if_else(Chromosome == 23, true = "X", false = as.character(Chromosome)))

    mega_assoc <- mega_p |>
        transmute(CHR, SNP=MarkerName, BP=Position, A1=EA, A2=NEA,
                  P=`P-value_association`, BETA=beta_0, SE=se_0, NMISS=Nsample,
                  CPID = str_glue("{CHR}:{BP}:{A2}:{A1}"))

    mega_ancestry <- mega_p |>
        transmute(CHR, SNP=MarkerName, BP=Position, A1=EA, A2=NEA,
                   P=`P-value_ancestry_het`, BETA=beta_0, SE=se_0, NMISS=Nsample,
                   CPID = str_glue("{CHR}:{BP}:{A2}:{A1}"))

    mega_residual <- mega_p |>
        transmute(CHR, SNP=MarkerName, BP=Position, A1=EA, A2=NEA,
                   P=`P-value_residual_het`, BETA=beta_0, SE=se_0, NMISS=Nsample,
                   CPID = str_glue("{CHR}:{BP}:{A2}:{A1}"))


    write_tsv(mega_p, "${analysis}.gz")
    write_csv(datasets, "${analysis}.csv")

    write_tsv(mega_assoc, "${analysis}.assoc")
    write_tsv(mega_ancestry, "${analysis}.ancestry.assoc")
    write_tsv(mega_residual, "${analysis}.residual.assoc")
    """

}


/* Query VCF to text input file 
*/
process FIXED_IN {
    tag "${dataset}"
    label 'tools'

    cpus = 1
    memory = 1.GB
    time = '30m'

    input:
    tuple val(details), val(dataset), path(vcf)

    output:
    tuple val(details), val(dataset), path("${dataset}.query.txt"), path("${dataset}.csv", includeInputs: true)

    script:
    """
    # bcf query to get required columns
    # remove variants with MAF < 0.5%, INFO < 0.4
    # rename missing variants names to CPID
    echo "#CHROM\tPOS\tID\tALT\tREF\tES\tSE\tLP\tSI\tAFCAS\tAFCON\tSS\tNC\tNE" > ${dataset}.query.txt
    bcftools norm -m- ${dataset}.vcf.gz |\
    bcftools query \
    -s ${dataset} \
    -i 'FORMAT/AF >= 0.005 & FORMAT/AF <= 0.995 & (FORMAT/SI >= 0.4 | FORMAT/SI = ".")' \
	  -f "%CHROM\\t%POS\\t%ID\\t%ALT\\t%REF\\t[%ES]\\t[%SE]\\t[%LP]\\t[%SI]\\t[%AFCAS]\\t[%AFCON]\\t[%SS]\\t[%NC]\\t[%NE]\\n" |\
    awk '{if(\$3 == ".") \$3 = \$1":"\$2; print \$0}' >> ${dataset}.query.txt
    """
/*
    columns of bcf query
     #CHROM
     POS
     ID
     ALT
     REF
     ES
     SE
     LP
     SI
     AFCAS
     AFCON
     SS
     NC
     NE
*/
}

// Fixed effects meta analysis
// Effect sizes weighted by inverse variance
// Sample information that is weighted by N
process FIXED {
    tag "${metaset.metaset}-${metaset.pheno}-${metaset.cluster}"
    label 'analysis'

    cpus = 6
    memory = 48.GB
    time = '30m'

    input:
    tuple val(metaset), val(details), val(datasets), path(gwas), path(csv)

    output:
    tuple val(metaset), val("${metaset.metaset}-fixed-${metaset.pheno}-${metaset.cluster}"), path("meta.tsv"), path(csv, includeInputs: true)

    script:
    """
    #!python3
    import os
    
    os.environ["POLARS_MAX_THREADS"] = "${task.cpus}"
    import polars as pl

    gwas_paths = "${gwas}"

    gwas = pl.concat(
        [pl.scan_csv(path, separator= "\\t", null_values = ".", schema_overrides = {"NE": pl.Float64, "SI": pl.Float64, "AFCAS": pl.Float64, "AFCON": pl.Float64}) for path in gwas_paths.split()],
        how = "vertical"
    )

    # set up weights
    weights = (
      gwas.select(
          pl.col("#CHROM").alias("CHROM"),
          pl.col("POS"),
          pl.col("ID"),
          pl.col("ALT"),
          pl.col("REF"),
          pl.col("ES").alias("BETA"),
          pl.col("SE"),
          pl.col("SI").alias("INFO"),
          pl.col("AFCAS"),
          pl.col("AFCON"),
          pl.col("SS").alias("N"),
          pl.col("NC").alias("NCAS"),
          pl.col("NE").alias("NEFF"),
          # calculate number of controls
          (pl.col("SS") - pl.col("NC")).alias("NCON")
      )
      .select(
          pl.col("*"),
          # weight effect by inverse variance
          (1 / pl.col("SE")**2).alias("W"),
      )
      .select(
          pl.col("*"),
          (pl.col("BETA") * pl.col("W")).alias("wBETA"),
          # weight INFO by N, AF by cases and controls
          (pl.col("INFO") * pl.col("N")).alias("wINFO"),
          (pl.col("AFCAS") * pl.col("NCAS")).alias("wAFCAS"),
          (pl.col("AFCON") * pl.col("NCON")).alias("wAFCON"),
          # replicate missingness over to a new N column
          pl.when(pl.col("INFO").is_null()).then(None).otherwise(pl.col("N")).alias("nINFO"),
          pl.when(pl.col("AFCAS").is_null()).then(None).otherwise(pl.col("NCAS")).alias("nAFCAS"),
          pl.when(pl.col("AFCON").is_null()).then(None).otherwise(pl.col("NCON")).alias("nAFCON")
      )
    )
      
    # meta-analyse
    meta = ( 
      weights.group_by("CHROM", "POS", "ID", "ALT", "REF")
      .agg(
          (pl.sum("wBETA") / pl.sum("W")).alias("BETA"),
          (1/pl.sum("W")).sqrt().alias("SE"),
          (pl.sum("wINFO") / pl.sum("nINFO")).alias("INFO"),
          (pl.sum("wAFCAS") / pl.sum("nAFCAS")).alias("AFCAS"),
          (pl.sum("wAFCON") / pl.sum("nAFCON")).alias("AFCON"),
          pl.sum("NCAS"),
          pl.sum("NCON"),
          pl.sum("NEFF"),
          pl.sum("N").alias("NTOT")
      )
    )
    
    # heterogeneity (sum of squared deviations)
    # remove if there is only one study
    het = (
      weights
      .join(meta, on = ["CHROM", "POS", "ID", "ALT", "REF"], how = "inner", suffix = "_bar")
      .select(
        pl.col("CHROM"),
        pl.col("POS"),
        pl.col("ID"),
        pl.col("ALT"),
        pl.col("REF"),
        ((pl.col("BETA") - pl.col("BETA_bar"))**2).alias("SQD"),
        pl.col("W")
      )
      .select(
        pl.col("*"),
        (pl.col("SQD") * pl.col("W")).alias("wSQD")
      )
      .group_by("CHROM", "POS", "ID", "ALT", "REF")
      .agg(
        pl.sum("wSQD").alias("Q"),
        pl.len().alias("studies")
      )
      .with_columns(
        pl.when(pl.col("studies") == 1)
        .then(None)
        .otherwise(pl.col("Q"))
        .alias("Q")
      )
    )
    
    meta_het = (
      meta.join(het, on = ["CHROM", "POS", "ID", "ALT", "REF"], how = "inner")
      .collect()
    )
    
    meta_het.write_csv("meta.tsv", separator = "\\t", null_value = "NA")
    """
}

// Fixed effects postprocessing
process FIXED_POST {
    tag "${analysis}"
    label 'analysis'

    publishDir "results/meta/${metaset.metaset}", pattern: "*.{gz,csv}", mode: 'copy'

    cpus = 2
    memory = 24.GB
    time = '30m'

    input:
    tuple val(metaset), val(analysis), path(meta), path(datasets)

    output:
    tuple val(metaset), val(analysis), path("${analysis}.gz"), path("${analysis}.csv"), path("*.assoc")

    script:
    """
    #!Rscript
	  library(dplyr)
	  library(readr)
    library(stringr)

    meta_col_types <- cols(CHROM = col_character(), POS = col_integer())
    meta <- read_table("${meta}", col_types = meta_col_types)
    
    # datasets info
    dataset_files <- str_split("${datasets}", pattern = " ")
    datasets <- bind_rows(lapply(dataset_files, read_csv)) |>
      arrange(cohort)

    # calculate p values
    meta_p <- meta |>
        mutate(CHISQ = BETA^2 / SE^2) |>
        mutate(P = pchisq(CHISQ, df = 1, lower.tail = F),
               QP = pchisq(Q, df = studies - 1, lower.tail = F)) |>
        select(CHROM, POS, ID, REF, ALT, studies,
               BETA, SE, CHISQ, P, Q, QP, INFO,
               AFCAS, AFCON, NCAS, NCON, NEFF, NTOT)


    assoc <- meta_p |>
        transmute(CHR = str_remove(CHROM, 'chr'), SNP = ID, BP = POS,
                  A1 = ALT, A2 = REF, P, BETA, ESS=NEFF,
                  CPID = str_glue("{CHR}:{BP}:{A2}:{A1}"))
    
    het <- meta_p |>
      filter(!is.na(QP)) |>
      transmute(CHR = str_remove(CHROM, 'chr'), SNP = ID, BP = POS,
                A1 = ALT, A2 = REF, P = QP, BETA, ESS=NEFF,
                CPID = str_glue("{CHR}:{BP}:{A2}:{A1}"))

    write_tsv(meta_p, "${analysis}.gz")
    write_csv(datasets, "${analysis}.csv")
    write_tsv(assoc, "${analysis}.assoc")
    write_tsv(het, "${analysis}.het.assoc")
    """
}

/* gwsig variants */
process GW {
    tag "${assoc.baseName}"

    cpus = 1
    memory = 1.GB
    time = '10m'

    input:
    tuple val(metaset), val(analysis), path(assoc)

    output:
    tuple val(metaset), val(analysis), path(assoc, includeInputs: true), env(GW)

    script:
    """
    GW=\$(cat ${assoc} | awk '\$6 <= 5e-8' | wc -l)
    """

}

/* manhattan plot */
process MANHATTAN {
	tag "${assoc}"
  label 'analysis'

	publishDir "results/meta/${metaset.metaset}", mode: 'copy'
	
	cpus = 1
	memory = 16.GB
	time = '30m'
	
	input:
	tuple val(metaset), val(analysis), path(assoc)
	
	output:
	path "*.png"
	
	script:
	"""
	#!Rscript
	library(dplyr)
	library(readr)
	library(fastman)
	
	assoc <- read_tsv("${assoc}", col_types = cols(CHR = col_character())) |> 
        mutate(CHR = case_match(CHR, "X" ~ 23, .default = as.numeric(CHR)))

    assoc_p <- assoc |> filter(P > 0)
	
	png("${assoc.baseName}.manhattan.png", width=10, height=6, units="in", res=300)
	fastman(assoc_p, chr = "CHR", bp = "BP", p = "P", maxP=NULL)
	dev.off()
	
	png("${assoc.baseName}.qq.png", width=6, height=6, units="in", res=300)
	fastqq(p1=pull(assoc_p, P), maxP=NULL)
	dev.off()
	"""
}

/* ref variants IDs as CPID */
process REF_CPID {
    tag "${ref}"

    cpus = 1
    memory = 8.GB
    time = '1h'

    input:
    tuple val(ref), path(pgen)

    output:
    tuple path("${ref}.pgen", includeInputs: true), path("${ref}.psam", includeInputs: true), path("cpid.pvar.zst")

    script:
    """
    plink2 \
    --make-just-pvar 'zs' \
    --pfile 'vzs' ${ref} \
    --set-all-var-ids @:#:\\\$r:\\\$a \
    --new-id-max-allele-len 500 error \
    --allow-extra-chr \
    --out cpid \
	  --threads ${task.cpus} \
	  --memory ${task.memory.bytes.intdiv(1000000)}
    """
}

process CLUMP {
    tag "${assoc.baseName}"

    errorStrategy 'ignore'

    cpus = 4
    memory = 16.GB
    time = '1h'

    input:
    tuple val(metaset), val(analysis), path(assoc), path(pgen), path(psam), path(pvar)

    output:
    tuple val(metaset), val(analysis), path("${assoc.baseName}.clumps"), path("${assoc.baseName}.log")

    script:
    if(metaset.containsKey("cluster")) {
      """
      plink2 \
      --clump ${assoc} \
      --clump-id-field CPID \
      --clump-p1 5e-7 \
      --clump-p2 1.0 \
      --clump-r2 0.1 \
      --clump-kb 1000 \
      --keep-cat-names ${metaset.cluster} \
      --keep-cat-pheno SuperPop \
      --keep-founders \
      --pgen ${pgen} \
      --psam ${psam} \
      --pvar ${pvar} \
      --allow-extra-chr \
      --out ${assoc.baseName} \
      --threads ${task.cpus} \
      --memory ${task.memory.bytes.intdiv(1000000)}
      """
    } else {
      """
      plink2 \
      --clump ${assoc} \
      --clump-id-field CPID \
      --clump-p1 5e-7 \
      --clump-p2 1.0 \
      --clump-r2 0.1 \
      --clump-kb 1000 \
      --keep-founders \
      --pgen ${pgen} \
      --psam ${psam} \
      --pvar ${pvar} \
      --allow-extra-chr \
      --out ${assoc.baseName} \
      --threads ${task.cpus} \
      --memory ${task.memory.bytes.intdiv(1000000)}
      """
    }
}

process POST {
    tag "${analysis}"
    label 'analysis'

    publishDir "results/meta/${metaset.metaset}", pattern: "*.tsv", mode: 'copy'

    cpus = 1
    memory = 31.GB
    time = '30m'

    input:
    tuple val(metaset), val(analysis), val(meta), val(datasets), path(assoc), path(clumps), path(logs)

    output:
    tuple val(metaset), val(analysis), path("${analysis}.clumps.tsv")

    script:
    """
    #!Rscript
    library(readr)
    library(dplyr)
    library(tidyr)
    library(stringr)
    library(plyranges)

    meta <- read_tsv("${meta}")

    clump_paths <- str_split("${clumps}", pattern = " ")[[1]]
    names(clump_paths) <- clump_paths

    clump_list <- lapply(clump_paths, read_table, col_types = cols("#CHROM" = col_character()))
    clumps <- bind_rows(clump_list, .id = "dataset") |>
      filter(P <= 5e-8)

    # deduce ranges from tagged SNPs
    clump_ranges <- clumps |>
      mutate(tagged = str_split(SP2, pattern = ",")) |>
      select(dataset, seqnames = `#CHROM`, POS, ID, tagged) |>
      unnest_longer(tagged) |>
      separate(tagged, into = c("chr", "pos", "ref", "alt")) |>
      mutate(pos = as.numeric(pos)) |>
      group_by(dataset, seqnames, POS, ID) |>
      summarize(start = min(pos), end = max(pos)) |>
      ungroup() |>
      mutate(start = if_else(is.na(start), true = POS, false = start),
             end = if_else(is.na(end), true = POS, false = end))

    # genomic ranges
    clump_gr <- as_granges(clump_ranges)
    # MHC region
    mhc <- tibble(seqnames = 6, start = 28510120, end = 33480577)
    mhc_gr <- as_granges(mhc)
    # group into loci
    locus_gr <- reduce_ranges(bind_ranges(clump_gr, mhc_gr))

    clump_loci <- find_overlaps(locus_gr, clump_gr) |>
      as_tibble() |>
      mutate(seqnames = as.character(seqnames)) |>
      select(-ID)

    if(str_detect("${meta}", "mrmega")) {

      # merge clump results with MR-MEGA output
      # keep most significant of each association for each locus
      meta_clumps <- meta |>
        mutate(seqnames = if_else(Chromosome == 23,
                                  true = "X",
                                  false = as.character(Chromosome))) |>
        inner_join(clump_loci, by = c("seqnames" = "seqnames", "Position" = "POS")) |>
        arrange(Chromosome, Position) |>
        group_by(Chromosome, start, end) |>
        mutate(Locus = cur_group_id()) |>
        filter((`P-value_association` == min(`P-value_association`) & 
                `P-value_association` <= 5e-8) | 
               (`P-value_ancestry_het` == min(`P-value_ancestry_het`) & 
                `P-value_ancestry_het` <= 5e-8) |
               (`P-value_residual_het` == min(`P-value_residual_het`) & 
                `P-value_residual_het` <= 5e-8)) |>
        ungroup() |>
        select(-seqnames, -width, -strand, -dataset) |>
        select(Locus, everything())

        write_tsv(meta_clumps, "${analysis}.clumps.tsv")
    } else {
      clump_loci_chr <- clump_loci |>
        mutate(CHROM = str_c("chr", seqnames))

      meta_clumps <- meta |>
        inner_join(clump_loci_chr, by = c("CHROM", "POS")) |>
        mutate(CHR = if_else(CHROM == "chrX", true = 23, false = as.numeric(str_remove(CHROM, "chr")))) |>
        arrange(CHR, POS) |>
        group_by(CHROM, start, end) |>
        mutate(LOCUS = cur_group_id()) |>
        filter(P == min(P)) |>
        ungroup() |>
        select(-CHR, -seqnames, -width, -strand, -dataset) |>
        select(LOCUS, everything())

        write_tsv(meta_clumps, "${analysis}.clumps.tsv")
    }
    """
}
