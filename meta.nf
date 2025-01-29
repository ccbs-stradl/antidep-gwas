nextflow.enable.dsl=2

/*
    Meta analyse GWASVCF files
*/

params.vcf = "vcf/gwas/GRCh38/*.{vcf.gz,vcf.gz.tbi,csv}"
params.metasets = "metasets/*.csv"
params.ref = "reference/all_hg38.{pgen,psam,pvar.zst}"
params.genes = "reference/glist_ensgid_hg38_v40_symbol_gene_names.txt"

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

    GENES_CH = Channel
        .fromPath(params.genes)
        
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
        
    /*
        MR-MEGA
        https://genomics.ut.ee/en/tools
    */
    
      // process for input into MR-MEGA
      // group by phenotype and metaset
      // conduct meta-analysis if there are at least 4 cohorts
      MEGA_IN_CH = MEGA_IN(VCF_IN_METASET_CH)
        .transpose(by: 4)
        .map { it -> [[metaset: it[4].metaset, pheno: it[0].pheno]].plus(it) }
        .groupTuple(by: 0)
        .filter { it -> it[1].size() >= 4 }
      
	
      MEGA_CH = MEGA(MEGA_IN_CH)
     
      MEGA_POST_CH = MEGA_POST(MEGA_CH)

      MEGA_ASSOC_CH = MEGA_POST_CH
          .map { it -> it[2] }
          .flatten()

    /*
        Fixed effects (polars)
    */
    
    // process for input into fixed-effects analysis
    // group by phenotype, cluster, and metaset
    // conduct meta-analysis if there are at least 2 cohorts
    FIXED_IN_CH = FIXED_IN(VCF_IN_METASET_CH)
      .transpose(by: 4)
      .map { it -> [[metaset: it[4].metaset, pheno: it[0].pheno, cluster: it[0].cluster]].plus(it) }
      .groupTuple(by: 0)
      .filter { it -> it[1].size() >= 2 }
      
    FIXED_CH = FIXED(FIXED_IN_CH)
    
    FIXED_POST_CH = FIXED_POST(FIXED_CH)

    FIXED_ASSOC_CH = FIXED_POST_CH
        .map { it -> it[2] }
         .flatten()

    /* 
        Post analysis
    */
    
    ASSOC_CH = MEGA_ASSOC_CH
         .concat(FIXED_ASSOC_CH)

    MANHATTAN(ASSOC_CH)

    // count number of GWsig variants
    // further processing of results with 1 or more loci
    ASSOC_GW_CH = GW(ASSOC_CH)
	     .filter { it -> it[1].toInteger() > 0 }
         .map { it -> it[0] }

    // match association and reference SNPs with CPID
    REF_CPID_CH = REF_CPID(REF_CH)

    ASSOC_REF_CH = ASSOC_GW_CH
         .combine(REF_CPID_CH)

    // clumping
    CLUMP_CH = CLUMP(ASSOC_REF_CH)
	
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
    tuple val(details), val(dataset), path(vcf), val(metasets)

    output:
    tuple val(details), val(dataset), path("${dataset}.txt.gz"), path("${dataset}.csv", includeInputs: true), val(metasets)

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
    tag "${meta.metaset}-${meta.pheno}"

    publishDir "meta", pattern: "*.log", mode: "copy"

    cpus = 2
    memory = 48.GB
    time = '3h'

    input:
    tuple val(meta), val(details), val(datasets), path(gwas), path(csv), val(metaset)

    output:
    tuple val("${meta.metaset}-mrmega-${meta.pheno}"), path("*.result"), path("*.log"), path(csv, includeInputs: true)

    script:
    """
    ls *.txt.gz > mrmega.in

    MR-MEGA --filelist mrmega.in \
    --pc ${[gwas.size() - 3, 3].min()} \
    --out ${meta.metaset}-mrmega-${meta.pheno}
    """
}

process MEGA_POST {
    tag "${result}"
    label 'analysis'

    publishDir "meta", pattern: "*.{gz,csv}", mode: 'copy'

    cpus = 2
    memory = 32.GB
    time = '30m'

    input:
    tuple val(meta), path(result), path(log), path(datasets)

    output:
    tuple path("${meta}.gz"), path("${meta}.csv"), path("*.assoc")

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
    arrange(Chromosome, Position)

    mega_assoc <- mega_p |>
        transmute(CHR=Chromosome, SNP=MarkerName, BP=Position, A1=EA, A2=NEA,
                  P=`P-value_association`, BETA=beta_0, SE=se_0, NMISS=Nsample,
                  CPID = str_glue("{CHR}:{BP}:{A2}:{A1}"))

    mega_ancestry <- mega_p |>
        transmute(CHR=Chromosome, SNP=MarkerName, BP=Position, A1=EA, A2=NEA,
                   P=`P-value_ancestry_het`, BETA=beta_0, SE=se_0, NMISS=Nsample,
                   CPID = str_glue("{CHR}:{BP}:{A2}:{A1}"))

    mega_residual <- mega_p |>
        transmute(CHR=Chromosome, SNP=MarkerName, BP=Position, A1=EA, A2=NEA,
                   P=`P-value_residual_het`, BETA=beta_0, SE=se_0, NMISS=Nsample,
                   CPID = str_glue("{CHR}:{BP}:{A2}:{A1}"))


    write_tsv(mega_p, "${meta}.gz")
    write_csv(datasets, "${meta}.csv")

    write_tsv(mega_assoc, "${meta}.assoc")
    write_tsv(mega_ancestry, "${meta}.ancestry.assoc")
    write_tsv(mega_residual, "${meta}.residual.assoc")
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
    tuple val(details), val(dataset), path(vcf), val(metasets)

    output:
    tuple val(details), val(dataset), path("${dataset}.query.txt"), path("${dataset}.csv", includeInputs: true), val(metasets)

    script:
    """
    # bcf query to get required columns
    # remove variants with MAF < 0.5%, INFO < 0.4
    # rename missing variants names to CPID
    bcftools norm -m- ${dataset}.vcf.gz |\
    bcftools query -HH \
    -i 'FORMAT/AF >= 0.005 & FORMAT/AF <= 0.995 & (FORMAT/SI >= 0.4 | FORMAT/SI = ".")' \
	  -f "%CHROM\\t%POS\\t%ID\\t%ALT\\t%REF\\t[%ES]\\t[%SE]\\t[%LP]\\t[%SI]\\t[%AFCAS]\\t[%AFCON]\\t[%SS]\\t[%NC]\\t[%NE]\\n" |\
    awk '{if(\$1 == ".") \$1 = \$13":"\$14; print \$0}' > ${dataset}.query.txt
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
    tag "${meta.metaset}-${meta.pheno}-${meta.cluster}"
    label 'analysis'

    cpus = 6
    memory = 48.GB
    time = '30m'

    input:
    tuple val(meta), val(details), val(datasets), path(gwas), path(csv), val(metaset)

    output:
    tuple val("${meta.metaset}-fixed-${meta.pheno}-${meta.cluster}"), path("meta.tsv"), path(csv, includeInputs: true)

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
    )
    
    meta_het = (
      meta.join(het, on = ["CHROM", "POS", "ID", "ALT", "REF"], how = "inner")
      .collect()
    )
    
    meta_het.write_csv("meta.tsv", separator = "\\t")
    """
}

// Fixed effects postprocessing
process FIXED_POST {
    tag "${dataset}"
    label 'analysis'

    publishDir "meta", pattern: "*.{gz,csv}", mode: 'copy'

    cpus = 2
    memory = 24.GB
    time = '30m'

    input:
    tuple val(dataset), path(meta), path(datasets)

    output:
    tuple path("${dataset}.gz"), path("${dataset}.csv"), path("*.assoc")

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
      transmute(CHR = str_remove(CHROM, 'chr'), SNP = ID, BP = POS,
                A1 = ALT, A2 = REF, P = QP, BETA, ESS=NEFF,
                CPID = str_glue("{CHR}:{BP}:{A2}:{A1}"))

    write_tsv(meta_p, "${dataset}.gz")
    write_csv(datasets, "${dataset}.csv")
    write_tsv(assoc, "${dataset}.assoc")
    write_tsv(het, "${dataset}.het.assoc")
    """
}

/* gwsig variants */
process GW {
    tag "${assoc.baseName}"

    cpus = 1
    memory = 1.GB
    time = '10m'

    input:
    path(assoc)

    output:
    tuple path(assoc, includeInputs: true), env(GW)

    script:
    """
    GW=\$(cat ${assoc} | awk '\$6 <= 5e-8' | wc -l)
    """

}

/* manhattan plot */
process MANHATTAN {
	tag "${assoc}"
        label 'analysis'

	publishDir 'meta', mode: 'copy'
	
	cpus = 1
	memory = 16.GB
	time = '30m'
	
	input:
	path(assoc)
	
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

    publishDir 'meta', mode: 'copy'

    //errorStrategy 'ignore'

    cpus = 4
    memory = 16.GB
    time = '1h'

    input:
    tuple path(assoc), path(pgen), path(psam), path(pvar)

    output:
    tuple path("${assoc.baseName}.clumps"), path("${assoc.baseName}.log")

    script:
    """
    plink2 \
    --clump ${assoc} \
    --clump-id-field CPID \
    --clump-p1 0.0001 \
    --clump-p2 1.0 \
    --clump-r2 0.1 \
    --clump-kb 1000 \
    --clump-allow-overlap \
    --pgen ${pgen} \
    --psam ${psam} \
    --pvar ${pvar} \
    --allow-extra-chr \
    --out ${assoc.baseName} \
	  --threads ${task.cpus} \
	  --memory ${task.memory.bytes.intdiv(1000000)}
    """
}


