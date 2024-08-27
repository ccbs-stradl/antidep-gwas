params.bfile = "gs://fc-aou-datasets-controlled/v7/microarray/plink_v7.1/arrays.{bed,bim,fam}"
params.flagged = "gs://fc-aou-datasets-controlled/v7/wgs/short_read/snpindel/aux/qc/flagged_samples.tsv"
params.ancestry = "gs://fc-aou-datasets-controlled/v7/wgs/short_read/snpindel/aux/ancestry/ancestry_preds.tsv"
params.relatedness = "gs://fc-aou-datasets-controlled/v7/wgs/short_read/snpindel/aux/relatedness/relatedness.tsv"
// high LD regions file from https://github.com/gabraham/flashpca/blob/master/exclusion_regions_hg19.txt
params.regions = "https://raw.githubusercontent.com/gabraham/flashpca/master/exclusion_regions_hg19.txt"
params.seed = 64074892

workflow {
    // genetic data
    BFILE_CH = Channel
        .fromFilePairs(params.bfile, size: 3) { it -> it.baseName }
    
    // filter file for cluster assignment of each sample
    ANCESTRY_CH = Channel
        .fromPath(params.ancestry)
    
    // parse cluster file to get names of each cluster
    CLUSTER_NAMES_CH = ANCESTRY_CH
        .splitCsv(sep: "\t", skip: 1, header: ['iid', 'cluster', 'probs', 'features', 'cluster_other'])
        .map { it -> it.cluster }
        .unique()
    
    // flagged samples to remove
    FLAGGED_CH = Channel
        .fromPath(params.flagged)
    
    // related samples
    RELATEDNESS_CH = Channel
        .fromPath(params.relatedness)
    
    // regions to exclude
    REGIONS_CH = Channel
        .fromPath(params.regions)
    
    // format keep/remove files
    SAMPLES_CH = SAMPLES(FLAGGED_CH, RELATEDNESS_CH, ANCESTRY_CH)

    // QC genotypes
    QC_CH = QC(BFILE_CH, REGIONS_CH)

    // extract and QC genotypes for each cluster
    CLUSTERS_QC_CH = CLUSTER(CLUSTER_NAMES_CH, QC_CH, SAMPLES_CH)

    // prune samples and SNPs
    CLUSTERS_PRUNE_CH = PRUNE(CLUSTERS_QC_CH)

    // calculate PCs and project
    CLUSTERS_PCS_CH = PCS(CLUSTERS_PRUNE_CH)
    CLUSTERS_PROJECT_CH = PROJECT(CLUSTERS_PCS_CH)
    
}

// process flagged samples and clusters into sample filter files
process SAMPLES {
    
    cpus = 1
    memory = 1.GB
    time = '10m'
    disk = 1.GB
    
    container = 'rocker/tidyverse:4.4.1'
    
    input:
    path(flagged)
    path(relatedness)
    path(ancestry)
    
    output:
    tuple path("remove.txt"), path("clusters.txt")
    
    script:
    """
    #!Rscript
    
    library(dplyr)
    library(readr)

    flagged_samples <- read_tsv("${flagged}")
    relatedness <- read_tsv("${relatedness}")
    ancestry <- read_tsv("${ancestry}")
    
    removes <- bind_rows(
        tibble(FID = 0, IID = pull(flagged_samples, s)),
        tibble(FID = 0, IID = pull(relatedness, i.s)),
        tibble(FID = 0, IID = pull(relatedness, j.s)),
        ) |>
        distinct()

    write_tsv(removes, "remove.txt")

    clusters <- ancestry |>
      transmute(FID = 0, IID = research_id, cluster = ancestry_pred)

    write_tsv(clusters, "clusters.txt")
    """
}


process QC {
    
    cpus = 1
    memory = 6.GB
    disk = 250.GB
    time = '2h'
    
    input:
    tuple val(bfile), path(bedbimfam)
    path(regions)
    
    output:
    tuple val("qc"), path("qc.{bed,bim,fam}")
    
    script:
    """
    plink2 \
    --make-bed \
    --bfile ${bfile} \
    --maf 0.05 \
    --geno 0.02 \
    --chr 1-22 \
    --exclude range ${regions} \
    --min-alleles 2 \
    --max-alleles 2 \
    --snps-only 'just-acgt' \
    --threads ${task.cpus} \
	--memory ${task.memory.bytes.intdiv(1000000)} \
    --out qc
    """
}

process CLUSTER {
    tag "${cluster}"
    
    cpus = 1
    memory = 6.GB
    disk = 150.GB
    time = '2h'
    
    input:
    each cluster
    tuple val(bfile), path(bedbimfam)
    tuple path(remove), path(clusters)
    
    output:
    tuple val(cluster), path("${cluster}.{bed,bim,fam}"), path(remove, includeInputs: true)
    
    script:
    """
    plink2 \
    --make-bed \
    --bfile ${bfile} \
    --keep-col-match ${clusters} ${cluster} \
    --keep-col-match-num 3 \
    --maf 0.05 \
    --geno 0.02 \
    --hwe 1e-10 \
    --threads ${task.cpus} \
	--memory ${task.memory.bytes.intdiv(1000000)} \
    --out ${cluster}
    """
}

// get unrelated random sample
process KING {
tag "${cluster}"
    
    cpus = 8
    memory = 32.GB
    disk = 150.GB
    time = '3h'
    
    input:
    each cluster
    tuple val(bfile), path(bedbimfam), path(keep), path(include)
    tuple path(remove), path(clusters)
    
    output:
    tuple val(cluster), val(bfile), path(bedbimfam), path(include, includeInputs: true), path("${cluster}.king.cutoff.in.id"), path(clusters, includeInputs: true)
    
    script:
    """
    plink2 \
    --king-cutoff 0.0442 \
    --bfile ${bfile} \
    --keep keep \
    --keep-col-match ${clusters} ${cluster} \
    --keep-col-match-num 3 \
    --thin-indiv-count 10000 \
    --threads ${task.cpus} \
	--memory ${task.memory.bytes.intdiv(1000000)} \
    --out ${cluster}
    """
}

// prune SNPs
process PRUNE {
tag "${cluster}"
    
    cpus = 8
    memory = 32.GB
    disk = 150.GB
    time = '1h'
    
    input:
    tuple val(cluster), path(bedbimfam), path(remove)
    
    output:
    tuple val(cluster), path(bedbimfam, includeInputs: true), path("${cluster}.id"), path("${cluster}.prune.in")
    
    script:
    """
    plink2 \
    --indep-pairwise 1000 80 0.1 \
    --write-samples \
    --bfile ${cluster} \
    --remove ${remove} \
    --thin-indiv-count 1000 \
    --threads ${task.cpus} \
    --memory ${task.memory.bytes.intdiv(1000000)} \
    --seed ${params.seed} \
    --out ${cluster}
    """
}

// calculate PCs in pruned samples
process PCS {
tag "${cluster}"
    
    cpus = 8
    memory = 32.GB
    disk = 150.GB
    time = '4h'
    
    input:
    tuple val(cluster), path(bedbimfam), path(keep), path(extract)
    
    output:
    tuple val(cluster), path(bedbimfam, includeInputs: true), path("${cluster}.acount"), path("${cluster}.eigenvec.allele")
    
    script:
    """
    plink2 \
    --freq counts \
    --pca allele-wts \
    --bfile ${cluster} \
    --keep ${keep} \
    --extract ${extract} \
    --out ${cluster} \
	--threads ${task.cpus} \
	--memory ${task.memory.bytes.intdiv(1000000)}
    """
}

// project cluster into PCs
process PROJECT {
tag "${cluster}"
    
    cpus = 8
    memory = 32.GB
    disk = 150.GB
    time = '4h'

    publishDir "pcs"
    
    input:
    tuple val(cluster), path(bedbimfam), path(acount), path(eigenvec)
    
    output:
    tuple val(cluster), path("${cluster}-projection.sscore")
    
    script:
    """
    plink2 --bfile ${cluster} \
    --read-freq ${acount} \
    --score ${eigenvec} 2 6 header-read no-mean-imputation \
        variance-standardize \
    --score-col-nums 7-16 \
    --out ${cluster}-projection \
    --threads ${task.cpus} \
    --memory ${task.memory.bytes.intdiv(1000000)}
    """
}