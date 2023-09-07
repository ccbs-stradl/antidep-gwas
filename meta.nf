nextflow.enable.dsl=2

/*
    Meta analyse GWASVCF files
*/

params.vcf = "vcf/*.vcf.gz"
params.ref = "reference/all_hg38.{pgen,psam,pvar.zst}"
params.genes = "reference/glist_ensgid_hg38_v40_symbol_gene_names.txt"

workflow {
    VCF_CH = Channel
        .fromPath(params.vcf) 

    REF_CH = Channel
        .fromFilePairs(params.ref, size: 3)

    GENES_CH = Channel
        .fromPath(params.genes)
    /*
        MR-MEGA
        https://genomics.ut.ee/en/tools
    */

    // key to phenotype
    MEGA_VCF_CH = VCF_CH
		.map { it -> [it.baseName.split("-")[1], it] }

    // conduct meta-analysis if there are at least 4 cohorts
    MEGA_IN_CH = MEGA_IN(MEGA_VCF_CH)
        .groupTuple()
		.filter { it -> it[1].size() >= 4 }
	
    MEGA_CH = MEGA(MEGA_IN_CH)
    
    MEGA_POST_CH = MEGA_POST(MEGA_CH)

    MEGA_ASSOC_CH = MEGA_POST_CH
        .map { it -> it[1] }
        .flatten()

    /*
        Fixed effects (PLINK)
        https://www.cog-genomics.org/plink/1.9/postproc#meta_analysis
    */

    // key to phenotype and cluster
    FIXED_VCF_CH = VCF_CH
        .map { it -> [ it.simpleName.split("-"), it ] }
        .map { it -> [ ["pheno": it[0][1], "cluster": it[0][2]], it[1] ] }

    // conduct meta-analysis if there are at least 2 cohorts
    FIXED_IN_CH = FIXED_IN(FIXED_VCF_CH)
        .groupTuple()
        .filter { it -> it[1].size() >= 2 }
    
    FIXED_CH = FIXED(FIXED_IN_CH)
    
    FIXED_POST_CH = FIXED_POST(FIXED_CH)

    FIXED_ASSOC_CH = FIXED_POST_CH
        .map { it -> it[1] }
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

    ASSOC_REF_CH = ASSOC_GW_CH
        .combine(REF_CH)
    
    ASSOC_BED_CH = REF_BED(ASSOC_REF_CH)
    
    ASSOC_GENES_CH = ASSOC_BED_CH
         .combine(GENES_CH)

    CLUMP_CH = CLUMP(ASSOC_GENES_CH)

    //MA_CH = MA(MEGA_ASSOC_CH)
	
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
    tag "${vcf}"

    cpus = 1
    memory = 1.GB
    time = '30m'

    input:
    tuple val(pheno), path(vcf)

    output:
    tuple val(pheno), path("${vcf.simpleName}.txt.gz")

    script:
    """
    # bcf query to get required columns
    # remove variants with MAF < 0.5%, INFO < 0.4
    # rename missing variants names to CPID
    bcftools query \
    -i 'FORMAT/AF >= 0.005 & FORMAT/AF <= 0.995 & (FORMAT/SI >= 0.4 | FORMAT/SI = ".")' \
	-f "%ID\\t%ALT\\t%REF\\t[%ES]\\t[%SE]\\t[%AF]\\t[%SS]\\t[%NC]\\t%CHROM\\t%POS\\n" \
	${vcf} | awk '{if(\$1 == ".") \$1 = \$9":"\$10; print \$0}' > ${vcf.simpleName}.query.txt

    # header
    echo "" | awk 'OFS = "\\t" {print "MARKERNAME", "EA", "NEA", "OR", "OR_95L", "OR_95U", "EAF", "N", "CHROMOSOME", "POSITION"}' > ${vcf.simpleName}.txt

    # reformat
    # remove "chr" from chromosome names
    # turn beta effect size back into OR
    # calculate N as effective sample size
    cat ${vcf.simpleName}.query.txt | awk \
    'OFS = "\\t" {gsub("chr", "", \$9); print \$1, \$2, \$3, exp(\$4), exp(\$4-1.96*\$5), exp(\$4+1.96*\$5), \$6, 4*\$8*(\$7-\$8)/(\$7), \$9, \$10}' >> ${vcf.simpleName}.txt

    # compress
    gzip ${vcf.simpleName}.txt
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
     9	%CHROM
    10	%POS
*/
}

process MEGA {
    tag "${pheno}"

    publishDir "meta", pattern: "*.log", mode: "copy"

    cpus = 2
    memory = 32.GB
    time = '3h'

    input:
    tuple val(pheno), path(gwas)

    output:
    tuple path("*.result"), path("*.log")

    script:
    """
    ls *.txt.gz > mrmega.in

    MR-MEGA --filelist mrmega.in \
    --pc ${[gwas.size() - 3, 3].min()} \
    --out mrmega-${pheno}
    """
}

process MEGA_POST {
    tag "${result}"

    publishDir "meta", pattern: "*.gz", mode: 'copy'

    cpus = 2
    memory = 32.GB
    time = '30m'

    input:
    tuple path(result), path(log)

    output:
    tuple path("${result}.gz"), path("*.assoc")

    script:
    """
    #!Rscript
	library(dplyr)
	library(readr)
	library(fastman)
	
	mega <- read_tsv("${result}")

    # recalculate P-values to allow for smaller values
    mega_p <- mega |>
        filter(!is.na(beta_0)) |>
        mutate(`P-value_association` = pchisq(chisq_association, ndf_association, lower.tail=FALSE),
               `P-value_ancestry_het` = pchisq(chisq_ancestry_het, ndf_ancestry_het, lower.tail=FALSE),
               `P-value_residual_het` = pchisq(chisq_residual_het, ndf_residual_het, lower.tail=FALSE)
        ) |>
    arrange(Chromosome, Position) |>
    mutate(Chromosome = if_else(Chromosome == 23, true = 'X', false = as.character(Chromosome)))

    mega_assoc <- mega_p |>
        transmute(CHR=Chromosome, SNP=MarkerName, BP=Position, A1=EA, A2=NEA,
                  P=`P-value_association`, BETA=beta_0, SE=se_0, NMISS=Nsample)
    mega_ancestry <- mega_p |>
        transmute(CHR=Chromosome, SNP=MarkerName, BP=Position, A1=EA, A2=NEA,
                   P=`P-value_ancestry_het`, BETA=beta_0, SE=se_0, NMISS=Nsample)
    mega_residual <- mega_p |>
        transmute(CHR=Chromosome, SNP=MarkerName, BP=Position, A1=EA, A2=NEA,
                   P=`P-value_residual_het`, BETA=beta_0, SE=se_0, NMISS=Nsample)


    write_tsv(mega_p, "${result}.gz")

    write_tsv(mega_assoc, "${result}.assoc")
    write_tsv(mega_ancestry, "${result}.ancestry.assoc")
    write_tsv(mega_residual, "${result}.residual.assoc")
    """

}


/* Query VCF to create PLINK .assoc input file 
    1) SNP – snp name
    2) A1 – effect allele
    3) A2 – non effect allele
    4) OR - odds ratio
    5) SE - standard error of log(OR)
    6) P - p-value
    7) ESS - effective sample size
    8) CHR  - chromosome of marker
    9) BP - position of marker
*/
process FIXED_IN {
    tag "${vcf}"

    cpus = 1
    memory = 1.GB
    time = '30m'

    input:
    tuple val(dataset), path(vcf)

    output:
    tuple val(dataset), path("${vcf.simpleName}.assoc")

    script:
    """
    # bcf query to get required columns
    # remove variants with MAF < 0.5%, INFO < 0.4
    # rename missing variants names to CPID
    bcftools query \
    -i 'FORMAT/AF >= 0.005 & FORMAT/AF <= 0.995 & (FORMAT/SI >= 0.4 | FORMAT/SI = ".")' \
	-f "%ID\\t%ALT\\t%REF\\t[%ES]\\t[%SE]\\t[%LP]\\t[%AF]\\t[%SS]\\t[%NC]\\t%CHROM\\t%POS\\n" \
	${vcf} | awk '{if(\$1 == ".") \$1 = \$10":"\$11; print \$0}' > ${vcf.simpleName}.query.txt

    # header
    echo "" | awk 'OFS = "\\t" {print "SNP", "A1", "A2", "OR", "SE", "P", "ESS", "CHR", "BP"}' > ${vcf.simpleName}.assoc

    # reformat
    # turn beta effect size back into OR
    # turn -log10(p) back into p
    # calculate N as effective sample size
    cat ${vcf.simpleName}.query.txt | awk \
    'OFS = "\\t" {print \$1, \$2, \$3, exp(\$4), \$5, 10^(-\$6), 4*\$9*(\$8-\$9)/(\$8), \$10, \$11}' >> ${vcf.simpleName}.assoc
    """
/*
    column numbers of bcf query
     1	%ID
     2	%ALT
     3	%REF
     4	[%ES]
     5	[%SE]
     6  [%LP]
     7	[%AF]
     8	[%SS]
     9	[%NC]
    10	%CHROM
    11	%POS
*/
}

process FIXED {
    tag "${dataset.pheno}-${dataset.cluster}"

    publishDir "meta", pattern: "*.log", mode: "copy"

    cpus = 4
    memory = 8.GB
    time = '1h'

    input:
    tuple val(dataset), path(gwas)

    output:
    tuple val(dataset), path("*.meta"), path("*.log")

    script:
    """
    plink \
    --meta-analysis ${gwas} \
    --output-chr 'M' \
    --out fixed-${dataset.pheno}-${dataset.cluster} \
	--threads ${task.cpus} \
	--memory ${task.memory.bytes.intdiv(1000000)}
    """
}

process FIXED_POST {
    tag "${meta}"

    publishDir "meta", pattern: "*.gz", mode: 'copy'

    cpus = 1
    memory = 16.GB
    time = '10m'

    input:
    tuple val(dataset), path(meta), path(log)

    output:
    tuple path("${meta}.gz"), path("*.assoc")

    script:
    """
    #!Rscript
	library(dplyr)
	library(readr)

    meta <- read_table("${meta}", col_types=cols(CHR=col_character()))

    assoc <- meta |>
        select(CHR, SNP, BP, A1, A2, P, OR)
    
    het <- meta |>
        filter(!is.na(Q)) |>
        select(CHR, SNP, BP, A1, A2, P=Q, OR)

    write_tsv(meta, "${meta}.gz")
    write_tsv(assoc, "${meta.baseName}.assoc")
    write_tsv(het, "${meta.baseName}.het.assoc")
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
	
	assoc <- read_tsv("${assoc}", col_types=cols(CHR = col_character()))

    assoc_p <- assoc |> filter(P > 0)
	
	png("${assoc.baseName}.manhattan.png", width=10, height=6, units="in", res=300)
	fastman(assoc_p, chr = "CHR", bp = "BP", p = "P", maxP=NULL)
	dev.off()
	
	png("${assoc.baseName}.qq.png", width=6, height=6, units="in", res=300)
	fastqq(p1=pull(assoc_p, P), maxP=NULL)
	dev.off()
	"""
}

process REF_BED {
    tag "${assoc.baseName}-${ref}"

    cpus = 1
    memory = 8.GB
    time = '10m'

    input:
    tuple path(assoc), val(ref), path(pgen)

    output:
    tuple path(assoc, includeInputs: true), val("ref"), path("ref.{bed,bim,fam}")

    script:
    """
    # chr/pos to extract
    cat ${assoc} | awk 'NR > 1 {print \$1, \$3, \$3}' > ${assoc.baseName}.bed1
    # rename CPIDs to rsID
    cat ${assoc} | awk 'NR > 1 {print \$1":"\$3":"\$5":"\$4, \$2}' > ${assoc.baseName}.names

    plink2 \
    --make-bed \
    --pfile 'vzs' ${ref} \
    --extract 'bed1' ${assoc.baseName}.bed1 \
    --set-all-var-ids @:#:\\\$r:\\\$a \
    --new-id-max-allele-len 500 error \
    --out ref-cpid \
    --allow-extra-chr \
	--threads ${task.cpus} \
	--memory ${task.memory.bytes.intdiv(1000000)}

    plink2 \
    --make-bed \
    --bfile ref-cpid \
    --update-name ${assoc.baseName}.names \
    --out ref \
	--threads ${task.cpus} \
	--memory ${task.memory.bytes.intdiv(1000000)}
    """
}

process CLUMP {
    tag "${assoc.baseName}"

    publishDir 'meta', mode: 'copy'

    //errorStrategy 'ignore'

    cpus = 8
    memory = 8.GB
    time = '1h'

    input:
    tuple path(assoc), val(ref), path(bed), path(genes)

    output:
    tuple path("${assoc.baseName}.clumped"), path("${assoc.baseName}.clumped.ranges"), path("${assoc.baseName}.log")

    script:
    """
    tail -n +2 ${genes} > ${genes.baseName}.bed1

    plink \
    --clump ${assoc} \
    --clump-range ${genes.baseName}.bed1 \
    --clump-range-border 1000 \
    --clump-p1 5e-8 \
    --clump-p2 5e-5 \
    --clump-r2 0.4 \
    --clump-kb 3000 \
    --bfile ${ref} \
    --out ${assoc.baseName} \
	--threads ${task.cpus} \
	--memory ${task.memory.bytes.intdiv(1000000)}
    """
}

process MA {
    tag "${assoc}"

    cpus = 1
    memory = 1.GB
    time = '10m'

    input: 
    path(assoc)

    output:
    path("${assoc.baseName}.ma")

    script:
    """
    cat ${assoc} | awk '{OFS = "\\t"; if(NR == 1) {print "SNP", "A1", "A2", "freq", "BETA", "SE", "P", "N"} else {print \$2, \$4, \$5, }}'
    """
}

// (CHR, SNP, BP, A1=EA, A2=NEA, BETA=beta_0, SE=se_0, NMISS=Nsample, P=`P-value_association`)