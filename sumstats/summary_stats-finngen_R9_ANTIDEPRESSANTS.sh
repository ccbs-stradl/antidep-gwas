#!/bin/sh

sumstats=$1
output=$2

    #  1  #chrom
    #  2  pos
    #  3  ref
    #  4  alt
    #  5  rsids
    #  6  nearest_genes
    #  7  pval
    #  8  mlogp
    #  9  beta
    # 10  sebeta
    # 11  af_alt
    # 12  af_alt_cases
    # 13  af_alt_controls

gunzip -c $sumstats | awk 'OFS = "\t" {if(NR == 1) {print "chr", "pos", "ea", "oa", "beta", "se", "pval", "ncase", "ncontrol", "snp", "eaf", "imp_info", "eaf_case", "eaf_control"} else if ($1 == 23) {print "chrX", $2, $4, $3, $9, $10, $7, 131176, 337396, $5, $11, NA, $12, $13} else {print "chr"$1, $2, $4, $3, $9, $10, $7, 131176, 337396, $5, $11, NA, $12, $13}}' > $output

# output columns
# chr
# pos
# ea
# oa
# beta
# se
# pval
# ncase
# ncontrol
# snp
# eaf
# imp_info
# eaf_case
# eaf_control