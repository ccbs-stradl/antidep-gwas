#!/bin/sh

sumstats=$1
output=$2

    #  1  v
    #  2  CHR
    #  3  POS
    #  4  SNPID
    #  5  Allele1
    #  6  Allele2
    #  7  AC_Allele2
    #  8  AF_Allele2
    #  9  imputationInfo
    # 10  N
    # 11  BETA
    # 12  SE
    # 13  Tstat
    # 14  p.value
    # 15  p.value.NA
    # 16  Is.SPA.converge
    # 17  varT
    # 18  varTstar
    # 19  AF.Cases
    # 20  AF.Controls

gunzip -c $sumstats | awk 'OFS = "\t" {if(NR == 1) {print "chr", "pos", "ea", "oa", "beta", "se", "pval", "ncase", "ncontrol", "snp", "eaf", "imp_info", "eaf_case", "eaf_control"} else {print "chr"$2, $3, $6, $5, $11, $12, $14, 3288, 175438, $4, $8, $9, $19, $20}}' > $output

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