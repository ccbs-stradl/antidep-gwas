# Documentation for figures, tables, data associated with the manuscript write up

Enables us to keep track of where figures and tables come from.

Maybe include an .Rmd which outputs text with key numeric values, eg. sample sizes.

## Figures

|        Figure Name        |                             Legend (_draft version_)                             |                         Location                         |         Script         |
| :-----------------------: | :------------------------------------------------------------------------------: | :------------------------------------------------------: | :--------------------: |
| Supplementary Figure XXXX |                           _chr3 - PLCL2, 1 SNP in CS_                            | `fineMapping/plots/region_plot_3:16751665:16962555.png`  | `fine_mapping_plots.R` |
| Supplementary Figure XXXX |                      _chr7 - no gene mapped, 2 SNPs in CS_                       | `fineMapping/plots/region_plot_7:41140716:41340716.png`  | `fine_mapping_plots.R` |
| Supplementary Figure XXXX |                           _chr17 - WNT3, 1 SNP in CS_                            | `fineMapping/plots/region_plot_17:44752612:44952612.png` | `fine_mapping_plots.R` |
| Supplementary Figure XXXX | _Plot showing SNPs in credible sets are only sig in EUR not AFR or SAS samples._ |        `fineMapping/plots/POST-HOC_PROB_POP.png`         | `fine_mapping_plots.R` |

## Tables

| Table Name |             Legend (_draft version_)       |        Location                |           Script         |
| :--------: | :----------------: | :---------------------------: | :------------------------: |
|  Table X   |     Summary of cohorts included in MR-MEGA meta-analysis. Each row contains sample size information for a meta-analysis of each anti-depressant phenotype. This includes the number of cases, controls and effective sample size (neff) are in each cohort: All of Us, Biobank Japan (BBJ), FinnGen, Generation Scotland (GenScot) and UK Biobank (UKB). Cells with NA mean that cohort was not included in that meta-analysis.  | `antidep-gwas/manuscript/tables/meta_analysis_cohort_summary_table_mrmega.csv` | `antidep-gwas/manuscript/scripts/meta_analysis_cohort_summary_table.R` |
|  Table X   |   Summary of cohorts included in fixed meta-analysis. Each row contains sample size information for a meta-analysis of each ancestry cluster and anti-depressant phenotype combination. This includes the number of cases, controls and effective sample size (neff) are in each cohort: All of Us, Biobank Japan (BBJ), FinnGen, Generation Scotland (GenScot) and UK Biobank (UKB). Cells with NA mean that cohort was not included in that meta-analysis.  | `antidep-gwas/manuscript/tables/meta_analysis_cohort_summary_table_fixed.csv`  | `antidep-gwas/manuscript/scripts/meta_analysis_cohort_summary_table.R` |
|  Table X   | Overlap of anti-depressant exposure associated genes with genes identified in MDD GWAS (Adams et al. 2025). The denominator is total number of genes identified by each method in the antidepressant exposure GWAS. The numerator is the number of those genes also itentified by that method in the MDD GWAS. The columns represent the anti-depressant meta-analysis or cross-ancestry if using a multi-ancestry method. |    `antidep-gwas/manuscript/tables/antidep_gwas_mdd_gwas_summary_table.csv`    |        `antidep-gwas/manuscript/scripts/MDD_GWAS_gene_overlap*`        |

|        Table Name        |                                                         Legend (_draft version_)                                                          |                            Location                             |                        Script                         |
| :----------------------: | :---------------------------------------------------------------------------------------------------------------------------------------: | :-------------------------------------------------------------: | :---------------------------------------------------: |
| Supplementary Table XXXX | Significant positionally mapped genes using mBAT-combo (Bonferroni corrected P \< 0.05), for all ancestries and antidepressant subgroups. |         `antidep-gwas/manuscript/tables/mBAT-combo.csv`         | `antidep-gwas/manuscript/scripts/tables_mBAT-combo.R` |
| Supplementary Table XXXX |         _Summary file output from susiex limited to significant credible set snps. Columns separated into respective ancestries._         | `antidep-gwas/manuscript/tables/susiex_significant_summary.csv` |   `antidep-gwas/manuscript/scripts/tables_susiex.R`   |
| Supplementary Table XXXX |      _Credible set file output from susiex limited to significant credible set snps. Columns separated into respective ancestries._       |   `antidep-gwas/manuscript/tables/susiex_significant_cs.csv`    |   `antidep-gwas/manuscript/scripts/tables_susiex.R`   |
