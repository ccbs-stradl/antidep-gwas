# Documentation for figures, tables, data associated with the manuscript write up

Enables us to keep track of where figures and tables come from.

Maybe include an .Rmd which outputs text with key numeric values, eg. sample sizes. 

## Figures

| Figure Name | Legend (*draft version*) | Location | Script |
| :------: | :------: | :----: | :----: |
| Supplementary Figure XXXX | *chr3 - PLCL2, 1 SNP in CS* | `fineMapping/plots/region_plot_3:16751665:16962555.png` | `fine_mapping_plots.R` |
| Supplementary Figure XXXX | *chr7 - no gene mapped, 2 SNPs in CS* | `fineMapping/plots/region_plot_7:41140716:41340716.png` | `fine_mapping_plots.R` |
| Supplementary Figure XXXX | *chr17 - WNT3, 1 SNP in CS* | `fineMapping/plots/region_plot_17:44752612:44952612.png` | `fine_mapping_plots.R` |
| Supplementary Figure XXXX | *Plot showing SNPs in credible sets are only sig in EUR not AFR or SAS samples.* | `fineMapping/plots/POST-HOC_PROB_POP.png` | `fine_mapping_plots.R` |

## Tables

| Table Name | Legend (*draft version*) | Location | Script |
| :------: | :-------: | :----: | :----: |
| Table X | Overlap of anti-depressant exposure associated genes with genes identified in MDD GWAS (Adams et al. 2025). The denominator is total number of genes identified by each method in the antidepressant exposure GWAS. The numerator is the number of those genes also itentified by that method in the MDD GWAS. The columns represent the anti-depressant meta-analysis or cross-ancestry if using a multi-ancestry method.

| Table Name | Legend (*draft version*) | Location | Script |
| :------: | :-------: | :----: | :----: |
| Supplementary Table XXXX | Significant positionally mapped genes using mBAT-combo (Bonferroni corrected P < 0.05), for all ancestries and antidepressant subgroups. | `antidep-gwas/manuscript/tables/mBAT-combo.csv`  | `antidep-gwas/manuscript/scripts/tables_mBAT-combo.R`
| Supplementary Table XXXX | *Summary file output from susiex limited to significant credible set snps. Columns separated into respective ancestries.* | `antidep-gwas/manuscript/tables/susiex_significant_summary.csv` |  `antidep-gwas/manuscript/scripts/tables_susiex.R` |
| Supplementary Table XXXX | *Credible set file output from susiex limited to significant credible set snps. Columns separated into respective ancestries.* | `antidep-gwas/manuscript/tables/susiex_significant_cs.csv` |  `antidep-gwas/manuscript/scripts/tables_susiex.R` |
