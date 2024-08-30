# Antidepressant exposure GWAS: All of Us

Genome regression and association of antidepressant exposure phenotypes in [All of Us](https://allofus.nih.gov/). Analysis conducted in the Researcher Workbench with dataset "All of Us Controlled Tier Dataset v7" under the project "[Genetics of antidepressant exposure](https://www.researchallofus.org/research-projects-directory/?searchBy=workspaceNameLike&directorySearch=Genetics+of+antidepressant+exposure)"

## Notebooks

- `AD exposure pheno.ipynb`: Create phenotype and covariate files based on cohorts and datasets constructed using the Workbench.
- `AD exposure sample sizes.ipynb`: Summary of GWAS sample sizes.

## Queries

- `queries`: SQL queries to define EHR, N06A, N06AA, and N06AB datasets.

## Workflows

- `pcs.nf`: Calculate principal components within each pre-defined genetic similarity cluster ("ancestry")
- `regenie.nf`: GWAS of each cluster using [regenie](https://rgcgithub.github.io/regenie/).