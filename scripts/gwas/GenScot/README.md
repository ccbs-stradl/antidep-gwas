# Antidepressent exposure GWAS: Generation Scotland

Genome regression and association of antidepressant exposure phenotypes in [Generation Scotland](https://genscot.ed.ac.uk).

## Setup

These scripts are primarily for documentation/reference purposes and filepath setup is not fully described. Derivation of ATC codes requires [NHSBSA dm+d](https://isd.digital.nhs.uk/trud/user/guest/group/0/pack/6/subpack/24/releases) dictionary files downloaded from [NHS TRUD](https://isd.digital.nhs.uk/trud/).

## Notebooks

- `ad-atc.Rmd`: Extracts ATC codes from TRUD download.
- `ad-exposure.Rmd`: AD exposure case/control phenotype from PIS data.

## Workflows

- `gs-topmed-regenie.nf`: Runs regenie on the TOPMed imputation.
