# GWAS Meta-analysis of Antidepressant Prescribing

This repository contains scripts, documentation and configuration files for conducting a genome-wide association study (GWAS) meta-analysis of antidepressant exposure.

## Repository Structure

-   `docs/`: Contains documentation for running the pipelines
    -   `workflow_documentation.md` step-by-step instructions for setting up the environment, downloading necessary files, and running the various Nextflow pipelines to reproduce the GWAS meta-analysis of antidepressant exposure and downstream analyses.
    -   `workflow_references.md` instructions to download required reference files for the workflows.
    -   `workflow_requirements.md` instructions for installing any requirements.
-   `workflows/` : Contains Nextflow scripts for running the pipelines.
-   `scripts/`: Contains R markdown scripts (and their outputs), jupyter notebooks, bash shell scripts for data pre-processeing and formatting/plotting outputs from the Nextflow scripts. See the README.md in this directory for more info on what each script does.
-   `inputs/`
    -   `datasets.csv`: table of input GWAS with cohort, phenotype, cluster, build and version information
    -   `metasets`: tables of cohort/versions used for each version of the meta-analysis
    -   `gwas.json`: formatting file for VCF conversion
-   `results/`: Contains outputs from the nextflow pipeline scripts. Sub-folders are made inside this directory for each nextflow pipeline. Most of these subfolders are not git commited, instead symlinks should be made to point to their contents (see `docs/workflow_documentation.md`)
-   `config/`: Contains nextflow configuration files
-   `env/`: Contains yaml files for setting up virtual conda environments.
-   `manuscript/` : Contains scripts and outputs used directly for the manuscript, eg. creating tables and pasting sentences with sample size numbers that can be copied into the main text of the manuscript.
-   `reference/` : Contains reference files needed in nextflow scripts, see `docs/workflow_references.md`

## Additional instructions for Edinburgh University colleagues

Results from Nextflow pipelines are on DataStore here: `GenScotDepression/data/AMBER/antidep-gwas`.

Symlinks to the directories `fineMapping/`, `maps/` `meta/` & `models/` need to be set up in the root of this repo for some scripts in `scripts` and `manuscript` to run.

## Contact

For questions or collaboration inquiries, please contact Mark Adams ([mark.adams\@ed.ac.uk](mailto:mark.adams@ed.ac.uk)).
