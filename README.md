# WGS Pre-Processing and GWAS Pipeline for Sepsis Survival
[![DOI](https://zenodo.org/badge/1457159910.svg)](https://doi.org/10.5281/zenodo.24280139)
This repository contains the upstream bash pipeline used to pre-process, genotype, filter, and annotate Whole-Genome Sequencing (WGS) data for 1,065 sepsis patients from three cohorts (APROCCHSS, RECORDS Observational, and TRIAL cohorts) aligned to GRCh38. It also includes the downstream PLINK quality control and PCA computation used to define the final n=1,030 analysis cohort.

## Pipeline Architecture

The pipeline relies on GATK4 Best Practices, `bcftools`, and `plink`, and is modularized into sequential steps:

### Phase 1: VCF Generation & Recalibration
* **`01_gvcf_prep.sh`**: Downloads raw `.g.vcf.gz` files, verifies MD5 checksums, normalizes multiallelic sites, and filters out ENCODE blacklisted genomic regions. Processes samples concurrently in batches to optimize compute time.
* **`02_joint_genotyping.sh`**: Iterates through precise genomic windows to import samples into a GenomicsDB workspace and performs joint genotyping across the cohort using `GenotypeGVCFs` in parallel batches.
* **`03_vqsr_filtering.sh`**: Applies ExcessHet hard-filtering and Variant Quality Score Recalibration (VQSR) using Broad Institute GRCh38 resource bundles for both SNPs and INDELs. 

### Phase 2: Annotation
* **`04_bcf_annotation.sh`**: Converts the recalibrated VCF to binary BCF format, filters for `PASS` variants, and splits the dataset by variant type. Annotates variants with dbSNP (rsIDs), CADD v1.7 scores, GnomAD v4.1 allele frequencies, and Ensembl consequences (`bcftools csq`).

### Phase 3: GWAS Quality Control & PCA
* **`05_plink_qc.sh`**: Converts the dataset to PLINK binary format and enforces strict GWAS quality control thresholds (`--geno 0.01`, `--hwe 0.0001`, `--mind 0.01`). Identifies and removes related or duplicated individuals using KING, and extracts the final cohorts.
* **`06_R_PCA_bim_harmonization.R`**: R script used to format `.bim` variant IDs for exact matching with the 1000 Genomes (1KGP) and HGDP reference databases.
* **`06_PCA.sh`**: Merges the study cohort with the 1KGP+HGDP reference, filters common variants (`--maf 0.05`), prunes for Linkage Disequilibrium (`--indep 50 5 2`), and computes the top Principal Components for ancestry inference and population stratification correction.

## Dependencies
* GATK 4.6.1+
* bcftools 1.19+ (required for advanced `csq` features)
* htslib (for `tabix` and `bgzip`)
* PLINK 1.9+
* KING
* GNU Utilities (`awk`, `sed`, `wget`, `tar`)
* R packages: `data.table`, `dplyr`, `tidyr`, `readr`, `magrittr`
