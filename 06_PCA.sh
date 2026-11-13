#!/bin/bash
# ==============================================================================
# Script: 06_PCA.sh
# Purpose: Merge the 1030 cohort with 1KGP+HGDP, prune for LD, and compute PCA.
# Note: You MUST run the R script to harmonize .bim IDs and generate 
#       'snptoinclude' before executing this bash script.
# ==============================================================================

# Directories and Inputs
COHORT_PREFIX="../records_1030_5%"
REF_PREFIX="/data/records_wgs/syphax_analyses/BCF/records_1065/PCA/1KG_HGDP/hgdp_tgp"
INTERSECT_LIST="snptoinclude"

# ------------------------------------------------------------------------------
# 1. Extract Overlapping Variants
# ------------------------------------------------------------------------------
echo "1. Extracting overlapping SNPs from both datasets..."

# Extract intersection from the 1030 WGS cohort
plink --bfile "${COHORT_PREFIX}" \
      --extract "${INTERSECT_LIST}" \
      --make-bed \
      --out wgs_for_1kg_merge

# Extract intersection from the 1KGP+HGDP reference database
plink --bfile "${REF_PREFIX}" \
      --extract "${INTERSECT_LIST}" \
      --make-bed \
      --out 1kg_for_wgs_merge

# ------------------------------------------------------------------------------
# 2. Merge Datasets and Apply Quality Filters
# ------------------------------------------------------------------------------
echo "2. Merging cohort with reference genome database..."

# Merge the two datasets. Because the .bim IDs were harmonized in R, 
# this step will cleanly align the variants without ID collision errors.
plink --bfile wgs_for_1kg_merge \
      --bmerge 1kg_for_wgs_merge.bed 1kg_for_wgs_merge.bim 1kg_for_wgs_merge.fam \
      --make-bed \
      --out data_merged

echo "3. Filtering merged database..."

# --maf 0.05: Keep only common variants (frequency > 5%). Rare variants 
#             do not provide stable information for broad population clustering.
# --geno 0.01: Remove variants missing in >1% of the combined population.
plink --bfile data_merged \
      --maf 0.05 \
      --geno 0.01 \
      --make-bed \
      --out data_merged2

# ------------------------------------------------------------------------------
# 4. Linkage Disequilibrium (LD) Pruning
# ------------------------------------------------------------------------------
echo "4. Pruning variants for Linkage Disequilibrium..."

# --indep 50 5 2: Looks at a window of 50 SNPs, shifts the window by 5 SNPs, 
# and calculates the correlation (r^2) between them. It removes one of a pair 
# of SNPs if their r^2 is greater than 0.2 (represented as Variance Inflation Factor 2).
# This prevents highly correlated genomic regions from artificially skewing the PCA.
# Outputs: ldprune.prune.in (variants to keep) and ldprune.prune.out (variants removed).
plink --bfile data_merged2 \
      --indep 50 5 2 \
      --out ldprune

# Extract only the independent SNPs into a final dataset for PCA
plink --bfile data_merged2 \
      --extract ldprune.prune.in \
      --make-bed \
      --out data4pca

# ------------------------------------------------------------------------------
# 5. Principal Component Analysis (PCA)
# ------------------------------------------------------------------------------
echo "5. Computing Principal Components..."

# Calculates the top 20 principal components and eigenvalues.
# Outputs: plink.eigenvec (patient coordinates) and plink.eigenval (variance explained).
plink --bfile data4pca \
      --pca \
      --out plink

echo "Module 06 Complete. PCA coordinates generated in 'plink.eigenvec'."
