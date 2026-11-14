#!/bin/bash
# ==============================================================================
# Script: 03_vqsr_filtering.sh
# Purpose: ExcessHet hard-filtering, VQSR calculation, and metric generation.
# ==============================================================================

# ------------------------------------------------------------------------------
# 0. Variables & Resources
# ------------------------------------------------------------------------------
# Directory containing the Broad Institute's gold-standard truth sets for GRCh38
RES="/data/Ressources/broad_bundle_hg38"
# The raw, merged cohort VCF from Module 02
INPUT="../03_recalibrated_vcfs/records_allchr_raw.vcf.gz"

# ------------------------------------------------------------------------------
# 1. Hard Filter Excess Heterozygosity & Create Sites-Only VCF
# ------------------------------------------------------------------------------
echo "1. Filtering ExcessHet and creating Sites-Only VCF..."

# VariantFiltration: Tags variants that drastically violate Hardy-Weinberg equilibrium 
# (ExcessHet > 54.69). These are usually mapping artifacts (paralogous regions).
gatk --java-options "-Xmx100g -Xms100g" VariantFiltration \
    -V ${INPUT} \
    --filter-expression "ExcessHet > 54.69" \
    --filter-name ExcessHet \
    -O records_allchr_excesshet.vcf.gz \
    2> Logs/ExcessHet.log

# MakeSitesOnlyVcf: Strips out all individual patient genotypes (the 1065 sample columns).
# It only keeps the variant sites (chr, pos, ref, alt). This drastically reduces the 
# file size and speeds up the complex math in the next VQSR step.
gatk MakeSitesOnlyVcf \
    -I records_allchr_excesshet.vcf.gz \
    -O records_allchr_sitesonly.vcf.gz \
    2> Logs/MakeSitesOnlyVcf.log

# ------------------------------------------------------------------------------
# 2. Build the VQSR Models (Running in parallel)
# ------------------------------------------------------------------------------
echo "2. Calculating VQSLOD tranches for INDELs and SNPs..."

# VariantRecalibrator builds a machine learning model (Gaussian mixture model) 
# that separates true variants from sequencing artifacts based on annotations 
# like Depth (DP), Mapping Quality (MQ), etc.

# A. Build the model for INDELs (--max-gaussians 4 is standard for INDELs)
gatk --java-options "-Xmx50g -Xms50g" VariantRecalibrator \
    -V records_allchr_sitesonly.vcf.gz \
    --trust-all-polymorphic -mode INDEL --max-gaussians 4 \
    -tranche 100.0 -tranche 99.95 -tranche 99.9 -tranche 99.5 -tranche 99.0 -tranche 95.0 -tranche 90.0 \
    -an FS -an ReadPosRankSum -an MQRankSum -an QD -an SOR -an DP \
    -resource:mills,known=false,training=true,truth=true,prior=12 ${RES}/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz \
    -resource:axiomPoly,known=false,training=true,truth=false,prior=10 ${RES}/Axiom_Exome_Plus.genotypes.all_populations.poly.hg38.vcf.gz \
    -resource:dbsnp,known=true,training=false,truth=false,prior=2 ${RES}/dbsnp_146.hg38.vcf.gz \
    -O records_indels.recal \
    --tranches-file records_indels.tranches \
    2> Logs/VariantRecalibrator_indel.log &

# B. Build the model for SNPs (--max-gaussians 6 is standard for SNPs)
gatk --java-options "-Xmx50g -Xms50g" VariantRecalibrator \
    -V records_allchr_sitesonly.vcf.gz \
    --trust-all-polymorphic -mode SNP --max-gaussians 6 \
    -tranche 100.0 -tranche 99.95 -tranche 99.9 -tranche 99.5 -tranche 99.0 -tranche 98.0 \
    -an QD -an MQRankSum -an ReadPosRankSum -an FS -an MQ -an SOR -an DP \
    -resource:hapmap,known=false,training=true,truth=true,prior=15 ${RES}/hapmap_3.3.hg38.vcf.gz \
    -resource:omni,known=false,training=true,truth=true,prior=12 ${RES}/1000G_omni2.5.hg38.vcf.gz \
    -resource:1000G,known=false,training=true,truth=false,prior=10 ${RES}/1000G_phase1.snps.high_confidence.hg38.vcf.gz \
    -resource:dbsnp,known=true,training=false,truth=false,prior=7 ${RES}/dbsnp_146.hg38.vcf.gz \
    -O records_snps.recal \
    --tranches-file records_snps.tranches \
    2> Logs/VariantRecalibrator_snp.log &

# Wait barrier: We absolutely MUST wait for both models to finish building 
# before we try to apply them to our data in the next step.
wait
echo "Tranches calculated successfully."

# ------------------------------------------------------------------------------
# 3. Apply the VQSR Models (Sequential execution)
# ------------------------------------------------------------------------------
echo "3. Applying VQSR sequentially..."

# ApplyVQSR evaluates each variant against the model and filters out those falling 
# below the 99.0% truth sensitivity threshold.

# First, apply the INDEL model to the original VCF (the one with the patient genotypes)
gatk --java-options "-Xmx50g -Xms50g" ApplyVQSR \
    -V records_allchr_excesshet.vcf.gz \
    --recal-file records_indels.recal --tranches-file records_indels.tranches \
    --truth-sensitivity-filter-level 99.0 --create-output-variant-index true -mode INDEL \
    -O records_indels.recal.vcf.gz \
    2> Logs/ApplyVQSR_indel.log

# as its input. It must run sequentially.
gatk --java-options "-Xmx50g -Xms50g" ApplyVQSR \
    -V records_indels.recal.vcf.gz \
    --recal-file records_snps.recal --tranches-file records_snps.tranches \
    --truth-sensitivity-filter-level 99.0 --create-output-variant-index true -mode SNP \
    -O records_indels_snps.recal.vcf.gz \
    2> Logs/ApplyVQSR_snp.log

# ------------------------------------------------------------------------------
# 4. Generate Quality Control Metrics
# ------------------------------------------------------------------------------
echo "4. Generating metrics..."

# CollectVariantCallingMetrics generates summary statistics (Ti/Tv ratio, dbSNP concordance)
gatk CollectVariantCallingMetrics \
    -I records_indels_snps.recal.vcf.gz \
    --DBSNP ${RES}/dbsnp_146.hg38.vcf.gz \
    -SD ${RES}/Homo_sapiens_assembly38.dict \
    -O records.variants.metrics \
    2> Logs/CollectVariantCallingMetrics.log &

# VariantEval produces highly detailed stratifications of the variant calls
gatk VariantEval \
    -R ${RES}/Homo_sapiens_assembly38.fasta \
    -eval records_indels_snps.recal.vcf.gz \
    -D ${RES}/dbsnp_146.hg38.vcf.gz \
    -O records_Eval.txt \
    2> Logs/VariantEval.log &

# Wait for both metric calculations to finish
wait
echo "Module 03 Complete."
