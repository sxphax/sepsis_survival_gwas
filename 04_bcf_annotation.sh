#!/bin/bash
# ==============================================================================
# Script: 04_bcf_annotation.sh
# Purpose: Convert to BCF, filter PASS variants, split by type, and annotate.
# Note: Ensure all database files (dbSNP, GnomAD, CADD, Ensembl) are pre-formatted.
# ==============================================================================

# Directories and Resources
IN_VCF="../03_recalibrated_vcfs/records_indels_snps.recal.vcf.gz"
OUT_DIR="."
ANN_DIR="annotation_files"
REF="/data/Ressources/bwa_reference/hs38DH.fa"

mkdir -p headers

# ------------------------------------------------------------------------------
# 1. Format Conversion and Base Filtering (Sequential)
# ------------------------------------------------------------------------------
echo "1. Converting VCF to BCF..."
# bcftools view -Ob: Converts the compressed VCF into binary BCF format. 
# BCF is significantly faster for downstream querying and requires less storage.
bcftools view "${IN_VCF}" -Ob -o "${OUT_DIR}/records_1065.bcf"
bcftools index "${OUT_DIR}/records_1065.bcf"

echo "2. Keeping only PASS variants..."
# -i 'FILTER="PASS"': Removes all variants that failed the VQSR truth thresholds in Step 3.
bcftools view -i 'FILTER="PASS"' "${OUT_DIR}/records_1065.bcf" -Ob -o "${OUT_DIR}/records_1065_pass.bcf"
bcftools index "${OUT_DIR}/records_1065_pass.bcf"

echo "3. Removing bloated INFO fields..."
# -x ^INFO/AC... : Strips out Allele Count, Allele Frequency, Allele Number, and Depth.
# We remove these because filtering out non-PASS variants makes the original cohort frequencies 
# mathematically incorrect. We will recalculate accurate frequencies later in PLINK.
bcftools annotate -x ^INFO/AC,INFO/AF,INFO/AN,INFO/DP "${OUT_DIR}/records_1065_pass.bcf" -Ob -o "${OUT_DIR}/records_1065_passx.bcf"
bcftools index "${OUT_DIR}/records_1065_passx.bcf"

# ------------------------------------------------------------------------------
# 2. Split into SNVs and INDELs
# ------------------------------------------------------------------------------
echo "4. Splitting into SNVs and non-SNVs..."
# We split the dataset because we will analyze exclusively SNVs

bcftools filter -i 'TYPE="snp"' "${OUT_DIR}/records_1065_passx.bcf" -Ob -o "${OUT_DIR}/records_1065_passx_snvs.bcf"
bcftools index "${OUT_DIR}/records_1065_passx_snvs.bcf"

bcftools filter -i 'TYPE!="snp"' "${OUT_DIR}/records_1065_passx.bcf" -Ob -o "${OUT_DIR}/records_1065_passx_notsnvs.bcf"
bcftools index "${OUT_DIR}/records_1065_passx_notsnvs.bcf"

# ------------------------------------------------------------------------------
# 3. Annotation Paths (Running SNV and INDEL branches in parallel)
# ------------------------------------------------------------------------------
echo "5. Annotating SNVs and INDELs (Parallel branches)..."

# Prepare headers for custom annotations
echo -e '##INFO=<ID=CADD,Number=.,Type=Float,Description="CADD v1.7 PHRED score">' > headers/cadd.hdr
echo -e '##INFO=<ID=gnomad41_genome_AF_grpmax,Number=.,Type=Float,Description="GnomADv4.1 PopMax AF">\n##INFO=<ID=gnomad41_genome_fafmax_faf95_max,Number=.,Type=Float,Description="GnomADv4.1 PopMax Filtered AF 95%CI">\n##INFO=<ID=gnomad41_genome_AF,Number=.,Type=Float,Description="GnomADv4.1 AF">\n##INFO=<ID=gnomad41_genome_faf95,Number=.,Type=Float,Description="GnomADv4.1 Filtered AF 95%CI">\n##INFO=<ID=gnomad41_genome_AF_nfe,Number=.,Type=Float,Description="GnomADv4.1 AF in NFE">' > headers/pop_af.hdr

# === PATH A: SNV ANNOTATION (Background Job) ===
(
    # A1. Add rsIDs from dbSNP
    bcftools annotate -a "${ANN_DIR}/dbSNPb156v2.vcf.gz" -c ID "${OUT_DIR}/records_1065_passx_snvs.bcf" -Ob -o "${OUT_DIR}/records_1065_snvs_rsid.bcf"
    bcftools index "${OUT_DIR}/records_1065_snvs_rsid.bcf"

    # A2. Add CADD Scores (Pathogenicity prediction)
    bcftools annotate -a "${ANN_DIR}/cadd1.7_snvs2.tsv.gz" -c CHROM,POS,-,REF,ALT,-,CADD -h headers/cadd.hdr -l CADD:append-missing "${OUT_DIR}/records_1065_snvs_rsid.bcf" -Ob -o "${OUT_DIR}/records_1065_snvs_rsid_cadd.bcf"
    bcftools index "${OUT_DIR}/records_1065_snvs_rsid_cadd.bcf"

    # A3. Add GnomAD v4.1 Allele Frequencies
    bcftools annotate -a "${ANN_DIR}/hg38_gnomad41_popmax_genome_nfe_sorted.tsv.gz" -h headers/pop_af.hdr -c CHROM,POS,-,REF,ALT,gnomad41_genome_AF_grpmax,gnomad41_genome_fafmax_faf95_max,gnomad41_genome_AF,gnomad41_genome_faf95,gnomad41_genome_AF_nfe -l gnomad41_genome_AF_grpmax:append-missing,gnomad41_genome_fafmax_faf95_max:append-missing,gnomad41_genome_AF:append-missing,gnomad41_genome_faf95:append-missing,gnomad41_genome_AF_nfe:append-missing "${OUT_DIR}/records_1065_snvs_rsid_cadd.bcf" -Ob -o "${OUT_DIR}/records_1065_snvs_rsid_cadd_AF.bcf"
    bcftools index "${OUT_DIR}/records_1065_snvs_rsid_cadd_AF.bcf"

    # A4. Add Consequences using Ensembl GFF3
    # bcftools csq predicts the functional consequence of the variant (e.g., missense, synonymous)
    bcftools csq -f "${REF}" -g "${ANN_DIR}/Homo_sapiens.GRCh38.111.gff3.gz" -Ob -o "${OUT_DIR}/records_1065_snvs_rsid_cadd_AF_csq.bcf" --verbose 2 --phase a "${OUT_DIR}/records_1065_snvs_rsid_cadd_AF.bcf"
    bcftools index "${OUT_DIR}/records_1065_snvs_rsid_cadd_AF_csq.bcf"
) &

# === PATH B: INDEL ANNOTATION (Background Job) ===
(
    # B1. Add rsIDs from dbSNP
    bcftools annotate -a "${ANN_DIR}/dbSNPb156v2.vcf.gz" -c ID "${OUT_DIR}/records_1065_passx_notsnvs.bcf" -Ob -o "${OUT_DIR}/records_1065_notsnvs_rsid.bcf"
    bcftools index "${OUT_DIR}/records_1065_notsnvs_rsid.bcf"

    # B2. Add GnomAD v4.1 Allele Frequencies
    bcftools annotate -a "${ANN_DIR}/hg38_gnomad41_popmax_genome_nfe_sorted.tsv.gz" -h headers/pop_af.hdr -c CHROM,POS,-,REF,ALT,gnomad41_genome_AF_grpmax,gnomad41_genome_fafmax_faf95_max,gnomad41_genome_AF,gnomad41_genome_faf95,gnomad41_genome_AF_nfe -l gnomad41_genome_AF_grpmax:append-missing,gnomad41_genome_fafmax_faf95_max:append-missing,gnomad41_genome_AF:append-missing,gnomad41_genome_faf95:append-missing,gnomad41_genome_AF_nfe:append-missing "${OUT_DIR}/records_1065_notsnvs_rsid.bcf" -Ob -o "${OUT_DIR}/records_1065_notsnvs_rsid_AF.bcf"
    bcftools index "${OUT_DIR}/records_1065_notsnvs_rsid_AF.bcf"
) &

# Wait for both SNV and INDEL branches to finish annotating
wait
echo "Module 04 Complete. Fully annotated BCFs are ready for PLINK conversion."
