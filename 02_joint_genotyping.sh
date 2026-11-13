#!/bin/bash
# ==============================================================================
# Script: 02_joint_genotyping.sh
# Purpose: Import cohort to GenomicsDB and GenotypeGVCFs using precise 
#          sub-chromosome windows, executed in parallel batches of 8.
# ==============================================================================

# ------------------------------------------------------------------------------
# 0. Define Variables
# ------------------------------------------------------------------------------
REF="/home/NGS/Ressources/bwa_reference/hs38DH.fa"
DB_DIR="../recordsDB"
MAP_FILE="${DB_DIR}/sample_map.txt"
OUT_RAW="raw_vcfs"

mkdir -p "${DB_DIR}" "${OUT_RAW}" TMP

# ------------------------------------------------------------------------------
# 1. Build the Sample Map
# ------------------------------------------------------------------------------
echo "Building sample map..."
rm -f "${MAP_FILE}"
for fn in ../01_gVCFs/norm/*.norm.g.vcf.gz; do
    id=$(echo "${fn}" | cut -d/ -f4 | cut -d. -f1)
    echo -e "${id}\t${fn}" >> "${MAP_FILE}"
done

# ------------------------------------------------------------------------------
# 2. Define Precise Genomic Windows (the 50Mb genomic windows)
# ------------------------------------------------------------------------------
# Each entry is formatted as "ID:INTERVAL". This allows us to keep the output
# files cleanly named (e.g., records_db_1_1) while processing exact coordinates.
CHUNKS=(
    "22:chr22"
    "1_1:chr1:1-50000000"
    "1_2:chr1:50000001-100000000"
    "1_3:chr1:100000001-150000000"
    "1_4:chr1:150000001-200000000"
    "1_5:chr1:200000001-248956422"
    "21:chr21"
    "20:chr20"
    "19:chr19"
    "18_1:chr18:1-50000000"
    "18_2:chr18:50000001-80373285"
    "17:chr17"
    "16:chr16"
    "15:chr15"
    "14_1:chr14:1-50000000"
    "14_2:chr14:50000001-107043718"
    "13_1:chr13:1-50000000"
    "13_2:chr13:50000001-114364328"
    "12_1:chr12:1-50000000"
    "12_2:chr12:50000001-100000000"
    "12_3:chr12:100000001-133275309"
    "11_1:chr11:1-50000000"
    "11_2:chr11:50000001-100000000"
    "11_3:chr11:100000001-135086622"
    "10_1:chr10:1-50000000"
    "10_2:chr10:50000001-100000000"
    "10_3:chr10:100000001-133797422"
    "9_1:chr9:1-50000000"
    "9_2:chr9:50000001-100000000"
    "9_3:chr9:100000001-138394717"
    "8_1:chr8:1-50000000"
    "8_2:chr8:50000001-100000000"
    "8_3:chr8:100000001-145138636"
    "7_1:chr7:1-50000000"
    "7_2:chr7:50000001-100000000"
    "7_3:chr7:100000001-159345973"
    "6_1:chr6:1-50000000"
    "6_2:chr6:50000001-100000000"
    "6_3:chr6:100000001-150000000"
    "6_4:chr6:150000001-170805979"
    "5_1:chr5:1-50000000"
    "5_2:chr5:50000001-100000000"
    "5_3:chr5:100000001-150000000"
    "5_4:chr5:150000001-181538259"
    "4_1:chr4:1-50000000"
    "4_2:chr4:50000001-100000000"
    "4_3:chr4:100000001-150000000"
    "4_4:chr4:150000001-190214555"
    "3_1:chr3:1-50000000"
    "3_2:chr3:50000001-100000000"
    "3_3:chr3:100000001-150000000"
    "3_4:chr3:150000001-198295559"
    "2_1:chr2:1-50000000"
    "2_2:chr2:50000001-100000000"
    "2_3:chr2:100000001-150000000"
    "2_4:chr2:150000001-200000000"
    "2_5:chr2:200000001-242193529"
    "X_1:chrX:1-50000000"
    "X_2:chrX:50000001-100000000"
    "X_3:chrX:100000001-156040895"
    "Y:chrY"
)

# ------------------------------------------------------------------------------
# 3. Import to DB in batches of 8
# ------------------------------------------------------------------------------
echo "Running GenomicsDBImport in batches of 8..."
count=0
for CHUNK in "${CHUNKS[@]}"; do
    # Extract the name (e.g., 1_1) and the interval (e.g., chr1:1-50000000)
    NAME="${CHUNK%%:*}"
    INTERVAL="${CHUNK#*:}"
    
    gatk --java-options "-Xmx4g -Xms4g" GenomicsDBImport \
        --genomicsdb-update-workspace-path "${DB_DIR}/records_db_${NAME}" \
        --batch-size 32 \
        -L "${INTERVAL}" \
        --sample-name-map "${MAP_FILE}" \
        --reader-threads 2 \
        --tmp-dir TMP \
        2>> "GenomicsDBImport_${NAME}.log" &
        
    ((count++))
    if (( count % 8 == 0 )); then
        wait
        echo "Finished 8 imports. Starting next batch..."
    fi
done
wait # Ensure the final batch finishes

# ------------------------------------------------------------------------------
# 4. Genotype in batches of 8
# ------------------------------------------------------------------------------
echo "Running GenotypeGVCFs in batches of 8..."
count=0
for CHUNK in "${CHUNKS[@]}"; do
    NAME="${CHUNK%%:*}"
    
    gatk --java-options "-Xmx4g -Xms4g -XX:ParallelGCThreads=2" GenotypeGVCFs \
        -R "${REF}" \
        -V "gendb://${DB_DIR}/records_db_${NAME}" \
        -O "${OUT_RAW}/records_${NAME}_raw.vcf.gz" \
        --tmp-dir TMP \
        2> "GenotypeGVCFs_${NAME}.log" &
        
    ((count++))
    if (( count % 8 == 0 )); then
        wait
        echo "Finished 8 genotypes. Starting next batch..."
    fi
done
wait

# ------------------------------------------------------------------------------
# 5. Merge all windows into one final VCF
# ------------------------------------------------------------------------------
echo "Merging chromosomes into final raw VCF..."
ls ${OUT_RAW}/*_raw.vcf.gz > raw_vcf.list

gatk --java-options "-Xmx100g -Xms100g" MergeVcfs \
    -I raw_vcf.list \
    -O ../03_recalibrated_vcfs/records_allchr_raw.vcf.gz \
    2> ../03_recalibrated_vcfs/Logs/MergeVcfs.log

echo "Module 02 Complete."
