#!/bin/bash
# ==============================================================================
# Script: 01_gvcf_prep.sh
# Purpose: Download, extract, verify, normalize, and index CNRGH g.vcf files.
# ==============================================================================

# ------------------------------------------------------------------------------
# 0. Define global variables to make the script easily adaptable
# ------------------------------------------------------------------------------
CNRGH_USER="RECORDS"
CNRGH_PASS="YOUR_PASSWORD"
# The BED file containing regions to exclude (centromeres, telomeres, dangerous regions)
BLACKLIST="/data/records_wgs/variants/03_recalibrated_vcfs/GRCh38_unified_blacklist.bed"
# Directory where the final normalized files will be stored
OUT_DIR="/data/records_wgs/variants/01_gVCFs/norm"

# Create the output directory
mkdir -p "${OUT_DIR}"

# ------------------------------------------------------------------------------
# 1. Download and Extract Archives
# ------------------------------------------------------------------------------
echo "Starting data download..."

# Loop through the specific delivery sets provided by CNRGH (e.g., sets 5 to 12)
for SET in 5 6 7 8 12; do
    echo "Downloading Set ${SET}..."
    # wget flags explained:
    # -m : mirror (download entire directory)
    # -c : continue getting a partially-downloaded file
    # --no-parent : don't ascend to the parent directory
    # -nH -nd : disable host directories and don't create a local directory hierarchy
    # -A .tar.gz : only accept files ending in .tar.gz
    wget -m -c --no-parent -nH -nd --cut-dirs=1 -A .tar.gz \
         --user="${CNRGH_USER}" --password="${CNRGH_PASS}" \
         "https://www.cnrgh.fr/data/RECORDS/variants/set${SET}/"
done

echo "Uncompressing archives..."
# gunzip uncompresses all .gz files in the current directory, turning them into .tar files
gunzip *.gz

# Loop through every .tar file to extract only the specific GATK HaplotypeCaller g.vcf files
for archive in *.tar; do 
    # -x : extract
    # -f : from file
    # --wildcards : allows the use of '*' to grab specific files inside the tarball
    tar -xf "${archive}" --wildcards "*gatk4HC.g.vcf.gz*"
done

# ------------------------------------------------------------------------------
# 2. Verify File Integrity (MD5 Checksums)
# ------------------------------------------------------------------------------
echo "Verifying MD5 checksums..."
# The original md5 files contained a path prefix (e.g., /something_hs38me/). 
# sed -i modifies the file in-place, removing everything up to '_hs38me/' so the paths match our current folder.
sed -i 's!/.*_hs38me/!!' *g.vcf.gz.md5

# md5sum -c reads the checksums and verifies the files. We redirect the output to a log file.
if md5sum -c *g.vcf.gz.md5 > records.gVCF.md5.log; then
    echo "SUCCESS: All MD5 checksums passed."
else
    echo "WARNING: Some MD5 checksums failed. Please check records.gVCF.md5.log"
fi

# Clean up the md5 files as they are no longer needed
rm *g.vcf.gz.md5

# ------------------------------------------------------------------------------
# 3. Normalize, Filter Blacklist, and Index (Parallel Processing)
# ------------------------------------------------------------------------------
echo "Starting normalization and indexing..."

# Initialize a counter to track how many background jobs we have launched
count=0

# Loop dynamically through every single g.vcf.gz file in the directory
for filename in *.g.vcf.gz; do
    
    # Extract the patient ID. 
    # cut -d"_" -f2 splits the filename by underscores and grabs the 2nd field (e.g., C001SOW)
    id=$(echo "${filename}" | cut -d"_" -f2)
    out_vcf="${OUT_DIR}/${id}.norm.g.vcf.gz"
    
    # Enclose the commands in parentheses (...) and add an ampersand '&' at the end.
    # This creates a "subshell" that runs in the background, allowing the main loop to continue.
    (
        # bcftools view -T ^${BLACKLIST} : Removes any variants falling inside the blacklist BED regions.
        # The output is piped '|' directly into bcftools norm.
        # bcftools norm -d none -m +any : Splits multiallelic sites into separate biallelic records.
        # -O z -o : Outputs a compressed vcf.gz to our target destination.
        bcftools view -T ^"${BLACKLIST}" "${filename}" | \
        bcftools norm -d none -m +any -O z -o "${out_vcf}"
        
        # Once normalization for this patient is done, immediately index the new file.
        # -f : force overwrite if an index exists
        # -p vcf : specifies the format
        tabix -f -p vcf "${out_vcf}"
    ) &
    
    # Increment the counter
    ((count++))
    
    # If the counter is evenly divisible by 8 (our chosen batch size):
    if (( count % 8 == 0 )); then
        # 'wait' pauses the script here until the 8 background jobs we just launched have completely finished.
        # This prevents overloading the server RAM/CPU.
        wait 
    fi
done

# Wait one final time to ensure any remaining jobs (if the total wasn't perfectly divisible by 8) finish.
wait 

echo "Module 01 Complete. Normalized files are in ${OUT_DIR}"
