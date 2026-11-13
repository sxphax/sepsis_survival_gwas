# ==============================================================================
# .bim ID Harmonization (Run BEFORE 06_PCA.sh)
# ==============================================================================
cat("Harmonizing .bim IDs for 1KG merge...\n")

bimKG <- fread("/data/records_wgs/syphax_analyses/BCF/records_1065/PCA/1KG_HGDP/hgdp_tgp.bim")
bim <- fread("/data/records_wgs/syphax_analyses/BCF/records_1065/PCA/records_1030_5%.bim")

# Standardize IDs to chr_pos_ref_alt
bim$V2 <- paste(bim$V1, bim$V4, bim$V5, bim$V6, sep="_")
bimKG$V2 <- paste(bimKG$V1, bimKG$V4, bimKG$V5, bimKG$V6, sep="_")

write_tsv(bim, "/data/records_wgs/syphax_analyses/BCF/records_1065/PCA/records_1030_5%.bim", col_names=FALSE)
write_tsv(bimKG, "/data/records_wgs/syphax_analyses/BCF/records_1065/PCA/1KG_HGDP/hgdp_tgp.bim", col_names=FALSE)
write_tsv(as.data.frame(intersect(bimKG$V2, bim$V2)), 
          "/data/records_wgs/syphax_analyses/BCF/records_1065/PCA/records1030_1KG_HGDP/snptoinclude", col_names=FALSE)
