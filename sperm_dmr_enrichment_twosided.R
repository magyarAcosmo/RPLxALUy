##############################################################################
# Sperm-specific DMR enrichment analysis — Two-Sided Version
##############################################################################

library(LOLA)
library(GenomicRanges)
library(rtracklayer)
library(dplyr)

# ── 1. User Inputs ──────────────────────────────────────────────────────────
dmr_file  <- "/Users/jordanmoore/Desktop/THEWAY/RPL/dmrs_2026-01-17_woAGE.bed"
univ_file <- "/Users/jordanmoore/Desktop/THEWAY/RPL/HMRs_2026-01-17.bed"
lola_parent_db <- "/Users/jordanmoore/Desktop/THEWAY/RPL/lola_db"
GENOME_OUT <- "hg38"

# ── 2. Load Data ─────────────────────────────────────────────────────────────
message("\n── Loading DMRs and Universe ──")
dmr_gr  <- rtracklayer::import(dmr_file,  format = "BED")
univ_gr <- rtracklayer::import(univ_file, format = "BED")
genome(dmr_gr)  <- GENOME_OUT
genome(univ_gr) <- GENOME_OUT

# Load LOLA DB
lola_db <- loadRegionDB(lola_parent_db, collections = "sperm_features")

# FIX: Correctly extract the regions and annotations from the LOLA object
region_list <- lola_db$regionGRL
anno_table  <- as.data.frame(lola_db$regionAnno)

# ── 3. Manual Two-Sided Fisher Test for Every Feature ───────────────────────
message("\n── Running Two-Sided Fisher Tests ──")

results_list <- list()

for (i in seq_along(region_list)) {
  # Get name from annotation table safely
  feature_name <- anno_table$antibody[i]
  if (is.null(feature_name) || is.na(feature_name) || feature_name == "") {
    feature_name <- anno_table$filename[i]
  }
  
  feat_gr <- region_list[[i]]
  
  # Contingency Table
  a <- length(subsetByOverlaps(dmr_gr, feat_gr))
  b <- length(subsetByOverlaps(univ_gr, feat_gr)) - a
  c <- length(dmr_gr) - a
  d <- length(univ_gr) - (a + b + c)
  
  fisher_res <- fisher.test(matrix(c(a, b, c, d), nrow = 2, byrow = TRUE), 
                            alternative = "two.sided")
  
  results_list[[i]] <- data.frame(
    antibody  = feature_name,
    oddsRatio = as.numeric(fisher_res$estimate),
    p_value   = as.numeric(fisher_res$p.value),
    support   = a
  )
}

# ── 4. Special Case: Bivalency (K4 + K27) ───────────────────────────────────
message("Calculating Bivalency...")

# Find the indices using case-insensitive matching
idx_k4  <- grep("H3K4me3",  anno_table$antibody, ignore.case = TRUE)[1]
idx_k27 <- grep("H3K27me3", anno_table$antibody, ignore.case = TRUE)[1]

if (!is.na(idx_k4) & !is.na(idx_k27)) {
  k4_gr  <- region_list[[idx_k4]]
  k27_gr <- region_list[[idx_k27]]
  
  # Define bivalent regions
  bivalent_gr <- subsetByOverlaps(k4_gr, k27_gr)
  
  # Fisher math for bivalency
  a_biv <- length(subsetByOverlaps(dmr_gr, bivalent_gr))
  univ_biv_gr <- subsetByOverlaps(subsetByOverlaps(univ_gr, k4_gr), k27_gr)
  b_biv <- length(univ_biv_gr) - a_biv
  c_biv <- length(dmr_gr) - a_biv
  d_biv <- length(univ_gr) - (a_biv + b_biv + c_biv)
  
  fisher_biv <- fisher.test(matrix(c(a_biv, b_biv, c_biv, d_biv), nrow = 2, byrow = TRUE), 
                            alternative = "two.sided")
  
  results_list[[length(results_list) + 1]] <- data.frame(
    antibody  = "Bivalent (K4+K27)",
    oddsRatio = as.numeric(fisher_biv$estimate),
    p_value   = as.numeric(fisher_biv$p.value),
    support   = a_biv
  )
} else {
  message("Warning: Could not find H3K4me3 or H3K27me3 for bivalency check.")
}

# ── 5. Combine and Adjust ───────────────────────────────────────────────────
final_table <- bind_rows(results_list) %>%
  mutate(q_value = p.adjust(p_value, method = "BH")) %>%
  mutate(status = case_when(
    p_value < 0.05 & oddsRatio > 1 ~ "Enriched",
    p_value < 0.05 & oddsRatio < 1 ~ "Depleted",
    TRUE ~ "NS"
  )) %>%
  arrange(p_value)

# ── 6. Permutation Testing (Two-Sided & Integrated) ─────────────────────────
library(regioneR)
NPERM <- 10000

message("\n── Running Permutations (Skipping SNPs) ──")

# Create empty columns to hold the permutation data
final_table$perm_zscore <- NA
final_table$perm_p_value <- NA

for (i in 1:nrow(final_table)) {
  feature_name <- final_table$antibody[i]
  
  # Skip SNPs (Permuting 15 million regions takes days)
  if (grepl("SNP", feature_name, ignore.case = TRUE)) {
    message("  Skipping permutations for ", feature_name, " (too large)")
    next
  }
  
  message("  Permuting ", feature_name, "...")
  
  # Figure out which GRanges object to use
  if (feature_name == "Bivalent (K4+K27)") {
    feat_gr <- bivalent_gr
  } else {
    idx <- which(anno_table$antibody == feature_name | anno_table$filename == feature_name)[1]
    feat_gr <- region_list[[idx]]
  }
  
  # Run the permutation
  # FIX: randomize.function changed to resampleRegions for custom universes!
  pt <- permTest(
    A = dmr_gr, 
    B = feat_gr, 
    universe = univ_gr,
    randomize.function = resampleRegions, 
    evaluate.function = numOverlaps, 
    ntimes = NPERM,
    alternative = "auto" 
  )
  
  # Extract the stats and add them to the table
  final_table$perm_zscore[i]  <- pt[[1]]$zscore
  final_table$perm_p_value[i] <- pt[[1]]$pval
  
  # Save the visual plot
  safe_name <- gsub("[^A-Za-z0-9_]", "_", feature_name)
  pdf(paste0("perm_dist_", safe_name, ".pdf"), width = 6, height = 4)
  plot(pt[[1]], main = paste("Permutation:", feature_name))
  dev.off()
}

# ── 7. Final Output ─────────────────────────────────────────────────────────
message("\n── Final Results with Permutations ──")
print(final_table)
write.csv(final_table, "sperm_enrichment_FINAL_with_perms_noAGE.csv", row.names = FALSE)
message("Success! Results saved to sperm_enrichment_FINAL_with_perms_noAGE.csv")