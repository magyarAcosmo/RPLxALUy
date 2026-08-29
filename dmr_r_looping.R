library(tidyverse)
library(rtracklayer)
library(GenomicRanges)
library(regioneR)
library(BSgenome.Hsapiens.UCSC.hg38)
# r looping reference from RLBase
library(RLHub)

# read my DMR data
dmr_df <- read.table("dmrs_2026-01-17_pval.bed", header = FALSE, sep ="\t")

# add column names to dmr_df
colnames(dmr_df) <- c("chr", "start", "end", "id", "pvalue", "strand")

# load my DMR data in GRanges
dmrs <- GRanges(
  seqnames = dmr_df$chr,
  ranges = IRanges(start = dmr_df$start, end = dmr_df$end),
  strand = dmr_df$strand,
  pval = dmr_df$pvalue
)

# meta has location column -- parse for start/end to build GRanges
rl_meta <- rlregions_meta()

# parse location "chr1:633760-634280:." string into separate columns
rl_meta_parsed <- rl_meta %>%
  separate(location, into = c("chr", "coords", "strand"), sep = ":") %>%
  separate(coords, into = c("start", "end"), sep = "-") %>%
  mutate(
    start = as.integer(start),
    end = as.integer(end),
    strand = ifelse(strand == ".", "*", strand) # replace . with * for GRanges format
  )

# convert to GRanges
rl_regions_gr <- GRanges(
  seqnames = rl_meta_parsed$chr,
  ranges = IRanges(start = rl_meta_parsed$start, end = rl_meta_parsed$end),
  strand = rl_meta_parsed$strand
)

# quick sanity check
rl_regions_gr
length(rl_regions_gr)

# overlap dmrs with R loop regions
hits <- findOverlaps(dmrs, rl_regions_gr)

# which dmrs overlap an R loop?
dmr_with_rloop <- dmrs[queryHits(hits)]
n_overlap <- length(unique(queryHits(hits)))

cat("DMRs overlapping R loops:", n_overlap, "of", length(dmrs), "\n")
cat("Percentage:", round(n_overlap / length(dmrs) * 100, 1), "%\n")

# permutation test for enrichment
perm_t <- permTest(
  A = dmrs,
  B = rl_regions_gr,
  randomize.function = randomizeRegions,
  evaluate.function = numOverlaps,
  genome = "hg38",
  ntimes = 1000,
  verbose = TRUE
)
capture.output(print(perm_t), file = "r_loop_perm_results.txt")
summary(perm_t)
plot(perm_t)

##############################
# Check: Ignore strand for overlap — often appropriate for DMR/R-loop analysis
hits_any <- findOverlaps(dmrs, rl_regions_gr, ignore.strand = TRUE)
n_overlap_any <- length(unique(queryHits(hits_any)))

cat("DMRs overlapping R loops (any strand):", n_overlap_any, "of", length(dmrs), "\n")

# Also rerun permutation test ignoring strand -- ~same pvalue & insignificant
perm_t2 <- permTest(
  A                  = dmrs,
  B                  = rl_regions_gr,
  randomize.function = randomizeRegions,
  evaluate.function  = numOverlaps,
  genome             = "hg38",
  ntimes             = 1000,
  verbose            = TRUE,
  ignore.strand      = TRUE
)

summary(perm_t2)
plot(perm_t2)
