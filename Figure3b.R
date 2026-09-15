# ============================================================================
# Figure 3d: DMR coverage across the AluY consensus sequence
# (colored by number of overlapping AluY elements per DMR)
# ============================================================================

# ----------------------------------------------------------------------------
# PRE-FIGURE PREP
# Produces: dmrs_2026-01-17_woAGE_effect.bed, AluY_elements_annot.txt
# ----------------------------------------------------------------------------

# DMR EFFECT
## match up noAGE DMRs with AGE DMRs in dmrs_2026-01-17_effect.bed
library(GenomicRanges)
library(rtracklayer)

# read both files
woAGE <- read.delim("DMRs_woAGE/dmrs_2026-01-17_woAGE.bed", header = FALSE,
                    col.names = c("chr", "start", "end"))

effect <- read.delim("DMRs_wAGE/dmrs_2026-01-17_effect.bed", header = FALSE,
                     col.names = c("chr", "start", "end", "x_counter", "effect_value", "strand"))

# convert to GRanges
gr_woAGE <- GRanges(seqnames = woAGE$chr,
                    ranges = IRanges(start = woAGE$start + 1, end = woAGE$end))
gr_effect <- GRanges(seqnames = effect$chr,
                     ranges = IRanges(start = effect$start + 1, end = effect$end),
                     effect_value = effect$effect_value)

# Find overlaps: for each woAGE DMR, find matching effect DMR
hits <- findOverlaps(gr_woAGE, gr_effect, type = "equal")   # exact coordinate match

# Build Output
df_out <- woAGE[queryHits(hits), ]
df_out$effect_value <- effect$effect_value[subjectHits(hits)]

# check if every woAGE DMR got a match
cat("woAGE DMRs:", nrow(woAGE), "\n")
cat("Matched DMRs:", nrow(df_out), "\n")

# write output file
write.table(df_out,
            file = "dmrs_2026-01-17_woAGE_effect.bed",
            sep = "\t",
            quote = FALSE,
            row.names = FALSE,
            col.names = FALSE)
######
#  Results: all 294 woAGE DMRs matched!
######

# AluY coordinates
## rmsk.txt.gz contains all AluY coordinates
## include all subfamilies that start with AluY

# read rmsk.txt file
rmsk <- read.delim("rmsk.txt", header = FALSE)

# filter: repName (col 11) starts with "AluY"
aluY <- rmsk[grepl("^AluY", rmsk[, 11]), ]

# check what subfamilies were included
table(aluY[, 11])   # 28 subfamilies

cat("Total AluY elements found:", nrow(aluY), "\n")   #151627

write.table(aluY,
            file = "rmsk_AluY.txt",
            sep = "\t",
            quote = FALSE,
            row.names = FALSE,
            col.names = FALSE)

# clean AluY elements to relevant columns
aluY_clean <- data.frame(
  chr = aluY[, 6],        # genoName
  start = aluY[, 7],      # genoStart (0-based, BED-style)
  end = aluY[, 8],        # genoEnd
  strand = aluY[, 10],    # strand
  repName = aluY[, 11],   # aluY subfamily
  repStart = aluY[, 14],  # start position within AluY consensus
  repEnd = aluY[, 15],    # end position wihtin AluY consensus
  repLeft = aluY[, 16]    # bases remaining to end of consensus
)

head(aluY_clean, 5)
head(aluY_clean$repLeft, 5)

# Normalize consensus coordinates by strand
## UCSC reference: https://genome.ucsc.edu/cgi-bin/hgTables?db=hg38&hgta_group=rep&hgta_track=rmsk&hgta_table=rmsk&hgta_doSchema=describe+table+schema
aluY_clean$consensus_start <- ifelse(aluY_clean$strand == "+", aluY_clean$repStart, aluY_clean$repLeft)
aluY_clean$consensus_end <- aluY_clean$repEnd   # always the middle/true end coordinate

# check: consensus_start always <= consensus_end
sum(aluY_clean$consensus_start > aluY_clean$consensus_end)

head(aluY_clean[, c("chr", "start", "end", "strand", "repStart", "repEnd", "repLeft",
                    "consensus_start", "consensus_end")], 5)

write.table(aluY_clean,
            file = "AluY_elements_annot.txt",
            sep = "\t",
            quote = FALSE,
            row.names = FALSE,

# ----------------------------------------------------------------------------
# FIGURE 3d CODE
# ----------------------------------------------------------------------------

library(ggplot2)
library(GenomicRanges)

# Read DMR file
dmrs <- read.delim("dmrs_2026-01-17_woAGE_effect.bed", header = FALSE,
                    col.names = c("chr","start","end","effect"))
gr_dmrs <- GRanges(seqnames = dmrs$chr, ranges = IRanges(start = dmrs$start + 1, end = dmrs$end),
                    effect = dmrs$effect)

# Read AluY file
aluY_clean <- read.delim("AluY_elements_annot.txt", header = TRUE)

# rebuild consensus_start/end, including all 151627 elements
aluY_clean$consensus_start <- ifelse(aluY_clean$strand == "+", aluY_clean$repStart, aluY_clean$repLeft)
aluY_clean$consensus_end   <- aluY_clean$repEnd

# swap instead of drop, for the 2 previously-excluded fragments
swap_idx <- aluY_clean$consensus_start > aluY_clean$consensus_end
sum(swap_idx)  # 2

tmp <- aluY_clean$consensus_start[swap_idx]
aluY_clean$consensus_start[swap_idx] <- aluY_clean$consensus_end[swap_idx]
aluY_clean$consensus_end[swap_idx]   <- tmp

nrow(aluY_clean)                                            # 151627
sum(aluY_clean$consensus_start > aluY_clean$consensus_end)  # 0

# build AluY GRanges
gr_alu <- GRanges(seqnames = aluY_clean$chr,
                   ranges = IRanges(start = aluY_clean$start + 1, end = aluY_clean$end),
                   strand = aluY_clean$strand,
                   repName = aluY_clean$repName,
                   consensus_start = aluY_clean$consensus_start,
                   consensus_end = aluY_clean$consensus_end)

# find overlaps
hits <- findOverlaps(gr_dmrs, gr_alu, type = "any")
overlap_ranges <- pintersect(gr_dmrs[queryHits(hits)], gr_alu[subjectHits(hits)])

alu_hit <- gr_alu[subjectHits(hits)]
g_start <- start(alu_hit); g_end <- end(alu_hit)
cons_start <- alu_hit$consensus_start; cons_end <- alu_hit$consensus_end
over_start <- start(overlap_ranges); over_end <- end(overlap_ranges)
is_minus <- as.character(strand(alu_hit)) == "-"

frac_start <- ifelse(is_minus, (g_end - over_end) / (g_end - g_start),
                     (over_start - g_start) / (g_end - g_start))
frac_end <- ifelse(is_minus, (g_end - over_start) / (g_end - g_start),
                   (over_end - g_start) / (g_end - g_start))

df_map <- data.frame(
  start     = cons_start + frac_start * (cons_end - cons_start),
  end       = cons_start + frac_end   * (cons_end - cons_start),
  dmr_chr   = as.character(seqnames(gr_dmrs))[queryHits(hits)],
  dmr_start = start(gr_dmrs)[queryHits(hits)],
  dmr_end   = end(gr_dmrs)[queryHits(hits)]
)

cat("Total overlap rows:", nrow(df_map), "\n") # 176

# count AluY overlaps per DMR
df_map$dmr_id <- paste(df_map$dmr_chr, df_map$dmr_start, df_map$dmr_end, sep = "_")
dmr_counts <- table(df_map$dmr_id)
df_map$overlap_count <- as.integer(dmr_counts[df_map$dmr_id])

table(dmr_counts)  # shows 1 / 2 / 3 breakdown (18 overlapping)

df_map$overlap_count <- factor(df_map$overlap_count, levels = c(1, 2, 3))

# keep only what's needed for plotting
df_map <- df_map[, c("start", "end", "overlap_count")]

# -----------------------------
# INPUTS
# -----------------------------
alu_len <- 300
cpg_pos <- c(4,8,10,20,48,53,57,64,78,98,109,138,142,150,
             154,174,198,206,213,230,238)
# -----------------------------
# PROCESS DATA
# -----------------------------
df_map$coverage_pct <- (df_map$end - df_map$start) / alu_len * 100
df_map <- df_map[order(df_map$coverage_pct), ]
df_map$y <- seq_len(nrow(df_map))

top_y  <- max(df_map$y)
top_y2 <- top_y + 5

# -----------------------------
# PLOT
# -----------------------------
ggplot() +
  geom_rect(aes(xmin = 0, xmax = 300,
                ymin = top_y2 - 1.5, ymax = top_y2 + 1.5),
            fill = "#ff4d4d", alpha = 0.25) +

  geom_point(data = data.frame(x = cpg_pos),
             aes(x = x, y = top_y2),
             color = "black", size = 1.5) +

  geom_segment(data = df_map,
               aes(x = start, xend = end, y = y, yend = y,
                   color = overlap_count),
               linewidth = 0.6) +

  scale_color_manual(
    values = c(`1` = "steelblue", `2` = "darkorange", `3` = "firebrick"),
    labels = c(`1` = "1 AluY element", `2` = "2 AluY elements", `3` = "3 AluY elements"),
    name = "DMR overlaps"
  ) +

  scale_x_continuous(position = "top",
                     breaks = seq(0, 300, 50),
                     name = "AluY position (bp)") +
  theme_classic() +
  theme(
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title.x.bottom = element_blank(),
    axis.text.x.bottom = element_blank(),
    axis.ticks.x.bottom = element_blank(),
    axis.line.y = element_blank(),

    axis.title.x = element_text(size = 12),
    axis.text.x  = element_text(size = 10),
    legend.title = element_text(size = 12),
    legend.text  = element_text(size = 11)
  )
