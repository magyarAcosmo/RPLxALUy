# ============================================================================
# Figure 3c: Length concordance chart
# (DMR length vs. length of its overlapping AluY element)
# ============================================================================

# ----------------------------------------------------------------------------
# PRE-FIGURE PREP
# Produces: dmrs_2026-01-17_woAGE_effect.bed, AluY_elements_annot.txt
# ----------------------------------------------------------------------------

# DMR EFFECT
## match up noAGE DMRs with AGE DMRs in dmrs_2026-01-17_effect.bed
library(GenomicRanges)
library(rtracklayer)
library(ggplot2)

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
hits <- findOverlaps(gr_woAGE, gr_effect, type = "equal")

# Build Output
df_out <- woAGE[queryHits(hits), ]
df_out$effect_value <- effect$effect_value[subjectHits(hits)]

# check if every woAGE DMR got a match
cat("woAGE DMRs:", nrow(woAGE), "\n")
cat("Matched DMRs:", nrow(df_out), "\n")

write.table(df_out,
            file = "dmrs_2026-01-17_woAGE_effect.bed",
            sep = "\t",
            quote = FALSE,
            row.names = FALSE,
            col.names = FALSE)

# AluY coordinates
## rmsk.txt.gz contains all AluY coordinates; include all subfamilies that start with AluY

# read rmsk.txt file
rmsk <- read.delim("rmsk.txt", header = FALSE)

# filter: repName (col 11) starts with "AluY"
aluY <- rmsk[grepl("^AluY", rmsk[, 11]), ]

table(aluY[, 11])   # 28 subfamilies

cat("Total AluY elements found:", nrow(aluY), "\n")   # 151627

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

# ----------------------------------------------------------------------------
# FIGURE 3c CODE
# ----------------------------------------------------------------------------

df <- data.frame(
  DMR_length = width(gr_dmrs)[queryHits(hits)],
  Alu_length = width(gr_alu)[queryHits(hits)]
)

nrow(df)  # 176

# medians for dashed lined
med_dmr <- median(df$DMR_length)
med_alu <- median(df$Alu_length)

ggplot(df) +

  geom_segment(aes(x = 1, xend = 2,
                   y = DMR_length, yend = Alu_length),
               color = "steelblue", linewidth = 0.5, alpha = 0.7) +
  geom_point(aes(x = 1, y = DMR_length),
             color = "#d2b48c", size = 2) +
  geom_point(aes(x = 2, y = Alu_length),
             color = "#5c3a21", size = 2) +

  # THICK MEDIAN LINES
  geom_hline(yintercept = med_dmr,
             linetype = "dashed",
             color = "#d2b48c",
             linewidth = 2) +
  geom_hline(yintercept = med_alu,
             linetype = "dashed",
             color = "#5c3a21",
             linewidth = 2) +
  scale_x_continuous(breaks = c(1,2),
                     labels = c("DMR length", "AluY length")) +
  labs(y = "DMR length")+
  scale_y_continuous(breaks = seq(0, max(df$DMR_length, df$Alu_length), 300)) +
  theme_classic()+
  theme(
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 12),
    axis.title.x = element_blank()
  )