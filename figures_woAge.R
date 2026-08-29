# Figures for DMRs without Age included
## patients = dmrs without age (dmrs_2026-01-17_woAGE.bed)

# Control: #2C7BB6
# Patients: #D7191C
# ------------------------------------------------------------------------------------------------------------------
# Figure prep
## 3d: DMR effect file, AluY file -> one row per DMR that overlaps an AluY element (where does DMR fall along consensus)

library(ggplot2)
library(GenomicRanges)

# Read DMR file w/ effect size
dmrs <- read.delim("dmrs_2026-01-17_woAGE_effect.bed", header = FALSE, col.names = c("chr", "start", "end", "effect"))
gr_dmrs <- GRanges(seqnames = dmrs$chr, ranges = IRanges(start = dmrs$start + 1, end = dmrs$end), effect = dmrs$effect)
aluY_clean <- read.delim("AluY_elements_annot.txt", header = TRUE)

# Build ALuY GRanges (genomic coords + consensus info as metadata)
gr_alu <- GRanges(seqnames = aluY_clean$chr,
                  ranges = IRanges(start = aluY_clean$start +1, end = aluY_clean$end),
                  strand = aluY_clean$strand,
                  repName = aluY_clean$repName,
                  consensus_start = aluY_clean$consensus_start,
                  consensus_end = aluY_clean$consensus_end)

# Find overlaps
hits <- findOverlaps(gr_dmrs, gr_alu, type = "any") # full and partial overlaps

# genomic intersection of each overlapping DMR-AluY pair
overlap_ranges <- pintersect(gr_dmrs[queryHits(hits)], gr_alu[subjectHits(hits)])

# map genomic overlap -> consensus coordinates
alu_hit <- gr_alu[subjectHits(hits)]
g_start <- start(alu_hit)
g_end <- end(alu_hit)
cons_start <- alu_hit$consensus_start
cons_end <- alu_hit$consensus_end
over_start <- start(overlap_ranges)
over_end <- end(overlap_ranges)
is_minus <- as.character(strand(alu_hit)) == "-"

# fraction along the genomic span (0 to 1), flipped for minus strand (genomic_end corresponds to consensus_start) -- how far into the element is the DMR?
frac_start <- ifelse(is_minus, (g_end - over_end) / (g_end - g_start),
                     (over_start - g_start) / (g_end - g_start))
frac_end <- ifelse(is_minus, (g_end - over_start) / (g_end - g_start),
                   (over_end - g_start) / (g_end - g_start))

df_map <- data.frame(
  start = cons_start + frac_start * (cons_end - cons_start),  # fraction onto consensus range
  end = cons_start + frac_end * (cons_end - cons_start),
  effect = gr_dmrs$effect[queryHits(hits)],
  repName = alu_hit$repName,
  dmr_chr = as.character(seqnames(gr_dmrs))[queryHits(hits)],
  dmr_start = start(gr_dmrs)[queryHits(hits)],
  dmr_end = end(gr_dmrs)[queryHits(hits)]
)

cat("Total DMRs:", length(gr_dmrs), "\n") # 294
cat("DMRs overlapping an ALuY element:", nrow(df_map), "\n") # 176 (59.86% of DMRs overlap with AluY elements)
sum(df_map$start > df_map$end)  # 0 - coordinate mapping marked for both strand directions!

head(df_map)

# count unique DMRs among 176 overlapping rows

# Use the original DMR coordinates as the unique identifier
dmr_id <- paste(df_map$dmr_chr, df_map$dmr_start, df_map$dmr_end, sep = "_")

n_unique_dmrs <- length(unique(dmr_id))
n_total_rows  <- nrow(df_map)

cat("Total overlap rows:", n_total_rows, "\n")  # 176 total DMRs overlapping with AluY elements
cat("Unique DMRs represented:", n_unique_dmrs, "\n")  # 158 unique DMRs
cat("DMRs overlapping >1 AluY element:", n_total_rows - n_unique_dmrs, "\n")  # 18 DMRs that overlap more than 1 AluY element

# Optional: see which DMRs (if any) hit multiple Alu elements
table(dmr_id)[table(dmr_id) > 1]


# Figure 3d code: DMRs coverage across AluY consensus
library(ggplot2)

## finalize df_map (what I just created with our DMR/AluY data)
df_map <- df_map[, c("start", "end")]   # keep only what original code expects

# -----------------------------
# INPUTS
# -----------------------------
alu_len <- 300
cpg_pos <- c(4,8,10,20,48,53,57,64,78,98,109,138,142,150,
             154,174,198,206,213,230,238)


# new processing to show DMRs that overlap with multiple AluY elements
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
sum(swap_idx)  # should be 2 -- yes

tmp <- aluY_clean$consensus_start[swap_idx]
aluY_clean$consensus_start[swap_idx] <- aluY_clean$consensus_end[swap_idx]
aluY_clean$consensus_end[swap_idx]   <- tmp

nrow(aluY_clean)                                            # should be 151627 -- yes
sum(aluY_clean$consensus_start > aluY_clean$consensus_end)  # should be 0 -- yes

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

cat("Total overlap rows:", nrow(df_map), "\n") #176

# count AluY overlaps per DMR
df_map$dmr_id <- paste(df_map$dmr_chr, df_map$dmr_start, df_map$dmr_end, sep = "_")
dmr_counts <- table(df_map$dmr_id)
df_map$overlap_count <- as.integer(dmr_counts[df_map$dmr_id])

table(dmr_counts)  # check: should show the 1 / 2 / 3 breakdown -- yes, 18 overlapping

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



# -----------------------------------------------------------------------------------
# Figure 3c: Length concordance chart – DMR length and its overlapping AluY length
## df: one row per DMR-AluY overlap pair w/ full length of each DMR & full length of AluY element (not just overlap/intersection region)
## reuse some data processing from Figure 3d

# build df for Figure 3c
# same 176 pairs as figure 3d
df <- data.frame(
  DMR_length = width(gr_dmrs)[queryHits(hits)],
  Alu_length = width(gr_alu)[queryHits(hits)]
)

nrow(df)  # 176

# medians for dashed lined
med_dmr <- median(df$DMR_length)
med_alu <- median(df$Alu_length)

med_dmr #407
med_alu #300

ggplot(df) +

  geom_segment(aes(x = 1, xend = 2,
                   y = DMR_length, yend = Alu_length),
               color = "steelblue", linewidth = 0.5, alpha = 0.7) +
  geom_point(aes(x = 1, y = DMR_length),
             color = "#d2b48c", size = 2) +
  geom_point(aes(x = 2, y = Alu_length),
             color = "#5c3a21", size = 2) +

  # THICK MEDIAN LINES (FIXED)
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
  labs(
    y = "DMR length"
  )+
  scale_y_continuous(breaks = seq(0, max(df$DMR_length, df$Alu_length), 300)) +
  theme_classic()+
  theme(
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 12),
    axis.title.x = element_blank(),

  )

# ------------------------------------------------------------------------------------------------------------------
# Figure 3b– CpG methylation at AluY consensus
ggplot(df_long, aes(x = Position, color = Group, fill = Group)) +

  # IQR ribbon
  geom_ribbon(
    aes(ymin = Donor_LQ, ymax = Donor_HQ),
    alpha = 0.2,
    color = NA
  ) +

  # mean line
  geom_line(
    aes(y = Donor_Mean, group = Group),
    linewidth = 1.2
  ) +

  # mean points
  geom_point(
    aes(y = Donor_Mean),
    size = 2.5
  ) +

  # CpG markers at bottom
  geom_point(
    data = distinct(df_long, Position),
    aes(x = Position, y = 0),
    inherit.aes = FALSE,
    color = "black",
    size = 2
  ) +

  scale_color_manual(values = c(
    "Fertile Donors" = "#1f77b4",
    "uRPL Patients" = "#d62728"
  )) +

  scale_fill_manual(values = c(
    "Fertile Donors" = "#1f77b4",
    "uRPL Patients" = "#d62728"
  )) +

  labs(
    x = "CpG positions at AluY consensus",
    y = "Mean DNA methylation"
  ) +

  theme_classic()

# ------------------------------------------------------------------------------------------------------------------
# Figure 2a : Diverging bar plot
# read DMR effect data
dmrs <- read.delim("dmrs_2026-01-17_woAGE_effect.bed", header = FALSE,
                   col.names = c("chr", "start", "end", "effect"))

# check effect size range to choose reasonable bins
summary(dmrs$effect)  #-0.21 to 0.17
hist(dmrs$effect, breaks = 30)

# bin effect size
bin_width <- 0.02

dmrs$bin <- cut(dmrs$effect, breaks = seq(floor(min(dmrs$effect) / bin_width) * bin_width,
                                                  ceiling(max(dmrs$effect) / bin_width) * bin_width,
                                          by = bin_width), include.lowest = TRUE)
# assign sign based on effect direction
dmrs$sign <- ifelse(dmrs$effect >= 0, "positive", "negative")

# count DMRs per bin, signed
df_counts <- aggregate(effect ~ bin + sign, data = dmrs, FUN = length)
names(df_counts)[3] <- "count"

# make negative-direction counts actually negative (for diverging bars)
df_counts$count_signed <- ifelse(df_counts$sign == "negative", -df_counts$count, df_counts$count)

# order bins properly on the y-axis (low to high effect size)
df_counts$bin <- factor(df_counts$bin, levels = levels(dmrs$bin))

head(df_counts) # correct column names, bins are properly ordered, distribution matches histogram
nrow(df_counts) # 20

max_x <- max(abs(df_counts$count_signed))

p <- ggplot(df_counts, aes(y = bin, x = count_signed, fill = sign)) +
  geom_bar(stat = "identity", width = 0.9) +
  scale_fill_manual(values = c("positive" = "#2C7BB6", "negative" = "#D7191C")) +
  scale_x_continuous(
    labels = abs,
    limits = c(-max_x, max_x)) +
  labs(
    y = "Effect Size",
    x = "Number of DMRs",
    title = "Effect Size Distribution of DMRs") +

  theme_classic() +
  theme(legend.position = "none") +
  geom_vline(xintercept = 0, color = "black")

# -----------------------------------------------------------------------
# Figure 2a & 2b combined -> Volcano plot
## x = methylation difference (hypo - 0 - hyper)
## y = p val 0.05 threshold
# red = hyper, blue = hypo

library(GenomicRanges)

# subset DMRs woAGE to pval
woAGE_effect <- read.delim("dmrs_2026-01-17_woAGE_effect.bed", header = FALSE,
                            col.names = c("chr", "start", "end", "effect"))

pval <- read.delim("DMRs_wAGE/dmrs_2026-01-17_pval.bed", header = FALSE,
                    col.names = c("chr", "start", "end", "name", "pvalue", "strand"))

# convert to GRanges (BED is 0-based half-open; GRanges is 1-based inclusive)
gr_woAGE <- GRanges(seqnames = woAGE_effect$chr,
                     ranges = IRanges(start = woAGE_effect$start + 1, end = woAGE_effect$end),
                     effect = woAGE_effect$effect)

gr_pval <- GRanges(seqnames = pval$chr,
                    ranges = IRanges(start = pval$start + 1, end = pval$end),
                    pvalue = pval$pvalue)

hits <- findOverlaps(gr_woAGE, gr_pval, type = "equal")

cat("woAGE DMRs:", length(gr_woAGE), "\n")  #294
cat("Matched to a p-value:", length(hits), "\n")  #294

df_volcano <- data.frame(
  chr    = as.character(seqnames(gr_woAGE))[queryHits(hits)],
  start  = start(gr_woAGE)[queryHits(hits)],
  end    = end(gr_woAGE)[queryHits(hits)],
  effect = gr_woAGE$effect[queryHits(hits)],
  pvalue = gr_pval$pvalue[subjectHits(hits)]
)

# direction label + -log10(p) FOR VOLCANO Y-AXIS
df_volcano$direction <- ifelse(df_volcano$effect >= 0, "hyper", "hypo")
df_volcano$neg_log10_p <- -log10(df_volcano$pvalue)

head(df_volcano)
nrow(df_volcano)  # 294: all matched

# volcano plot
library(ggplot2)

sig_threshold <- 2

ggplot(df_volcano, aes(x = effect, y = neg_log10_p, color = direction)) +
  geom_point(alpha = 0.7, size = 2) +
  scale_color_manual(values = c("hyper" = "#D7191C", "hypo" = "#2C7BB6"),
                     labels = c("hyper" = "Hypermethylated", hypo = "Hypomethylated"),
                      name = "Methylation Status") +
  geom_vline(xintercept = 0, color = "black") +
  geom_hline(yintercept = sig_threshold, linetype = "dashed", color = "grey40") +
  scale_x_continuous(
    limits = c(-1, 1),
    breaks = seq(-1, 1, by = 0.1),
    labels = function(x) ifelse(round(x * 10) %% 2 == 0, sprintf("%.1f", x), ""))+
  labs(
    x = "Methylation Difference",
    y = expression(paste("-log"[10], " (p-value)"))
  ) +
  theme_classic() +
  theme(axis.title = element_text(size = 12),
        legend.position = "inside",
        legend.position.inside = c(0.95, 0.95),
        legend.justification = c("right", "top"),
        legend.background = element_rect(fill ="white", color = "grey70"))

# -----------------------------------------------------------------------
# Figure 1b: Coverage – reads per site
library(ggplot2)
library(ggdist)

ggplot(df_plot, aes(y = reads)) +

    # DENSITY (Control)
    stat_halfeye(
        data = subset(df_plot, group == "Control"),
        aes(x = 1),
        side = "left",
        adjust = 0.5,
        width = 0.45,
        .width = 0,
        point_interval = NULL,
        fill = "#2C7BB6",
        alpha = 0.6
    ) +

    # DENSITY (uRPL)
    stat_halfeye(
        data = subset(df_plot, group == "uRPL"),
        aes(x = 2),
        side = "right",
        adjust = 0.5,
        width = 0.45,
        .width = 0,
        point_interval = NULL,
        fill = "#D7191C",
        alpha = 0.6
    ) +

    # BOX (OUTSIDE LEFT)
    geom_boxplot(
        data = subset(df_plot, group == "Control"),
        aes(x = 0.35, y = reads),
        width = 0.10,
        outlier.shape = NA,
        fill = "#2C7BB6",
        alpha = 0.85
    ) +

    # BOX (OUTSIDE RIGHT)
    geom_boxplot(
        data = subset(df_plot, group == "uRPL"),
        aes(x = 2.65, y = reads),
        width = 0.10,
        outlier.shape = NA,
        fill = "#D7191C",
        alpha = 0.85
    ) +

    # POINTS (OUTSIDE + MATCH COLOR)
    geom_point(
        data = subset(df_plot, group == "Control"),
        aes(x = 0.25, y = reads),
        color = "#2C7BB6",
        size = 2,
        alpha = 0.8,
        position = position_jitter(width = 0.02)
    ) +

    geom_point(
        data = subset(df_plot, group == "uRPL"),
        aes(x = 2.75, y = reads),
        color = "#D7191C",
        size = 2,
        alpha = 0.8,
        position = position_jitter(width = 0.02)
    ) +

    # AXIS
    scale_x_continuous(
        breaks = c(1, 2),
        labels = c("Control", "uRPL"),
        limits = c(0, 3)
    ) +

    theme_classic(base_size = 14) +

    theme(
        legend.position = "none",
        plot.title = element_text(hjust = 0.5, face = "bold")
    ) +

    labs(
        x = NULL,
        y = "Number of reads per site",
        title = "Read Distribution Across Samples"
    )

