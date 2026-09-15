# ============================================================================
# Figure 2b: Volcano plot (methylation difference vs. -log10 p-value)
# ============================================================================
# hyper = red (#D7191C), hypo = blue (#2C7BB6)
# ----------------------------------------------------------------------------
# PRE-FIGURE
# Produces: dmrs_2026-01-17_woAGE_effect.bed
# NOTE: this figure also reads DMRs_wAGE/dmrs_2026-01-17_pval.bed directly --
# that file is a raw input, not something previously generates.
# ----------------------------------------------------------------------------

# DMR EFFECT
## match up woAGE DMRs with AGE DMRs in dmrs_2026-01-17_effect.bed
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

# ----------------------------------------------------------------------------
# FIGURE 2b CODE
# ----------------------------------------------------------------------------

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
