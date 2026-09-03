#!/usr/bin/env Rscript

# clusterProfiler-inspired dotplot for the sperm enrichment table.
# The table is already an enrichment summary, so this script emulates the
# familiar clusterProfiler/enrichplot visual grammar with ggplot2.
#
# Usage:
#   Rscript sperm_enrichment_clusterprofiler_style.R [input_csv] [output_prefix]

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(scales)
})

args <- commandArgs(trailingOnly = TRUE)

input_csv <- if (length(args) >= 1) {
  args[[1]]
} else {
  "/Users/jordanmoore/Desktop/THEWAY/RPL/sperm_enrichment_FINAL_with_perms_noAGE.csv"
}

output_prefix <- if (length(args) >= 2) {
  args[[2]]
} else {
  "/Users/jordanmoore/Documents/Codex/2026-07-13/ca/outputs/sperm_enrichment_clusterprofiler_style"
}

dir.create(dirname(output_prefix), recursive = TRUE, showWarnings = FALSE)

df <- read.csv(input_csv, stringsAsFactors = FALSE, check.names = FALSE) %>%
  mutate(
    Description = gsub("_", " ", antibody),
    neg_log10_q = -log10(q_value),
    significant = q_value < 0.05
  ) %>%
  arrange(oddsRatio) %>%
  mutate(Description = factor(Description, levels = Description))

q_score_breaks <- c(1, 5, 10, 15, 20)
q_score_breaks <- q_score_breaks[
  q_score_breaks <= ceiling(max(df$neg_log10_q, na.rm = TRUE))
]

p <- ggplot(
  df,
  aes(
    x = oddsRatio,
    y = Description,
    size = support,
    fill = neg_log10_q
  )
) +
  geom_vline(
    xintercept = 1,
    linewidth = 0.45,
    linetype = "dashed",
    color = "grey55"
  ) +
  geom_point(
    shape = 21,
    color = "grey25",
    stroke = 0.25,
    alpha = 0.96
  ) +
  geom_point(
    data = df %>% filter(!significant),
    shape = 21,
    fill = NA,
    color = "black",
    stroke = 1.05,
    show.legend = FALSE
  ) +
  geom_text(
    data = df %>% filter(!significant),
    aes(x = oddsRatio * 1.15, y = Description, label = "NS"),
    inherit.aes = FALSE,
    hjust = 0,
    size = 3.2,
    color = "grey15"
  ) +
  scale_x_log10(
    breaks = c(0.125, 0.25, 0.5, 1, 2, 4, 8, 16),
    labels = c("0.125", "0.25", "0.5", "1", "2", "4", "8", "16"),
    limits = c(0.10, 28),
    expand = expansion(mult = c(0.01, 0.04))
  ) +
  scale_size_continuous(
    range = c(3.2, 13.5),
    limits = c(0, 250),
    breaks = c(2, 10, 50, 100, 250),
    name = "Support"
  ) +
  scale_fill_viridis_c(
    option = "viridis",
    direction = 1,
    limits = c(0, ceiling(max(df$neg_log10_q, na.rm = TRUE))),
    breaks = q_score_breaks,
    name = expression(-log[10]("FDR q"))
  ) +
  labs(
    title = "Sperm enrichment across chromatin annotations",
    subtitle = "clusterProfiler-style dotplot using odds ratio as the enrichment axis",
    x = "Odds ratio",
    y = NULL,
    caption = "Size = support. Fill = -log10(FDR q). Black outline/NS marks q >= 0.05. Dashed line = OR 1."
  ) +
  guides(
    fill = guide_colorbar(order = 1, barheight = unit(54, "pt")),
    size = guide_legend(
      order = 2,
      override.aes = list(shape = 21, fill = "grey70", color = "grey25", alpha = 1)
    )
  ) +
  theme_classic(base_family = "Helvetica", base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", size = 15, color = "grey12"),
    plot.subtitle = element_text(size = 10, color = "grey35", margin = margin(b = 12)),
    plot.caption = element_text(size = 8.5, color = "grey45", hjust = 0, margin = margin(t = 10)),
    axis.text.y = element_text(size = 10.5, color = "grey15"),
    axis.text.x = element_text(size = 9.5, color = "grey25"),
    axis.title.x = element_text(size = 10.5, color = "grey15", margin = margin(t = 8)),
    axis.line = element_line(color = "grey30", linewidth = 0.35),
    axis.ticks = element_line(color = "grey30", linewidth = 0.35),
    legend.title = element_text(face = "bold", size = 9.2),
    legend.text = element_text(size = 8.8),
    legend.position = "right",
    plot.margin = margin(12, 18, 10, 10)
  )

if (requireNamespace("ragg", quietly = TRUE)) {
  ggsave(
    paste0(output_prefix, ".png"),
    plot = p,
    device = ragg::agg_png,
    width = 7.2,
    height = 4.6,
    units = "in",
    dpi = 450,
    bg = "white"
  )
} else {
  ggsave(
    paste0(output_prefix, ".png"),
    plot = p,
    width = 7.2,
    height = 4.6,
    units = "in",
    dpi = 450,
    bg = "white"
  )
}

ggsave(
  paste0(output_prefix, ".pdf"),
  plot = p,
  width = 7.2,
  height = 4.6,
  units = "in",
  bg = "white"
)

message("Saved: ", paste0(output_prefix, ".png"))
message("Saved: ", paste0(output_prefix, ".pdf"))
