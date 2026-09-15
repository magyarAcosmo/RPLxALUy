# ============================================================================
# Figure 1a – Fraction of sites covered
# ============================================================================
# Control: #2C7BB6
# Patients: #D7191C
# ----------------------------------------------------------------------------
# FIGURE 1b PLOTTING CODE -- requires df2
# ----------------------------------------------------------------------------

library(ggplot2)
library(ggdist)

ggplot(df2, aes(y = frac)) +

    # =========================
    # DENSITY (CENTER)
    # =========================

    stat_halfeye(
        data = subset(df2, group == "Control"),
        aes(x = 1),
        side = "left",
        adjust = 0.6,
        width = 0.35,
        .width = 0,
        point_interval = NULL,
        fill = "#2C7BB6",
        alpha = 0.6
    ) +

    stat_halfeye(
        data = subset(df2, group == "uRPL"),
        aes(x = 1.6),
        side = "right",
        adjust = 0.6,
        width = 0.35,
        .width = 0,
        point_interval = NULL,
        fill = "#D7191C",
        alpha = 0.6
    ) +

    # =========================
    # BOX (OUTER LANE)
    # =========================

    geom_boxplot(
        data = subset(df2, group == "Control"),
        aes(x = 0.55, y = frac),
        width = 0.08,
        outlier.shape = NA,
        fill = "#2C7BB6",
        alpha = 0.85
    ) +

    geom_boxplot(
        data = subset(df2, group == "uRPL"),
        aes(x = 2.05, y = frac),
        width = 0.08,
        outlier.shape = NA,
        fill = "#D7191C",
        alpha = 0.85
    ) +

    # =========================
    # POINTS (OUTER EDGE, COLOR MATCH)
    # =========================

    geom_point(
        data = subset(df2, group == "Control"),
        aes(x = 0.40, y = frac),
        color = "#2C7BB6",
        size = 2,
        alpha = 0.8,
        position = position_jitter(width = 0.005)
    ) +

    geom_point(
        data = subset(df2, group == "uRPL"),
        aes(x = 2.20, y = frac),
        color = "#D7191C",
        size = 2,
        alpha = 0.8,
        position = position_jitter(width = 0.005)
    ) +

    # =========================
    # AXIS
    # =========================

    scale_x_continuous(
        breaks = c(1, 1.6),
        labels = c("Control", "uRPL"),
        limits = c(0.2, 2.4)
    ) +
    theme_classic(base_size = 14) +
    theme(
       legend.position = "none",
        plot.title = element_text(hjust = 0.5, face = "bold")
    ) +

    labs(
        x = NULL,
        y = "Fraction of sites covered",
        title = "Genome Coverage Distribution Across Samples"
    )

