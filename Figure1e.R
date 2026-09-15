# ============================================================================
# Figure 1e: Mean methylation per sample
# ============================================================================
# Control: #2C7BB6
# Patients: #D7191C
# ----------------------------------------------------------------------------
# PRE-FIGURE PREP - requires df
# ----------------------------------------------------------------------------


library(ggplot2)
library(ggdist)

ggplot(df, aes(y = meth)) +

    # =========================
    # DENSITY (CENTER)
    # =========================

    stat_halfeye(
        data = subset(df, group == "Control"),
        aes(x = 1),
        side = "left",
        adjust = 0.5,
        width = 0.40,
        .width = 0,
        point_interval = NULL,
        fill = "#2C7BB6",
        alpha = 0.6
    ) +

    stat_halfeye(
        data = subset(df, group == "uRPL"),
        aes(x = 1.6),
        side = "right",
        adjust = 0.5,
        width = 0.40,
        .width = 0,
        point_interval = NULL,
        fill = "#D7191C",
        alpha = 0.6
    ) +

    # =========================
    # BOX (PUSHED OUTSIDE KDE)
    # =========================

    geom_boxplot(
        data = subset(df, group == "Control"),
        aes(x = 0.55, y = meth),
        width = 0.08,
        outlier.shape = NA,
        fill = "#2C7BB6",
        alpha = 0.85
    ) +

    geom_boxplot(
        data = subset(df, group == "uRPL"),
        aes(x = 2.05, y = meth),
        width = 0.08,
        outlier.shape = NA,
        fill = "#D7191C",
        alpha = 0.85
    ) +

    # =========================
    # POINTS (FURTHER OUTSIDE BOX)
    # =========================

    geom_point(
        data = subset(df, group == "Control"),
        aes(x = 0.40, y = meth),
        color = "#2C7BB6",
        size = 2,
        alpha = 0.8,
        position = position_jitter(width = 0.01)
    ) +

    geom_point(
        data = subset(df, group == "uRPL"),
        aes(x = 2.20, y = meth),
        color = "#D7191C",
        size = 2,
        alpha = 0.8,
        position = position_jitter(width = 0.01)
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
        y = "Mean methylation per sample",
        title = "Sperm DNA Methylation Distribution"
    )