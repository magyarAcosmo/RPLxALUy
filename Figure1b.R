# ============================================================================
# Figure 1b: Read coverage per site (Control vs. uRPL)
# ============================================================================
#
# INCOMPLETE: this figure's code depends on a data frame called `df_plot`
# (with columns group ["Control"/"uRPL"] and reads), but the pre-figure prep
# script that builds df_plot was not found in either figures_woAge.R or
# pre_figure.R. Please locate/supply that prep code (likely reading a
# coverage/reads-per-site file and reshaping it by group), and it can be
# inserted above the plotting code below.
#
# ----------------------------------------------------------------------------
# FIGURE 1b PLOTTING CODE -- requires df_plot
# ----------------------------------------------------------------------------

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

