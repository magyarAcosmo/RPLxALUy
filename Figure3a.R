# ============================================================================
# Figure 3b: CpG methylation across the AluY consensus
# ============================================================================
#
# INCOMPLETE: this figure's code depends on a data frame called `df_long`
# (with columns Position, Group, Donor_Mean, Donor_LQ, Donor_HQ), but the
# pre-figure prep script that builds df_long was not found in either
# figures_woAge.R or pre_figure.R. Please locate/supply that prep code
# (likely a separate script that summarizes per-CpG methylation by group),
# and it can be inserted above the plotting code below.
#
# ----------------------------------------------------------------------------
# FIGURE 3b PLOTTING CODE -- requires df_long
# ----------------------------------------------------------------------------

library(ggplot2)

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
