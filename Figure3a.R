# ============================================================================
# Figure 3b: CpG methylation across the AluY consensus
# ============================================================================
# df_long = 3a_methylation_by_position_in_AluY.csv
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
