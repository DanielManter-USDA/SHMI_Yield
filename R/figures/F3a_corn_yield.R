################################################################################
# F3a_corn_yield.R  (11 Nov 2026)   Fig. 3a: corn yield against SHMI within sites
# Reads:  plotdata_yield_main.rds (B2_SHMI_yield_main.R), element fig3a
# Writes: Outputs/figures/F3a_corn_yield.png
# Points: one per management system at a site (or per SHMI level, as B2 was run),
# mean adjusted yield +/- SE; line and band: the within-site slope at 150 kg N/ha and
# the average system age.
################################################################################
source(here::here("R", "figures", "fig_common.R"))
d <- read_pd("yield_main")$fig3a

p <- ggplot() +
  geom_ribbon(data = d$line, aes(SHMI, ymin = lower, ymax = upper), fill = "darkgreen", alpha = 0.15) +
  geom_line(data = d$line, aes(SHMI, fit), colour = "darkgreen", linewidth = 1.1) +
  geom_errorbar(data = d$points, aes(SHMI, ymin = y - se, ymax = y + se), width = 0, colour = "grey35") +
  geom_point(data = d$points, aes(SHMI, y, fill = site), shape = 21, size = 2.6, colour = "black") +
  scale_fill_manual(values = d$site_cols, name = NULL) +
  corner_labels(lab_slope(d$slope, d$lower, d$upper), lab_p(d$p)) +
  labs(x = XLAB_SHMI, y = YLAB_CORN_ADJ) +
  theme_fig

save_fig(p + panel_size(3.6, 2.6), "F3a_corn_yield.png", width = 6.8, height = 4.6)
