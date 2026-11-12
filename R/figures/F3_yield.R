################################################################################
# F3_yield.R  (11 Nov 2026)   Fig. 3: crop yield and SHMI within sites
# Reads:  plotdata_yield_main.rds (B2, fig3a), plotdata_soybean.rds (B2b),
#         plotdata_stability.rds (B3)
# Writes: Outputs/figures/F3_yield.png (and each panel alone: F3a_corn.png, F3b_soybean.png,
#         F3c_years.png)
#   a  corn yield against SHMI: one point per management system at a site (mean adjusted
#      yield +/- SE); line and band: the within-site slope at 150 kg N/ha and average age
#   b  the same for soybean
#   c  the corn slope in poor, typical and good years by three definitions (county yield,
#      June-August SPEI, crop-insurance causes of loss)
# Options (set before sourcing):
#   FIG3_ONE_LEGEND  TRUE (default): one legend of all experiments, under a and b (FALSE: each panel its own)
################################################################################
source(here::here("R", "figures", "fig_common.R"))
if (!exists("FIG3_ONE_LEGEND")) FIG3_ONE_LEGEND <- TRUE
corn <- read_pd("yield_main")$fig3a
soy  <- read_pd("soybean")
yrs  <- read_pd("stability")

# ---- a, b: yield against SHMI ----------------------------------------------------------
sites_all <- unique(c(as.character(corn$points$site), as.character(soy$points$site)))
sites_all <- names(corn$site_cols)[names(corn$site_cols) %in% sites_all]   # legend in B2's site order
yield_panel <- function(d, col, ylab) {
  ggplot() +
    geom_ribbon(data = d$line, aes(SHMI, ymin = lower, ymax = upper), fill = col, alpha = 0.15) +
    geom_line(data = d$line, aes(SHMI, fit), colour = col, linewidth = 1.1) +
    geom_errorbar(data = d$points, aes(SHMI, ymin = y - se, ymax = y + se), width = 0, colour = "grey35") +
    geom_point(data = d$points, aes(SHMI, y, fill = site), shape = 21, size = 2.6, colour = "black") +
    scale_fill_manual(values = corn$site_cols, limits = if (FIG3_ONE_LEGEND) sites_all else NULL, drop = FALSE, name = NULL,
                      guide = guide_legend(nrow = if (FIG3_ONE_LEGEND) 2 else NULL,
                                           override.aes = list(shape = 21, size = 3, colour = "black"))) +
    corner_labels(lab_slope(d$slope, d$lower, d$upper), lab_p(d$p)) +
    labs(x = XLAB_SHMI, y = ylab) +
    theme_fig
}
pa <- yield_panel(corn, "darkgreen", YLAB_CORN_ADJ)
pb <- yield_panel(soy, "#2e7d32", paste0("Soybean yield adjusted for site, year\nand other treatments (", KG_HA, ")"))
# one legend for a and b: all experiments, from panel a only, in two rows under a and b
if (FIG3_ONE_LEGEND) pb <- pb + guides(fill = "none")

# ---- c: the slope in poor, typical and good years ---------------------------------------
pc <- ggplot(yrs$years, aes(cond, b, colour = cond)) +
  geom_hline(yintercept = 0, colour = "grey60", linetype = "dashed") +
  geom_pointrange(aes(ymin = lower, ymax = upper), size = 0.4) +
  geom_text(aes(y = upper, label = lab), vjust = -0.6, size = 2.3, colour = "grey30") +
  facet_grid(. ~ head, scales = "free_x", space = "free_x") +
  scale_colour_manual(values = yrs$colours, guide = "none") +
  labs(x = NULL, y = paste0("Corn yield change per 10 SHMI\n(", KG_HA, ", 95% CI)")) +
  theme_fig + theme(strip.background = element_blank(), strip.text = element_text(face = "bold", size = 8),
                    axis.text.x = element_text(angle = 30, hjust = 1), panel.grid.major.x = element_blank())

# ---- each panel alone (for checking) and the figure -----------------------------------------
save_fig(pa + panel_size(3.6, 2.6), "F3a_corn.png", width = 6.8, height = 4.6)
save_fig(pb + panel_size(3.6, 2.6), "F3b_soybean.png", width = 6.8, height = 4.6)
save_fig(pc, "F3c_years.png", width = 7.2, height = 3.8)
top <- if (FIG3_ONE_LEGEND) {
  (pa | pb) + plot_layout(guides = "collect") & theme(legend.position = "bottom", legend.text = element_text(size = 8))
} else (pa | pb)
bottom <- (plot_spacer() | pc | plot_spacer()) + plot_layout(widths = c(0.1, 0.8, 0.1))   # c centred at ~80% width, as before
fig <- (top / bottom) + plot_layout(heights = c(1, 0.9)) + plot_annotation(tag_levels = "a") & TAG_THEME
save_fig(fig, "F3_yield.png", width = 12, height = 9.2)
