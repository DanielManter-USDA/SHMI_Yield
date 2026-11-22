################################################################################
# F3_yield.R  (17 Nov 2026)   Fig. 3: SHMI predicts crop yield in independent experiments
# Reads (through yield_panels.R): plotdata_yield_main.rds (B2, fig3a and fig4c),
#         plotdata_soybean.rds (B2b)
# Writes: Outputs/figures/F3_yield.png
#   a  corn and soybean yield against SHMI within sites (pooled slopes; one legend)
#   b  each experiment's own corn slope at 150 kg N/ha and the average age (95% CI; filled:
#      interval excludes zero), ordered by slope. Dashed green line: the pooled slope of a.
#      Experiments whose SHMI range came almost entirely from tillage (Inverse Disturbance
#      share of the within-site SHMI variance >= TILLAGE_ONLY, B2's own rule) are drawn in
#      orange and labelled "tillage only". Intervals beyond the y range are cut (arrows).
################################################################################
source(here::here("R", "figures", "fig_common.R"), local = TRUE)   # local: keeps the pipeline's globals intact
if (!exists("TILLAGE_ONLY")) TILLAGE_ONLY <- 0.8
e <- new.env(parent = globalenv()); e$FIG_PANELS_ONLY <- TRUE
sys.source(here::here("R", "figures", "yield_panels.R"), envir = e)

# ---- a: corn | soybean ----------------------------------------------------------------
pa <- (enlarge(e$pa) | enlarge(e$pb)) + plot_layout(guides = "collect") &
  theme(legend.position = "bottom", legend.text = element_text(size = 8 * FONT_SCALE))

# ---- b: each experiment ---------------------------------------------------------------
c4 <- read_pd("yield_main")$fig4c
COL_T <- "#b35806"; COL_O <- "darkgreen"
dd <- c4$sites %>%
  mutate(sig = slope_adj_lower > 0 | slope_adj_upper < 0,
         till = !is.na(share_dist) & share_dist >= TILLAGE_ONLY,
         col = ifelse(till, COL_T, COL_O),
         xlab = sprintf("%s\n(%.0f–%.0f)%s", site, SHMI_min, SHMI_max, ifelse(till, "\ntillage only", ""))) %>%
  arrange(slope_adj) %>% mutate(i = row_number())
r <- range(c(dd$slope_adj, 0, c4$global_slope), na.rm = TRUE)
YL <- if (exists("FIG3B_YLIM")) FIG3B_YLIM else c(floor(r[1] / 100) * 100 - 300, ceiling(r[2] / 100) * 100 + 300)
dd <- dd %>% mutate(lo_c = pmax(slope_adj_lower, YL[1]), hi_c = pmin(slope_adj_upper, YL[2]),
                    cut_lo = slope_adj_lower < YL[1], cut_hi = slope_adj_upper > YL[2])
pb <- ggplot(dd, aes(i, slope_adj)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey60") +
  geom_hline(yintercept = c4$global_slope, colour = COL_O, linetype = "dashed") +
  geom_linerange(aes(ymin = lo_c, ymax = hi_c, colour = col), linewidth = 0.7) +
  geom_point(data = filter(dd, cut_hi), aes(y = YL[2], colour = col), shape = 17, size = 2) +
  geom_point(data = filter(dd, cut_lo), aes(y = YL[1], colour = col), shape = 25, size = 2, fill = "grey40") +
  geom_point(aes(colour = col, fill = ifelse(sig, col, "white")), shape = 21, size = 3, stroke = 1) +
  scale_colour_identity() + scale_fill_identity() +
  scale_x_continuous(breaks = dd$i, labels = dd$xlab, expand = expansion(add = 0.5)) +
  coord_cartesian(ylim = YL, clip = "off") +
  labs(x = "Experiment (SHMI range)", y = YLAB_SLOPE) +
  theme_fig + theme(panel.grid.major.x = element_blank(), panel.grid.minor = element_blank(),
                    axis.text.x = element_text(size = 7.5, lineheight = 0.9))

fig <- wrap_elements(full = pa) / enlarge(pb) + plot_layout(heights = c(1.15, 0.85))
save_fig(tag_fig(fig), "F3_yield.png", width = 12, height = 10.5)
cut <- dd$site[dd$cut_lo | dd$cut_hi]
message(sprintf("Fig. 3b: y range %.0f to %.0f; intervals cut: %s; tillage only: %s", YL[1], YL[2],
                if (length(cut)) paste(cut, collapse = ", ") else "none",
                if (any(dd$till)) paste(dd$site[dd$till], collapse = ", ") else "none"))
