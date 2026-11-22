################################################################################
# moderator_panels.R  (was F4_yield_moderators.R; renamed 12 Nov 2026)
# Builds pa (SHMI x N), pb (SHMI x age) and pc (each experiment's slope) for Figs 4 and 5.
# Sourced by F4_yield.R and F5_outcome.R with FIG_PANELS_ONLY = TRUE; run alone, it still
# writes its own figure (Outputs/figures/moderator_panels.png).
# Reads:  plotdata_yield_main.rds (B2_SHMI_yield_main.R), elements fig4ab and fig4c
# Writes: Outputs/figures/F4_yield_moderators.png
#   a  the within-site slope against N applied (95% band; points at 50, 100, 150, 200)
#   b  the slope against time under management (points at 1, 10, 30 years)
#   c  each experiment's own slope at 150 kg N/ha and the average age, experiments
#      ordered by the share of their SHMI spread that came from crop diversity; filled
#      if its 95% CI excludes zero; dashed green line: the slope of Fig. 3a
# Options (set before sourcing):
#   FIG4C_COLOUR  "diversity" (default): colour by the diversity share, disturbance share
#                 printed under each experiment; "disturbance": the reverse
#   FIG4_SAME_Y   TRUE (default): one y range for a-c, intervals beyond it cut (arrows)
#   FIG4_YLIM     c(lo, hi) to set that range
#   FIG4_ONE_ROW  FALSE (default): a and b above c; TRUE: a | b | c
################################################################################
source(here::here("R", "figures", "fig_common.R"), local = TRUE)   # local: keeps the pipeline's globals intact
pd <- read_pd("yield_main")
if (!exists("FIG4C_COLOUR")) FIG4C_COLOUR <- "diversity"
if (!exists("FIG4_SAME_Y"))  FIG4_SAME_Y  <- TRUE
if (!exists("FIG4_ONE_ROW")) FIG4_ONE_ROW <- FALSE
if (!exists("FIG_PANELS_ONLY")) FIG_PANELS_ONLY <- FALSE
c4 <- pd$fig4c

# ---- shared y range ---------------------------------------------------------------
YL <- if (exists("FIG4_YLIM")) FIG4_YLIM else {
  r <- range(c(c4$bands$lo, c4$bands$hi, c4$sites$slope_adj, c4$groups$lo, c4$groups$hi, 0), na.rm = TRUE)
  c(floor(r[1] / 100) * 100 - 50, ceiling(r[2] / 100) * 100 + 50)
}
same_y <- if (FIG4_SAME_Y) coord_cartesian(ylim = YL) else NULL

# ---- a, b: the slope against N and against age --------------------------------------
slope_panel <- function(m, xlab, p_lab, col, label_right = FALSE) {
  ggplot(m$curve, aes(x, b)) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey60") +
    geom_ribbon(aes(ymin = lower, ymax = upper), fill = col, alpha = 0.18) +
    geom_line(colour = col, linewidth = 1.2) +
    geom_point(data = m$marks, size = 2.6, colour = col) +
    geom_text(data = m$marks, aes(label = lab), vjust = -1.1, size = 3) +
    geom_rug(data = tibble(x = m$rug), aes(x = x), inherit.aes = FALSE, sides = "b", alpha = 0.25, colour = "grey30") +
    annotate("text", x = if (label_right) Inf else -Inf, y = Inf, hjust = if (label_right) 1.05 else -0.05,
             vjust = 1.4, size = 2.8, label = paste0("Yield change per 10 SHMI, 95% band; ", lab_p(m$p, p_lab))) +
    labs(x = xlab, y = YLAB_SLOPE) + theme_fig
}
pa <- slope_panel(pd$fig4ab$N, "N applied (kg N/ha)\nfertilizer + available manure N + legume credit",
                  "SHMI × N", "#b30000", label_right = TRUE)
pb <- slope_panel(pd$fig4ab$age, "Time under management (years)", "SHMI × age", "#238b45")

# ---- c: each experiment, ordered by diversity share -------------------------------------
col_div <- FIG4C_COLOUR == "diversity"
dd <- c4$sites %>%
  arrange(share_div, slope_adj) %>% mutate(i = row_number()) %>%
  mutate(lo_c = pmax(slope_adj_lower, YL[1]), hi_c = pmin(slope_adj_upper, YL[2]),
         cut_lo = slope_adj_lower < YL[1], cut_hi = slope_adj_upper > YL[2],
         ccol = pmax(if (col_div) share_div else share_dist, 0),
         xlab = if (col_div) sprintf("%s\nDist. %.0f%%", as.character(lab), 100 * pmax(share_dist, 0))
                else sprintf("%s\n%.0f%%", as.character(lab), 100 * share_div))
DCOL <- if (col_div) c("#E5F5E0", "#006D2C") else c("#56B4E9", "#D55E00")
CLIM <- if (col_div) c(0, max(0.5, ceiling(10 * max(dd$ccol)) / 10)) else c(0, 1)
pc <- ggplot(dd, aes(i, slope_adj)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey60") +
  geom_hline(yintercept = c4$global_slope, colour = "darkgreen", linetype = "dashed") +
  geom_linerange(aes(ymin = lo_c, ymax = hi_c, colour = ccol), linewidth = 0.6) +
  geom_point(data = filter(dd, cut_hi), aes(y = YL[2], colour = ccol), shape = 17, size = 1.8) +
  geom_point(data = filter(dd, cut_lo), aes(y = YL[1], colour = ccol), shape = 25, size = 1.8, fill = "grey40") +
  geom_point(aes(colour = ccol, fill = ifelse(sig, ccol, NA)), shape = 21, size = 2.8, stroke = 1) +
  scale_x_continuous(breaks = dd$i, labels = dd$xlab, expand = expansion(add = 0.5)) +
  scale_colour_gradient(low = DCOL[1], high = DCOL[2], limits = CLIM, labels = function(v) paste0(round(100 * v), "%"),
                        name = if (col_div) "Diversity\nshare" else "Disturbance\nshare") +
  scale_fill_gradient(low = DCOL[1], high = DCOL[2], limits = CLIM, na.value = "white", guide = "none") +
  labs(x = if (col_div) "Experiment (SHMI range; share of its SHMI spread from Inverse Disturbance), ordered by diversity share"
           else "Experiment (SHMI range), ordered by the diversity share of its SHMI spread", y = YLAB_SLOPE) +
  theme_fig + theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
                    panel.grid.major.x = element_blank(), panel.grid.minor = element_blank(),
                    legend.position = c(0.99, 0.02), legend.justification = c(1, 0),
                    legend.key.height = unit(8, "pt"), legend.title = element_text(size = 7),
                    legend.text = element_text(size = 6.5), legend.background = element_blank())
if (FIG4_SAME_Y) pc <- pc + coord_cartesian(ylim = YL, clip = "off")

# ---- assemble -----------------------------------------------------------------------
if (!FIG_PANELS_ONLY) {
  if (FIG4_ONE_ROW) {
    fig <- wrap_plots(pa + same_y, pb + same_y + labs(y = NULL), pc + labs(y = NULL), nrow = 1, widths = c(1, 1, 1.35))
    w <- 13; h <- 4.6
  } else {
    fig <- (pa + same_y | pb + same_y + labs(y = NULL)) / pc + plot_layout(heights = c(1, 1.15))
    w <- 9; h <- 8.4
  }
  fig <- fig + plot_annotation(tag_levels = "a") & TAG_THEME
  save_fig(fig, "moderator_panels.png", width = w, height = h)
  cut <- dd$site[dd$cut_lo | dd$cut_hi]
  message(sprintf("Fig. 4: y range %.0f to %.0f; intervals cut: %s", YL[1], YL[2], if (length(cut)) paste(cut, collapse = ", ") else "none"))
}
