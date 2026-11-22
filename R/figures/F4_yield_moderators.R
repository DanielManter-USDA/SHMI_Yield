################################################################################
# F4_yield_moderators.R  (17 Nov 2026)   Fig. 4: when the yield benefit is larger, and how steady it is
# Reads (through moderator_panels.R and yield_panels.R): plotdata_yield_main.rds (B2, fig4ab),
#         plotdata_stability.rds (B3)
# Writes: Outputs/figures/F4_yield_moderators.png
#   a  the within-site corn slope against fertilizer N rate (95% band)
#   b  the slope against time under management
#   c  the slope in poor, typical and good years (county yield, June-August SPEI,
#      crop-insurance causes of loss)
# a and b share one y range (their bands only; the experiments are in Fig. 3b).
################################################################################
source(here::here("R", "figures", "fig_common.R"), local = TRUE)   # local: keeps the pipeline's globals intact
panels <- function(script) {
  e <- new.env(parent = globalenv()); e$FIG_PANELS_ONLY <- TRUE; e$FIG4_SAME_Y <- FALSE
  sys.source(here::here("R", "figures", script), envir = e); e
}
m <- panels("moderator_panels.R")
y <- panels("yield_panels.R")
ab <- read_pd("yield_main")$fig4ab
r <- range(c(ab$N$curve$lower, ab$N$curve$upper, ab$age$curve$lower, ab$age$curve$upper, 0), na.rm = TRUE)
YL <- if (exists("FIG4_YLIM")) FIG4_YLIM else c(floor(r[1] / 100) * 100 - 50, ceiling(r[2] / 100) * 100 + 100)
same_y <- coord_cartesian(ylim = YL)
top <- enlarge(m$pa + same_y) | enlarge(m$pb + same_y + labs(y = NULL))
fig <- top / enlarge(y$pc) + plot_layout(heights = c(1, 0.9))
save_fig(tag_fig(fig), "F4_yield_moderators.png", width = 12, height = 9.5)
