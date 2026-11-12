################################################################################
# fig_common.R  (11 Nov 2026)
# Shared look and helpers for the figure scripts in R/figures/. Each figure script
# reads plot data written by an analysis script (Outputs/.../plotdata_*.rds) and only
# draws, so a figure can be changed and redrawn in seconds without rerunning the
# analysis. Change fonts, sizes and colours here once for every figure.
################################################################################
suppressPackageStartupMessages({ library(dplyr); library(ggplot2); library(patchwork) })

FIG_DIR <- here::here("Outputs", "figures")          # where the figure scripts write
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)

# where each analysis writes its plot data (the analysis folders as run for the paper)
PLOTDATA <- list(
  yield_main = here::here("Outputs", "yield_main_default_total_agcros", "plotdata_yield_main.rds"),   # B2
  soybean    = here::here("Outputs", "yield_soybean", "plotdata_soybean.rds"),                          # B2b
  stability  = here::here("Outputs", "yield_stability_default_total_agcros", "plotdata_stability.rds"))  # B3

read_pd <- function(name) {
  f <- PLOTDATA[[name]]
  if (is.null(f) || !file.exists(f))
    stop("No plot data for '", name, "' (", f, "). Run the analysis script that writes it first.")
  readRDS(f)
}

# ---- the look -----------------------------------------------------------------
theme_fig <- theme_bw(base_size = 10) +
  theme(plot.caption = element_text(size = 7.5, colour = "grey30", hjust = 0), plot.caption.position = "plot")
TAG_THEME <- theme(plot.tag = element_text(face = "bold", size = 14))

# ---- labels -------------------------------------------------------------------
KG_HA     <- "kg ha⁻¹"
XLAB_SHMI <- "Soil Health Management Index (SHMI)"
YLAB_CORN_ADJ <- paste0("Corn yield adjusted for site, year,\nN and other treatments (", KG_HA, ")")
YLAB_SLOPE    <- paste0("Corn yield change\nper 10 SHMI (", KG_HA, ")")
lab_slope <- function(est, lo, hi, what = "Yield change per 10 SHMI")
  sprintf("%s: %+.0f kg/ha (95%% CI %.0f%s%.0f)", what, est, lo, if (lo < 0) " to " else "–", hi)
lab_p <- function(p, what = NULL)
  paste0(if (!is.null(what)) paste0(what, ": "), ifelse(p < 0.001, "p < 0.001", sprintf("p = %.3f", p)))
corner_labels <- function(top, bottom) list(
  annotate("text", x = -Inf, y = Inf,  hjust = -0.1, vjust = 1.1, size = 2.8, label = top),
  annotate("text", x = Inf,  y = -Inf, hjust = 1.1,  vjust = -1,  size = 2.8, label = bottom))

# fixed panel size, as in the analysis scripts' figures (ggh4x)
panel_size <- function(w_in, h_in) ggh4x::force_panelsizes(rows = unit(h_in, "in"), cols = unit(w_in, "in"))

save_fig <- function(p, file, width, height, dpi = 300) {
  ggsave(file.path(FIG_DIR, file), p, width = width, height = height, dpi = dpi, bg = "white")
  message("  wrote ", file.path("Outputs", "figures", file))
  invisible(file.path(FIG_DIR, file))
}
