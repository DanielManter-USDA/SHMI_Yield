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

# 12 Nov 2026: panels saved as ggplot objects by the soil-health scripts (Outputs/panels and the
# folders below). A figure script can restyle them (theme, labels, scales, layout); changing
# what a panel plots still means rerunning the script that saves it.
PANEL_FILES <- list(
  path          = here::here("Outputs", "panels", "panel_path.rds"),                     # C1
  additivity    = here::here("Outputs", "panels", "panel_additivity.rds"),               # C1 (from C1d)
  label_vs_shmi = here::here("Outputs", "panels", "panel_label_vs_SHMI.rds"),            # C1g
  label_within  = here::here("Outputs", "panels", "panel_label_within.rds"),             # C1g
  shmi_vs_sh    = here::here("Outputs", "panels", "panel_SHMI_vs_SH.rds"),               # C1h (list: all, between, within)
  shmi_vs_shape = here::here("Outputs", "panels", "panel_SHMI_vs_SHAPE.rds"),            # C1h (same)
  records       = here::here("Outputs", "panels", "panel_records.rds"),                  # C1c (list: a, b)
  data_budget   = here::here("Outputs", "Practice_additivity", "Fig_data_budget_panels.rds"),  # C1f (list)
  parts         = here::here("Outputs", "Manuscript_work", "Fig2a_panel.rds"),           # A4b
  meta_reg      = here::here("Outputs", "panels", "panel_meta_reg.rds"),                 # C1i
  meta_data     = here::here("Outputs", "Pair_meta", "plotdata_pair_meta.rds"))          # C1i (pooled, metareg, pairs)
read_panel <- function(name) {
  f <- PANEL_FILES[[name]]
  if (is.null(f) || !file.exists(f))
    stop("No saved panel '", name, "' (", f, "). Run the script that saves it first.")
  readRDS(f)
}
# result tables read directly by a figure script
TABLE_FILES <- list(
  data_budget = here::here("Outputs", "Practice_additivity", "data_budget.csv"),     # C1f (Fig. 3c)
  sem       = here::here("Outputs", "sem_noMX_standardized.csv"),                 # A8 (Part A)
  benchmark = here::here("Outputs", "Benchmark", "benchmark_summary.csv"),        # A4a (Part A)
  yield_parts = here::here("Outputs", "yield_parts", "yield_parts.csv"),          # B2d
  label_vs_shmi = here::here("Outputs", "Practice_additivity", "label_vs_SHMI.csv"))   # C1g
read_tab <- function(name) {
  f <- TABLE_FILES[[name]]
  if (is.null(f) || !file.exists(f)) stop("No table '", name, "' (", f, "). Run the script that writes it first.")
  read.csv(f, stringsAsFactors = FALSE, check.names = FALSE)
}
# assemble and save a multi-panel figure with bold letters a, b, c ...
tag_fig <- function(p) p + plot_annotation(tag_levels = "a") & TAG_THEME

# ---- text size (13 Nov 2026) ---------------------------------------------------------------
# Figures are drawn ~12.5 in wide and printed at ~7 in (180 mm), i.e. at ~57%. FONT_SCALE
# enlarges every text element and in-panel label so text prints at ~7-8 pt. enlarge() works
# on a ggplot or a patchwork (recursively); apply it to panels before wrap_elements().
if (!exists("FONT_SCALE")) FONT_SCALE <- 1.4
enlarge <- function(p, k = FONT_SCALE) {
  if (inherits(p, "patchwork") && length(p$patches$plots))
    p$patches$plots <- lapply(p$patches$plots, enlarge, k = k)
  if (!inherits(p, "ggplot")) return(p)
  for (i in seq_along(p$layers)) {
    g <- p$layers[[i]]$geom
    if (inherits(g, c("GeomText", "GeomLabel"))) {
      s <- p$layers[[i]]$aes_params$size
      p$layers[[i]]$aes_params$size <- (if (is.null(s)) 3.87 else s) * k
    }
  }
  th <- p$theme
  for (nm in names(th)) {
    e <- th[[nm]]
    if (inherits(e, "element_text") && is.numeric(e$size) && !inherits(e$size, "rel")) th[[nm]]$size <- e$size * k
  }
  p$theme <- th
  if (is.null(p$theme$text)) p <- p + theme(text = element_text(size = 11 * k))
  p
}

read_pd <- function(name) {
  f <- PLOTDATA[[name]]
  if (is.null(f) || !file.exists(f))
    stop("No plot data for '", name, "' (", f, "). Run the analysis script that writes it first.")
  readRDS(f)
}

# ---- the look -----------------------------------------------------------------
theme_fig <- theme_bw(base_size = 10) +
  theme(plot.caption = element_text(size = 7.5, colour = "grey30", hjust = 0), plot.caption.position = "plot")
TAG_THEME <- theme(plot.tag = element_text(face = "bold", size = 20))

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
