################################################################################
# B7_Assemble_manuscript_figures.R
#
# Builds the multi-panel manuscript figures from the panel PNGs the analysis
# scripts write, adding a / b / c labels:
#   Fig. 3     a  <MAIN>/Fig1_within_sites.png         (B2_SHMI_yield_main.R)
#              b  <MAIN>/Fig2_SHMI_x_N.png
#              c  <MAIN>/Fig3_SHMI_x_age.png
#   Fig. 4     a  <STAB>/FigS1_sensitivity.png         (B3_SHMI_yield_stability.R)
#              b  <STAB>/FigS3_drought_summary.png
#              c  <STAB>/FigS2_downside.png
#   ED Fig. 4  a  <MAIN>/Fig_SHMI_ranges.png           (B2_SHMI_yield_main.R, section 7b)
#              b  <MAIN>/EDFig7b_forest.png
# Layout: Fig. 3, a across the top, b and c side by side below; Fig. 4, a (centred)
# above b and c (full width); ED Fig. 4, a above b (sites in the same order).
# <MAIN> = Outputs/yield_main_<RUN>, <STAB> = Outputs/yield_stability_<RUN>.
# Output: Outputs/manuscript_figures_<RUN>/Fig3.png, Fig4.png, ED_Fig4_sites.png
# (Before 25 Oct 2026 Fig. 4 was the forest plot plus the worst years, and the year-
# condition panels were an Extended Data figure.)
#
# RUN: after B2_SHMI_yield_main.R and B3_SHMI_yield_stability.R. Set RUN_TAG first to assemble another run's
# figures (e.g. RUN_TAG <- "default_total" for the 7-site run).
# Needs the magick package (install.packages("magick")).
################################################################################
library(magick)

setwd(here::here())                 # the project root (SHMI_Yield.Rproj)

# Folders written by B2_SHMI_yield_main.R and B3_SHMI_yield_stability.R ("_agcros" = with the AgCROS sites)
RUN  <- if (exists("RUN_TAG")) RUN_TAG else "default_total_agcros"
MAIN <- file.path("Outputs", paste0("yield_main_", RUN))
STAB <- file.path("Outputs", paste0("yield_stability_", RUN))
OUTM <- file.path("Outputs", paste0("manuscript_figures_", RUN))
cat("Panels from", MAIN, "and", STAB, "\n")

# every panel must exist before anything is written (a missing panel usually means
# the stability script was run on another site set, or not yet)
PANELS <- c(file.path(MAIN, c("Fig1_within_sites.png", "Fig2_SHMI_x_N.png", "Fig3_SHMI_x_age.png",
                              "Fig_SHMI_ranges.png", "EDFig7b_forest.png")),
            file.path(STAB, c("FigS2_downside.png", "FigS1_sensitivity.png", "FigS3_drought_summary.png")))
miss <- PANELS[!file.exists(PANELS)]
if (length(miss)) stop("Missing panels (run B2_SHMI_yield_main.R and B3_SHMI_yield_stability.R first):\n  ",
                       paste(miss, collapse = "\n  "))
# stability panels made before the latest main run may come from an older version of it
old <- PANELS[grepl(STAB, PANELS, fixed = TRUE) & file.mtime(PANELS) < max(file.mtime(PANELS[1:5]))]
if (length(old)) warning("Stability panels are older than the main-run panels; rerun B3_SHMI_yield_stability.R?\n  ",
                         paste(old, collapse = "\n  "))

dir.create(OUTM, recursive = TRUE, showWarnings = FALSE)
WIDTH <- 3600                                   # final width in pixels (~12 in at 300 dpi)

rd <- function(f) image_trim(image_read(f))
lab <- function(img, letter) {
  img <- image_border(img, "white", "40x60")
  image_annotate(img, letter, size = 90, weight = 700, gravity = "northwest", location = "+10+0")
}
fit_w <- function(img, w) image_scale(img, as.character(w))
fit_h <- function(img, h) image_scale(img, paste0("x", h))
pad_to <- function(img, w, h) image_extent(img, paste0(w, "x", h), gravity = "center", color = "white")
save_fig <- function(img, name) {
  image_write(image_background(img, "white"), file.path(OUTM, name), density = 300)
  cat("Saved", name, "\n")
}

# ---- Fig. 3 ---------------------------------------------------------------------------
pa <- lab(rd(file.path(MAIN, "Fig1_within_sites.png")), "a")
pb <- lab(rd(file.path(MAIN, "Fig2_SHMI_x_N.png")), "b")
pc <- lab(rd(file.path(MAIN, "Fig3_SHMI_x_age.png")), "c")
top <- fit_w(pa, round(WIDTH * 0.62)); top <- pad_to(top, WIDTH, image_info(top)$height)
half <- WIDTH / 2
bc <- list(fit_w(pb, half), fit_w(pc, half))
h <- max(sapply(bc, function(x) image_info(x)$height))
bottom <- image_append(image_join(lapply(bc, pad_to, w = half, h = h)))
save_fig(image_append(c(top, bottom), stack = TRUE), "Fig3.png")

# ---- Fig. 4: the SHMI yield benefit in poor, dry and good years ---------------------
# a  the slope against the year's conditions (county yield deviation)
# b  the slope at SPEI -1.5 / 0 / +1.5: June-August and May-July, all and rainfed sites
# c  worst years: relative to the site-year mean and to the plot's own mean
pa <- lab(rd(file.path(STAB, "FigS1_sensitivity.png")), "a")
pb <- lab(rd(file.path(STAB, "FigS3_drought_summary.png")), "b")
pc <- lab(rd(file.path(STAB, "FigS2_downside.png")), "c")
top <- fit_w(pa, round(WIDTH * 0.62)); top <- pad_to(top, WIDTH, image_info(top)$height)
save_fig(image_append(c(top, fit_w(pb, WIDTH), fit_w(pc, WIDTH)), stack = TRUE), "Fig4.png")

# ---- Extended Data Fig. 4: the experiments ------------------------------------------
# a  each experiment's stretch of the SHMI scale; b  each experiment's slope (forest plot)
pa <- lab(rd(file.path(MAIN, "Fig_SHMI_ranges.png")), "a")
pb <- lab(rd(file.path(MAIN, "EDFig7b_forest.png")), "b")
w_a <- round(WIDTH * 0.85); w_b <- round(WIDTH * 0.85)
top    <- pad_to(fit_w(pa, w_a), WIDTH, image_info(fit_w(pa, w_a))$height)
bottom <- pad_to(fit_w(pb, w_b), WIDTH, image_info(fit_w(pb, w_b))$height)
save_fig(image_append(c(top, bottom), stack = TRUE), "ED_Fig4_sites.png")

cat("All figures in", normalizePath(OUTM), "\n")
