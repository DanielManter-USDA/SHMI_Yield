################################################################################
# Assemble_manuscript_figures.R
#
# Builds the multi-panel manuscript figures from the panel PNGs the analysis
# scripts write, adding a / b / c labels:
#   Fig. 3     a  <MAIN>/Fig1_within_sites.png         (SHMI_yield_main.R)
#              b  <MAIN>/Fig2_SHMI_x_N.png
#              c  <MAIN>/Fig3_SHMI_x_age.png
#   Fig. 4     a  <MAIN>/Fig4_forest.png               (SHMI_yield_main.R)
#              b  <STAB>/FigS2_downside.png            (SHMI_yield_stability.R)
#   ED Fig. 8  a  <STAB>/FigS1_sensitivity.png         (SHMI_yield_stability.R)
#              b  <STAB>/FigS3_drought_summary.png
# Layout: Fig. 3, a across the top, b and c side by side below; Fig. 4, a on the
# left, b on the right (top-aligned); ED Fig. 8, a (centred) above b (full width).
# <MAIN> = Outputs/yield_main_<RUN>, <STAB> = Outputs/yield_stability_<RUN>.
# Output: Outputs/manuscript_figures_<RUN>/Fig3.png, Fig4.png, ED_Fig8_stability.png
#
# RUN: after SHMI_yield_main.R and SHMI_yield_stability.R. Set RUN_TAG first to assemble another run's
# figures (e.g. RUN_TAG <- "default_total" for the 7-site run).
# Needs the magick package (install.packages("magick")).
################################################################################
library(magick)

setwd(here::here())                 # the project root (SHMI_Yield.Rproj)

# Folders written by SHMI_yield_main.R and SHMI_yield_stability.R ("_agcros" = with the AgCROS sites)
RUN  <- if (exists("RUN_TAG")) RUN_TAG else "default_total_agcros"
MAIN <- file.path("Outputs", paste0("yield_main_", RUN))
STAB <- file.path("Outputs", paste0("yield_stability_", RUN))
OUTM <- file.path("Outputs", paste0("manuscript_figures_", RUN))
cat("Panels from", MAIN, "and", STAB, "\n")

# every panel must exist before anything is written (a missing panel usually means
# the stability script was run on another site set, or not yet)
PANELS <- c(file.path(MAIN, c("Fig1_within_sites.png", "Fig2_SHMI_x_N.png", "Fig3_SHMI_x_age.png",
                              "Fig4_forest.png")),
            file.path(STAB, c("FigS2_downside.png", "FigS1_sensitivity.png", "FigS3_drought_summary.png")))
miss <- PANELS[!file.exists(PANELS)]
if (length(miss)) stop("Missing panels (run SHMI_yield_main.R and SHMI_yield_stability.R first):\n  ",
                       paste(miss, collapse = "\n  "))
# stability panels made before the latest main run may come from an older version of it
old <- PANELS[grepl(STAB, PANELS, fixed = TRUE) & file.mtime(PANELS) < max(file.mtime(PANELS[1:4]))]
if (length(old)) warning("Stability panels are older than the main-run panels; rerun SHMI_yield_stability.R?\n  ",
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

# ---- Fig. 4 ---------------------------------------------------------------------------
# a forest plot (left, full height); b SHMI benefit in the worst years (right, top-aligned)
pa <- lab(rd(file.path(MAIN, "Fig4_forest.png")), "a")
pb <- lab(rd(file.path(STAB, "FigS2_downside.png")), "b")
left_w <- round(WIDTH * 0.52); right_w <- WIDTH - left_w
left  <- fit_w(pa, left_w)
right <- fit_w(pb, right_w)
h <- max(image_info(left)$height, image_info(right)$height)
left  <- image_extent(left,  paste0(left_w,  "x", h), gravity = "north", color = "white")
right <- image_extent(right, paste0(right_w, "x", h), gravity = "north", color = "white")
save_fig(image_append(c(left, right)), "Fig4.png")

# ---- Extended Data Fig. 8: yield stability ------------------------------------------
# a  SHMI benefit against the year's conditions (county yield deviation)
# b  SHMI benefit at SPEI -1.5 / 0 / +1.5: June-August and May-July, all and rainfed sites
pa <- lab(rd(file.path(STAB, "FigS1_sensitivity.png")), "a")
pb <- lab(rd(file.path(STAB, "FigS3_drought_summary.png")), "b")
top <- fit_w(pa, round(WIDTH * 0.62)); top <- pad_to(top, WIDTH, image_info(top)$height)
save_fig(image_append(c(top, fit_w(pb, WIDTH)), stack = TRUE), "ED_Fig8_stability.png")

cat("All figures in", normalizePath(OUTM), "\n")
