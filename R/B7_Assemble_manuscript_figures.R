################################################################################
# B7_Assemble_manuscript_figures.R
#
# Builds the multi-panel manuscript figures from the panel PNGs the analysis
# scripts write, adding a / b / c labels:
#   Fig. 3     a  <MAIN>/Fig1_within_sites.png         (B2_SHMI_yield_main.R)
#              b  Outputs/yield_record_length/Fig_record_length.png  (B2c: record vs current season)
#              c  <MAIN>/Fig2_SHMI_x_N.png
#              d  <MAIN>/Fig3_SHMI_x_age.png
#   Fig. 4     a  <STAB>/FigS1_sensitivity.png         (B3_SHMI_yield_stability.R)
#              b  <STAB>/FigS3_drought_summary.png
#              c  <STAB>/FigS2_downside.png
#   ED Fig. 4  a  <MAIN>/Fig_SHMI_ranges_slopes.png    (B2_SHMI_yield_main.R, section 7b: range and own slope)
#              b  Outputs/yield_parts/Fig_yield_parts.png   (B2d: SHMI against its parts for yield)
# Layout: Fig. 3, a and b side by side on top, c and d below; Fig. 4, a (centred)
# above b and c (full width); ED Fig. 4, a above b.
# (Before 28 Oct 2026, Fig. 3 had no record-length panel (it was Extended Data Fig. 6)
# and Extended Data Fig. 4b was the forest plot of each experiment's slope, still written
# by B2 as EDFig7b_forest.png.) Without B2c or B2d output, Fig. 3 is drawn without b
# and ED Fig. 4 without b, with a message.
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
REC   <- file.path("Outputs", "yield_record_length", "Fig_record_length.png")   # B2c (optional here)
PARTS <- file.path("Outputs", "yield_parts", "Fig_yield_parts.png")            # B2d (optional here)
PANELS <- c(file.path(MAIN, c("Fig1_within_sites.png", "Fig2_SHMI_x_N.png", "Fig3_SHMI_x_age.png",
                              "Fig_SHMI_ranges.png")),
            file.path(STAB, c("FigS2_downside.png", "FigS1_sensitivity.png", "FigS3_drought_summary.png")))
miss <- PANELS[!file.exists(PANELS)]
if (length(miss)) stop("Missing panels (run B2_SHMI_yield_main.R and B3_SHMI_yield_stability.R first):\n  ",
                       paste(miss, collapse = "\n  "))
# stability panels made before the latest main run may come from an older version of it
old <- PANELS[grepl(STAB, PANELS, fixed = TRUE) & file.mtime(PANELS) < max(file.mtime(PANELS[1:4]))]
if (length(old)) warning("Stability panels are older than the main-run panels; rerun B3_SHMI_yield_stability.R?\n  ",
                         paste(old, collapse = "\n  "))

dir.create(OUTM, recursive = TRUE, showWarnings = FALSE)
WIDTH <- 3600                                   # final width in pixels (~12 in at 300 dpi)

rd <- function(f) image_trim(image_read(f))
# panel letter in a white strip of height `top` above the panel (letter size in pixels)
lab <- function(img, letter, size = 90, top = 60) {
  img <- image_border(img, "white", "40x0")
  w <- image_info(img)$width
  img <- image_append(c(image_blank(w, top, "white"), img, image_blank(w, 60, "white")), stack = TRUE)
  image_annotate(img, letter, size = size, weight = 700, gravity = "northwest", location = "+10+0")
}
fit_w <- function(img, w) image_scale(img, as.character(w))
fit_h <- function(img, h) image_scale(img, paste0("x", h))
pad_to <- function(img, w, h) image_extent(img, paste0(w, "x", h), gravity = "center", color = "white")
save_fig <- function(img, name) {
  image_write(image_background(img, "white"), file.path(OUTM, name), density = 300)
  cat("Saved", name, "\n")
}

# ---- Fig. 3 ---------------------------------------------------------------------------
# a the within-site relationship and b record vs current season side by side (same
# height), c SHMI x N and d SHMI x age below; smaller letters, set above each panel so
# they clear the axis titles
# panels are scaled first and lettered after, so every letter has the same size
L3 <- function(img, letter) lab(img, letter, size = 65, top = 85)
half <- WIDTH / 2
row2 <- function(files, letters) {         # two panels side by side, each half the width
  x <- mapply(function(f, l) L3(fit_w(rd(f), half - 80), l), files, letters, SIMPLIFY = FALSE)
  h <- max(sapply(x, function(i) image_info(i)$height))
  image_append(image_join(lapply(x, pad_to, w = half, h = h)))
}
if (file.exists(REC)) {
  # top row: a and b at a common height, together the full width
  ra <- rd(file.path(MAIN, "Fig1_within_sites.png")); rb <- rd(REC)
  hh <- max(image_info(ra)$height, image_info(rb)$height)
  ra <- fit_h(ra, hh); rb <- fit_h(rb, hh)
  k  <- (WIDTH - 160) / (image_info(ra)$width + image_info(rb)$width)
  top <- image_append(c(L3(fit_w(ra, round(k * image_info(ra)$width)), "a"),
                        L3(fit_w(rb, round(k * image_info(rb)$width)), "b")))
  bottom <- row2(file.path(MAIN, c("Fig2_SHMI_x_N.png", "Fig3_SHMI_x_age.png")), c("c", "d"))
} else {
  message("Fig. 3 drawn without the record-length panel: ", REC, " not found (run B2c_SHMI_record_length.R)")
  top <- L3(fit_w(rd(file.path(MAIN, "Fig1_within_sites.png")), round(WIDTH * 0.62)), "a")
  top <- pad_to(top, WIDTH, image_info(top)$height)
  bottom <- row2(file.path(MAIN, c("Fig2_SHMI_x_N.png", "Fig3_SHMI_x_age.png")), c("b", "c"))
}
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

# ---- Extended Data Fig. 4: the experiments, and SHMI against its parts for yield -----
# a each experiment's stretch of the SHMI scale; b SHMI against each part alone, equal
# weights and the current season (B2d)
RS <- file.path(MAIN, "Fig_SHMI_ranges_slopes.png")    # ranges beside each experiment's slope (B2, since 30 Oct)
pa <- lab(rd(if (file.exists(RS)) RS else file.path(MAIN, "Fig_SHMI_ranges.png")), "a")
w_a <- round(WIDTH * 0.95)
top <- pad_to(fit_w(pa, w_a), WIDTH, image_info(fit_w(pa, w_a))$height)
if (file.exists(PARTS)) {
  pb <- lab(rd(PARTS), "b")
  bottom <- pad_to(fit_w(pb, w_a), WIDTH, image_info(fit_w(pb, w_a))$height)   # same width as a, so similar type size
  save_fig(image_append(c(top, bottom), stack = TRUE), "ED_Fig4_sites.png")
} else {
  message("Extended Data Fig. 4 drawn without panel b: ", PARTS, " not found (run B2d_SHMI_yield_parts.R)")
  save_fig(top, "ED_Fig4_sites.png")
}

cat("All figures in", normalizePath(OUTM), "\n")
