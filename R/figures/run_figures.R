################################################################################
# run_figures.R  (11 Nov 2026)
# Redraws every figure that has a figure script, from the saved plot data only
# (seconds; no analysis is rerun). Run after the analyses, or on its own after
# changing a figure script or fig_common.R.
################################################################################
FIGS <- c("F3_yield.R", "F4_yield_moderators.R")
for (f in FIGS) {
  message("Figure: ", f)
  sys.source(here::here("R", "figures", f), envir = new.env(parent = globalenv()))
}
