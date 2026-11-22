################################################################################
# run_figures.R  (11 Nov 2026; four-figure layout 21 Nov 2026)
# Redraws every figure from the saved plot data and panels only (seconds; no analysis is
# rerun). yield_panels.R and moderator_panels.R build panels for F3 and F4.
################################################################################
FIGS <- c("F1_soil_health_model.R", "F2_system_score.R", "F3_yield.R", "F4_yield_moderators.R")
for (f in FIGS) {
  message("Figure: ", f)
  tryCatch(sys.source(here::here("R", "figures", f), envir = new.env(parent = globalenv())),
           error = function(e) message("  not drawn: ", conditionMessage(e)))
}
