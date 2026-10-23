################################################################################
# Run_calibration_scenarios.R
#
# Calibrates SHMI separately for each dataset in datasets.R, each against the
# SH target fitted to that same dataset:
#   all        all units
#   noMX       without the Mexican sites
#   noGZ       without perennial-only systems
#   noMX_noGZ  US annual croplands (PRIMARY_DATASET)
# plus each candidate primary dataset (CANDIDATE_DATASETS) against its
# GAM-adjusted target (sensitivity).
#
# Every run also predicts SHMI for ALL units with its calibrated weights
# (<scenario>_floorXXX_all_unit_predictions.csv), so each dataset's weights
# can be checked on the units it excluded (Summarize_results.R).
#
# All runs: calendar seasons, gamma Diversity (cap 8), Organic Inputs by the
# "either" rule (amendment or animals in a year), Diversity share DIV_FLOOR.
# Iterations are raised by Run_pipeline.R via options().
#
# Inputs  (Outputs/): shmi_inputs_v120.rds, sem_<dataset>_scores.rds,
#                     sem_gam_<candidate>_scores.rds
# Outputs (Outputs/Calibration/): per-scenario calibration files,
#                     scenario_comparison_floorXXX.rds and .xlsx
################################################################################

library(dplyr)
library(tidyr)
library(purrr)
library(tibble)
library(here)

source(here::here("R", "datasets.R"))
BASE      <- here::here("Outputs")
SCRIPT    <- here::here("R", "SHMI_calibration_v3.R")
OUT_DIR   <- file.path(BASE, "Calibration")
DIV_FLOOR <- 0.10
DIVERSITY <- "gamma"        # total diversity, cap 8 (SHMI 1.2.0 decision)
MAXDIV_GRID <- c(6, 7, 8, 9, 10, 12)   # Diversity caps profiled in each run (Section 4g)
B_SENS    <- getOption("shmi.B_sens", 500)   # bootstrap draws for the GAM sensitivity run

# One calibration per dataset (its own SH target holds only its units),
# plus each candidate dataset with its GAM-adjusted target
SCENARIOS <- c(
  lapply(setNames(names(DATASETS), names(DATASETS)), function(k)
    list(scores = file.path(BASE, paste0("sem_", k, "_scores.rds")), subset = function(d) d)),
  setNames(lapply(CANDIDATE_DATASETS, function(k)
             list(scores = file.path(BASE, paste0("sem_gam_", k, "_scores.rds")),
                  subset = function(d) d, B_boot = B_SENS)),
           paste0("gam_", CANDIDATE_DATASETS))
)

tag <- function(nm) sprintf("%s_floor%03d_", nm, round(100 * DIV_FLOOR))

# ---- Run ---------------------------------------------------------------------
for (nm in names(SCENARIOS)) {
  CFG_OVERRIDE <- c(SCENARIOS[[nm]],
                    list(div_floor = DIV_FLOOR, diversity = DIVERSITY, cover_seasons = "calendar",
         maxdiv_grid = MAXDIV_GRID,
                         out_dir = OUT_DIR, out_prefix = tag(nm)))
  local({
    CFG_OVERRIDE <- CFG_OVERRIDE
    source(SCRIPT, local = TRUE)
  })
}

# ---- Compare -----------------------------------------------------------------
res <- imap(SCENARIOS, ~ readRDS(file.path(OUT_DIR, paste0(tag(.y), "results.rds"))))

# Size, selected hill, and out-of-sample performance
size_perf <- imap_dfr(res, function(r, nm) {
  m <- r$metrics %>% group_by(model) %>%
    summarise(R2 = mean(R2), RMSE = mean(RMSE), .groups = "drop")
  tibble(scenario = nm, units = r$n_units, sites = r$n_sites, hill = r$hill,
         R2_equal_weights = m$R2[m$model == "M0"], R2_estimated = m$R2[m$model == "M1"],
         RMSE_equal_weights = m$RMSE[m$model == "M0"], RMSE_estimated = m$RMSE[m$model == "M1"])
})

# Weights with bootstrap intervals: long format (for Excel) and wide text (for the console)
weights_long <- imap_dfr(res, ~ mutate(.x$weight_summary, scenario = .y)) %>%
  select(scenario, weight, estimate, boot_median, lo95, hi95, p_zero)
weights_wide <- weights_long %>%
  mutate(value = sprintf("%.3f [%.3f, %.3f]", estimate, lo95, hi95)) %>%
  select(weight, scenario, value) %>%
  pivot_wider(names_from = scenario, values_from = value)

# Diversity floor and warm-share cost profiles for every scenario
floor_all <- imap_dfr(res, ~ mutate(.x$floor_table, scenario = .y, .before = 1))
warm_all  <- imap_dfr(res, ~ mutate(.x$warm_table,  scenario = .y, .before = 1))
maxdiv_all <- imap_dfr(res, ~ if (is.null(.x$maxdiv_table)) NULL else
                          mutate(.x$maxdiv_table, scenario = .y, .before = 1))
floor_primary <- res[[PRIMARY_DATASET]]$floor_table
warm_primary  <- res[[PRIMARY_DATASET]]$warm_table

cat("\n=== Size and out-of-sample performance ===\n");      print(size_perf, width = Inf)
cat("\n=== Weights, estimate [95% site-bootstrap interval] ===\n"); print(weights_wide, width = Inf)
cat("\n=== Cost of a fixed Diversity share, every scenario ===\n"); print(floor_all, n = Inf, width = Inf)
cat("\n=== Cost of a fixed warm-season share, every scenario (NA = estimated) ===\n"); print(warm_all, n = Inf, width = Inf)

# ---- Export ------------------------------------------------------------------
stem <- file.path(OUT_DIR, sprintf("scenario_comparison_floor%03d", round(100 * DIV_FLOOR)))
saveRDS(list(size_perf = size_perf, weights = weights_long, floor_primary = floor_primary,
             warm_primary = warm_primary, floor_all = floor_all, warm_all = warm_all,
             maxdiv_all = maxdiv_all,
             div_floor = DIV_FLOOR), paste0(stem, ".rds"))
writexl::write_xlsx(list(performance   = size_perf,
                         weights       = weights_long,
                         weights_wide  = weights_wide,
                         floor_all     = floor_all,
                         warm_all      = warm_all,
                         maxdiv_all    = if (nrow(maxdiv_all)) maxdiv_all else tibble(note = "not run")),
                    paste0(stem, ".xlsx"))
cat("\nSaved:", paste0(stem, c(".rds", ".xlsx")), sep = "\n  ")
