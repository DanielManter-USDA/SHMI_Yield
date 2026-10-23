################################################################################
# Run_within_calibration_scenarios.R  (v4)
#
# Reruns SHMI_calibration_plus.R for a grid of
#   what the weights are fitted to: site averages only (between sites), all
#                  variation (the original calibration), or within-site contrasts
#                  only (add 0.5 to LAMBDAS for a half-within version)
#   fixed_div      Diversity's fixed share: 0.10, 0.15, 0.20, 0.25
# and collects weights and cross-validated performance in one table and figure.
# Each scenario writes its usual outputs to Outputs/Calib_scenarios/L<lambda>_D<div>/
# (B_D<div>/ for the between-site fits).
# Scenarios already run are skipped. Grazing and amendment gaps stay at 0.25 and 0.60.
#
# Needs the version of SHMI_calibration_plus.R with within_lambda and between_only. Set
# CALIB_SCRIPT if it is not in R/.
################################################################################
suppressPackageStartupMessages({ library(dplyr); library(tidyr); library(ggplot2); library(here) })

CALIB_SCRIPT <- if (exists("CALIB_SCRIPT")) CALIB_SCRIPT else here::here("R", "SHMI_calibration_plus.R")
OUT_ROOT     <- here::here("Outputs", "Calib_scenarios")
LAMBDAS      <- c(0, 1)
BETWEEN      <- TRUE           # also fit to site averages only
DIVS         <- c(0.10, 0.15, 0.20, 0.25)
stopifnot(file.exists(CALIB_SCRIPT))
if (!any(grepl("between_only", readLines(CALIB_SCRIPT))))
  stop("CALIB_SCRIPT is an older SHMI_calibration_plus.R without between_only: ", CALIB_SCRIPT)
dir.create(OUT_ROOT, recursive = TRUE, showWarnings = FALSE)

scen <- expand.grid(within_lambda = LAMBDAS, fixed_div = DIVS)
scen$between_only <- FALSE
if (isTRUE(BETWEEN)) scen <- rbind(data.frame(within_lambda = 0, fixed_div = DIVS, between_only = TRUE), scen)
scen$dir <- file.path(OUT_ROOT, ifelse(scen$between_only, sprintf("B_D%.2f", scen$fixed_div),
                                       sprintf("L%.1f_D%.2f", scen$within_lambda, scen$fixed_div)))
scen$fit <- ifelse(scen$between_only, "between", as.character(scen$within_lambda))

for (i in seq_len(nrow(scen))) {
  if (file.exists(file.path(scen$dir[i], "final_results.rds"))) {
    message("Skipping ", basename(scen$dir[i]), " (already run; delete the folder to rerun)"); next
  }
  cat(sprintf("\n========== Scenario %d of %d: %s, fixed_div = %.2f ==========\n",
              i, nrow(scen), basename(scen$dir[i]), scen$fixed_div[i]))
  CFG_OVERRIDE <- list(within_lambda = scen$within_lambda[i], between_only = scen$between_only[i],
                       fixed_div = scen$fixed_div[i], out_dir = scen$dir[i])
  source(CALIB_SCRIPT, local = new.env())
  rm(CFG_OVERRIDE)
}

# ------------------------------------------------------------------------------
# Collect
# ------------------------------------------------------------------------------
collect <- function(i) {
  res <- readRDS(file.path(scen$dir[i], "final_results.rds"))
  w <- res$weights %>% filter(weight %in% c("cover", "diversity", "disturbance", "growing_share"))
  gw <- function(nm, col = "final") w[[col]][w$weight == nm]
  cv <- res$us_cv %>% filter(score == "SHMI")
  mx <- res$test %>% filter(score == "SHMI", grepl("^Mexico", set))
  pick <- function(nm) if (nm %in% names(cv)) cv[[nm]] else NA_real_
  mxv <- function(set, col) { v <- mx[[col]][mx$set == set]; if (length(v)) v[1] else NA_real_ }
  tibble(fit = scen$fit[i], within_lambda = scen$within_lambda[i], fixed_div = scen$fixed_div[i],
         cover = gw("cover"), cover_lo = gw("cover", "lo95"), cover_hi = gw("cover", "hi95"),
         diversity = gw("diversity"),
         disturbance = gw("disturbance"), disturbance_lo = gw("disturbance", "lo95"), disturbance_hi = gw("disturbance", "hi95"),
         growing_share = gw("growing_share"), growing_share_lo = gw("growing_share", "lo95"),
         growing_share_hi = gw("growing_share", "hi95"),
         R2 = pick("R2_median"), R2_lo = pick("R2_lo95"), R2_hi = pick("R2_hi95"),
         between_r = pick("between_r_median"), within_r = pick("within_r_median"),
         sites_ranked_correctly = pick("within_share_positive_median"),
         within_median_rho = pick("within_median_rho_median"),
         mexico_r_US_SEM = mxv("Mexico: SH from US SEM (independent test)", "r"),
         mexico_r_all_SEM = mxv("Mexico: SH from all-data SEM (check)", "r"),
         mexico_within_rho = mxv("Mexico: SH from US SEM (independent test)", "within_median_rho"))
}
FIT_LEVELS <- c(between = "Between sites", `0` = "All variation", `0.5` = "Half within", `1` = "Within sites")
tab <- bind_rows(lapply(seq_len(nrow(scen)), collect)) %>%
  mutate(fit_label = factor(FIT_LEVELS[fit], FIT_LEVELS[FIT_LEVELS %in% FIT_LEVELS[unique(fit)]])) %>%
  arrange(fit_label, fixed_div)
write.csv(tab, file.path(OUT_ROOT, "scenario_summary.csv"), row.names = FALSE)

cat("\nWeights (median [95% range over 1,000 fits]) and cross-validated performance by scenario\n")
ci <- function(v, lo, hi) sprintf("%.2f [%.2f, %.2f]", v, lo, hi)
print(tab %>% transmute(fit = fit_label, fixed_div,
                        Cover = ci(cover, cover_lo, cover_hi), InvDist = ci(disturbance, disturbance_lo, disturbance_hi),
                        `Season share` = ci(growing_share, growing_share_lo, growing_share_hi),
                        R2 = round(R2, 3), between_r = round(between_r, 3), within_r = round(within_r, 3),
                        ranked = round(sites_ranked_correctly, 3), mexico_r = round(mexico_r_US_SEM, 3)) %>%
        as.data.frame(), row.names = FALSE)

# ------------------------------------------------------------------------------
# Figure: weights and performance by scenario
# ------------------------------------------------------------------------------
div_col <- setNames(c("#c6dbef", "#6baed6", "#2171b5", "#08306b")[seq_along(DIVS)], sprintf("%.2f", DIVS))
pw <- tab %>% select(fit_label, fixed_div, Cover = cover, Diversity = diversity, InvDist = disturbance) %>%
  pivot_longer(c(Cover, Diversity, InvDist), names_to = "pillar", values_to = "weight") %>%
  mutate(fit = fit_label, div = sprintf("Diversity fixed %.2f", fixed_div)) %>%
  ggplot(aes(fit, weight, fill = pillar)) + geom_col(width = 0.7) +
  geom_text(aes(label = sprintf("%.2f", weight)), position = position_stack(vjust = 0.5), size = 3, colour = "white") +
  facet_wrap(~ div, nrow = 1) +
  scale_fill_manual(values = c(Cover = "#2e7d32", Diversity = "#f2a541", InvDist = "#1b6ca8"), name = NULL) +
  labs(x = "Weights fitted to", y = "Weight", title = "(a) Pillar weights") +
  theme_bw(base_size = 10) + theme(axis.text.x = element_text(angle = 30, hjust = 1))
pp <- tab %>% select(fit_label, fixed_div, `Overall R²` = R2, `Between-site r` = between_r,
                     `Within-site r` = within_r, `Sites ranked correctly` = sites_ranked_correctly,
                     `Mexico r` = mexico_r_US_SEM) %>%
  pivot_longer(-c(fit_label, fixed_div), names_to = "metric") %>%
  mutate(fit = fit_label, metric = factor(metric, unique(metric)),
         div = sprintf("%.2f", fixed_div)) %>%
  ggplot(aes(fit, value, colour = div, group = div)) +
  geom_line() + geom_point(size = 2) + facet_wrap(~ metric, scales = "free_y", nrow = 1) +
  scale_colour_manual(values = div_col, name = "Diversity fixed at") +
  labs(x = "Weights fitted to", y = NULL, title = "(b) Cross-validated performance (median of 100 repeats)") +
  theme_bw(base_size = 10) + theme(axis.text.x = element_text(angle = 30, hjust = 1), legend.position = "bottom")
p <- if (requireNamespace("patchwork", quietly = TRUE)) patchwork::wrap_plots(pw, pp, ncol = 1, heights = c(1, 1)) else pp
ggsave(file.path(OUT_ROOT, "Fig_calibration_scenarios.png"), p, width = 13, height = 8, dpi = 300)
cat("\nSaved scenario_summary.csv and Fig_calibration_scenarios.png in", OUT_ROOT, "\n")
