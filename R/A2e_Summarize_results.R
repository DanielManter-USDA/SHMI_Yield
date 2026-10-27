################################################################################
# A2e_Summarize_results.R
#
# Collects the results of the SHMI 1.2.0 calibration pipeline into one Excel
# workbook (one sheet per table) and a set of figures (PNG).
#
# Run after:
#   1. A1b_Prepare_SEM.R
#   2. A1c_soil_health_SEM.R
#   3. A2d_SHMI_calibration_final.R
# Sections whose input files are missing are skipped with a message. The
# weight-structure tests and calibration scenarios of the SHMI development build
# (1.2.x) are not part of this release, so their sections are skipped.
#
# Design: each dataset in datasets.R (all, noMX, noGZ, noMX_noGZ) goes through
# adjustment -> SEM -> weight tests -> calibration on its own; the tables
# compare the datasets side by side. PRIMARY_DATASET is reported as the main
# result; each dataset's weights are also checked on the units it excluded.
#
# Output: Outputs/Summary/SHMI_summary_tables_floorXXX.xlsx
#         Outputs/Summary/F*.png
################################################################################

library(dplyr)
library(tidyr)
library(purrr)
library(tibble)
library(ggplot2)
library(here)

OUT   <- function(...) here::here("Outputs", ...)
source(here::here("R", "datasets.R"))
DS    <- names(DATASETS)
FLOOR <- 0.10       # which calibration pass to summarize (Diversity floor used in the run)
SUM   <- OUT("Summary")
dir.create(SUM, showWarnings = FALSE, recursive = TRUE)

TAG      <- sprintf("floor%03d", round(100 * FLOOR))
CAL_DIR  <- OUT("Calibration")
TABLES   <- list()                      # collected for the Excel workbook
theme_set(theme_bw(base_size = 11))

save_fig <- function(p, name, w = 7, h = 4.5) {
  f <- file.path(SUM, paste0(name, ".png"))
  ggsave(f, p, width = w, height = h, dpi = 300)
  cat("  figure:", basename(f), "\n")
}
have <- function(f) {
  ok <- file.exists(f)
  if (!ok) message("Missing: ", f, " - section skipped")
  ok
}


# ==============================================================================
# 1. Data and scope
# ==============================================================================
cat("\n== 1. Data and scope ==\n")
if (have(OUT("scope_v120.rds"))) {
  scope <- readRDS(OUT("scope_v120.rds"))
  TABLES$T1_systems <- scope$systems %>%
    count(system, in_scope, name = "units") %>% arrange(desc(units))
  print(TABLES$T1_systems)
}
if (have(OUT("data_checks_for_review.csv"))) {
  TABLES$T2_data_checks <- read.csv(OUT("data_checks_for_review.csv")) %>%
    count(type, in_scope, name = "records")
}


# ==============================================================================
# 2. Soil-health target (SEM)
# ==============================================================================
cat("\n== 2. SEM target ==\n")
if (have(OUT("sem_summary_all_runs.rds"))) {
  sem_sum <- readRDS(OUT("sem_summary_all_runs.rds"))
  TABLES$T3_sem_runs      <- sem_sum$runs
  TABLES$T4_sem_agreement <- sem_sum$agreement
  print(sem_sum$runs, width = Inf)

  p <- sem_sum$runs %>%
    select(run, cfi.robust, rmsea.robust, srmr) %>%
    pivot_longer(-run, names_to = "index", values_to = "value") %>%
    ggplot(aes(run, value, fill = run)) +
    geom_col(show.legend = FALSE) +
    facet_wrap(~ index, scales = "free_y") +
    labs(x = NULL, y = NULL, title = "SEM fit by run") +
    theme(axis.text.x = element_text(angle = 30, hjust = 1))
  save_fig(p, "F1_sem_fit_by_run", w = 8, h = 3.5)
}
if (exists("sem_sum")) {
  if (!is.null(sem_sum$invariance) && nrow(sem_sum$invariance)) TABLES$T5_invariance_system <- sem_sum$invariance
  if (!is.null(sem_sum$adjustment) && nrow(sem_sum$adjustment)) TABLES$T7_adjustment_summary <- sem_sum$adjustment
}


# ==============================================================================
# 3. Weight structure tests
# ==============================================================================
cat("\n== 3. Weight structure tests ==\n")
if (have(OUT("weight_tests_results.rds"))) {
  wt <- readRDS(OUT("weight_tests_results.rds"))
  TABLES$T8_weight_tests_wald  <- wt$comparison$A
  TABLES$T9_weight_tests_cv    <- wt$comparison$B
  TABLES$T10_mimic_wald        <- wt$comparison$D_wald
  TABLES$T11_direct_effects    <- wt$comparison$D_direct
  print(wt$comparison$A)

  free_w <- imap_dfr(wt$results, ~ mutate(.x$C, scenario = .y))
  TABLES$T12_free_weights <- free_w

  p <- free_w %>%
    filter(!is.na(estimate)) %>%
    mutate(group = ifelse(grepl("^C_|warm", weight), "Seasonal cover", "Organic inputs"),
           weight = sub("^C_|^O_", "", weight)) %>%
    ggplot(aes(weight, estimate, colour = scenario)) +
    geom_pointrange(aes(ymin = lo95, ymax = hi95), position = position_dodge(width = 0.6)) +
    facet_wrap(~ group, scales = "free_x") +
    labs(x = NULL, y = "Weight (95% site-bootstrap interval)",
         title = "Free seasonal and organic-input weights by scenario")
  save_fig(p, "F2_weight_tests_free_weights", w = 9, h = 4.5)

  p <- wt$comparison$A %>%
    pivot_longer(-H0, names_to = "scenario", values_to = "p") %>%
    ggplot(aes(scenario, H0, fill = p < 0.05)) +
    geom_tile(colour = "white") +
    geom_text(aes(label = signif(p, 2)), size = 3.5) +
    scale_fill_manual(values = c(`TRUE` = "#f4a582", `FALSE` = "#d1e5f0"),
                      labels = c(`TRUE` = "rejected (p < 0.05)", `FALSE` = "not rejected"), name = NULL) +
    labs(x = NULL, y = NULL, title = "Weight-structure hypotheses: wild-bootstrap p") +
    theme(axis.text.x = element_text(angle = 30, hjust = 1))
  save_fig(p, "F3_weight_tests_pvalues", w = 7, h = 3.5)
}


# ==============================================================================
# 4. Calibration scenarios
# ==============================================================================
cat("\n== 4. Calibration (", TAG, ") ==\n")
cmp_file <- file.path(CAL_DIR, paste0("scenario_comparison_", TAG, ".rds"))
if (have(cmp_file)) {
  cmp <- readRDS(cmp_file)
  TABLES$T13_calibration_performance <- cmp$size_perf
  TABLES$T14_calibration_weights     <- cmp$weights
  TABLES$T15_diversity_floor_cost    <- if (!is.null(cmp$floor_all)) cmp$floor_all else cmp$floor_primary
  print(cmp$size_perf, width = Inf)

  p <- cmp$size_perf %>%
    select(scenario, `Equal weights` = R2_equal_weights, `Estimated weights` = R2_estimated) %>%
    pivot_longer(-scenario, names_to = "weights", values_to = "R2") %>%
    mutate(scenario = factor(scenario, levels = cmp$size_perf$scenario)) %>%
    ggplot(aes(scenario, R2, fill = weights)) +
    geom_col(position = "dodge") +
    labs(x = NULL, y = "Out-of-sample R\u00b2", fill = NULL,
         title = "Out-of-sample R\u00b2 of SHMI for soil health")
  save_fig(p, "F4_calibration_R2")

  p <- cmp$weights %>%
    filter(weight != "warm_share") %>%
    mutate(weight = recode(weight, w_cover = "Cover", w_diversity = "Diversity",
                           w_invdist = "Inverse disturbance", w_orginput = "Organic inputs")) %>%
    ggplot(aes(weight, estimate, colour = scenario)) +
    geom_pointrange(aes(ymin = lo95, ymax = hi95), position = position_dodge(width = 0.6)) +
    labs(x = NULL, y = "Pillar weight (95% site-bootstrap interval)",
         title = "Calibrated pillar weights by scenario")
  save_fig(p, "F5_calibration_weights")

  p <- cmp$weights %>%
    filter(weight == "warm_share") %>%
    ggplot(aes(scenario, estimate)) +
    geom_pointrange(aes(ymin = lo95, ymax = hi95)) +
    geom_hline(yintercept = 0.5, linetype = 2) +
    labs(x = NULL, y = "Warm-season share of Cover weight",
         title = "Warm-season share by scenario (dashed = equal seasons)")
  save_fig(p, "F6_warm_share")

  # Diversity floor: pillar weights at each fixed share, one panel per scenario
  if (!is.null(cmp$floor_all)) {
    p <- cmp$floor_all %>%
      filter(scenario %in% DS) %>%
      select(scenario, floor, Cover = w_cover, Diversity = w_diversity,
             `Inverse disturbance` = w_invdist, `Organic inputs` = w_orginput) %>%
      pivot_longer(-c(scenario, floor), names_to = "pillar", values_to = "weight") %>%
      mutate(scenario = factor(scenario, DS)) %>%
      ggplot(aes(floor, weight, colour = pillar)) +
      geom_line() + geom_point() +
      facet_wrap(~ scenario) +
      labs(x = "Fixed Diversity share", y = "Pillar weight", colour = NULL,
           title = "Pillar weights when the Diversity share is fixed, by dataset")
    save_fig(p, "F7_diversity_floor", w = 8, h = 5.5)
  }

  # Warm-season share: out-of-sample R2 at each fixed share, one line per
  # scenario; the estimated share is marked with a cross
  if (!is.null(cmp$warm_all)) {
    TABLES$T15b_warm_share_cost <- cmp$warm_all
    wt  <- cmp$warm_all %>% filter(!is.na(warm_fixed))
    est <- cmp$warm_all %>% filter(is.na(warm_fixed))
    p <- ggplot(wt, aes(warm_fixed, R2, colour = scenario)) +
      geom_line() + geom_point(size = 2) +
      geom_point(data = est, aes(x = warm_share, y = R2), shape = 4, size = 4, stroke = 1.2) +
      labs(x = "Fixed warm-season share of Cover (summer + fall)",
           y = "Out-of-sample R\u00b2", colour = NULL,
           title = "Cost of fixing the seasonal split of Cover, by scenario",
           caption = "x = warm share estimated from the data")
    save_fig(p, "F7b_warm_share_cost", w = 8, h = 5)
  }

}


# ==============================================================================
# 4c. Diversity cap profile (Section 4g of each calibration run)
#     A cap is admissible only if some units reach Diversity = 100; the share at
#     100 is compared with Cover and Inverse Disturbance.
# ==============================================================================
cat("\n== 4c. Diversity cap profile ==\n")
if (exists("cmp") && !is.null(cmp$maxdiv_all) && nrow(cmp$maxdiv_all) &&
    "max_div" %in% names(cmp$maxdiv_all)) {
  TABLES$T35_diversity_cap <- cmp$maxdiv_all
  print(TABLES$T35_diversity_cap, n = Inf, width = Inf)
  md <- cmp$maxdiv_all %>% filter(scenario %in% DS) %>% mutate(scenario = factor(scenario, DS))
  refs <- md %>% group_by(scenario) %>%
    summarise(Cover = first(cover_share_at_100), `Inverse disturbance` = first(invdist_share_at_100),
              .groups = "drop") %>%
    pivot_longer(-scenario, names_to = "pillar", values_to = "share")
  p1 <- ggplot(md, aes(max_div, R2, colour = scenario)) + geom_line() + geom_point(size = 2) +
    labs(x = "Diversity cap (effective species scoring 100)", y = "Out-of-sample R\u00b2",
         colour = NULL, title = "Prediction")
  p2 <- ggplot(md, aes(max_div, share_at_100, colour = scenario)) + geom_line() +
    geom_point(aes(shape = spans_0_100), size = 2.5) +
    scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 1), name = "Reaches 100") +
    geom_hline(data = refs, aes(yintercept = share, linetype = pillar), colour = "grey40") +
    scale_y_continuous(labels = scales::percent) +
    labs(x = "Diversity cap (effective species scoring 100)", y = "Units at Diversity = 100",
         colour = NULL, linetype = "Share at 100 for", title = "Spread of the scale")
  if (requireNamespace("patchwork", quietly = TRUE)) {
    p <- patchwork::wrap_plots(p1, p2, ncol = 2) +
      patchwork::plot_annotation(title = "Choosing the Diversity cap")
    save_fig(p, "F18_diversity_cap", w = 11, h = 4.5)
  } else {
    save_fig(p1, "F18a_diversity_cap_R2", w = 6, h = 4.5)
    save_fig(p2, "F18b_diversity_cap_spread", w = 6, h = 4.5)
  }
}

# ==============================================================================
# 5. Ranking treatments within sites (out-of-fold predictions), every dataset;
#    observed vs predicted plots for the primary dataset
# ==============================================================================
cat("\n== 5. Prediction and within-site ranking ==\n")
within_stats <- function(oof) {
  within <- oof %>% group_by(sitecode) %>% filter(n() >= 3) %>%
    mutate(y_dev = y - mean(y), p_dev = M1 - mean(M1)) %>% ungroup()
  site_rho <- within %>% group_by(sitecode) %>% filter(sd(M1) > 0, sd(y) > 0) %>%
    summarise(rho = cor(y, M1, method = "spearman"),
              rho_equal = if (sd(M0) > 0) cor(y, M0, method = "spearman") else NA_real_,
              n = n(), .groups = "drop")
  no_sep <- within %>% group_by(sitecode) %>% summarise(s = sd(M1), .groups = "drop") %>% filter(s == 0)
  n_nz <- sum(site_rho$rho != 0)
  list(summary = tibble(
         sites_3plus = n_distinct(within$sitecode), sites_not_separable = nrow(no_sep),
         pooled_within_r = cor(within$y_dev, within$p_dev),
         median_rho = median(site_rho$rho), share_positive = mean(site_rho$rho > 0),
         sign_test_p = if (n_nz) binom.test(sum(site_rho$rho > 0), n_nz)$p.value else NA_real_,
         median_rho_5plus = median(site_rho$rho[site_rho$n >= 5]),
         median_rho_equal = median(site_rho$rho_equal, na.rm = TRUE),
         share_calibrated_better = mean(site_rho$rho > site_rho$rho_equal, na.rm = TRUE),
         share_calibrated_worse  = mean(site_rho$rho < site_rho$rho_equal, na.rm = TRUE)),
       by_site = site_rho)
}
ws_all <- list()
for (k in DS) {
  f <- file.path(CAL_DIR, paste0(k, "_", TAG, "_oof_predictions.csv"))
  if (file.exists(f)) ws_all[[k]] <- within_stats(read.csv(f))
}
if (length(ws_all)) {
  TABLES$T16_within_site <- imap_dfr(ws_all, ~ mutate(.x$summary, dataset = .y, .before = 1))
  TABLES$T17_within_site_by_site <- imap_dfr(ws_all, ~ mutate(.x$by_site, dataset = .y, .before = 1))
  print(TABLES$T16_within_site, width = Inf)

  p <- TABLES$T17_within_site_by_site %>%
    mutate(dataset = factor(dataset, DS)) %>%
    ggplot(aes(rho)) +
    geom_histogram(breaks = seq(-1, 1, 0.25), fill = "grey60", colour = "white") +
    geom_vline(data = TABLES$T16_within_site %>% mutate(dataset = factor(dataset, DS)),
               aes(xintercept = median_rho), colour = "red") +
    facet_wrap(~ dataset) +
    labs(x = "Spearman rho between SHMI prediction and SH within a site", y = "Sites",
         title = "Ranking of treatments within sites, by dataset (red = median)")
  save_fig(p, "F9_within_site_rho", w = 8, h = 5.5)
}
oof_file <- file.path(CAL_DIR, paste0(PRIMARY_DATASET, "_", TAG, "_oof_predictions.csv"))
if (file.exists(oof_file)) {
  oof <- read.csv(oof_file)
  p <- ggplot(oof, aes(M1, y)) +
    geom_point(aes(colour = system), alpha = 0.6) +
    geom_smooth(method = "lm", formula = y ~ x, se = TRUE, colour = "black") +
    labs(x = "Predicted SH (out-of-fold, estimated weights)", y = "Observed SH (0-100)", colour = NULL,
         title = sprintf("%s: observed vs out-of-fold predicted SH (r = %.2f)",
                         ds_label(PRIMARY_DATASET), cor(oof$y, oof$M1)))
  save_fig(p, "F8_observed_vs_predicted", w = 6.5, h = 5)
}

# ==============================================================================
# 5b. Between-site and within-site relationships (out-of-fold predictions)
#     Between: each site's mean SH against its mean predicted SH (one point per
#     site). Within: each unit's SH and prediction as deviations from its site
#     mean (sites with >= 2 units), pooled across sites.
# ==============================================================================
cat("\n== 5b. Between- and within-site relationships ==\n")
bw <- list()
for (k in DS) {
  f <- file.path(CAL_DIR, paste0(k, "_", TAG, "_oof_predictions.csv"))
  if (!file.exists(f)) next
  o <- read.csv(f)
  site_means <- o %>% group_by(sitecode) %>%
    summarise(y = mean(y), pred = mean(M1), n = n(), .groups = "drop")
  within <- o %>% group_by(sitecode) %>% filter(n() >= 2) %>%
    mutate(y = y - mean(y), pred = M1 - mean(M1)) %>% ungroup()
  bw[[k]] <- bind_rows(
    site_means %>% transmute(dataset = k, level = "Between sites (site means)", y, pred),
    within %>% transmute(dataset = k, level = "Within sites (deviations from site mean)", y, pred))
}
if (length(bw)) {
  bwd <- bind_rows(bw) %>% mutate(dataset = factor(dataset, DS))
  TABLES$T34_between_within <- bwd %>%
    group_by(dataset, level) %>%
    summarise(points = n(), r = cor(y, pred),
              slope = unname(coef(lm(y ~ pred))[2]),
              r2 = cor(y, pred)^2, .groups = "drop")
  print(TABLES$T34_between_within, width = Inf)
  lab <- TABLES$T34_between_within %>%
    mutate(label = sprintf("r = %.2f\nslope = %.2f\nn = %d", r, slope, points))
  p <- ggplot(bwd, aes(pred, y)) +
    geom_hline(yintercept = 0, colour = "grey85") + geom_vline(xintercept = 0, colour = "grey85") +
    geom_point(alpha = 0.45, size = 1.4) +
    geom_smooth(method = "lm", formula = y ~ x, colour = "black", se = TRUE) +
    geom_text(data = lab, aes(x = -Inf, y = Inf, label = label),
              hjust = -0.1, vjust = 1.2, size = 3.2, inherit.aes = FALSE) +
    facet_grid(dataset ~ level, scales = "free") +
    labs(x = "Predicted soil health (out-of-fold SHMI prediction)", y = "Measured soil health (SH)",
         title = "How SHMI tracks soil health between sites and within sites")
  save_fig(p, "F17_between_within_sites", w = 9, h = 7)
}

# ==============================================================================
# 6. Transfer: each dataset's calibrated weights applied to all units, split
#    into the units it was calibrated on (in-sample) and the groups it excluded.
#    Correlations and within-site rankings are scale-free, so SH comes from the
#    all-units target for every group.
# ==============================================================================
cat("\n== 6. Transfer to units outside each calibration set ==\n")
transfer_stats <- function(d) {
  ws <- d %>% group_by(sitecode) %>% filter(n() >= 3, sd(pred) > 0, sd(SH) > 0) %>%
    summarise(rho = cor(SH, pred, method = "spearman"), .groups = "drop")
  n_pos <- sum(ws$rho > 0); n_nz <- sum(ws$rho != 0)
  tibble(units = nrow(d), sites = n_distinct(d$sitecode),
         pearson = cor(d$SH, d$pred), spearman = cor(d$SH, d$pred, method = "spearman"),
         sites_ranked = nrow(ws),
         within_median_rho = if (nrow(ws)) median(ws$rho) else NA_real_,
         within_share_pos  = if (nrow(ws)) mean(ws$rho > 0) else NA_real_,
         sign_test_p = if (n_nz > 0) binom.test(n_pos, n_nz)$p.value else NA_real_)
}
tr <- list(); pa_cand <- list()
for (k in DS) {
  f <- file.path(CAL_DIR, paste0(k, "_", TAG, "_all_unit_predictions.csv"))
  if (!file.exists(f)) next
  pa <- read.csv(f) %>%
    mutate(SH = SH_all, system = coalesce(system, "unknown"),
           group = case_when(in_calibration ~ "calibration set (in-sample)",
                             system == "perennial" & region == "US" ~ "excluded: US perennial",
                             system == "perennial" ~ "excluded: Mexico perennial",
                             region == "Mexico" ~ "excluded: Mexico cropland",
                             TRUE ~ "excluded: other")) %>%
    filter(!is.na(SH), !is.na(pred))
  tr[[k]] <- pa %>% group_by(group) %>% group_modify(~ transfer_stats(.x)) %>% ungroup() %>%
    mutate(dataset = k, .before = 1)
  if (k %in% CANDIDATE_DATASETS) pa_cand[[k]] <- pa
}
if (length(tr)) {
  TABLES$T18_transfer <- bind_rows(tr)
  print(TABLES$T18_transfer, n = Inf, width = Inf)
}
for (k in names(pa_cand)) {
  p <- ggplot(pa_cand[[k]], aes(pred, SH)) +
    geom_point(alpha = 0.5) +
    geom_smooth(method = "lm", formula = y ~ x, colour = "black") +
    facet_wrap(~ group, scales = "free") +
    labs(x = paste0("SHMI prediction (weights calibrated on: ", ds_label(k), ")"),
         y = "Soil health (SH, all-units target)",
         title = paste0("Weights from '", k, "' applied inside and outside its calibration set"))
  save_fig(p, paste0("F10_transfer_", k), w = 9, h = 6)
}

# ==============================================================================
# 6b. Pooling test (criterion 3 of the decision rule): on the units of each
#     narrower dataset, do the all-data weights predict SH worse out-of-sample
#     than that dataset's own weights? Both sets of predictions are out-of-fold;
#     both are scored against the same SH (the all-units target) and against the
#     narrower dataset's own target, with scale-free measures.
# ==============================================================================
cat("\n== 6b. Pooling test: all-data weights vs each dataset's own weights ==\n")
oof_of <- function(k) {
  f <- file.path(CAL_DIR, paste0(k, "_", TAG, "_oof_predictions.csv"))
  if (file.exists(f)) read.csv(f) else NULL
}
oof_all <- oof_of("all")
sh_all  <- if (file.exists(OUT("sem_all_scores.rds")))
  readRDS(OUT("sem_all_scores.rds")) %>% select(MGT_combo, SH_all = SH) else NULL
if (!is.null(oof_all) && !is.null(sh_all)) {
  score_preds <- function(d, pred, y) {
    ws <- d %>% group_by(sitecode) %>% filter(n() >= 3, sd(.data[[pred]]) > 0, sd(.data[[y]]) > 0) %>%
      summarise(rho = cor(.data[[y]], .data[[pred]], method = "spearman"), .groups = "drop")
    list(stats = c(r2 = cor(d[[y]], d[[pred]])^2,
                   spearman = cor(d[[y]], d[[pred]], method = "spearman"),
                   within_median_rho = median(ws$rho), within_share_pos = mean(ws$rho > 0)),
         sites = ws)
  }
  pool <- list(); pool_sites <- list()
  for (k in setdiff(DS, "all")) {
    own <- oof_of(k); if (is.null(own)) next
    d <- own %>% select(MGT_combo, sitecode, y_own = y, pred_own = M1) %>%
      inner_join(oof_all %>% select(MGT_combo, pred_all = M1), by = "MGT_combo") %>%
      inner_join(sh_all, by = "MGT_combo")
    for (tg in c("SH_all", "y_own")) {
      a <- score_preds(d, "pred_all", tg); o <- score_preds(d, "pred_own", tg)
      pool[[paste(k, tg)]] <- tibble(
        dataset = k, units = nrow(d), sites = n_distinct(d$sitecode),
        scored_against = ifelse(tg == "SH_all", "all-units SH", "dataset's own SH"),
        measure = names(a$stats), all_data_weights = a$stats, own_weights = o$stats,
        difference_own_minus_all = o$stats - a$stats)
      if (tg == "SH_all") {
        j <- inner_join(a$sites, o$sites, by = "sitecode", suffix = c("_all", "_own"))
        pool_sites[[k]] <- tibble(dataset = k, sites_compared = nrow(j),
                                  share_own_better = mean(j$rho_own > j$rho_all),
                                  share_all_better = mean(j$rho_all > j$rho_own),
                                  share_tied = mean(j$rho_all == j$rho_own))
      }
    }
  }
  if (length(pool)) {
    TABLES$T23_pooling_test <- bind_rows(pool)
    TABLES$T24_pooling_by_site <- bind_rows(pool_sites)
    print(TABLES$T23_pooling_test, n = Inf, width = Inf)
    print(TABLES$T24_pooling_by_site, width = Inf)

    p <- TABLES$T23_pooling_test %>%
      filter(scored_against == "all-units SH", measure %in% c("r2", "within_median_rho")) %>%
      pivot_longer(c(all_data_weights, own_weights), names_to = "weights", values_to = "value") %>%
      mutate(weights = recode(weights, all_data_weights = "all-data weights",
                              own_weights = "dataset's own weights"),
             measure = recode(measure, r2 = "Out-of-sample R\u00b2 (all-units SH)",
                              within_median_rho = "Median within-site rho")) %>%
      ggplot(aes(dataset, value, fill = weights)) +
      geom_col(position = "dodge") +
      facet_wrap(~ measure, scales = "free_y") +
      labs(x = "Units of dataset", y = NULL, fill = NULL,
           title = "Pooling test: do all-data weights predict each subset as well as its own weights?")
    save_fig(p, "F12_pooling_test", w = 9, h = 4.5)
  }
}

# ==============================================================================
# 7. Information available for each SHMI component, by dataset
# ==============================================================================
cat("\n== 7. Information per component ==\n")
# Scores as in A2d_SHMI_calibration_final.R (shmi_components, STIR matched to T-DISC);
# Diversity's rotation part = 4-year-span Diversity minus per-calendar-year Diversity.
if (have(OUT("shmi_inputs_v120.rds")) && have(OUT("sem_all_scores.rds")) &&
    have(OUT("site_climate_monthly.rds"))) tryCatch({
  suppressPackageStartupMessages(library(SHMI))
  inp <- readRDS(OUT("shmi_inputs_v120.rds"))
  units <- readRDS(OUT("sem_all_scores.rds")) %>% select(MGT_combo, sitecode, system)
  irr <- read.csv(here::here("Inputs", "napeshm", "treatments.csv"), stringsAsFactors = FALSE) %>%
    transmute(MGT_combo = paste0("NAPESHM_", sitecode, "_", treatmentid),
              irrigated = !is.na(eu_irrigation) & !(trimws(tolower(eu_irrigation)) %in% c("none", "", "n_a")))
  climate <- inp$mgt %>%
    transmute(MGT_combo, sitecode = sub("^NAPESHM_(.*)_[^_]+$", "\\1", MGT_combo)) %>%
    left_join(readRDS(OUT("site_climate_monthly.rds")), by = "sitecode") %>%
    left_join(irr, by = "MGT_combo") %>%
    mutate(irrigated = coalesce(irrigated, FALSE)) %>% select(-sitecode)
  comp_all <- shmi_components(inp, climate, dist_meth = "STIR") %>% select(MGT_combo, InvDist)
  div4 <- compute_diversity(inp$crop, inp$rot_bounds, window = 4) %>% select(MGT_combo, D4 = Diversity)
  div1 <- compute_diversity(inp$crop, inp$rot_bounds, window = 1) %>% select(MGT_combo, D1 = Diversity)
  org  <- compute_orginput_components(inp$rot_bounds, inp$amend, inp$animal) %>%
    select(MGT_combo, O_amend = p_amend, O_animal = p_animal)
  comp_all <- comp_all %>% left_join(div4, by = "MGT_combo") %>% left_join(div1, by = "MGT_combo") %>%
    left_join(org, by = "MGT_combo") %>%
    mutate(across(c(D4, D1, O_amend, O_animal), ~ coalesce(.x, 0)),
           D_rotation = pmax(D4 - D1, 0))
  info <- function(x, site, name) tibble(component = name, units_nonzero = sum(x > 0),
    sites_nonzero = n_distinct(site[x > 0]),
    sites_varying = sum(tapply(x, site, function(v) diff(range(v)) > 0)))
  info_tabs <- list()
  for (k in DS) {
    ids <- ds_filter(units, k)$MGT_combo
    comp <- comp_all %>% filter(MGT_combo %in% ids) %>%
      mutate(site = sub("^NAPESHM_([A-Za-z0-9]+)_.*$", "\\1", MGT_combo))
    info_tabs[[k]] <- bind_rows(
      info(100 - comp$InvDist,  comp$site, "Inverse disturbance (any tillage)"),
      info(comp$D4,             comp$site, "Diversity: more than one species in a 4-year span"),
      info(comp$D_rotation,     comp$site, "Diversity: added by rotation (4-year span minus per year)"),
      info(comp$D1,             comp$site, "Diversity: more than one species within a year"),
      info(comp$O_amend,        comp$site, "Organic inputs: amendments"),
      info(comp$O_animal,       comp$site, "Organic inputs: animals")) %>% mutate(dataset = k, .before = 1)
  }
  TABLES$T19_component_information <- bind_rows(info_tabs)
  print(TABLES$T19_component_information, n = Inf)
}, error = function(e) message("Section 7 skipped: ", conditionMessage(e)))


# ==============================================================================
# 8. Cover seasons x Diversity definition tests (if run)
# ==============================================================================
cat("\n== 8. Cover x Diversity tests ==\n")
TEST_DIR <- file.path(CAL_DIR, "tests"); FU_DIR <- file.path(CAL_DIR, "followup")
if (file.exists(file.path(TEST_DIR, "cover_diversity_comparison.xlsx"))) {
  TABLES$T20_cover_diversity_perf <- readxl::read_excel(file.path(TEST_DIR, "cover_diversity_comparison.xlsx"),
                                                        sheet = "performance")
}
if (file.exists(file.path(FU_DIR, "followup_comparison.xlsx"))) {
  TABLES$T21_followup_perf <- readxl::read_excel(file.path(FU_DIR, "followup_comparison.xlsx"),
                                                 sheet = "performance")
}

# Paired comparison of two test cases across CV repeats (same folds). The SE
# reflects fold variation only, not sampling of sites (see Nadeau & Bengio 2003).
paired <- function(tag_a, tag_b, dir_a, dir_b = dir_a) {
  fa <- file.path(dir_a, paste0(tag_a, "results.rds")); fb <- file.path(dir_b, paste0(tag_b, "results.rds"))
  if (!file.exists(fa) || !file.exists(fb)) return(NULL)
  a <- readRDS(fa)$metrics %>% filter(model == "M1")
  b <- readRDS(fb)$metrics %>% filter(model == "M1")
  d <- inner_join(a, b, by = "rep_id", suffix = c("_a", "_b")) %>% mutate(dR2 = R2_a - R2_b)
  tibble(case_a = tag_a, case_b = tag_b, mean_dR2 = mean(d$dR2),
         se = sd(d$dR2) / sqrt(nrow(d)), share_a_better = mean(d$dR2 > 0))
}
pairs_list <- list(
  list("calendar_combined_floor010_", "calendar_alpha_floor010_", TEST_DIR),
  list("calendar_combined_floor010_", "calendar_beta_floor010_",  TEST_DIR),
  list("calendar_combined_floor010_", "calendar_gamma_floor010_", TEST_DIR),
  list("calendar_combined_floor010_", "climate_combined_floor010_", TEST_DIR),
  list("noMX_calendar_combined_t05_floor010_", "noMX_calendar_alpha_t05_floor010_", FU_DIR),
  list("noMX_calendar_combined_t05_floor010_", "noMX_calendar_beta_t05_floor010_",  FU_DIR),
  list("noMX_calendar_combined_t05_floor010_", "noMX_calendar_gamma_t05_floor010_", FU_DIR))
pc <- bind_rows(lapply(pairs_list, function(x) paired(x[[1]], x[[2]], x[[3]])))
if (nrow(pc)) { TABLES$T22_paired_comparisons <- pc; print(pc, width = Inf) }

if (!is.null(TABLES$T20_cover_diversity_perf)) {
  p <- TABLES$T20_cover_diversity_perf %>%
    mutate(floor = factor(div_floor)) %>%
    ggplot(aes(diversity, R2_estimated, fill = floor)) +
    geom_col(position = "dodge") +
    facet_wrap(~ cover_seasons) +
    labs(x = "Diversity definition", y = "Out-of-sample R\u00b2", fill = "Diversity floor",
         title = "Cover seasons x Diversity definition (in scope, with Mexico)")
  save_fig(p, "F11_cover_diversity_tests", w = 8, h = 4.5)
}


# ==============================================================================
# 8b. Spatial vs temporal Diversity (Run_diversity_split_tests.R, if run)
#     gamma (total diversity) vs the split w x Alpha + (1 - w) x Beta_norm, per
#     dataset; profile of out-of-sample R2 over the split w_alpha
# ==============================================================================
cat("\n== 8b. Diversity: total (gamma) vs spatial/temporal split ==\n")
SPLIT_RUNS <- c(combined_norm = "diversity_split", gamma_split = "diversity_gamma_split")
split_perf <- list(); split_prof <- list()
for (st in names(SPLIT_RUNS)) {
  f <- file.path(CAL_DIR, SPLIT_RUNS[[st]], "diversity_split_comparison.xlsx")
  if (!file.exists(f)) next
  split_perf[[st]] <- readxl::read_excel(f, sheet = "performance") %>% mutate(split_type = st, .before = 1)
  split_prof[[st]] <- readxl::read_excel(f, sheet = "split_profile") %>% mutate(split_type = st, .before = 1)
  if (st == "combined_norm") TABLES$T27_record_cycles <- readxl::read_excel(f, sheet = "record_cycles")
}
if (length(split_perf)) {
  TABLES$T25_diversity_split_perf    <- bind_rows(split_perf) %>% distinct()
  TABLES$T26_diversity_split_profile <- bind_rows(split_prof)
  print(TABLES$T25_diversity_split_perf, n = Inf, width = Inf)

  gam <- TABLES$T25_diversity_split_perf %>%
    filter(diversity == "gamma", div_floor > 0) %>%
    distinct(dataset, .keep_all = TRUE) %>%
    transmute(dataset = factor(dataset, DS), R2_gamma = R2_estimated)
  prof <- TABLES$T26_diversity_split_profile %>%
    mutate(dataset = factor(dataset, DS),
           split_type = recode(split_type, gamma_split = "gamma's parts (raw beta)",
                               combined_norm = "Alpha + Beta_norm"))
  if (nrow(prof)) {
    p <- ggplot(prof, aes(w_alpha, R2, colour = split_type)) +
      geom_hline(data = gam, aes(yintercept = R2_gamma), linetype = 2, colour = "grey40") +
      geom_line() + geom_point(size = 2) +
      facet_wrap(~ dataset, scales = "free_y") +
      scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1),
                         labels = c("0\nrotation\nonly", "0.25", "0.5", "0.75", "1\nwithin-year\nonly")) +
      labs(x = "Share of Diversity on within-year (spatial) diversity",
           y = "Out-of-sample R\u00b2", colour = "Split",
           title = "Splitting Diversity into spatial and temporal parts, by dataset",
           caption = "Dashed line = total diversity (gamma); for gamma's parts, w = 0.5 is gamma. Diversity share fixed at 0.10.") +
      theme(legend.position = "bottom")
    save_fig(p, "F14_diversity_split_profile", w = 8, h = 6)
  }
}


# ==============================================================================
# 8c. Climate seasons with a moisture criterion (Run_climate_moisture_tests.R,
#     if run): calendar vs climate (>= 5 C) vs climate (>= 5 C and not dry)
# ==============================================================================
cat("\n== 8c. Cover seasons: calendar vs climate (temperature, moisture) ==\n")
CM_FILE <- file.path(CAL_DIR, "climate_moist", "climate_moisture_comparison.xlsx")
if (file.exists(CM_FILE)) {
  TABLES$T28_climate_moisture <- readxl::read_excel(CM_FILE, sheet = "comparison")
  print(TABLES$T28_climate_moisture, n = Inf, width = Inf)
  p <- TABLES$T28_climate_moisture %>%
    mutate(dataset = factor(dataset, DS)) %>%
    ggplot(aes(dataset, R2_estimated, fill = seasons)) +
    geom_col(position = "dodge") +
    labs(x = NULL, y = "Out-of-sample R\u00b2", fill = "Cover seasons",
         title = "Calendar vs climate-based Cover seasons, by dataset")
  save_fig(p, "F15_climate_moisture", w = 8, h = 4.5)
}


# ==============================================================================
# 8d. SHMI v1.1.0 vs the current v1.2.0 candidate (Compare_v110_v120.R, if run)
# ==============================================================================
cat("\n== 8d. v1.1.0 vs v1.2.0 candidate ==\n")
VC_FILE <- file.path(SUM, "v110_vs_v120_summary.csv")
if (file.exists(VC_FILE)) {
  TABLES$T29_v110_vs_v120 <- read.csv(VC_FILE)
  vs <- file.path(SUM, "v110_vs_v120_scores.csv")
  if (file.exists(vs)) {
    sc <- read.csv(vs)
    TABLES$T30_v110_vs_v120_by_group <- sc %>%
      group_by(index, group) %>%
      summarise(units = n(), pearson = cor(v110, v120), mean_change = mean(v120 - v110),
                units_moved_gt10 = sum(abs(v120 - v110) > 10), .groups = "drop")
    # the units that move most, per index (for checking rule edge cases)
    TABLES$T31_v110_vs_v120_largest_changes <- sc %>%
      mutate(change = v120 - v110) %>%
      group_by(index) %>% slice_max(abs(change), n = 10, with_ties = FALSE) %>% ungroup() %>%
      arrange(index, desc(abs(change)))
  }
  print(TABLES$T29_v110_vs_v120)
}

# ==============================================================================
# 8e. Organic Inputs folded into Diversity (Run_org_merge_tests.R, if run)
# ==============================================================================
cat("\n== 8e. Organic Inputs as effective species in Diversity ==\n")
OM_FILE <- file.path(CAL_DIR, "org_merge", "org_merge_comparison.xlsx")
if (file.exists(OM_FILE)) {
  TABLES$T32_org_merge <- readxl::read_excel(OM_FILE, sheet = "org_merge_profile")
  # Requirement: like Cover and InvDist, Diversity must span 0-100, i.e. some
  # units must reach the cap. A cap no unit reaches is too high.
  if ("div_max" %in% names(TABLES$T32_org_merge)) {
    TABLES$T32_org_merge <- TABLES$T32_org_merge %>%
      mutate(spans_0_100 = div_max >= 99.95)
  }
  print(TABLES$T32_org_merge, n = Inf, width = Inf)
  om <- TABLES$T32_org_merge %>% mutate(dataset = factor(dataset, DS))
  if (!"hill" %in% names(om)) om$hill <- 1
  if (!"mode" %in% names(om)) om$mode <- ifelse(is.na(om$x_amend), "reference", "days")
  four  <- om %>% filter(mode == "reference") %>% select(dataset, diversity_share, hill, R2_four = R2)
  three <- om %>% filter(mode != "reference") %>%
    left_join(four, by = c("dataset", "diversity_share", "hill"))
  for (md in unique(three$mode)) {
    dd <- three %>% filter(mode == md) %>%
      mutate(animal = factor(x_animal), hill_lab = paste0("Hill q = ", hill))
    unit_lab <- if (md == "days") "days a plant species grows" else "Diversity points added"
    p <- ggplot(dd, aes(x_amend, R2, colour = animal)) +
      geom_hline(aes(yintercept = R2_four), colour = "grey30", linetype = 2) +
      geom_line() + geom_point(aes(shape = spans_0_100), size = 2) +
      scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 1), name = "Reaches 100") +
      facet_grid(dataset + diversity_share ~ hill_lab, scales = "free_y") +
      labs(x = paste0("Each amendment year: ", unit_lab),
           colour = paste0("Each animal year:\n", unit_lab), y = "Out-of-sample R\u00b2",
           title = paste0("Organic inputs folded into Diversity (", md, " version)"),
           caption = "Dashed line: four pillars with Organic Inputs separate (same Hill order).")
    save_fig(p, paste0("F16_org_merge_", md), w = 11, h = 8)
  }
  TABLES$T33_org_merge_spread <- om %>%
    filter(diversity_share == first(diversity_share)) %>%
    select(dataset, hill, mode, x_amend, x_animal, div_max, div_p99, share_at_cap, share_zero, spans_0_100)
}

# ==============================================================================
# 9. Export all tables
# ==============================================================================
xlsx <- file.path(SUM, paste0("SHMI_summary_tables_", TAG, ".xlsx"))
writexl::write_xlsx(TABLES, xlsx)
cat("\nTables:", xlsx, "\n")
cat("Sheets:", paste(names(TABLES), collapse = ", "), "\n")
