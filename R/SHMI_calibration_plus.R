################################################################################
# SHMI_calibration_plus.R
#
# Final calibration of SHMI (three pillars plus the organic amendment bonus)
# and SHMI_base (the three pillars alone). Same machinery as SHMI_calibration_final.R (the four-pillar
# calibration, kept as the comparison in Outputs/Final).
#
# Structure
#   Cover/Roots  living cover in the active and dormant season, by site climate;
#                active-season share estimated
#   Diversity    plant species per 4-year rotation window (sliding, averaged),
#                100 x (richness - 1) / 7, plus grazing credit:
#                Diversity + graze_gap x p_animal x (100 - Diversity)
#                (p_animal = share of years with grazing animals);
#                share of the pillar weights fixed at fixed_div
#   InvDist      inverse disturbance (STIR matched to T-DISC, as in the final calibration)
#
#   SHMI_base = w_C Cover/Roots + w_D Diversity + w_I InvDist, weights sum to 1
#   SHMI      = SHMI_base + amend_gap x p_amend x (100 - SHMI_base)
#            (p_amend = share of years with a manure or compost amendment);
#            no score exceeds 100
#
# Grazing and amendment gap shares were chosen in SHMI_orginput_redesign_test.R
# (gap versions: R2 flat for amend_gap 0.5-0.8 and graze_gap 0.25-0.5).
#
# 1. Calibration on US units: 10-fold leave-sites-out CV repeated 100 times
#    (1,000 fits), sites weighted equally. The core weights are fitted by NNLS
#    with the amendment term in place (SHMI), so Cover/Roots does not absorb
#    the amendment effect; the same weights give SHMI_base. Both scores are mapped to
#    soil health by their own site-weighted line on the training sites.
#    Final weights = medians across fits (Cover/Roots and InvDist renormalized to
#    sum to 1 - fixed_div); intervals = 2.5-97.5 percentiles.
# 2. Independent test on the Mexican units, for both scores.
#
# within_lambda sets what the weights are fitted to: 0 = all variation (the
# original calibration), 1 = differences between treatments within a site only
# (each column and soil health minus lambda x its site mean), values between
# mix the two. between_only = TRUE fits the weights to site averages only (each
# column and soil health replaced by its site mean). The line mapping SHMI to soil
# health always uses raw soil health.
#
# Outputs: Outputs/Final_plus/ final_weights.csv, final_results.xlsx,
#          shmi_scores.csv, us_oof_predictions.csv, final_results.rds, F1-F4
################################################################################

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(purrr); library(tibble); library(ggplot2); library(here); library(SHMI)
})
BASE <- here::here("Outputs")
CFG <- list(
  inputs      = file.path(BASE, "shmi_inputs_v120.rds"),
  us_scores   = file.path(BASE, "sem_noMX_scores.rds"),
  mx_scores   = file.path(BASE, "sem_noMX_projected_MX.rds"),
  all_scores  = file.path(BASE, "sem_all_scores.rds"),
  climate     = file.path(BASE, "site_climate_monthly.rds"),
  treatments  = here::here("Inputs", "CSV Files", "treatments.csv"),
  out_dir     = file.path(BASE, "Final_plus"),
  div_window  = 4,           # Diversity window (years)
  fixed_div   = 0.15,        # Diversity's share of the three pillar weights
  graze_gap   = 0.25,        # share of (100 - Diversity) closed by grazing every year
  amend_gap   = 0.60,        # share of (100 - SHMI_base) closed by an amendment every year (SHMI)
  within_lambda = 0,         # 0 = fit weights to all variation; 1 = within-site contrasts only
  between_only  = FALSE,     # TRUE = fit weights to site averages only (overrides within_lambda)
  fixed_weights = NULL,      # e.g. c(cover = 1/3, diversity = 1/3, disturbance = 1/3, growing_share = 0.5):
                             # use these weights in every fold instead of fitting (set fixed_div to the
                             # same diversity value); only the line to soil health is fitted
  test_mexico = TRUE,
  K           = 10,
  R           = 100,
  seed        = 2026         # same folds as SHMI_calibration_final.R
)
if (exists("CFG_OVERRIDE")) CFG <- utils::modifyList(CFG, CFG_OVERRIDE)
dir.create(CFG$out_dir, recursive = TRUE, showWarnings = FALSE)
OUTF <- function(x) file.path(CFG$out_dir, x)

# ==============================================================================
# Components for every unit
# ==============================================================================
inp <- readRDS(CFG$inputs)
irr <- read.csv(CFG$treatments, stringsAsFactors = FALSE) %>%
  transmute(MGT_combo = paste0("NAPESHM_", sitecode, "_", treatmentid),
            irrigated = !is.na(eu_irrigation) & !(trimws(tolower(eu_irrigation)) %in% c("none", "", "n_a")))
climate <- inp$mgt %>%
  transmute(MGT_combo, sitecode = sub("^NAPESHM_(.*)_[^_]+$", "\\1", MGT_combo)) %>%
  left_join(readRDS(CFG$climate), by = "sitecode") %>%
  left_join(irr, by = "MGT_combo") %>%
  mutate(irrigated = coalesce(irrigated, FALSE)) %>% select(-sitecode)

options(SHMI.diversity_window = CFG$div_window)
stopifnot("window" %in% names(formals(SHMI::compute_diversity)))
org_parts <- getExportedValue("SHMI", "compute_orginput_components")
oc <- org_parts(inp$rot_bounds, inp$amend, inp$animal) %>% select(MGT_combo, p_amend, p_animal)

# Plant-only Diversity: the grazing credit is added below with CFG$graze_gap.
# Since the SHMI/SHMI_base package update, shmi_components() returns Diversity
# with the grazing credit already added and plant-only Diversity as
# Diversity_plants; older versions return plant-only Diversity as Diversity.
X <- shmi_components(inp, climate, dist_meth = "STIR")      # NAPESHM records STIR
if (!"Diversity_plants" %in% names(X)) X$Diversity_plants <- X$Diversity
X <- X %>%
  transmute(MGT_combo, NonGrow = Cover_nongrowing, Grow = Cover_growing,
            Diversity_plants, Richness, InvDist) %>%
  left_join(oc, by = "MGT_combo") %>%
  mutate(across(c(p_amend, p_animal), ~ coalesce(.x, 0)),
         Diversity = Diversity_plants + CFG$graze_gap * p_animal * (100 - Diversity_plants))
OTHER <- c("NonGrow", "Grow", "InvDist")

us <- readRDS(CFG$us_scores) %>% transmute(MGT_combo, sitecode = as.character(sitecode), system, SH) %>%
  inner_join(X, by = "MGT_combo") %>% filter(!is.na(SH))
if (isTRUE(CFG$test_mexico)) {
  mx <- readRDS(CFG$mx_scores) %>% transmute(MGT_combo, sitecode = as.character(sitecode), system, SH) %>%
    inner_join(X, by = "MGT_combo") %>% filter(!is.na(SH))
  mx_all <- readRDS(CFG$all_scores) %>% transmute(MGT_combo, SH_all = SH) %>% filter(MGT_combo %in% mx$MGT_combo)
  mx <- mx %>% left_join(mx_all, by = "MGT_combo")
} else {
  mx <- us[0, ] %>% mutate(SH_all = numeric(0))
}
cat("US units:", nrow(us), "at", length(unique(us$sitecode)), "sites | Mexican units:", nrow(mx),
    "at", length(unique(mx$sitecode)), "sites\n")
cat("US units with amendments:", sum(us$p_amend > 0), "| with grazing:", sum(us$p_animal > 0), "\n")

# ==============================================================================
# Helpers
# ==============================================================================
site_w <- function(s) 1 / as.numeric(table(s)[s])
wmean  <- function(x, w) sum(w * x) / sum(w)
# SHMI = base (1 - c p) + 100 c p with base = sum_o w_o P_o + s Div. With
# k = s / (1 - s), NNLS on (P_o + k Div)(1 - c p) + (1 + k) 100 c p fits b_o >= 0;
# B = (1 + k) sum(b_o) is the soil-health points per SHMI point.
fit_plus <- function(D, w) {
  if (!is.null(CFG$fixed_weights)) {                 # fixed weights: nothing fitted, only the SH line
    fw <- CFG$fixed_weights[c("cover", "diversity", "disturbance", "growing_share")]
    return(setNames(as.numeric(fw), names(fw)))
  }
  k <- CFG$fixed_div / (1 - CFG$fixed_div); cp <- CFG$amend_gap * D$p_amend
  Xm <- (as.matrix(D[, OTHER]) + k * D$Diversity) * (1 - cp) + (1 + k) * 100 * cp
  yv <- D$SH
  if (isTRUE(CFG$between_only)) {                   # keep only site means (between-site fit)
    Xm <- apply(Xm, 2, function(v) stats::ave(v, D$sitecode))
    yv <- stats::ave(yv, D$sitecode)
  } else if (CFG$within_lambda > 0) {               # remove lambda x site means
    Xm <- Xm - CFG$within_lambda * apply(Xm, 2, function(v) stats::ave(v, D$sitecode))
    yv <- yv - CFG$within_lambda * stats::ave(yv, D$sitecode)
  }
  xm <- colSums(Xm * w) / sum(w); ym <- wmean(yv, w)
  f <- nnls::nnls(sweep(Xm, 2, xm) * sqrt(w), (yv - ym) * sqrt(w))
  b <- setNames(f$x, OTHER); B <- (1 + k) * sum(b)
  cover <- b[["NonGrow"]] + b[["Grow"]]
  wts <- c(cover = cover / B, diversity = CFG$fixed_div, disturbance = b[["InvDist"]] / B,
           growing_share = if (cover > 0) b[["Grow"]] / cover else NA_real_)
  wts
}
shmi_core <- function(D, w) {                       # SHMI_base
  cover <- (1 - w[["growing_share"]]) * D$NonGrow + w[["growing_share"]] * D$Grow
  w[["cover"]] * cover + w[["diversity"]] * D$Diversity + w[["disturbance"]] * D$InvDist
}
shmi_plus <- function(D, w) { s <- shmi_core(D, w); s + CFG$amend_gap * D$p_amend * (100 - s) }   # SHMI
line_fit <- function(s, y, w) {                      # site-weighted line y ~ s
  ms <- wmean(s, w); my <- wmean(y, w)
  b <- sum(w * (s - ms) * (y - my)) / sum(w * (s - ms)^2); c(a = my - b * ms, b = b)
}
within_site <- function(obs, pred, site, min_n = 3) {
  tibble(obs, pred, site) %>% group_by(site) %>% filter(n() >= min_n, sd(pred) > 0, sd(obs) > 0) %>%
    summarise(rho = cor(obs, pred, method = "spearman"), .groups = "drop")
}
r2 <- function(obs, pred) 1 - sum((obs - pred)^2) / sum((obs - mean(obs))^2)
ccc <- function(x, y) { n <- length(x); sx <- sqrt(sum((x - mean(x))^2) / n); sy <- sqrt(sum((y - mean(y))^2) / n)
                        2 * cor(x, y) * sx * sy / (sx^2 + sy^2 + (mean(x) - mean(y))^2) }
perf_of <- function(y, pred, site, r, score) {
  ws <- within_site(y, pred, site)
  so <- tapply(y, site, mean); sp <- tapply(pred, site, mean)
  dy <- y - ave(y, site); dp <- pred - ave(pred, site)
  tibble(repeat_id = r, score = score, R2 = r2(y, pred), RMSE = sqrt(mean((y - pred)^2)),
         between_r = cor(so, sp), within_r = cor(dy, dp),
         within_median_rho = median(ws$rho), within_share_positive = mean(ws$rho > 0))
}

# ==============================================================================
# 1. Calibration on US units: K-fold leave-sites-out, repeated R times
# ==============================================================================
sites <- unique(us$sitecode); y <- us$SH; wus <- site_w(us$sitecode)
fits <- list(); perf <- list(); oof_plus <- numeric(nrow(us)); oof_core <- numeric(nrow(us))
for (r in seq_len(CFG$R)) {
  set.seed(CFG$seed + r)
  fold_of_site <- sample(rep_len(seq_len(CFG$K), length(sites))); fold <- fold_of_site[match(us$sitecode, sites)]
  p_plus <- numeric(nrow(us)); p_core <- numeric(nrow(us))
  for (k in seq_len(CFG$K)) {
    tr <- fold != k; w <- site_w(us$sitecode[tr]); Dtr <- us[tr, ]; Dte <- us[!tr, ]
    wt <- fit_plus(Dtr, w)
    lp <- line_fit(shmi_plus(Dtr, wt), Dtr$SH, w); p_plus[!tr] <- lp[["a"]] + lp[["b"]] * shmi_plus(Dte, wt)
    lc <- line_fit(shmi_core(Dtr, wt), Dtr$SH, w); p_core[!tr] <- lc[["a"]] + lc[["b"]] * shmi_core(Dte, wt)
    fits[[length(fits) + 1]] <- c(repeat_id = r, fold = k, wt)
  }
  oof_plus <- oof_plus + p_plus; oof_core <- oof_core + p_core
  perf[[length(perf) + 1]] <- perf_of(y, p_plus, us$sitecode, r, "SHMI")
  perf[[length(perf) + 1]] <- perf_of(y, p_core, us$sitecode, r, "SHMI_base")
  if (r %% 10 == 0) cat("repeat", r, "of", CFG$R, "\n")
}
fit_tab  <- as_tibble(do.call(rbind, fits))
perf_tab <- bind_rows(perf)

WN <- c("cover", "diversity", "disturbance", "growing_share")
weight_summary <- map_dfr(WN, function(wn) {
  v <- fit_tab[[wn]]
  tibble(weight = wn, median = median(v, na.rm = TRUE), lo95 = quantile(v, 0.025, na.rm = TRUE),
         hi95 = quantile(v, 0.975, na.rm = TRUE), mean = mean(v, na.rm = TRUE),
         share_zero = mean(v <= 1e-10, na.rm = TRUE))
})
w_final <- setNames(weight_summary$median, weight_summary$weight)
oth <- c("cover", "disturbance")
w_final[oth] <- w_final[oth] / sum(w_final[oth]) * (1 - CFG$fixed_div)
w_final[["diversity"]] <- CFG$fixed_div
weight_summary <- weight_summary %>% mutate(final = unname(w_final[weight])) %>%
  bind_rows(tibble(weight = c("graze_gap", "amend_gap"), final = c(CFG$graze_gap, CFG$amend_gap)))
cat("\nFinal weights (medians of", nrow(fit_tab), "fits):\n"); print(weight_summary, width = Inf)
us_cv <- perf_tab %>% group_by(score) %>%
  summarise(across(-repeat_id, list(median = median, lo95 = ~ quantile(.x, 0.025), hi95 = ~ quantile(.x, 0.975))),
            .groups = "drop")

oof <- us %>% transmute(MGT_combo, sitecode, system, SH, pred_oof = oof_plus / CFG$R, pred_oof_base = oof_core / CFG$R) %>%
  group_by(sitecode) %>% mutate(err = SH - pred_oof, site_err = mean(err), within_err = err - site_err) %>% ungroup()
error_split <- tibble(
  score = "SHMI", RMSE_total = sqrt(mean(oof$err^2)),
  share_between_sites = sum(oof$site_err^2) / sum(oof$err^2),
  share_within_sites  = sum(oof$within_err^2) / sum(oof$err^2),
  RMSE_within_only = sqrt(mean(oof$within_err^2)), RMSE_between_only = sqrt(mean(oof$site_err^2)))
w_full <- fit_plus(us, wus)

# ==============================================================================
# 2. Final scores; SH-scale calibration on US; independent Mexico test
# ==============================================================================
add_scores <- function(D) D %>% mutate(
  Cover = (1 - w_final[["growing_share"]]) * NonGrow + w_final[["growing_share"]] * Grow,
  SHMI_base = shmi_core(., w_final), SHMI = shmi_plus(., w_final), Amend_bonus = SHMI - SHMI_base)
all_units <- add_scores(X); us <- add_scores(us); mx <- add_scores(mx)
cal_core <- stats::lm(SH ~ SHMI_base, data = us, weights = wus)
cal_plus <- stats::lm(SH ~ SHMI, data = us, weights = wus)
us$pred <- predict(cal_plus, us); us$pred_base <- predict(cal_core, us)
mx$pred <- predict(cal_plus, mx); mx$pred_base <- predict(cal_core, mx)

test_stats <- function(D, obs_col, score_col, pred_col, label) {
  obs <- D[[obs_col]]; ok <- !is.na(obs); s <- D[[score_col]][ok]; p <- D[[pred_col]][ok]; o <- obs[ok]
  ws <- within_site(o, s, D$sitecode[ok])
  sm <- tibble(o, p, site = D$sitecode[ok]) %>% group_by(site) %>% summarise(o = mean(o), p = mean(p), .groups = "drop")
  tibble(set = label, score = score_col,
         units = sum(ok), sites = length(unique(D$sitecode[ok])),
         r = cor(o, s), R2 = r2(o, p), RMSE = sqrt(mean((o - p)^2)), bias = mean(p - o),
         score_bias = mean(s - o), CCC_score = ccc(s, o), between_site_r = cor(sm$o, sm$p),
         sites_ranked = nrow(ws), within_median_rho = median(ws$rho), within_share_positive = mean(ws$rho > 0),
         sign_test_p = if (nrow(ws)) binom.test(sum(ws$rho > 0), sum(ws$rho != 0))$p.value else NA_real_)
}
test_table <- bind_rows(lapply(list(c("SHMI", "pred"), c("SHMI_base", "pred_base")), function(sc) bind_rows(
  test_stats(us, "SH", sc[1], sc[2], "US units, final weights (in-sample)"),
  if (nrow(mx)) test_stats(mx, "SH", sc[1], sc[2], "Mexico: SH from US SEM (independent test)"),
  if (nrow(mx)) test_stats(mx, "SH_all", sc[1], sc[2], "Mexico: SH from all-data SEM (check)"))))
cat("\nPerformance:\n"); print(test_table, width = Inf)
cat("\nUS cross-validated performance (median and 95% range over", CFG$R, "repeats):\n"); print(us_cv, width = Inf)
cat("\nUS units at 100: SHMI_base", sum(us$SHMI_base >= 99.95), "| SHMI", sum(us$SHMI >= 99.95),
    "| max SHMI", round(max(all_units$SHMI), 2), "\n")

# ==============================================================================
# Save
# ==============================================================================
calib <- bind_rows(tibble(score = "SHMI_base", term = names(coef(cal_core)), estimate = unname(coef(cal_core))),
                   tibble(score = "SHMI", term = names(coef(cal_plus)), estimate = unname(coef(cal_plus))))
write.csv(weight_summary, OUTF("final_weights.csv"), row.names = FALSE)
write.csv(all_units %>% select(MGT_combo, NonGrow, Grow, Cover, Diversity_plants, Diversity, Richness, InvDist,
                               p_amend, p_animal, SHMI_base, Amend_bonus, SHMI),
          OUTF("shmi_scores.csv"), row.names = FALSE)
write.csv(oof, OUTF("us_oof_predictions.csv"), row.names = FALSE)
writexl::write_xlsx(list(final_weights = weight_summary,
                         full_fit_weights = tibble(weight = names(w_full), estimate = unname(w_full)),
                         us_cv_performance = us_cv, error_split = error_split, test = test_table,
                         calibration = calib, fold_weights = fit_tab, repeats = perf_tab,
                         mexico_units = mx %>% select(MGT_combo, sitecode, SH, SH_all, SHMI_base, SHMI, pred_base, pred)),
                    OUTF("final_results.xlsx"))
saveRDS(list(config = CFG, weights = weight_summary, w_final = w_final, w_full = w_full, fits = fit_tab,
             us_cv = us_cv, error_split = error_split, oof = oof, repeats = perf_tab, test = test_table,
             calibration = list(SHMI = cal_plus, SHMI_base = cal_core), us = us, mx = mx, scores = all_units),
        OUTF("final_results.rds"))

# ==============================================================================
# Figures
# ==============================================================================
theme_set(theme_bw(base_size = 11))
lab <- c(cover = "Cover/Roots", diversity = "Diversity\n(fixed)", disturbance = "Inverse\nDisturbance")
wl <- fit_tab %>% select(cover, diversity, disturbance) %>% pivot_longer(everything(), names_to = "weight") %>%
  mutate(weight = factor(lab[weight], lab))
wf <- weight_summary %>% filter(weight %in% names(lab)) %>% mutate(weight = factor(lab[weight], lab))
p1a <- ggplot(wl, aes(weight, value)) +
  geom_violin(data = filter(wl, weight != lab[["diversity"]]), fill = "grey85", colour = NA, scale = "width") +
  geom_boxplot(data = filter(wl, weight != lab[["diversity"]]), width = 0.15, outlier.size = 0.4) +
  geom_errorbar(data = wf, aes(y = final, ymin = lo95, ymax = hi95), width = 0.12, colour = "#d1495b") +
  geom_point(data = wf, aes(y = final), colour = "#d1495b", size = 3) +
  geom_text(data = wf, aes(y = final, label = sprintf("%.2f", final)), nudge_x = 0.28, colour = "#d1495b", size = 3.6) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(x = NULL, y = "Weight across 1,000 fits", title = "(a) Pillar weights (sum to 1)")
gs <- weight_summary %>% filter(weight == "growing_share")
p1b <- ggplot(fit_tab, aes("Active-season share\nin Cover/Roots", growing_share)) +
  geom_violin(fill = "#cfe8c8", colour = NA) + geom_boxplot(width = 0.15, outlier.size = 0.4) +
  geom_errorbar(data = gs, aes(y = final, ymin = lo95, ymax = hi95), width = 0.12, colour = "#d1495b") +
  geom_point(data = gs, aes(y = final), colour = "#d1495b", size = 3) +
  geom_text(data = gs, aes(y = final, label = sprintf("%.2f", final)), nudge_x = 0.28, colour = "#d1495b", size = 3.6) +
  scale_y_continuous(limits = c(0.5, 1)) +
  labs(x = NULL, y = "Share", title = "(b) Within Cover/Roots")
p1 <- if (requireNamespace("patchwork", quietly = TRUE)) {
  patchwork::wrap_plots(p1a, p1b, widths = c(3, 1.2)) +
    patchwork::plot_annotation(title = sprintf(
      "SHMI weights: medians of 1,000 fits and 95%% ranges (red). Grazing gap %.2f in Diversity; amendment gap %.2f in SHMI",
      CFG$graze_gap, CFG$amend_gap))
} else p1a
ggsave(OUTF("F1_weights.png"), p1, width = 10, height = 4.5, dpi = 200)

p2 <- perf_tab %>% select(score, R2, within_share_positive) %>% pivot_longer(-score) %>%
  mutate(name = recode(name, R2 = "Out-of-fold R²", within_share_positive = "Sites ranked correctly")) %>%
  ggplot(aes(value, fill = score)) + geom_histogram(bins = 25, colour = "white", position = "identity", alpha = 0.6) +
  facet_wrap(~ name, scales = "free") + scale_fill_manual(values = c(SHMI_base = "#7a7a7a", SHMI = "#1b6ca8"), name = NULL) +
  labs(x = NULL, y = "Repeats", title = "US cross-validated performance over 100 repeats")
ggsave(OUTF("F2_us_cv.png"), p2, width = 9, height = 3.8, dpi = 200)

p3 <- bind_rows(us %>% transmute(set = "US (calibration)", score = "SHMI", SH, pred),
                us %>% transmute(set = "US (calibration)", score = "SHMI_base", SH, pred = pred_base),
                mx %>% transmute(set = "Mexico (independent test)", score = "SHMI", SH, pred),
                mx %>% transmute(set = "Mexico (independent test)", score = "SHMI_base", SH, pred = pred_base)) %>%
  mutate(score = factor(score, c("SHMI", "SHMI_base"))) %>%
  ggplot(aes(pred, SH)) + geom_point(alpha = 0.5, size = 1.2) + geom_abline(linetype = 2, colour = "grey40") +
  geom_smooth(method = "lm", formula = y ~ x, colour = "#d1495b") + facet_grid(score ~ set) +
  labs(x = "Predicted soil health (US calibration)", y = "Soil health (SH)", title = "Observed vs predicted")
ggsave(OUTF("F3_observed_vs_predicted.png"), p3, width = 9, height = 7, dpi = 200)

ws_all <- bind_rows(within_site(us$SH, us$SHMI, us$sitecode) %>% mutate(set = "US"),
                    within_site(mx$SH, mx$SHMI, mx$sitecode) %>% mutate(set = "Mexico"))
p4 <- ggplot(ws_all, aes(rho, fill = set)) + geom_histogram(breaks = seq(-1, 1, 0.25), position = "dodge", colour = "white") +
  geom_vline(xintercept = 0) + scale_fill_manual(values = c(US = "#1b6ca8", Mexico = "#d1495b"), name = NULL) +
  labs(x = "Within-site rank correlation (SHMI vs SH)", y = "Sites", title = "Ranking treatments within sites")
ggsave(OUTF("F4_within_site.png"), p4, width = 7, height = 4, dpi = 200)

cat("\nSaved in", CFG$out_dir, "\n")
cat("\nFor the package: Cover/Roots", round(w_final[["cover"]], 3), "| Diversity", round(w_final[["diversity"]], 3),
    "| Inverse Disturbance", round(w_final[["disturbance"]], 3), "| active-season share", round(w_final[["growing_share"]], 3),
    "| grazing gap", CFG$graze_gap, "| amendment gap", CFG$amend_gap, "\n")
