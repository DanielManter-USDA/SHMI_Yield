################################################################################
# A2d_SHMI_calibration_final.R
#
# Final SHMI calibration.
# Structure:
#   Cover       living cover in the growing and non-growing season, by site
#               climate (month mean temperature >= 5 C and, unless irrigated,
#               precipitation >= 2 x temperature); growing-season share estimated
#   Diversity   average annual plant species richness, 100 x (richness - 1) / 7
#   Disturbance inverse disturbance from STIR
#   Organic     share of years with an organic amendment or animals
#
# 1. Calibration on the US units only (noMX; SEM target fitted on US units):
#    10-fold leave-sites-out cross-validation repeated 100 times (1,000 fits).
#    Each fit: non-negative least squares (sites weighted equally); pillar
#    weights normalized to sum to 1. Final weights = medians across fits,
#    renormalized; intervals = 2.5-97.5 percentiles.
# 2. Independent test on the Mexican units: their SH comes from the US SEM
#    applied to them (A2c_Project_SH.R); weights and the SH-scale calibration come
#    from the US units only. The all-data SH is reported as a check.
#
# Outputs: Outputs/Final/ final_weights.csv, final_results.xlsx,
#          shmi_final_scores.csv, final_results.rds and figures F1-F4.
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
  treatments  = here::here("Inputs", "napeshm", "treatments.csv"),
  out_dir     = file.path(BASE, "Final"),
  test_mexico = TRUE,        # FALSE: skip the independent Mexico test
  K           = 10,
  R           = 100,
  seed        = 2026
)
if (exists("CFG_OVERRIDE")) CFG <- utils::modifyList(CFG, CFG_OVERRIDE)
dir.create(CFG$out_dir, recursive = TRUE, showWarnings = FALSE)
OUTF <- function(x) file.path(CFG$out_dir, x)

# ==============================================================================
# Components for every unit
# ==============================================================================
inp <- readRDS(CFG$inputs)
# monthly climate for each unit: its site's normals, plus whether it is irrigated
irr <- read.csv(CFG$treatments, stringsAsFactors = FALSE) %>%
  transmute(MGT_combo = paste0("NAPESHM_", sitecode, "_", treatmentid),
            irrigated = !is.na(eu_irrigation) & !(trimws(tolower(eu_irrigation)) %in% c("none", "", "n_a")))
climate <- inp$mgt %>%
  transmute(MGT_combo, sitecode = sub("^NAPESHM_(.*)_[^_]+$", "\\1", MGT_combo)) %>%
  left_join(readRDS(CFG$climate), by = "sitecode") %>%
  left_join(irr, by = "MGT_combo") %>%
  mutate(irrigated = coalesce(irrigated, FALSE)) %>% select(-sitecode)

# the four sub-indices exactly as build_shmi() computes them; Cover is kept in
# its two seasons so the growing-season share can be estimated
X <- shmi_components(inp, climate, dist_meth = "STIR") %>%   # NAPESHM records STIR
  transmute(MGT_combo, NonGrow = Cover_nongrowing, Grow = Cover_growing,
            Diversity, Richness, InvDist, Organic = OrgInput)
COLS <- c("NonGrow", "Grow", "Diversity", "InvDist", "Organic")

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

# ==============================================================================
# Helpers
# ==============================================================================
site_w <- function(s) 1 / as.numeric(table(s)[s])
fit_nnls <- function(D, yy, ww) {
  Xm <- as.matrix(D[, COLS]); xm <- colSums(Xm * ww) / sum(ww); ym <- sum(yy * ww) / sum(ww)
  f <- nnls::nnls(sweep(Xm, 2, xm) * sqrt(ww), (yy - ym) * sqrt(ww))
  b <- stats::setNames(f$x, COLS); list(a = ym - sum(b * xm), b = b)
}
predict_fit <- function(f, D) f$a + as.vector(as.matrix(D[, COLS]) %*% f$b)
weights_of <- function(b) {
  cover <- b[["NonGrow"]] + b[["Grow"]]
  cc <- c(cover = cover, diversity = b[["Diversity"]], disturbance = b[["InvDist"]], organic = b[["Organic"]])
  tot <- sum(cc)
  c(if (tot > 0) cc / tot else cc * NA, growing_share = if (cover > 0) b[["Grow"]] / cover else NA_real_)
}
shmi_score <- function(D, w) {
  cover <- (1 - w[["growing_share"]]) * D$NonGrow + w[["growing_share"]] * D$Grow
  w[["cover"]] * cover + w[["diversity"]] * D$Diversity + w[["disturbance"]] * D$InvDist + w[["organic"]] * D$Organic
}
within_site <- function(obs, pred, site, min_n = 3) {
  tibble(obs, pred, site) %>% group_by(site) %>% filter(n() >= min_n, sd(pred) > 0, sd(obs) > 0) %>%
    summarise(rho = cor(obs, pred, method = "spearman"), .groups = "drop")
}
r2 <- function(obs, pred) 1 - sum((obs - pred)^2) / sum((obs - mean(obs))^2)

# ==============================================================================
# 1. Calibration on US units: K-fold leave-sites-out, repeated R times
# ==============================================================================
sites <- unique(us$sitecode); y <- us$SH; wus <- site_w(us$sitecode)
fits <- list(); perf <- list(); oof_sum <- numeric(nrow(us))
for (r in seq_len(CFG$R)) {
  set.seed(CFG$seed + r)
  fold_of_site <- sample(rep_len(seq_len(CFG$K), length(sites))); fold <- fold_of_site[match(us$sitecode, sites)]
  pred <- numeric(nrow(us))
  for (k in seq_len(CFG$K)) {
    tr <- fold != k
    f <- fit_nnls(us[tr, ], y[tr], site_w(us$sitecode[tr]))
    pred[!tr] <- predict_fit(f, us[!tr, ])
    fits[[length(fits) + 1]] <- c(repeat_id = r, fold = k, weights_of(f$b))
  }
  oof_sum <- oof_sum + pred
  ws <- within_site(y, pred, us$sitecode)
  perf[[r]] <- tibble(repeat_id = r, R2 = r2(y, pred), RMSE = sqrt(mean((y - pred)^2)),
                      within_median_rho = median(ws$rho), within_share_positive = mean(ws$rho > 0))
}
fit_tab  <- as_tibble(do.call(rbind, fits))
perf_tab <- bind_rows(perf)
cat("\nFits:", nrow(fit_tab), "\n")

WN <- c("cover", "diversity", "disturbance", "organic", "growing_share")
weight_summary <- map_dfr(WN, function(w) {
  v <- fit_tab[[w]]
  tibble(weight = w, median = median(v, na.rm = TRUE), lo95 = quantile(v, 0.025, na.rm = TRUE),
         hi95 = quantile(v, 0.975, na.rm = TRUE), mean = mean(v, na.rm = TRUE),
         share_zero = mean(v <= 1e-10, na.rm = TRUE))
})
w_final <- setNames(weight_summary$median, weight_summary$weight)
pill <- c("cover", "diversity", "disturbance", "organic")
w_final[pill] <- w_final[pill] / sum(w_final[pill])               # renormalize medians to sum to 1
weight_summary <- weight_summary %>% mutate(final = unname(w_final[weight]))
cat("\nFinal weights (medians of", nrow(fit_tab), "fits, renormalized):\n"); print(weight_summary, width = Inf)
us_cv <- perf_tab %>% summarise(across(-repeat_id, list(median = median, lo95 = ~ quantile(.x, 0.025),
                                                        hi95 = ~ quantile(.x, 0.975))))

# out-of-fold predictions (averaged over repeats) and the split of their error
# into a between-site part (whole sites too high or too low) and a within-site part
oof <- us %>% transmute(MGT_combo, sitecode, system, SH, pred_oof = oof_sum / CFG$R) %>%
  group_by(sitecode) %>% mutate(err = SH - pred_oof, site_err = mean(err), within_err = err - site_err) %>% ungroup()
error_split <- tibble(
  RMSE_total = sqrt(mean(oof$err^2)),
  share_between_sites = sum(oof$site_err^2) / sum(oof$err^2),
  share_within_sites  = sum(oof$within_err^2) / sum(oof$err^2),
  RMSE_within_only = sqrt(mean(oof$within_err^2)),
  RMSE_between_only = sqrt(mean(oof$site_err^2)))
cat("\nOut-of-fold error split:\n"); print(error_split, width = Inf)

# single full-data fit, for comparison
w_full <- weights_of(fit_nnls(us, y, wus)$b)

# ==============================================================================
# 2. Final SHMI for all units; SH-scale calibration on US; independent Mexico test
# ==============================================================================
all_units <- X %>% mutate(SHMI = shmi_score(., w_final))
us$SHMI <- shmi_score(us, w_final); mx$SHMI <- shmi_score(mx, w_final)
cal <- stats::lm(SH ~ SHMI, data = us, weights = wus)          # SH scale from US units only
us$pred <- predict(cal, us); mx$pred <- predict(cal, mx)

test_stats <- function(D, obs_col, label) {
  obs <- D[[obs_col]]; ok <- !is.na(obs)
  ws <- within_site(obs[ok], D$SHMI[ok], D$sitecode[ok])
  sm <- tibble(o = obs[ok], p = D$pred[ok], s = D$sitecode[ok]) %>% group_by(s) %>%
    summarise(o = mean(o), p = mean(p), .groups = "drop")
  tibble(set = label, units = sum(ok), sites = length(unique(D$sitecode[ok])),
         r = cor(obs[ok], D$SHMI[ok]), R2 = r2(obs[ok], D$pred[ok]),
         RMSE = sqrt(mean((obs[ok] - D$pred[ok])^2)), bias = mean(D$pred[ok] - obs[ok]),
         between_site_r = cor(sm$o, sm$p),
         sites_ranked = nrow(ws), within_median_rho = median(ws$rho), within_share_positive = mean(ws$rho > 0),
         sign_test_p = if (nrow(ws)) binom.test(sum(ws$rho > 0), sum(ws$rho != 0))$p.value else NA_real_)
}
test_table <- bind_rows(
  test_stats(us, "SH", "US units, final weights (in-sample)"),
  if (nrow(mx)) test_stats(mx, "SH", "Mexico: SH from US SEM (independent test)"),
  if (nrow(mx)) test_stats(mx, "SH_all", "Mexico: SH from all-data SEM (check)"))
cat("\nPerformance:\n"); print(test_table, width = Inf)
cat("\nUS cross-validated performance (median and 95% range over", CFG$R, "repeats):\n"); print(us_cv, width = Inf)

# ==============================================================================
# Save
# ==============================================================================
write.csv(weight_summary, OUTF("final_weights.csv"), row.names = FALSE)
write.csv(all_units, OUTF("shmi_final_scores.csv"), row.names = FALSE)
write.csv(oof, OUTF("us_oof_predictions.csv"), row.names = FALSE)
writexl::write_xlsx(list(final_weights = weight_summary,
                         full_fit_weights = tibble(weight = names(w_full), estimate = unname(w_full)),
                         us_cv_performance = us_cv, error_split = error_split, test = test_table,
                         calibration = tibble(term = names(coef(cal)), estimate = unname(coef(cal))),
                         fold_weights = fit_tab, repeats = perf_tab,
                         mexico_units = mx %>% select(MGT_combo, sitecode, SH, SH_all, SHMI, pred)),
                    OUTF("final_results.xlsx"))
saveRDS(list(config = CFG, weights = weight_summary, w_final = w_final, w_full = w_full, fits = fit_tab,
             us_cv = us_cv, error_split = error_split, oof = oof,
             repeats = perf_tab, test = test_table, calibration = cal,
             us = us, mx = mx, scores = all_units), OUTF("final_results.rds"))

# ==============================================================================
# Figures
# ==============================================================================
theme_set(theme_bw(base_size = 11))
lab <- c(cover = "Cover", diversity = "Diversity", disturbance = "Disturbance", organic = "Organic", growing_share = "Growing-season\nshare (in Cover)")
p1 <- fit_tab %>% select(all_of(WN)) %>% pivot_longer(everything(), names_to = "weight") %>%
  mutate(weight = factor(lab[weight], lab)) %>%
  ggplot(aes(weight, value)) + geom_violin(fill = "grey85", colour = NA, scale = "width") +
  geom_boxplot(width = 0.15, outlier.size = 0.4) +
  geom_point(data = weight_summary %>% mutate(weight = factor(lab[weight], lab)), aes(y = final), colour = "#d1495b", size = 3) +
  labs(x = NULL, y = "Weight across 1,000 fits", title = "Final SHMI weights (red: final value = median, renormalized)")
ggsave(OUTF("F1_weights.png"), p1, width = 9, height = 4.5, dpi = 200)

p2 <- perf_tab %>% select(R2, within_share_positive) %>% pivot_longer(everything()) %>%
  mutate(name = recode(name, R2 = "Out-of-fold R\u00b2", within_share_positive = "Sites ranked correctly")) %>%
  ggplot(aes(value)) + geom_histogram(bins = 25, fill = "#1b6ca8", colour = "white") + facet_wrap(~ name, scales = "free") +
  labs(x = NULL, y = "Repeats", title = "US cross-validated performance over 100 repeats")
ggsave(OUTF("F2_us_cv.png"), p2, width = 9, height = 3.8, dpi = 200)

p3 <- bind_rows(us %>% transmute(set = "US (calibration)", SH, pred, sitecode),
                mx %>% transmute(set = "Mexico (independent test)", SH, pred, sitecode)) %>%
  ggplot(aes(pred, SH)) + geom_point(alpha = 0.5, size = 1.4) + geom_abline(linetype = 2, colour = "grey40") +
  geom_smooth(method = "lm", formula = y ~ x, colour = "#d1495b") + facet_wrap(~ set) +
  labs(x = "Predicted soil health (from SHMI, US calibration)", y = "Soil health (SH)", title = "Observed vs predicted")
ggsave(OUTF("F3_observed_vs_predicted.png"), p3, width = 9, height = 4.5, dpi = 200)

ws_all <- bind_rows(within_site(us$SH, us$SHMI, us$sitecode) %>% mutate(set = "US"),
                    within_site(mx$SH, mx$SHMI, mx$sitecode) %>% mutate(set = "Mexico"))
p4 <- ggplot(ws_all, aes(rho, fill = set)) + geom_histogram(breaks = seq(-1, 1, 0.25), position = "dodge", colour = "white") +
  geom_vline(xintercept = 0) + scale_fill_manual(values = c(US = "#1b6ca8", Mexico = "#d1495b"), name = NULL) +
  labs(x = "Within-site rank correlation (SHMI vs SH)", y = "Sites", title = "Ranking treatments within sites")
ggsave(OUTF("F4_within_site.png"), p4, width = 7, height = 4, dpi = 200)
cat("\nSaved final calibration in", CFG$out_dir, "\n")
