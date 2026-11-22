################################################################################
# A4a_SHMI_benchmark.R  (v4, 16 Nov 2026: SHMI = published weights held fixed; refitted SHMI is a sensitivity row)
#
# Two questions, both on US units with exactly the calibration's machinery
# (A3a_SHMI_calibration_plus.R: sites weighted equally, 10-fold leave-sites-out CV
# repeated 100 times with the same seed, so the folds are identical to the ones
# behind the published OOF R2 of 0.26; anything with weights is refitted inside
# each fold and mapped to soil health by a site-weighted line):
#
# 1. Does SHMI beat its own four components, and what does each add?
#      SHMI          SHMI (refitted each fold; reproduces the paper) and SHMI_base
#      Drop-one      SHMI without Cover/Roots, Diversity, InvDist or amendments,
#                    and SHMI with equal weights (no calibration)
#      Component     Cover/Roots (season split fitted in the fold), Diversity
#                    (with grazing credit), InvDist, amendments (share of years)
#      Unconstrained weighted OLS and random forest on the same inputs, shown as
#                    a reference only (set include_blackbox = FALSE to drop)
#    Every candidate is scored overall, BETWEEN sites (site means of the OOF
#    predictions vs site-mean soil health) and WITHIN sites (deviations from the
#    site mean), plus the share of sites whose treatments are ranked correctly.
#
# 2. Why does the weighting matter so little?
#    A weight surface: every combination of Cover/Roots, Diversity and InvDist
#    weights on a grid (summing to 1; active-season share and amendment gap held
#    at their calibrated values), each scored out of fold with fixed weights.
#    Shows how flat the surface is, and whether between-site and within-site
#    performance prefer different weights.
#
# Uncertainty: across repeats (fold-assignment noise) and a site-cluster
# bootstrap of the repeat-averaged OOF predictions (sampling noise across sites)
# for the difference from SHMI in overall R2, between-site r and within-site r.
# Each candidate is also fitted on all US units and applied to Mexico.
#
# Input   Outputs/Final_plus/final_results.rds (from A3a_SHMI_calibration_plus.R)
# Outputs Outputs/Benchmark/
#   benchmark_summary.csv    one row per candidate: medians, 95% ranges, paired
#                            differences from SHMI, bootstrap intervals
#   benchmark_repeats.csv    every repeat x candidate
#   benchmark_paired.csv     paired differences from SHMI per repeat
#   benchmark_bootstrap.csv  site-cluster bootstrap of differences from SHMI
#   benchmark_mexico.csv     fitted on all US units, applied to Mexico
#   weight_grid.csv          weight surface (median over grid repeats)
#   benchmark_results.xlsx   all of the above (if writexl is installed)
#   Fig_benchmark.png        candidates x four metrics
#   Fig_weight_surface.png   weight surface, one panel per metric
################################################################################

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(tibble); library(ggplot2); library(here)
})
if (!requireNamespace("nnls", quietly = TRUE)) stop("Install nnls: install.packages('nnls')")

CFG_B <- list(
  final            = here::here("Outputs", "Final_plus", "final_results.rds"),
  out_dir          = here::here("Outputs", "Benchmark"),
  drop_systems     = character(0),  # e.g. c("reference") for cropland only (changes the folds)
  R                = NULL,          # NULL = same repeats as the calibration (100)
  boot_B           = 2000,
  include_blackbox = TRUE,
  rf_trees         = 500,
  grid_step        = 0.05,          # weight-surface resolution
  grid_R           = 20             # repeats for the weight surface (first grid_R calibration fold sets)
)
if (exists("CFG_B_OVERRIDE")) CFG_B <- utils::modifyList(CFG_B, CFG_B_OVERRIDE)
dir.create(CFG_B$out_dir, recursive = TRUE, showWarnings = FALSE)
OUTB <- function(x) file.path(CFG_B$out_dir, x)

# ==============================================================================
# Data
# ==============================================================================
res <- readRDS(CFG_B$final)
CFG <- res$config                               # K, R, seed, fixed_div, graze_gap, amend_gap
R_rep <- if (is.null(CFG_B$R)) CFG$R else CFG_B$R
W_FINAL <- res$w_final                          # cover, diversity, disturbance, growing_share
us <- as_tibble(res$us); mx <- as_tibble(res$mx)
need <- c("MGT_combo", "sitecode", "system", "SH", "NonGrow", "Grow", "Diversity_plants",
          "Diversity", "Richness", "InvDist", "p_amend", "p_animal")
stopifnot(all(need %in% names(us)))
if (length(CFG_B$drop_systems)) {
  us <- us %>% filter(!system %in% CFG_B$drop_systems)
  mx <- mx %>% filter(!system %in% CFG_B$drop_systems)
  message("Dropped systems: ", paste(CFG_B$drop_systems, collapse = ", "),
          " -- folds now differ from the calibration")
}
add_derived <- function(D) D %>% mutate(
  Cover_equal = 0.5 * NonGrow + 0.5 * Grow,
  Cover_cal   = (1 - W_FINAL[["growing_share"]]) * NonGrow + W_FINAL[["growing_share"]] * Grow)
us <- add_derived(us); mx <- add_derived(mx)
cat("US units:", nrow(us), "at", n_distinct(us$sitecode), "sites | Mexican units:", nrow(mx), "\n")

# ==============================================================================
# Helpers (same definitions as A3a_SHMI_calibration_plus.R)
# ==============================================================================
site_w <- function(s) 1 / as.numeric(table(s)[s])
wmean  <- function(x, w) sum(w * x) / sum(w)
r2     <- function(obs, pred) 1 - sum((obs - pred)^2) / sum((obs - mean(obs))^2)
safe_r <- function(a, b) if (length(a) > 2 && sd(a) > 0 && sd(b) > 0) cor(a, b) else NA_real_
line_fit <- function(s, y, w) {
  ms <- wmean(s, w); my <- wmean(y, w)
  b <- sum(w * (s - ms) * (y - my)) / sum(w * (s - ms)^2); c(a = my - b * ms, b = b)
}
OTHER <- c("NonGrow", "Grow", "InvDist")

# Generalized SHMI fit: NNLS on the parts in 'use', Diversity's share fixed at
# 'div', amendment gap 'gap'. With the defaults this is exactly fit_plus() from
# the calibration. Returns weights on the 0-100 scale, or NULL if the fit collapses.
fit_w <- function(D, w, use = OTHER, div = CFG$fixed_div, gap = CFG$amend_gap) {
  k <- if (div > 0) div / (1 - div) else 0
  cp <- gap * D$p_amend
  Xm <- (as.matrix(D[, use, drop = FALSE]) + k * D$Diversity) * (1 - cp) + (1 + k) * 100 * cp
  xm <- colSums(Xm * w) / sum(w); ym <- wmean(D$SH, w)
  f <- nnls::nnls(sweep(Xm, 2, xm) * sqrt(w), (D$SH - ym) * sqrt(w))
  b <- setNames(numeric(length(OTHER)), OTHER); b[use] <- f$x
  B <- (1 + k) * sum(b)
  if (B <= 1e-12) return(NULL)
  b / B
}
score_w <- function(D, b, div = CFG$fixed_div, gap = CFG$amend_gap) {
  base <- as.numeric(as.matrix(D[, OTHER]) %*% b[OTHER]) + div * D$Diversity
  base + gap * D$p_amend * (100 - base)
}

# A "model" is function(Dtr, w) returning a predictor function(D).
const_model <- function(Dtr, w) { m <- wmean(Dtr$SH, w); function(D) rep(m, nrow(D)) }
line_model <- function(score_fun) function(Dtr, w) {
  s <- score_fun(Dtr)
  if (sd(s) == 0) return(const_model(Dtr, w))
  l <- line_fit(s, Dtr$SH, w)
  function(D) l[["a"]] + l[["b"]] * score_fun(D)
}
var_model <- function(v) line_model(function(D) D[[v]])
shmi_model <- function(use = OTHER, div = CFG$fixed_div, gap = CFG$amend_gap, score_gap = gap) function(Dtr, w) {
  b <- fit_w(Dtr, w, use, div, gap)
  if (is.null(b)) return(const_model(Dtr, w))
  line_model(function(D) score_w(D, b, div, score_gap))(Dtr, w)
}
cover_fold_model <- function(Dtr, w) {          # Cover/Roots with the season split fitted in the fold
  b <- fit_w(Dtr, w)
  gs <- if (is.null(b) || sum(b[c("NonGrow", "Grow")]) <= 0) 0.5 else b[["Grow"]] / sum(b[c("NonGrow", "Grow")])
  line_model(function(D) (1 - gs) * D$NonGrow + gs * D$Grow)(Dtr, w)
}
equal_model <- line_model(function(D) {        # no calibration at all: thirds, equal seasons
  base <- (D$Cover_equal + D$Diversity + D$InvDist) / 3
  base + CFG$amend_gap * D$p_amend * (100 - base)
})
fixed_model <- line_model(function(D) {        # published weights held fixed (what a user applies);
  base <- W_FINAL[["cover"]] * D$Cover_cal +   # only the line to soil health is fitted in the fold,
          W_FINAL[["diversity"]] * D$Diversity + # so this row is slightly optimistic
          W_FINAL[["disturbance"]] * D$InvDist
  base + CFG$amend_gap * D$p_amend * (100 - base)
})
BB_VARS <- c("NonGrow", "Grow", "Diversity_plants", "Richness", "InvDist", "p_amend", "p_animal")
ols_model <- function(Dtr, w) {
  Dtr$.w <- w
  f <- stats::lm(stats::reformulate(BB_VARS, "SH"), data = Dtr, weights = .w)
  function(D) as.numeric(stats::predict(f, newdata = D))
}
rf_model <- function(Dtr, w) {
  f <- ranger::ranger(stats::reformulate(BB_VARS, "SH"), data = as.data.frame(Dtr[, c("SH", BB_VARS)]),
                      case.weights = w, num.trees = CFG_B$rf_trees, seed = 1, num.threads = 1)
  function(D) stats::predict(f, data = as.data.frame(D[, BB_VARS]))$predictions
}

# ==============================================================================
# Candidates (SHMI must stay first: every paired comparison is against it)
# ==============================================================================
# 16 Nov 2026: SHMI and its reduced versions use the PUBLISHED (calibrated) weights, held fixed,
# so the comparison with the parts is of the index as users apply it. Only the line from score
# to soil health is fitted in each fold. Dropping a part rescales the remaining published
# weights to sum to one. The version with weights refitted in each fold is kept as a sensitivity
# row (it shows the advantage is not from fitting the weights to the same data).
fixed_base <- function(D, use = c("cover", "diversity", "disturbance")) {
  X <- list(cover = D$Cover_cal, diversity = D$Diversity, disturbance = D$InvDist)
  w <- unlist(W_FINAL[use]); Reduce(`+`, Map(function(v, k) v * X[[k]], w / sum(w), use))
}
fixed_drop <- function(use, amend = TRUE) line_model(function(D) {
  base <- fixed_base(D, use)
  if (amend) base + CFG$amend_gap * D$p_amend * (100 - base) else base
})
cand <- list(
  list("SHMI",                          "SHMI",       fixed_model),
  list("SHMI, weights refitted in each fold", "SHMI", shmi_model()),
  list("SHMI_base (no amendment term)", "SHMI",       fixed_drop(c("cover", "diversity", "disturbance"), amend = FALSE)),
  list("SHMI without Cover/Roots",      "Drop-one",   fixed_drop(c("diversity", "disturbance"))),
  list("SHMI without Diversity",        "Drop-one",   fixed_drop(c("cover", "disturbance"))),
  list("SHMI without InvDist",          "Drop-one",   fixed_drop(c("cover", "diversity"))),
  list("SHMI, equal weights",           "Drop-one",   equal_model),
  list("Cover/Roots",                   "Component",  var_model("Cover_cal")),   # published season split
  list("Diversity",                     "Component",  var_model("Diversity")),
  list("InvDist",                       "Component",  var_model("InvDist")),
  list("Amendments",                    "Component",  var_model("p_amend"))
)
if (isTRUE(CFG_B$include_blackbox)) {
  cand <- c(cand, list(list("OLS on all inputs", "Unconstrained reference", ols_model)))
  if (requireNamespace("ranger", quietly = TRUE)) {
    cand <- c(cand, list(list("Random forest on all inputs", "Unconstrained reference", rf_model)))
  } else message("ranger not installed: random forest skipped")
}
M <- tibble(model = vapply(cand, `[[`, "", 1), group = vapply(cand, `[[`, "", 2),
            fun = lapply(cand, `[[`, 3))
nM <- nrow(M)

# ==============================================================================
# Performance: overall, between sites, within sites, ranking
# ==============================================================================
rankable <- us %>% group_by(sitecode) %>% summarise(ok = n() >= 3 & sd(SH) > 0, .groups = "drop")
N_RANKABLE <- sum(rankable$ok)
PERF_NAMES <- c("n", "R2", "RMSE", "r", "between_R2", "between_r", "within_R2", "within_r",
                "sites_ranked", "within_median_rho", "within_share_positive", "sites_correct_of_all")

perf <- function(obs, pred, site) {
  ok <- !is.na(pred) & !is.na(obs); o <- obs[ok]; p <- pred[ok]; s <- site[ok]
  out <- setNames(rep(NA_real_, length(PERF_NAMES)), PERF_NAMES); out[["n"]] <- length(o)
  if (length(o) < 3) return(out)
  mo <- tapply(o, s, mean); mp <- tapply(p, s, mean)              # site means
  do <- o - ave(o, s); dp <- p - ave(p, s)                        # deviations from site means
  rho <- vapply(split(seq_along(o), s), function(i) {
    if (length(i) < 3 || sd(o[i]) == 0 || sd(p[i]) == 0) NA_real_ else cor(o[i], p[i], method = "spearman")
  }, numeric(1))
  rho <- rho[!is.na(rho)]
  out[c("R2", "RMSE", "r")] <- c(r2(o, p), sqrt(mean((o - p)^2)), safe_r(o, p))
  out[c("between_R2", "between_r")] <- c(if (length(mo) > 2) r2(mo, mp) else NA_real_, safe_r(mo, mp))
  out[c("within_R2", "within_r")] <- c(if (sum(do^2) > 0) 1 - sum((do - dp)^2) / sum(do^2) else NA_real_,
                                        safe_r(do, dp))
  out[["sites_ranked"]] <- length(rho)
  if (length(rho)) out[c("within_median_rho", "within_share_positive")] <- c(median(rho), mean(rho > 0))
  out[["sites_correct_of_all"]] <- sum(rho > 0) / N_RANKABLE     # unrankable sites count as not correct
  out
}

fold_sets <- function(r, sites) {               # identical to the calibration's fold assignment
  set.seed(CFG$seed + r)
  sample(rep_len(seq_len(CFG$K), length(sites)))
}

# ==============================================================================
# 1. Candidates: cross-validation on identical folds
# ==============================================================================
sites <- unique(us$sitecode); y <- us$SH
reps <- vector("list", R_rep * nM); paired <- vector("list", R_rep * nM); ii <- 0
P_sum <- matrix(0, nrow(us), nM)
t0 <- Sys.time()
for (r in seq_len(R_rep)) {
  fold <- fold_sets(r, sites)[match(us$sitecode, sites)]
  P <- matrix(NA_real_, nrow(us), nM)
  for (k in seq_len(CFG$K)) {
    tr <- fold != k; Dtr <- us[tr, ]; Dte <- us[!tr, ]; w <- site_w(Dtr$sitecode)
    for (m in seq_len(nM)) P[!tr, m] <- M$fun[[m]](Dtr, w)(Dte)
  }
  P_sum <- P_sum + P
  p1 <- perf(y, P[, 1], us$sitecode)
  for (m in seq_len(nM)) {
    ii <- ii + 1
    pm <- if (m == 1) p1 else perf(y, P[, m], us$sitecode)
    reps[[ii]] <- c(repeat_id = r, model_id = m, pm)
    paired[[ii]] <- c(repeat_id = r, model_id = m, (pm - p1)[c("R2", "between_r", "within_r", "sites_correct_of_all")])
  }
  if (r %% 10 == 0) cat(sprintf("repeat %d of %d (%.1f min)\n", r, R_rep, as.numeric(difftime(Sys.time(), t0, units = "mins"))))
}
to_tab <- function(L) as_tibble(as.data.frame(do.call(rbind, L))) %>%
  mutate(model = M$model[model_id], .after = repeat_id) %>% select(-model_id)
rep_tab <- to_tab(reps)
pair_tab <- to_tab(paired) %>% rename_with(~ paste0("d_", .x), c(R2, between_r, within_r, sites_correct_of_all))

q <- function(x, p) unname(quantile(x, p, na.rm = TRUE))
summ <- rep_tab %>% group_by(model) %>%
  summarise(across(c(R2, between_r, within_r, sites_correct_of_all),
                   list(lo = ~ q(.x, 0.025), hi = ~ q(.x, 0.975)), .names = "{.col}_{.fn}"),
            across(c(n, R2, RMSE, r, between_R2, between_r, within_R2, within_r, sites_ranked,
                     within_median_rho, within_share_positive, sites_correct_of_all), ~ median(.x, na.rm = TRUE)),
            .groups = "drop")
pair_summ <- pair_tab %>% group_by(model) %>%
  summarise(share_repeats_SHMI_better_R2 = mean(d_R2 < 0),
            share_repeats_SHMI_better_within = mean(d_within_r < 0, na.rm = TRUE),
            across(starts_with("d_"), ~ median(.x, na.rm = TRUE)), .groups = "drop")

pub <- res$us_cv$R2_median[res$us_cv$score == "SHMI"]
now <- summ$R2[summ$model == "SHMI, weights refitted in each fold"]   # the calibration refits in each fold
if (length(pub) && length(CFG_B$drop_systems) == 0 && R_rep == CFG$R) {
  cat(sprintf("\nCheck: refitted SHMI OOF R2 here %.4f vs calibration %.4f %s\n", now, pub,
              if (abs(now - pub) < 1e-6) "(identical)" else "(DIFFERENT -- check inputs)"))
}

if ("SHMI" %in% names(us)) {                  # fixed-weight score should reproduce the published SHMI
  s_fixed <- environment(fixed_model)$score_fun(us)
  cat(sprintf("Fixed-weight score vs published SHMI: max abs difference %.2e\n", max(abs(s_fixed - us$SHMI))))
}

# ==============================================================================
# 2. Site-cluster bootstrap of the repeat-averaged OOF predictions
# ==============================================================================
Pbar <- P_sum / R_rep
rows_by_site <- split(seq_len(nrow(us)), us$sitecode)
BOOT_M <- c("R2", "between_r", "within_r")
set.seed(CFG$seed)
boot <- array(NA_real_, c(CFG_B$boot_B, nM, length(BOOT_M)))
for (b in seq_len(CFG_B$boot_B)) {
  pick <- sample(names(rows_by_site), replace = TRUE)
  idx <- unlist(rows_by_site[pick], use.names = FALSE)
  clus <- rep(seq_along(pick), lengths(rows_by_site[pick]))     # a site drawn twice = two clusters
  p1 <- perf(y[idx], Pbar[idx, 1], clus)[BOOT_M]
  for (m in seq_len(nM)) boot[b, m, ] <- perf(y[idx], Pbar[idx, m], clus)[BOOT_M] - p1
}
avg_perf <- t(sapply(seq_len(nM), function(m) perf(y, Pbar[, m], us$sitecode)[BOOT_M]))
boot_tab <- bind_rows(lapply(seq_along(BOOT_M), function(j) tibble(
  model = M$model, metric = BOOT_M[j], value_avg_oof = avg_perf[, j],
  diff_from_SHMI = apply(boot[, , j, drop = FALSE], 2, median, na.rm = TRUE),
  diff_lo = apply(boot[, , j, drop = FALSE], 2, q, 0.025),
  diff_hi = apply(boot[, , j, drop = FALSE], 2, q, 0.975),
  p_not_worse = apply(boot[, , j, drop = FALSE], 2, function(v) mean(v >= 0, na.rm = TRUE)))))
boot_wide <- boot_tab %>% select(model, metric, diff_from_SHMI, diff_lo, diff_hi) %>%
  pivot_wider(names_from = metric, values_from = c(diff_from_SHMI, diff_lo, diff_hi), names_glue = "boot_{metric}_{.value}")

# ==============================================================================
# 3. Fitted on all US units, applied to Mexico
# ==============================================================================
mx_tab <- if (nrow(mx)) bind_rows(lapply(seq_len(nM), function(m) {
  pm <- M$fun[[m]](us, site_w(us$sitecode))(mx)
  bind_rows(lapply(intersect(c("SH", "SH_all"), names(mx)), function(oc) {
    ok <- !is.na(mx[[oc]])
    pf <- perf(mx[[oc]], pm, mx$sitecode)
    tibble(model = M$model[m], target = if (oc == "SH") "SH from US SEM" else "SH from all-data SEM",
           bias = mean(pm[ok] - mx[[oc]][ok]), !!!as.list(pf[c("n", "R2", "RMSE", "r", "between_r", "within_r",
                                                                 "sites_ranked", "within_median_rho", "within_share_positive")]))
  }))
})) else tibble()

# ==============================================================================
# 4. Weight surface: fixed weights on a grid, scored out of fold
# ==============================================================================
st <- CFG_B$grid_step
grid <- expand.grid(cover = seq(0, 1, by = st), diversity = seq(0, 1, by = st)) %>%
  filter(cover + diversity <= 1 + 1e-9) %>% mutate(disturbance = pmax(0, 1 - cover - diversity))
grid_score <- function(D, g) {
  base <- g$cover * D$Cover_cal + g$diversity * D$Diversity + g$disturbance * D$InvDist
  base + CFG$amend_gap * D$p_amend * (100 - base)
}
GRID_M <- c("R2", "between_r", "within_r", "sites_correct_of_all")
gR <- min(CFG_B$grid_R, R_rep)
folds_g <- lapply(seq_len(gR), function(r) fold_sets(r, sites)[match(us$sitecode, sites)])
cat("\nWeight surface:", nrow(grid), "weight sets x", gR, "repeats\n")
grid_res <- bind_rows(lapply(seq_len(nrow(grid)), function(i) {
  g <- grid[i, ]; s <- grid_score(us, g)
  vals <- sapply(folds_g, function(fold) {
    p <- numeric(nrow(us))
    for (k in seq_len(CFG$K)) {
      tr <- fold != k
      if (sd(s[tr]) == 0) { p[!tr] <- wmean(y[tr], site_w(us$sitecode[tr])); next }
      l <- line_fit(s[tr], y[tr], site_w(us$sitecode[tr])); p[!tr] <- l[["a"]] + l[["b"]] * s[!tr]
    }
    perf(y, p, us$sitecode)[GRID_M]
  })
  tibble(g, !!!as.list(apply(vals, 1, median, na.rm = TRUE)))
}))
ref_pts <- tibble(label = c("Calibrated", "Equal thirds"),
                  cover = c(W_FINAL[["cover"]], 1/3), diversity = c(W_FINAL[["diversity"]], 1/3))
grid_ref <- function(cv, dv) {                   # nearest grid point
  grid_res[which.min((grid_res$cover - cv)^2 + (grid_res$diversity - dv)^2), ]
}
tol <- c(R2 = 0.01, between_r = 0.02, within_r = 0.02, sites_correct_of_all = 0.02)
grid_summary <- bind_rows(lapply(GRID_M, function(mm) {
  v <- grid_res[[mm]]; top <- grid_res[which.max(v), ]     # not 'best': tibble() would see the new column
  cal <- grid_ref(ref_pts$cover[1], ref_pts$diversity[1]); eq <- grid_ref(1/3, 1/3)
  tibble(metric = mm, best = max(v, na.rm = TRUE), best_cover = top$cover, best_diversity = top$diversity,
         best_disturbance = top$disturbance, at_calibrated = cal[[mm]], at_equal_thirds = eq[[mm]],
         tolerance = tol[[mm]], share_of_grid_within_tol = mean(v >= max(v, na.rm = TRUE) - tol[[mm]], na.rm = TRUE))
}))

# ==============================================================================
# Save and print
# ==============================================================================
summary_tab <- M %>% select(model, group) %>%
  left_join(summ, by = "model") %>% left_join(pair_summ, by = "model") %>% left_join(boot_wide, by = "model") %>%
  mutate(model = factor(model, M$model)) %>% arrange(model) %>% mutate(model = as.character(model))
write.csv(summary_tab, OUTB("benchmark_summary.csv"), row.names = FALSE)
write.csv(rep_tab, OUTB("benchmark_repeats.csv"), row.names = FALSE)
write.csv(pair_tab, OUTB("benchmark_paired.csv"), row.names = FALSE)
write.csv(boot_tab, OUTB("benchmark_bootstrap.csv"), row.names = FALSE)
write.csv(mx_tab, OUTB("benchmark_mexico.csv"), row.names = FALSE)
write.csv(grid_res, OUTB("weight_grid.csv"), row.names = FALSE)
if (requireNamespace("writexl", quietly = TRUE)) {
  writexl::write_xlsx(list(summary = summary_tab, bootstrap = boot_tab, mexico = mx_tab,
                           weight_grid_summary = grid_summary, weight_grid = grid_res,
                           paired_repeats = pair_tab, repeats = rep_tab), OUTB("benchmark_results.xlsx"))
}

fmt_ci <- function(v, lo, hi) sprintf("%.2f [%.2f, %.2f]", v, lo, hi)
fmt_d  <- function(m, v, lo, hi) ifelse(m == "SHMI", "--", sprintf("%+.3f [%+.3f, %+.3f]", v, lo, hi))
cat("\nOut-of-fold performance on US units (median [95% range] over", R_rep, "repeats; identical folds)\n")
summary_tab %>% transmute(
  model,
  `Overall R2`     = fmt_ci(R2, R2_lo, R2_hi),
  `Between r`      = fmt_ci(between_r, between_r_lo, between_r_hi),
  `Within r`       = fmt_ci(within_r, within_r_lo, within_r_hi),
  `Sites correct`  = sprintf("%.0f%% of %d", 100 * sites_correct_of_all, N_RANKABLE)) %>%
  as.data.frame() %>% print(row.names = FALSE)
cat("\nDifference from SHMI, site-cluster bootstrap (median [95% CI])\n")
summary_tab %>% transmute(
  model,
  `d overall R2` = fmt_d(model, boot_R2_diff_from_SHMI, boot_R2_diff_lo, boot_R2_diff_hi),
  `d between r`  = fmt_d(model, boot_between_r_diff_from_SHMI, boot_between_r_diff_lo, boot_between_r_diff_hi),
  `d within r`   = fmt_d(model, boot_within_r_diff_from_SHMI, boot_within_r_diff_lo, boot_within_r_diff_hi)) %>%
  as.data.frame() %>% print(row.names = FALSE)
cat("\nWeight surface (", gR, " repeats): best, calibrated and equal-thirds values; share of grid within tolerance of best\n", sep = "")
print(grid_summary %>% mutate(across(where(is.numeric), ~ round(.x, 3))) %>% as.data.frame(), row.names = FALSE)
if (nrow(mx_tab)) {
  cat("\nMexico (fitted on all US units)\n")
  print(mx_tab %>% mutate(across(where(is.numeric), ~ round(.x, 2))) %>% as.data.frame(), row.names = FALSE)
}

# ==============================================================================
# Figures
# ==============================================================================
METRIC_LAB <- c(R2 = "Overall out-of-fold R²", between_r = "Between sites: r of site means",
                within_r = "Within sites: r of deviations from site mean",
                sites_correct_of_all = "Sites with treatments ranked correctly")
grp_col <- c("SHMI" = "#1b6ca8", "Drop-one" = "#7fb3d5", "Component" = "#8a8a8a",
             "Unconstrained reference" = "#e6a141")
grp_col <- grp_col[names(grp_col) %in% M$group]
fig_d <- rep_tab %>% left_join(M %>% select(model, group), by = "model") %>%
  select(model, group, all_of(names(METRIC_LAB))) %>%
  pivot_longer(-c(model, group), names_to = "metric") %>%
  group_by(model, group, metric) %>%
  summarise(med = median(value, na.rm = TRUE), lo = q(value, 0.025), hi = q(value, 0.975), .groups = "drop") %>%
  mutate(model = factor(model, rev(M$model)), group = factor(group, names(grp_col)),
         metric = factor(METRIC_LAB[metric], METRIC_LAB))
ref <- fig_d %>% filter(model == "SHMI") %>% select(metric, ref = med)
pf <- ggplot(fig_d, aes(med, model, colour = group)) +
  geom_vline(data = ref, aes(xintercept = ref), linetype = 2, colour = "#1b6ca8") +
  geom_linerange(aes(xmin = lo, xmax = hi)) + geom_point(size = 2.2) +
  facet_wrap(~ metric, nrow = 1, scales = "free_x", labeller = label_wrap_gen(24)) +
  scale_colour_manual(values = grp_col, name = NULL) +
  labs(x = sprintf("Median and 95%% range over %d repeats of %d-fold leave-sites-out CV", R_rep, CFG$K), y = NULL,
       title = "SHMI against its components and reduced versions (US units, out of fold)") +
  theme_bw(base_size = 10) + theme(legend.position = "bottom", panel.grid.minor = element_blank())
ggsave(OUTB("Fig_benchmark.png"), pf, width = 13, height = 0.3 * nM + 2.4, dpi = 300)

surf <- grid_res %>% pivot_longer(all_of(GRID_M), names_to = "metric") %>%
  group_by(metric) %>% mutate(near_best = value >= max(value, na.rm = TRUE) - tol[metric[1]]) %>% ungroup() %>%
  mutate(metric = factor(METRIC_LAB[metric], METRIC_LAB))
best_pts <- surf %>% group_by(metric) %>% slice_max(value, n = 1, with_ties = FALSE) %>% ungroup()
ps <- ggplot(surf, aes(cover, diversity)) +
  geom_tile(aes(fill = value), width = st, height = st) +
  geom_tile(data = filter(surf, near_best), width = st, height = st, fill = NA, colour = "white", linewidth = 0.25) +
  geom_point(data = best_pts, shape = 8, size = 2.5, colour = "white") +
  geom_point(data = ref_pts, aes(shape = label), size = 2.8, colour = "#d1495b", fill = "#d1495b") +
  facet_wrap(~ metric, nrow = 1, labeller = label_wrap_gen(24)) +
  scale_fill_viridis_c(name = NULL) + scale_shape_manual(values = c(Calibrated = 21, `Equal thirds` = 24), name = NULL) +
  coord_equal() +
  labs(x = "Cover/Roots weight", y = "Diversity weight",
       title = "Out-of-fold performance across SHMI weights (InvDist = 1 - Cover - Diversity)",
       subtitle = sprintf(paste0("Fixed weights, %d repeats; active-season share %.2f and amendment gap %.2f held at calibrated values.\n",
                                 "White outline: within %s of the best (R² 0.01; r and ranking 0.02). Star: best."),
                          gR, W_FINAL[["growing_share"]], CFG$amend_gap, "tolerance")) +
  theme_bw(base_size = 10) + theme(legend.position = "bottom", panel.grid = element_blank())
ggsave(OUTB("Fig_weight_surface.png"), ps, width = 15, height = 5.2, dpi = 300)

cat("\nSaved in", CFG_B$out_dir, "\n")
