################################################################################
# SHMI_weight_tests.R
#
# Tests of the seasonal cover weights and the amendment/animal weights.
#
# Because Cover and OrgInput are linear in their components, every weighting
# hypothesis is a linear restriction on a regression of SH on
#   C_winter, C_spring, C_summer, C_fall, Diversity, InvDist, O_amend, O_animal
#   equal seasonal weights    <=> b_winter = b_spring = b_summer = b_fall
#   official seasonal weights <=> b_s / w_s(official) equal for all s
#   warm/cool two groups      <=> b_winter = b_spring and b_summer = b_fall
#   equal amend/animal        <=> b_amend = b_animal
# Pillar weights stay free in every test, so the seasonal and organic-input
# questions are isolated.
#
#   Part A  inferential: cluster-robust Wald tests, wild cluster bootstrap p
#   Part B  predictive : repeated leave-site-out CV of restricted vs free
#   Part C  uncertainty: site-cluster bootstrap of the free weights
#   Part D  latent     : MIMIC model - weights estimated inside the SEM
#                        (no factor-score step), plus a test of whether
#                        management acts on SH only through one SH factor
#
# Every part is run for each scenario in SCENARIOS (current workbook, SHMI
# 1.2.0 inputs, and the SEM targets from SH_SEM_refinement_clean.R):
#   one scenario per dataset in datasets.R (all, noMX, noGZ, noMX_noGZ), each
#   with the SH target fitted to that same dataset
# The final section compares the scenarios side by side. A weight structure
# is general if its conclusions hold across scenarios; if they change when a
# region is dropped, the weight mainly reflects that region.
################################################################################

library(SHMI)
library(dplyr)
library(tidyr)
library(purrr)
library(tibble)
library(nnls)
library(lavaan)
library(here)

OUT <- function(...) here::here("Outputs", ...)
source(here::here("R", "datasets.R"))

CFG <- list(
  # Nonlinear settings, fixed a priori for these tests. hill is selected
  # later, inside nested CV in SHMI_calibration_v2.R; here it is held at the
  # a priori value. max_stir = 317.5 puts annual STIR 80 (NRCS 345 ceiling
  # for reduced tillage) at the top of the reduced-tillage classes.
  nonlinear = list(dist_meth = "STIR", max_stir = 317.5, ti_rep = "max",
                   hill = 1, max_div = 8),
  target        = "SH",
  complete_only = FALSE,   # TRUE: only rows with all SH indicators observed
  B_wild        = 9999,    # wild cluster bootstrap draws (Part A)
  B_boot        = 2000,    # site bootstrap draws (Part C)
  R_cv          = 20,      # CV repeats (Part B)
  K_cv          = 5,       # CV folds (Part B)
  seed          = 2026,
  out           = OUT("weight_tests_results.rds"),
  out_xlsx      = OUT("weight_tests_comparison.xlsx")
)

# SH targets (from SH_SEM_refinement_clean.R) and SHMI inputs (from Prepare_SEM.R)
TARGETS <- lapply(setNames(names(DATASETS), names(DATASETS)), function(k)
  list(scores = OUT(paste0("sem_", k, "_scores.rds")),
       model  = OUT(paste0("sem_", k, "_model.rds")),
       report = OUT(paste0("sem_", k, "_report.rds"))))
SHMI_INPUTS <- OUT("shmi_inputs_v120.rds")

# Each dataset against its own SH target (the target already holds only that
# dataset's units)
SCENARIOS <- lapply(setNames(names(DATASETS), names(DATASETS)), function(k)
  list(target = k, subset = function(d) d))

W_OFF_SEASON <- c(C_winter = 0.1259, C_spring = 0.1260, C_summer = 0.3755, C_fall = 0.3726)
W_OFF_ORG    <- c(O_amend = 0.6615, O_animal = 0.3385)
C_cols <- names(W_OFF_SEASON)
O_cols <- names(W_OFF_ORG)
PRED   <- c(C_cols, "Diversity", "InvDist", O_cols)

# Hypotheses as lists of restriction rows (named coefficient vectors).
# A hypothesis is skipped in a scenario if it involves a component with no
# variation there (e.g. no animals in the subset).
HYPOTHESES <- list(
  season_equal    = list(c(C_winter = 1, C_spring = -1),
                         c(C_spring = 1, C_summer = -1),
                         c(C_summer = 1, C_fall   = -1)),
  season_official = lapply(1:3, function(i)
                      setNames(c(1 / W_OFF_SEASON[[i]], -1 / W_OFF_SEASON[[i + 1]]),
                               C_cols[i:(i + 1)])),
  season_2group   = list(c(C_winter = 1, C_spring = -1),
                         c(C_summer = 1, C_fall   = -1)),
  org_equal       = list(c(O_amend = 1, O_animal = -1)),
  org_official    = list(setNames(c(1 / W_OFF_ORG[[1]], -1 / W_OFF_ORG[[2]]), O_cols))
)

# ==============================================================================
# 0. Data (loaded once; scenarios are subsets)
# ==============================================================================
shmi_inputs <- readRDS(SHMI_INPUTS)
comp <- compute_shmi_components(shmi_inputs, settings = CFG$nonlinear)

# For each target: analysis data (SH + components), SEM syntax, SEM data
TARG <- lapply(TARGETS, function(t) {
  d <- readRDS(t$scores) %>%
    transmute(MGT_combo, system, sitecode = as.character(sitecode), n_missing,
              y = .data[[CFG$target]]) %>%
    inner_join(comp, by = "MGT_combo") %>%
    filter(!is.na(y))
  if (CFG$complete_only) d <- filter(d, n_missing == 0)
  list(d = d,
       sem_spec = readRDS(t$model),
       sem_data = readRDS(t$report)$prepared_data %>%
         select(-any_of(c("Cover", "Diversity", "InvDist", "OrgInput", "SHMI"))))
})


# ==============================================================================
# Helpers
# ==============================================================================

# Components with any variation in this data set
active_components <- function(d) PRED[vapply(d[PRED], function(x) sd(x) > 0, logical(1))]

# NNLS with free intercept (centering is exact for LS with an intercept)
fit_nnls <- function(Xm, y) {
  xm <- colMeans(Xm); ym <- mean(y)
  f  <- nnls::nnls(sweep(Xm, 2, xm), y - ym)
  list(a = ym - sum(f$x * xm), b = f$x)
}
pred_lin <- function(f, Xm) f$a + as.vector(Xm %*% f$b)

# Scenario description: size, where low-summer units come from, animals
describe_data <- function(d) {
  list(
    n_rows   = nrow(d),
    n_sites  = length(unique(d$sitecode)),
    cor      = round(cor(d[, C_cols]), 2),
    sd       = round(sapply(d[, C_cols], sd), 1),
    low_summer_n     = sum(d$C_summer < 50),
    low_summer_sites = sort(table(d$sitecode[d$C_summer < 50]), decreasing = TRUE),
    animals  = table(d$system, animals = d$O_animal > 0)
  )
}


# ==============================================================================
# Part A. Cluster-robust Wald tests with wild cluster restricted bootstrap
# ==============================================================================
part_A <- function(d, active, B) {
  X    <- cbind(`(Intercept)` = 1, as.matrix(d[, active]))
  yv   <- d$y
  n    <- nrow(X); k <- ncol(X)
  XtXi <- solve(crossprod(X))
  cl   <- as.integer(factor(d$sitecode))
  G    <- max(cl)

  vcov_cr1 <- function(e) {
    S   <- rowsum(X * e, cl)                               # G x k cluster scores
    adj <- (G / (G - 1)) * ((n - 1) / (n - k))
    adj * XtXi %*% crossprod(S) %*% XtXi
  }
  wald <- function(b, V, R) {
    Rb <- R %*% b
    as.numeric(t(Rb) %*% solve(R %*% V %*% t(R), Rb))
  }
  row_R <- function(coefs) {
    r <- setNames(numeric(k), colnames(X)); r[names(coefs)] <- coefs; r
  }

  b_hat <- as.vector(XtXi %*% crossprod(X, yv)); names(b_hat) <- colnames(X)
  V_hat <- vcov_cr1(as.vector(yv - X %*% b_hat))

  wcr_test <- function(R) {
    A     <- XtXi %*% t(R) %*% solve(R %*% XtXi %*% t(R))
    b_r   <- as.vector(b_hat - A %*% (R %*% b_hat))       # restricted estimate
    fit_r <- as.vector(X %*% b_r); e_r <- yv - fit_r
    W_obs <- wald(b_hat, V_hat, R)
    W_star <- replicate(B, {
      v  <- sample(c(-1, 1), G, replace = TRUE)[cl]       # Rademacher per site
      ys <- fit_r + e_r * v
      bs <- as.vector(XtXi %*% crossprod(X, ys))
      wald(bs, vcov_cr1(ys - as.vector(X %*% bs)), R)
    })
    tibble(W = W_obs, df = nrow(R),
           p_asymptotic = pchisq(W_obs, nrow(R), lower.tail = FALSE),
           p_wild_boot  = (1 + sum(W_star >= W_obs)) / (B + 1))
  }

  usable <- keep(HYPOTHESES, function(h) all(unlist(lapply(h, names)) %in% active))
  tests  <- imap_dfr(usable, function(h, nm) {
    R <- do.call(rbind, lapply(h, row_R))
    wcr_test(R) %>% mutate(H0 = nm, .before = 1)
  })

  list(coef  = tibble(term = names(b_hat), estimate = b_hat,
                      se = sqrt(diag(V_hat)), z = b_hat / sqrt(diag(V_hat))),
       tests = tests,
       G     = G)
}


# ==============================================================================
# Part B. Predictive comparison of weight structures (pillar weights free)
# ==============================================================================
SEASON_DESIGN <- list(
  equal    = function(d) tibble(Cov = rowMeans(d[, C_cols])),
  official = function(d) tibble(Cov = as.vector(as.matrix(d[, C_cols]) %*% W_OFF_SEASON)),
  twogroup = function(d) tibble(Cool = (d$C_winter + d$C_spring) / 2,
                                Warm = (d$C_summer + d$C_fall) / 2),
  free     = function(d) d[, C_cols]
)
ORG_DESIGN <- list(
  equal    = function(d) tibble(Org = rowMeans(d[, O_cols])),
  official = function(d) tibble(Org = as.vector(as.matrix(d[, O_cols]) %*% W_OFF_ORG)),
  free     = function(d) d[, O_cols]
)
SPECS <- expand.grid(season = names(SEASON_DESIGN), org = names(ORG_DESIGN),
                     stringsAsFactors = FALSE)

part_B <- function(d, R_cv, K_cv, seed) {
  yv   <- d$y
  site <- d$sitecode
  Xs   <- map(seq_len(nrow(SPECS)), function(j)
    as.matrix(bind_cols(SEASON_DESIGN[[SPECS$season[j]]](d),
                        d[, c("Diversity", "InvDist")],
                        ORG_DESIGN[[SPECS$org[j]]](d))))

  fold_by_site <- function(s) {
    u <- unique(site); set.seed(s)
    unname(setNames(sample(rep_len(seq_len(K_cv), length(u))), u)[site])
  }

  raw <- map_dfr(seq_len(R_cv), function(r) {
    fold <- fold_by_site(seed + r)
    map_dfr(seq_len(nrow(SPECS)), function(j) {
      p <- numeric(nrow(d))
      for (kk in seq_len(K_cv)) {
        tr <- fold != kk
        p[!tr] <- pred_lin(fit_nnls(Xs[[j]][tr, , drop = FALSE], yv[tr]),
                           Xs[[j]][!tr, , drop = FALSE])
      }
      tibble(rep_id = r, season = SPECS$season[j], org = SPECS$org[j],
             RMSE = sqrt(mean((yv - p)^2)),
             R2   = 1 - sum((yv - p)^2) / sum((yv - mean(yv))^2))
    })
  })

  ref <- raw %>% filter(season == "equal", org == "equal") %>% select(rep_id, RMSE_ref = RMSE)
  summary <- raw %>% left_join(ref, by = "rep_id") %>%
    group_by(season, org) %>%
    summarise(RMSE = mean(RMSE), R2 = mean(R2),
              gain = mean(RMSE_ref - RMSE), win_rate = mean(RMSE_ref - RMSE > 0),
              .groups = "drop") %>%
    arrange(desc(gain))

  list(raw = raw, summary = summary)
}


# ==============================================================================
# Part C. Site-cluster bootstrap of the free weights
# ==============================================================================
part_C <- function(d, active, B) {
  yv      <- d$y
  X_free  <- as.matrix(d[, PRED])
  rows_of <- split(seq_len(nrow(d)), d$sitecode)
  org_ok  <- all(O_cols %in% active)      # amend/animal split identifiable?

  weights_from_b <- function(b) {
    b  <- setNames(b, PRED)
    bs <- b[C_cols]; bo <- b[O_cols]
    c(if (sum(bs) > 0) bs / sum(bs) else setNames(rep(NA, 4), C_cols),
      if (org_ok && sum(bo) > 0) bo / sum(bo) else setNames(rep(NA, 2), O_cols),
      warm_share = if (sum(bs) > 0) sum(bs[c("C_summer", "C_fall")]) / sum(bs) else NA)
  }

  boot <- t(replicate(B, {
    idx <- unlist(rows_of[sample(names(rows_of), replace = TRUE)], use.names = FALSE)
    weights_from_b(fit_nnls(X_free[idx, ], yv[idx])$b)
  }))
  est <- weights_from_b(fit_nnls(X_free, yv)$b)

  tibble(weight   = colnames(boot),
         estimate = est,
         lo95     = apply(boot, 2, quantile, .025, na.rm = TRUE),
         hi95     = apply(boot, 2, quantile, .975, na.rm = TRUE),
         p_zero   = colMeans(boot <= 1e-10, na.rm = TRUE),
         official = c(W_OFF_SEASON, W_OFF_ORG, sum(W_OFF_SEASON[3:4])))
}


# ==============================================================================
# Part D. MIMIC: weights estimated inside the SEM
# ==============================================================================
LABELS <- c(C_winter = "bw", C_spring = "bsp", C_summer = "bsu", C_fall = "bf",
            Diversity = "bd", InvDist = "bi", O_amend = "ba", O_animal = "bn")

WALD_D <- list(
  season_equal    = list(need = c("bw", "bsp", "bsu", "bf"),
                         con  = "bw == bsp; bsp == bsu; bsu == bf"),
  season_official = list(need = c("bw", "bsp", "bsu", "bf"),
                         con  = with(as.list(W_OFF_SEASON / W_OFF_SEASON[["C_winter"]]),
                                     sprintf("bsp == %.6f*bw; bsu == %.6f*bw; bf == %.6f*bw",
                                             C_spring, C_summer, C_fall))),
  season_2group   = list(need = c("bw", "bsp", "bsu", "bf"),
                         con  = "bw == bsp; bsu == bf"),
  org_equal       = list(need = c("ba", "bn"), con = "ba == bn"),
  org_official    = list(need = c("ba", "bn"),
                         con  = sprintf("ba == %.6f*bn",
                                        W_OFF_ORG[["O_amend"]] / W_OFF_ORG[["O_animal"]]))
)

part_D <- function(d, active, sem_spec, sem_all) {
  # Components rescaled to 0-1 for numerical stability; weights are scale-free
  sem_dat <- sem_all %>%
    filter(MGT_combo %in% d$MGT_combo) %>%
    inner_join(comp %>% mutate(across(all_of(PRED), ~ .x / 100)), by = "MGT_combo")

  labs   <- LABELS[active]
  struct <- paste0("SH ~ ", paste0(labs, "*", active, collapse = " + "))
  if (all(C_cols %in% active)) {
    struct <- paste(struct,
      "w_winter := bw  / (bw + bsp + bsu + bf)",
      "w_spring := bsp / (bw + bsp + bsu + bf)",
      "w_summer := bsu / (bw + bsp + bsu + bf)",
      "w_fall   := bf  / (bw + bsp + bsu + bf)",
      "warm_share := (bsu + bf) / (bw + bsp + bsu + bf)", sep = "\n")
  }
  if (all(O_cols %in% active)) struct <- paste(struct, "w_amend := ba / (ba + bn)", sep = "\n")

  mimic <- tryCatch(
    sem(paste(sem_spec$syntax, struct, sep = "\n"), data = sem_dat,
        std.lv = TRUE, missing = "fiml", estimator = "MLR", cluster = "sitecode"),
    error = function(e) { message("MIMIC failed: ", conditionMessage(e)); NULL })
  if (is.null(mimic)) return(NULL)

  estimates <- parameterEstimates(mimic) %>%
    filter(op %in% c("~", ":="), lhs %in% c("SH", "w_winter", "w_spring", "w_summer",
                                            "w_fall", "warm_share", "w_amend")) %>%
    select(lhs, op, rhs, label, est, se, pvalue, ci.lower, ci.upper)

  wald <- imap_dfr(keep(WALD_D, ~ all(.x$need %in% labs)), function(h, nm) {
    w <- tryCatch(lavTestWald(mimic, constraints = h$con), error = function(e) NULL)
    if (is.null(w)) return(tibble(H0 = nm))
    tibble(H0 = nm, stat = w$stat, df = w$df, p = w$p.value)
  })

  # Direct paths from components to first-order factors (score tests).
  # Many significant paths would mean one scalar SHMI is too coarse.
  first_order <- setdiff(lavNames(mimic, "lv"), "SH")
  add_paths   <- as.vector(outer(first_order, active, paste, sep = " ~ "))
  score <- map_dfr(add_paths, function(p) {
    s <- tryCatch(lavTestScore(mimic, add = p)$test, error = function(e) NULL)
    if (is.null(s)) return(tibble(path = p))
    tibble(path = p, X2 = s$X2, df = s$df, p = s$p.value)
  }) %>% mutate(p_bonf = pmin(1, p * n())) %>% arrange(p)

  list(fit = mimic, estimates = estimates, wald = wald, score = score)
}


# ==============================================================================
# Run all parts for one scenario
# ==============================================================================
run_scenario <- function(name, d, target) {
  cat("\n\n", strrep("=", 78), "\n SCENARIO: ", name, "\n", strrep("=", 78), "\n", sep = "")
  set.seed(CFG$seed)

  desc   <- describe_data(d)
  active <- active_components(d)
  cat("Rows:", desc$n_rows, " Sites:", desc$n_sites,
      " Low-summer units (C_summer < 50):", desc$low_summer_n, "\n")
  if (desc$n_sites < 30) message("Few clusters: rely on wild bootstrap p-values.")
  if (length(setdiff(PRED, active)))
    cat("Components without variation (dropped):", setdiff(PRED, active), "\n")
  cat("\nSeasonal correlations:\n"); print(desc$cor)
  cat("\nSeasonal SDs:\n");          print(desc$sd)
  cat("\nLand use x animals present:\n"); print(desc$animals)

  cat("\n--- Part A: Wald tests (small p = reject H0) ---\n")
  A <- part_A(d, active, CFG$B_wild)
  print(A$coef, digits = 3)
  print(A$tests, digits = 3)

  cat("\n--- Part B: leave-site-out CV; gain = RMSE(equal,equal) - RMSE(spec) ---\n")
  B <- part_B(d, CFG$R_cv, CFG$K_cv, CFG$seed)
  print(B$summary, digits = 3)

  cat("\n--- Part C: free weights, site-bootstrap 95% intervals ---\n")
  C <- part_C(d, active, CFG$B_boot)
  print(C, digits = 3)

  cat("\n--- Part D: MIMIC ---\n")
  D <- part_D(d, active, TARG[[target]]$sem_spec, TARG[[target]]$sem_data)
  if (!is.null(D)) {
    print(D$estimates, digits = 3)
    print(D$wald, digits = 3)
    cat("\nDirect-effect score tests, Bonferroni p < .10:\n")
    print(filter(D$score, p_bonf < 0.10), n = Inf, digits = 3)
  }

  list(desc = desc, active = active, A = A, B = B, C = C,
       D = if (!is.null(D)) D[c("estimates", "wald", "score")])
}

RES <- imap(SCENARIOS, function(sc, nm)
  run_scenario(nm, sc$subset(TARG[[sc$target]]$d), sc$target))


# ==============================================================================
# Comparison across scenarios
# ==============================================================================
cat("\n\n", strrep("=", 78), "\n COMPARISON ACROSS SCENARIOS\n", strrep("=", 78), "\n", sep = "")

fmt_ci <- function(est, lo, hi) sprintf("%.3f [%.3f, %.3f]", est, lo, hi)

# Size of each scenario
cmp_size <- imap_dfr(RES, function(r, nm) {
  an <- r$desc$animals
  tibble(scenario = nm, rows = r$desc$n_rows, sites = r$desc$n_sites,
         low_summer = r$desc$low_summer_n, sd_summer = r$desc$sd[["C_summer"]],
         animal_units = if ("TRUE" %in% colnames(an)) sum(an[, "TRUE"]) else 0L)
})
cat("\nScenario sizes:\n"); print(cmp_size)

# Part A: wild-bootstrap p-values for each hypothesis
cmp_A <- imap_dfr(RES, ~ mutate(.x$A$tests, scenario = .y)) %>%
  select(H0, scenario, p_wild_boot) %>%
  pivot_wider(names_from = scenario, values_from = p_wild_boot)
cat("\nPart A - wild bootstrap p (small = reject):\n"); print(cmp_A, digits = 3)

# Part B: CV RMSE gain for every structure, and the best structure
cmp_B <- imap_dfr(RES, ~ mutate(.x$B$summary, scenario = .y)) %>%
  mutate(spec = paste0("season=", season, ", org=", org)) %>%
  select(spec, scenario, gain) %>%
  pivot_wider(names_from = scenario, values_from = gain) %>%
  arrange(desc(.data[[PRIMARY_DATASET]]))
cat("\nPart B - CV RMSE gain over (equal, equal):\n"); print(cmp_B, n = Inf, digits = 3)

cmp_B_best <- imap_dfr(RES, ~ .x$B$summary %>% slice_max(gain, n = 1) %>%
                         mutate(scenario = .y)) %>%
  select(scenario, season, org, gain, win_rate)
cat("\nPart B - best structure per scenario:\n"); print(cmp_B_best, digits = 3)

# Part C: free weights with bootstrap intervals
cmp_C <- imap_dfr(RES, ~ mutate(.x$C, scenario = .y)) %>%
  mutate(value = ifelse(is.na(estimate), "not identified", fmt_ci(estimate, lo95, hi95))) %>%
  select(weight, scenario, value) %>%
  pivot_wider(names_from = scenario, values_from = value)
cat("\nPart C - free weights, estimate [95% site-bootstrap CI]:\n"); print(cmp_C, width = Inf)

# Part D: MIMIC implied weights and Wald tests
cmp_D_w <- imap_dfr(RES, function(r, nm) {
  if (is.null(r$D)) return(tibble())
  r$D$estimates %>% filter(op == ":=") %>%
    transmute(weight = lhs, scenario = nm, value = fmt_ci(est, ci.lower, ci.upper))
}) %>% pivot_wider(names_from = scenario, values_from = value)
cat("\nPart D - MIMIC implied weights, estimate [95% CI]:\n"); print(cmp_D_w, width = Inf)

cmp_D_p <- imap_dfr(RES, function(r, nm) {
  if (is.null(r$D)) return(tibble())
  mutate(r$D$wald, scenario = nm)
}) %>% select(H0, scenario, p) %>%
  pivot_wider(names_from = scenario, values_from = p)
cat("\nPart D - MIMIC Wald p (small = reject):\n"); print(cmp_D_p, digits = 3)

cmp_D_direct <- imap_dfr(RES, function(r, nm) {
  if (is.null(r$D)) return(tibble())
  tibble(scenario = nm,
         direct_paths_bonf_05 = sum(r$D$score$p_bonf < 0.05, na.rm = TRUE),
         tested = sum(!is.na(r$D$score$p)))
})
cat("\nPart D - significant direct effects on first-order factors:\n"); print(cmp_D_direct)


# ==============================================================================
# Save
# ==============================================================================
saveRDS(list(config = CFG, scenarios = names(SCENARIOS), results = RES,
             comparison = list(size = cmp_size, A = cmp_A, B = cmp_B, B_best = cmp_B_best,
                               C = cmp_C, D_weights = cmp_D_w, D_wald = cmp_D_p,
                               D_direct = cmp_D_direct),
             shmi_version = as.character(packageVersion("SHMI")),
             session = sessionInfo()),
        CFG$out)
cat("\nSaved:", CFG$out, "\n")

# ---- Excel export of the comparison tables ----------------------------------
writexl::write_xlsx(list(scenario_sizes = cmp_size, wald_tests = cmp_A, cv_gain = cmp_B,
                         cv_best = cmp_B_best, free_weights = cmp_C,
                         mimic_weights = cmp_D_w, mimic_wald = cmp_D_p,
                         direct_effects = cmp_D_direct),
                    CFG$out_xlsx)
cat("Saved:", CFG$out_xlsx, "\n")
