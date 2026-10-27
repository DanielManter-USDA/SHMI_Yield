################################################################################
# A1c_soil_health_SEM.R
#
# Purpose
#   Choose the measurement model for soil health (SH) that will serve as the
#   calibration target for SHMI.
#
# Firewall rule
#   Nothing in this script looks at management data (Cover, Diversity,
#   InvDist, OrgInput, SHMI). The SH model is chosen on soil data alone, so the
#   target is independent of the index it will later be used to calibrate.
#
# What this script does
#   1. Prepare indicators: outlier removal, climate/texture adjustment,
#      z-scoring. A log-scale version is also built, for the record only.
#   1b. Compare three climate/texture adjustments of the indicators (OLS,
#      mixed model with a site random intercept, penalized GAM with a site
#      random intercept) and test each for climate/texture signal left in the
#      site means. The pre-specified OLS adjustment is kept unless the rule in
#      Step 1b says otherwise.
#   2. Compare the pre-specified candidate models on the full data.
#   3. ESEM diagnostic: which cross-loadings does the CFA force to zero?
#   4. Compare the original model with theory-based alternatives on
#      repeated split halves of sites (fit on unseen sites).
#   4b. Estimability: if the second-order level runs into a boundary on
#      smaller samples, test an equality constraint on the Bio and C/N
#      second-order loadings and adopt it if it costs no fit.
#   5. Evaluate the final model: local misfit, reliability, factor
#      determinacy, measurement invariance by cropping system, and leave-site-out
#      stability of the factor scores.
#   6. Export SH scores and the model for the weight tests and calibration.
#
# One run per dataset in datasets.R (all, noMX, noGZ, noMX_noGZ), each through
# Steps 1-6 on its own units (outlier screening, climate/texture adjustment and
# SEM all fitted within the dataset), plus GAM-adjusted runs of the candidate
# primary datasets (CANDIDATE_DATASETS) as a sensitivity check. A summary compares fit, determinacy,
# stability, invariance and the agreement of SH scores between datasets.
#
# Input : Outputs/semDF_prepped.rds from A1b_Prepare_SEM.R (all units, with
#         `in_scope`; before outlier removal and adjustment).
# Output: Outputs/sem_<dataset>_{scores,model,report}.rds for each dataset,
#         sem_gam_<candidate>_*.rds, sem_summary_all_runs.rds and .xlsx
################################################################################

library(lavaan)
library(dplyr)
library(purrr)
library(tibble)
library(tidyr)
library(here)
library(mgcv)

# ==============================================================================
# 0. Settings
# ==============================================================================
OUT <- function(...) here::here("Outputs", ...)
source(here::here("R", "datasets.R"))

CFG <- list(
  input      = OUT("semDF_prepped.rds"),   # from A1b_Prepare_SEM.R (all units, with in_scope)

  # Runs, all executed in one pass. Each fits the chosen SEM to one dataset
  # with one climate/texture adjustment and exports its own files:
  #   Outputs/<prefix><dataset>_{scores,model,report}.rds
  # Adjustments:
  #   lm_continuous    : OLS, indicator ~ MAT * clay+silt + MAP * clay+silt (pre-specified)
  #   mixed_continuous : as lm_continuous, plus a random intercept for site
  #   gam_mixed        : penalized smooth climate/texture terms + site random intercept
  # RUNS are built from datasets.R below (one per dataset + a GAM check)


  seed       = 2026,
  # TRUE: skip model selection (Steps 2-4b) and fit the specification chosen
  # on the full NAPESHM data (A with equal Bio/CN second-order loadings).
  # FALSE: rerun model selection, in the primary dataset's run only.
  fixed_model = TRUE,
  n_splits   = 20,     # repeated split-half comparisons (Steps 4 and 4b)
  bic_win    = 2/3,    # share of split-halves an alternative must win (Step 4)
  eq_alpha   = 0.05,   # Step 4b: adopt the equality constraint if p > eq_alpha
  K_scores   = 10,     # Step 5e: leave-site-out folds for score stability
  inv_min_n  = 30      # Step 5d: minimum units per cropping system in the invariance test
)
CFG$adjust <- "lm_continuous"   # default for prep_indicators()

# One run per dataset (pre-specified linear adjustment, with the adjustment
# comparison), plus each candidate primary dataset with the GAM adjustment
CFG$RUNS <- c(
  lapply(setNames(names(DATASETS), names(DATASETS)), function(k)
    list(dataset = k, adjust = "lm_continuous", compare_adjust = TRUE, prefix = "sem_")),
  setNames(lapply(CANDIDATE_DATASETS, function(k)
             list(dataset = k, adjust = "gam_mixed", compare_adjust = FALSE, prefix = "sem_gam_")),
           paste0("gam_", CANDIDATE_DATASETS)),
  # Side-by-side test: linear adjustment plus SHAPE soil groups, US sites only
  list(shape_noMX = list(dataset = "noMX", adjust = "lm_shape", compare_adjust = TRUE,
                         compare_methods = c("lm_continuous", "lm_shape"), prefix = "sem_shape_")))
# Run only some of the runs: define SEM_RUNS_ONLY (names of CFG$RUNS) before sourcing.
# By default the soil-group run is left out of the full pipeline.
CFG$RUNS <- if (exists("SEM_RUNS_ONLY")) CFG$RUNS[SEM_RUNS_ONLY] else CFG$RUNS[setdiff(names(CFG$RUNS), "shape_noMX")]

# The 11 soil health indicators used in the SEM
IND <- c("Cmin", "Nmin", "ARYL", "NAG",      # biological activity
         "SOC", "TN", "POXC", "ACE",         # carbon and nitrogen pools
         "SSI", "WAS", "Ksat")               # physical condition


# ==============================================================================
# 1. Indicator preparation
# ==============================================================================

# Sample skewness (used only to decide which indicators the log variant logs)
skew <- function(x) {
  x <- x[is.finite(x)]
  mean((x - mean(x))^3) / sd(x)^3
}

# Indicators logged in the log-scale variant: skewness > 1 and all values > 0.
# (This variant was pre-specified, fit worse, and is kept only for reporting.)
select_log_vars <- function(df, ind = IND) {
  ind[vapply(ind, function(v) {
    x <- df[[v]]; s <- skew(x)
    is.finite(s) && s > 1 && all(x[!is.na(x)] > 0)
  }, logical(1))]
}

# Extreme outliers (beyond 3 x IQR of all units) set to NA. Indicators with
# fewer than 20 observations are left untouched.
replace_outliers <- function(x) {
  if (sum(!is.na(x)) < 20) return(x)
  q   <- quantile(x, c(0.25, 0.75), na.rm = TRUE)
  iqr <- diff(q)
  x[x < q[1] - 3 * iqr | x > q[2] + 3 * iqr] <- NA
  x
}

# Remove climate and texture effects: the fitted climate/texture part of a
# model is subtracted from each indicator. Site random intercepts (mixed
# methods) are NOT removed: they are site deviations the climate and texture
# variables do not explain, including management. Terms from mgcv's
# predict(type = "terms") are centred, so the indicator keeps its level.
#   lm_*       : ordinary least squares; each treatment counts equally, so sites
#                with many treatments dominate the climate coefficients
#   mixed_*    : a site random intercept makes climate effects between-site
#                estimates and stops sites with many treatments dominating
#   gam_mixed  : penalized smooths (shrink towards linear unless the data
#                support curvature) plus the site random intercept
ADJ_METHODS <- c("lm_continuous", "mixed_continuous", "gam_mixed")
ADJ_VARS <- list(lm_shape         = c("Clay_Silt", "site_mean_temp", "site_precip", "shape_grp3"),
                 lm_continuous    = c("Clay_Silt", "site_mean_temp", "site_precip"),
                 mixed_continuous = c("Clay_Silt", "site_mean_temp", "site_precip", "sitecode"),
                 gam_mixed        = c("Clay_Silt", "site_mean_temp", "site_precip", "sitecode"))
ADJ_RHS <- list(
  # lm_continuous + SHAPE soil taxonomy group (S2-3 / S4 / S5), US sites only
  lm_shape         = "site_mean_temp * Clay_Silt + site_precip * Clay_Silt + shape_grp3",
  lm_continuous    = "site_mean_temp * Clay_Silt + site_precip * Clay_Silt",
  mixed_continuous = "site_mean_temp * Clay_Silt + site_precip * Clay_Silt + s(site_f, bs = 're')",
  gam_mixed        = paste("s(site_mean_temp, k = 5) + s(site_precip, k = 5) + s(Clay_Silt, k = 5) +",
                           "ti(site_mean_temp, Clay_Silt, k = 4) + ti(site_precip, Clay_Silt, k = 4) +",
                           "s(site_f, bs = 're')")
)

adjust_indicator <- function(df, var, method = "lm_continuous", return_model = FALSE) {
  keep <- complete.cases(df[, c(var, ADJ_VARS[[method]])])
  d    <- df[keep, ]
  d$site_f <- factor(d$sitecode)
  m <- mgcv::gam(as.formula(paste(var, "~", ADJ_RHS[[method]])), data = d, method = "REML")
  if (return_model) return(m)
  tt  <- stats::predict(m, type = "terms")
  env <- rowSums(tt[, !grepl("site_f", colnames(tt)), drop = FALSE])   # climate + texture only
  out <- rep(NA_real_, nrow(df))
  out[keep] <- d[[var]] - env
  out
}

# Diagnostics for one adjusted indicator, at the SITE level (site means, so
# each site counts once and climate, constant within a site, is not
# pseudo-replicated):
#   between_removed  share of between-site variance removed by the adjustment
#   within_kept      share of within-site variance kept (should be close to 1:
#                    within-site differences are mostly management)
#   left_smooth_*    remaining curved dependence of site means on MAT, MAP and
#                    clay+silt (penalized GAM): deviance explained and smallest p
#   left_other_*     remaining dependence on climate variables NOT used by the
#                    adjustments (moisture index, climatic moisture deficit,
#                    GDD10, elevation): R2 and F-test p
adjustment_diagnostics <- function(raw, adjusted, var) {
  d <- tibble(site = raw$sitecode, y_raw = raw[[var]], y_adj = adjusted,
              MAT = raw$site_mean_temp, MAP = raw$site_precip, CS = raw$Clay_Silt,
              mi = raw$mi, cmd = raw$hargreave_cmd,
              gdd = raw$gdd10, elev = raw$site_elevation) %>%
    filter(!is.na(y_adj))
  d <- d %>% group_by(site) %>% mutate(dev_raw = y_raw - mean(y_raw), dev_adj = y_adj - mean(y_adj)) %>% ungroup()
  within_kept <- var(d$dev_adj) / var(d$dev_raw)
  s <- d %>% group_by(site) %>%
    summarise(y_raw = mean(y_raw), y_adj = mean(y_adj), MAT = first(MAT), MAP = first(MAP),
              CS = mean(CS, na.rm = TRUE), mi = first(mi), cmd = first(cmd),
              gdd = first(gdd), elev = first(elev), .groups = "drop")
  between_removed <- 1 - var(s$y_adj) / var(s$y_raw)
  g  <- mgcv::gam(y_adj ~ s(MAT, k = 5) + s(MAP, k = 5) + s(CS, k = 5), data = s, method = "REML")
  gs <- summary(g)
  o  <- lm(y_adj ~ mi + cmd + gdd + elev, data = s)
  os <- summary(o)
  f  <- os$fstatistic
  tibble(between_removed = round(between_removed, 3),
         within_kept     = round(within_kept, 3),
         left_smooth_dev = round(gs$dev.expl, 3),
         left_smooth_p   = signif(min(gs$s.table[, "p-value"]), 3),
         left_other_r2   = round(os$r.squared, 3),
         left_other_p    = signif(stats::pf(f[1], f[2], f[3], lower.tail = FALSE), 3))
}

# Full preparation pipeline:
#   (optional log) -> outliers -> climate/texture adjustment -> z-score.
# Z-scoring is required: without it the indicators' variances differ by ~10^7
# (e.g. POXC vs TN), which makes lavaan's robust test statistic numerically
# singular. ML estimation is scale-invariant, so fit indices, standardized
# loadings and factor-score correlations are unchanged by it.
prep_indicators <- function(df, log_vars = character(0), ind = IND, method = CFG$adjust) {
  for (v in log_vars) df[[v]] <- log10(df[[v]])
  df <- df %>%
    mutate(across(all_of(ind), replace_outliers))
  for (v in ind) df[[v]] <- adjust_indicator(df, v, method)
  df %>%
    mutate(across(all_of(ind), ~ as.numeric(scale(.x))),
           sitecode = as.character(sitecode))
}


# ==============================================================================
# 2. Model syntax
# ==============================================================================

# Building blocks
MEAS <- list(
  bio      = "SH.Bio  =~ Cmin + Nmin + ARYL + NAG",
  cn       = "SH.CN   =~ SOC + TN + POXC + ACE",
  phys     = "SH.Phys =~ SSI + WAS + Ksat",
  phys_noK = "SH.Phys =~ SSI + WAS",
  hi       = "SH =~ SH.Bio + SH.CN + SH.Phys"      # second-order soil health factor
)

# Residual covariances in the original model
RC_BASE <- c("SOC ~~ TN",      # C:N stoichiometry
             "SOC ~~ WAS",     # organic matter as aggregate binding agent
             "SOC ~~ POXC",    # POXC is a fraction of SOC
             "POXC ~~ ACE")    # both labile organic pools

# Extra residual covariances tested in model D
RC_THEORY <- c("Cmin ~~ POXC",   # labile C supplies respiration substrate
               "Nmin ~~ ACE")    # protein N supplies mineralizable N

mk <- function(...) paste(c(...), collapse = "\n")

# Step 2 candidates. `data` says which prepared data set each model uses.
MODELS_STEP2 <- list(
  A_original   = list(data = "orig", syntax = mk(MEAS$bio, MEAS$cn, MEAS$phys,     MEAS$hi, RC_BASE)),
  B_log        = list(data = "log",  syntax = mk(MEAS$bio, MEAS$cn, MEAS$phys,     MEAS$hi, RC_BASE)),
  C_log_noKsat = list(data = "log",  syntax = mk(MEAS$bio, MEAS$cn, MEAS$phys_noK, MEAS$hi, RC_BASE)),
  D_log_theory = list(data = "log",  syntax = mk(MEAS$bio, MEAS$cn, MEAS$phys,     MEAS$hi, RC_BASE, RC_THEORY)),
  # Correlated first-order factors. Statistically equivalent to B_log: with
  # three first-order factors the second-order level is just-identified, so
  # global fit cannot test whether a single SH factor is adequate.
  E_log_corr3  = list(data = "log",  syntax = mk(MEAS$bio, MEAS$cn, MEAS$phys,              RC_BASE))
)

# Step 4 alternatives, justified by the Cornell CASH grouping, which treats
# ACE protein and active carbon (POXC) as biological indicators. When ACE or
# POXC share a factor with each other, SOC ~~ POXC and POXC ~~ ACE are dropped
# to avoid counting the shared variance twice.
MODELS_STEP4 <- list(
  A_original     = MODELS_STEP2$A_original$syntax,
  G_ACEbio       = mk("SH.Bio  =~ Cmin + Nmin + ARYL + NAG + ACE",
                      "SH.CN   =~ SOC + TN + POXC",
                      MEAS$phys, MEAS$hi, "SOC ~~ TN", "SOC ~~ WAS"),
  H_ACEbio_POXCx = mk("SH.Bio  =~ Cmin + Nmin + ARYL + NAG + ACE + POXC",
                      "SH.CN   =~ SOC + TN + POXC",
                      MEAS$phys, MEAS$hi, "SOC ~~ TN", "SOC ~~ WAS")
)


# Step 4b: original model with equal second-order loadings for Bio and C/N.
# With three first-order factors the second-order loadings are fixed by the
# three factor correlations, e.g. loading(Bio)^2 = r(Bio,CN) * r(Bio,Phys) /
# r(CN,Phys). When Bio and C/N correlate very highly this can exceed 1 in a
# subsample, the Bio disturbance is pushed to zero, and the fit fails to
# converge. With the constraint, loading(Bio)^2 = loading(CN)^2 = r(Bio,CN),
# which stays below 1, and the second-order level gains 1 df.
EQ_SECOND_ORDER <- "SH =~ a*SH.Bio + a*SH.CN + SH.Phys"
MODEL_A_EQ <- mk(MEAS$bio, MEAS$cn, MEAS$phys, EQ_SECOND_ORDER, RC_BASE)


# ==============================================================================
# 3. Fitting helpers
# ==============================================================================

# Fit one model with the settings used throughout:
#   MLR       robust SEs and test statistic (indicators are non-normal)
#   fiml      full-information ML for missing indicator values
#   cluster   site-level clustering (treatments nested in sites)
# If the fit does not converge, it is retried once with bounds = "standard",
# which keeps variances non-negative; the result records whether that
# happened. Returns NULL if both attempts fail.
fit_sem <- function(syntax, data, ...) {
  try_fit <- function(bounds) {
    tryCatch(
      suppressWarnings(sem(syntax, data = data, std.lv = TRUE, missing = "fiml",
                           estimator = "MLR", cluster = "sitecode",
                           bounds = bounds, control = list(iter.max = 20000), ...)),
      error = function(e) NULL)
  }
  ok <- function(f) !is.null(f) && lavInspect(f, "converged")

  f <- try_fit("none")
  if (ok(f)) { attr(f, "bounded") <- FALSE; return(f) }
  f <- try_fit("standard")
  if (ok(f)) { attr(f, "bounded") <- TRUE;  return(f) }
  NULL
}

FIT_NAMES <- c("chisq.scaled", "df", "cfi.robust", "tli.robust",
               "rmsea.robust", "rmsea.ci.lower.robust", "rmsea.ci.upper.robust",
               "srmr", "bic")

# Fit measures as a one-row tibble; NA with converged = FALSE if the fit failed
fit_row <- function(fit) {
  if (is.null(fit)) {
    return(tibble(converged = FALSE, bounded = NA,
                  !!!setNames(as.list(rep(NA_real_, length(FIT_NAMES))), FIT_NAMES)))
  }
  tibble(converged = TRUE, bounded = attr(fit, "bounded"),
         !!!as.list(round(fitMeasures(fit, FIT_NAMES), 3)))
}

# Heywood check: any negative residual variance of an observed indicator
has_heywood <- function(fit) {
  th <- diag(lavInspect(fit, "est")$theta)
  any(th < 0)
}


# ==============================================================================
# Step 2. Compare candidate models on the full data
# ==============================================================================
compare_full <- function(models, data_list) {
  imap_dfr(models, function(m, nm) {
    fit <- fit_sem(m$syntax, data_list[[m$data]])
    tibble(model = nm, data = m$data, fit_row(fit),
           heywood = if (is.null(fit)) NA else has_heywood(fit))
  })
}


# ==============================================================================
# Step 3. ESEM diagnostic
#   Exploratory 3-factor model (geomin rotation) with the residual covariances
#   that have a clear mechanistic basis. Loadings with |value| > .30 outside
#   the CFA pattern point to cross-loadings the CFA forces to zero. Treat the
#   result as a direction, not an estimate, and distrust it if it has a
#   Heywood case.
# ==============================================================================
run_esem <- function(data, rescov = c("SOC ~~ TN", "SOC ~~ WAS"), ind = IND) {
  syntax <- mk(paste0('efa("b")*F1 + efa("b")*F2 + efa("b")*F3 =~ ',
                      paste(ind, collapse = " + ")),
               rescov)
  fit <- tryCatch(
    suppressWarnings(sem(syntax, data = data, missing = "fiml", estimator = "MLR",
                         rotation = "geomin", cluster = "sitecode")),
    error = function(e) { message("ESEM failed: ", conditionMessage(e)); NULL })
  if (is.null(fit)) return(NULL)
  list(fit      = fit,
       loadings = round(lavInspect(fit, "std")$lambda, 2),
       heywood  = has_heywood(fit),
       rmsea    = unname(fitMeasures(fit, "rmsea.robust")))
}


# ==============================================================================
# Step 4. Repeated split-half comparison
#   Sites are split at random into two halves; every model is fit to each
#   half separately. Repeating this n_splits times shows which model fits
#   unseen sites best, and how reliably each model converges. BIC is
#   comparable across models within a split-half because all models are fit
#   to the same rows.
# ==============================================================================
compare_split_half <- function(models, data, n_splits, seed) {
  sites <- unique(data$sitecode)
  map_dfr(seq_len(n_splits), function(r) {
    set.seed(seed + r)
    half_of_site <- setNames(sample(rep(1:2, length.out = length(sites))), sites)
    half_id      <- half_of_site[data$sitecode]
    map_dfr(1:2, function(h) {
      d_h <- data[half_id == h, ]
      imap_dfr(models, function(syntax, nm) {
        tibble(split = r, half = h, model = nm, fit_row(fit_sem(syntax, d_h)))
      })
    })
  })
}

# Summaries of the split-half results
summarise_split_half <- function(cv) {
  by_model <- cv %>%
    group_by(model) %>%
    summarise(converged    = mean(converged),
              bounded      = mean(bounded, na.rm = TRUE),
              cfi.robust   = mean(cfi.robust,   na.rm = TRUE),
              rmsea.robust = mean(rmsea.robust, na.rm = TRUE),
              srmr         = mean(srmr,         na.rm = TRUE),
              .groups = "drop")

  # Head-to-head: among split-halves where every model converged, which one
  # has the lowest BIC?
  bic_wins <- cv %>%
    group_by(split, half) %>%
    filter(all(converged)) %>%
    slice_min(bic, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    count(model, name = "bic_wins") %>%
    mutate(share = bic_wins / sum(bic_wins))

  list(by_model = by_model, bic_wins = bic_wins)
}

# Decision rule, fixed before looking at the results:
#   keep the original model unless an alternative (a) wins BIC in at least
#   `bic_win` of comparable split-halves, (b) converges at least as often,
#   and (c) is no worse on mean robust CFI, RMSEA and SRMR.
choose_model <- function(summ, reference = "A_original", bic_win = 2/3) {
  ref <- summ$by_model %>% filter(model == reference)
  candidates <- summ$by_model %>%
    filter(model != reference) %>%
    left_join(summ$bic_wins, by = "model") %>%
    mutate(share = coalesce(share, 0),
           passes = share >= bic_win &
                    converged    >= ref$converged &
                    cfi.robust   >= ref$cfi.robust &
                    rmsea.robust <= ref$rmsea.robust &
                    srmr         <= ref$srmr)
  winner <- candidates %>% filter(passes) %>% slice_max(share, n = 1, with_ties = FALSE)
  list(chosen = if (nrow(winner) == 1) winner$model else reference,
       table  = candidates)
}



# ==============================================================================
# Step 4b. Estimability of the second-order level
#   Scaled chi-square difference test of the equality constraint on the full
#   data, plus convergence of both versions on the same split halves.
# ==============================================================================
test_equal_second_order <- function(syntax_free, syntax_eq, data, n_splits, seed) {
  f_free <- fit_sem(syntax_free, data)
  f_eq   <- fit_sem(syntax_eq,   data)
  lrt    <- lavTestLRT(f_free, f_eq)
  conv   <- compare_split_half(list(free = syntax_free, equal = syntax_eq),
                               data, n_splits, seed) %>%
    group_by(model) %>%
    summarise(converged = mean(converged), bounded = mean(bounded, na.rm = TRUE),
              .groups = "drop")
  list(lrt = lrt,
       p   = lrt[["Pr(>Chisq)"]][2],
       fit = bind_rows(free = fit_row(f_free), equal = fit_row(f_eq), .id = "version"),
       convergence = conv)
}


# ==============================================================================
# Step 5. Evaluation of the final model
# ==============================================================================

# Residual correlations larger than `cut` (local misfit)
local_misfit <- function(fit, cut = 0.10) {
  rc <- residuals(fit, type = "cor")$cov
  rc[upper.tri(rc, diag = TRUE)] <- NA
  big <- which(abs(rc) > cut, arr.ind = TRUE)
  tibble(var1 = rownames(rc)[big[, 1]], var2 = colnames(rc)[big[, 2]],
         r = round(rc[big], 3)) %>%
    arrange(desc(abs(r)))
}

# Composite reliability (omega), including the second-order SH factor.
# Needs semTools; returns NULL if it is not installed or fails.
reliability <- function(fit) {
  if (!requireNamespace("semTools", quietly = TRUE)) {
    message("semTools not installed: reliability skipped.")
    return(NULL)
  }
  tryCatch(semTools::compRelSEM(fit, higher = "SH"), error = function(e) {
    message("compRelSEM failed: ", conditionMessage(e)); NULL })
}

# Factor determinacy: correlation between each latent variable and its best
# linear predictor from the indicators (complete data). Values below about
# .80 mean the factor scores are a poor stand-in for the latent variable.
factor_determinacy <- function(fit) {
  L     <- lavInspect(fit, "est")$lambda
  Phi   <- lavInspect(fit, "cov.lv")[colnames(L), colnames(L)]
  Sigma <- lavInspect(fit, "cov.ov")[rownames(L), rownames(L)]
  C     <- Phi %*% t(L)                                   # cov(latent, indicators)
  round(sqrt(diag(C %*% solve(Sigma, t(C))) / diag(Phi)), 3)
}

# Multi-group syntax: a shared label such as a* would also force equality
# ACROSS groups, so give each group its own copy of every label (a1, a2, ...).
# Cross-group equality is then controlled only by group.equal.
multigroup_syntax <- function(syntax, n_groups) {
  labs <- unique(regmatches(syntax, gregexpr("\\b[a-zA-Z][a-zA-Z0-9_]*(?=\\*)",
                                             syntax, perl = TRUE))[[1]])
  for (l in labs) {
    vec <- paste0("c(", paste0(l, seq_len(n_groups), collapse = ","), ")")
    syntax <- gsub(paste0("\\b", l, "\\*"), paste0(vec, "*"), syntax, perl = TRUE)
  }
  syntax
}

# Measurement invariance across cropping systems (`system` from
# A1b_Prepare_SEM.R: annual, fallow_reference, perennial): configural (same structure),
# metric (equal loadings), scalar (equal loadings and intercepts).
# SH is pooled across systems when SHMI is calibrated, which presumes at
# least scalar invariance. Common heuristic for non-invariance between
# successive levels: drop in robust CFI > .01 or rise in robust RMSEA > .015.
invariance_by_group <- function(syntax, data, group = "system") {
  data <- data[!is.na(data[[group]]), ]
  n_g  <- length(unique(data[[group]]))
  syn  <- multigroup_syntax(syntax, n_g)
  levels_eq <- list(configural = character(0),
                    metric     = "loadings",
                    scalar     = c("loadings", "intercepts"))
  fits <- lapply(levels_eq, function(eq) {
    tryCatch(suppressWarnings(
      sem(syn, data = data, std.lv = TRUE, missing = "fiml", estimator = "MLR",
          cluster = "sitecode", group = group, group.equal = eq)),
      error = function(e) { message("Invariance fit failed: ", conditionMessage(e)); NULL })
  })
  ok <- vapply(fits, function(f) !is.null(f) && lavInspect(f, "converged"), logical(1))
  if (!all(ok)) return(list(fits = fits, ok = ok))
  list(fits  = fits, ok = ok,
       sizes = table(data[[group]]),
       lrt   = lavTestLRT(fits$configural, fits$metric, fits$scalar),
       fit   = bind_rows(lapply(fits, fit_row), .id = "level"))
}

# Leave-site-out stability of the factor scores: refit the model without a
# fold of sites, predict scores for the held-out sites, and correlate them
# with the full-data scores. Refits use plain ML (point estimates are the
# same as MLR; clustering affects only SEs and tests) and start from the
# full-data estimates, which helps the optimizer without changing the answer.
score_stability <- function(syntax, data, fit_full, K, seed, ind = IND) {
  lv_full <- lavPredict(fit_full, type = "lv")
  full <- tibble(row = lavInspect(fit_full, "case.idx"), as_tibble(lv_full))
  lv_names <- colnames(lv_full)

  sites <- unique(data$sitecode)
  set.seed(seed)
  fold_of_site <- setNames(sample(rep_len(seq_len(K), length(sites))), sites)
  fold <- fold_of_site[data$sitecode]
  has_data <- rowSums(!is.na(data[, ind])) > 0

  oos <- map_dfr(seq_len(K), function(k) {
    tr <- which(fold != k)
    te <- which(fold == k & has_data)
    f <- tryCatch(suppressWarnings(
      sem(syntax, data = data[tr, ], std.lv = TRUE, missing = "fiml",
          estimator = "ML", start = fit_full)), error = function(e) NULL)
    if (is.null(f) || !lavInspect(f, "converged")) return(tibble(fold = k, converged = FALSE))
    p <- tryCatch(lavPredict(f, newdata = data[te, ], type = "lv"), error = function(e) NULL)
    if (is.null(p) || nrow(p) != length(te)) return(tibble(fold = k, converged = TRUE))
    tibble(fold = k, converged = TRUE, row = te, as_tibble(p))
  })

  conv <- oos %>% distinct(fold, converged)
  paired <- oos %>% filter(!is.na(row)) %>%
    select(row, fold, all_of(lv_names)) %>%
    inner_join(full, by = "row", suffix = c("_oos", "_full"))

  stab <- map_dfr(lv_names, function(v) {
    a <- paired[[paste0(v, "_full")]]; b <- paired[[paste0(v, "_oos")]]
    tibble(latent = v, n = sum(complete.cases(a, b)),
           pearson  = round(cor(a, b, use = "complete.obs"), 3),
           spearman = round(cor(a, b, use = "complete.obs", method = "spearman"), 3))
  })
  list(summary = stab, folds_converged = conv, paired = paired)
}


# ==============================================================================
# Step 6. Export
#   Regression factor scores from the final model, raw and min-max scaled to
#   0-100. n_missing (number of unobserved indicators per row) is exported
#   because FIML scores are shrunk more for rows with more missing indicators;
#   the calibration can be repeated on complete rows as a sensitivity check.
# ==============================================================================
export_scores <- function(fit, data, ind = IND) {
  lv  <- lavPredict(fit, type = "lv")
  idx <- lavInspect(fit, "case.idx")
  minmax <- function(x) 100 * (x - min(x, na.rm = TRUE)) / diff(range(x, na.rm = TRUE))

  scores <- data %>%
    select(MGT_combo, system, sitecode) %>%
    mutate(n_missing = rowSums(is.na(data[, ind])))
  for (nm in colnames(lv)) {
    raw <- rep(NA_real_, nrow(data)); raw[idx] <- lv[, nm]
    scores[[paste0(nm, "_raw")]] <- raw
    scores[[nm]] <- minmax(raw)
  }
  scores
}

# ==============================================================================
# Step 6b. Projection object: everything needed to score NEW samples with this
#   run's model, without refitting anything on them: the outlier bounds, the
#   fitted climate/texture adjustment models, the standardization constants of
#   the adjusted indicators, the final SEM, and the 0-100 rescaling constants.
#   Used by A2c_Project_SH.R to score the Mexican samples with the US model.
# ==============================================================================
make_projection <- function(dat_raw, fit, method, ind = IND) {
  bounds <- lapply(ind, function(v) {
    x <- dat_raw[[v]]
    if (sum(!is.na(x)) < 20) return(c(-Inf, Inf))
    q <- quantile(x, c(0.25, 0.75), na.rm = TRUE); iqr <- diff(q)
    unname(c(q[1] - 3 * iqr, q[2] + 3 * iqr))
  })
  names(bounds) <- ind
  d <- dat_raw %>% mutate(across(all_of(ind), replace_outliers))
  models <- setNames(lapply(ind, function(v) adjust_indicator(d, v, method, return_model = TRUE)), ind)
  adj <- d; for (v in ind) adj[[v]] <- adjust_indicator(d, v, method)
  lv  <- lavPredict(fit, type = "lv")
  rng <- apply(lv, 2, range, na.rm = TRUE)
  list(indicators = ind, method = method, adj_vars = ADJ_VARS[[method]], bounds = bounds, models = models,
       center = vapply(ind, function(v) mean(adj[[v]], na.rm = TRUE), numeric(1)),
       sd     = vapply(ind, function(v) sd(adj[[v]], na.rm = TRUE), numeric(1)),
       fit = fit, lv_min = rng[1, ], lv_range = rng[2, ] - rng[1, ])
}

# ==============================================================================
# Step 1b. Compare the climate/texture adjustments
#   (a) per indicator: variance removed between sites, within-site variance
#       kept, and climate/texture signal LEFT OVER in site means (curved terms,
#       and climate variables not used in the adjustment)
#   (b) fit and determinacy of the final model under each adjustment
#   (c) agreement of each method's SH scores with the pre-specified method
# Decision rule (fixed in advance): keep the pre-specified lm_continuous
# adjustment unless its residuals show remaining climate/texture signal
# (left_smooth_p or left_other_p < 0.05/11, Bonferroni over indicators) for
# several indicators that a more robust method removes while keeping
# within-site variance and fitting at least as well. Among adequate methods,
# prefer the simplest.
# ==============================================================================
compare_adjustments <- function(dat_raw, syntax, ind = IND, methods = ADJ_METHODS) {
  # outliers are screened once, identically for all methods
  screened <- dat_raw %>%
    mutate(across(all_of(ind), replace_outliers))

  per_ind <- map_dfr(ind, function(v) {
    raw_diag <- adjustment_diagnostics(screened, screened[[v]], v) %>%
      mutate(indicator = v, method = "none (raw)", .before = 1)
    bind_rows(raw_diag, map_dfr(methods, function(meth) {
      adj <- tryCatch(adjust_indicator(screened, v, meth), error = function(e) NULL)
      if (is.null(adj)) return(tibble(indicator = v, method = meth))
      adjustment_diagnostics(screened, adj, v) %>% mutate(indicator = v, method = meth, .before = 1)
    }))
  })

  fits <- map(setNames(methods, methods), function(meth) {
    d <- tryCatch(prep_indicators(dat_raw, method = meth), error = function(e) NULL)
    if (is.null(d)) return(list(fit = NULL))
    list(data = d, fit = fit_sem(syntax, d))
  })

  model_tab <- imap_dfr(fits, function(x, meth) {
    if (is.null(x$fit)) return(tibble(method = meth, converged = FALSE))
    det <- factor_determinacy(x$fit)
    tibble(method = meth, fit_row(x$fit), det_SH = det[["SH"]], det_min = min(det))
  })

  sh <- map(fits, function(x) {
    if (is.null(x$fit)) return(NULL)
    tibble(row = lavInspect(x$fit, "case.idx"), SH = lavPredict(x$fit, type = "lv")[, "SH"])
  })
  ref <- sh[["lm_continuous"]]
  agree <- imap_dfr(sh, function(x, meth) {
    if (is.null(x) || is.null(ref)) return(tibble(method = meth))
    j <- inner_join(ref, x, by = "row", suffix = c("_ref", "_m"))
    tibble(method = meth, n = nrow(j),
           pearson  = round(cor(j$SH_ref, j$SH_m), 4),
           spearman = round(cor(j$SH_ref, j$SH_m, method = "spearman"), 4))
  })

  list(per_indicator = per_ind, model = model_tab, sh_agreement = agree)
}


# ==============================================================================
# One run: prepare -> (1b) -> (2-4b) -> fit -> evaluate (5) -> export (6)
# ==============================================================================
run_sem <- function(name, run) {
  cat("\n\n", strrep("#", 78), "\n# RUN: ", name, "  (dataset = ", run$dataset, ": ", ds_label(run$dataset),
      ", adjustment = ", run$adjust, ")\n", strrep("#", 78), "\n", sep = "")
  files <- list(scores = OUT(paste0(run$prefix, run$dataset, "_scores.rds")),
                model  = OUT(paste0(run$prefix, run$dataset, "_model.rds")),
                report = OUT(paste0(run$prefix, run$dataset, "_report.rds")))

  # ---- Step 1: data ----------------------------------------------------------
  dat_raw <- readRDS(CFG$input)
  dat_raw <- ds_filter(dat_raw, run$dataset)
  dat_raw <- dat_raw %>% mutate(sitecode = as.character(sitecode))
  needed <- unique(c(IND, "MGT_combo", "sitecode", "system", "in_scope",
                     "Clay_Silt", "site_mean_temp", "site_precip", ADJ_VARS[[run$adjust]]))
  missing_cols <- setdiff(needed, names(dat_raw))
  if (length(missing_cols)) stop("semDF_prepped.rds lacks: ", paste(missing_cols, collapse = ", "))
  cat(nrow(dat_raw), "units at", length(unique(dat_raw$sitecode)), "sites\n")
  print(table(system = dat_raw$system))

  LOG_VARS <- select_log_vars(dat_raw)
  DATA <- list(orig = prep_indicators(dat_raw, method = run$adjust),
               log  = prep_indicators(dat_raw, LOG_VARS, method = run$adjust))

  # ---- Step 1b: compare the climate/texture adjustments ------------------------
  step1b <- NULL
  if (isTRUE(run$compare_adjust)) {
    cat("\n=== Step 1b: climate/texture adjustments (model A, equal loadings) ===\n")
    step1b <- compare_adjustments(dat_raw, MODEL_A_EQ,
                                  methods = if (!is.null(run$compare_methods)) run$compare_methods else ADJ_METHODS)
    cat("\n(a) Per indicator (site-level diagnostics):\n")
    print(step1b$per_indicator, n = Inf, width = Inf)
    cat("\nSummary by method (indicators with remaining signal at p < 0.05/11):\n")
    print(step1b$per_indicator %>%
            group_by(method) %>%
            summarise(between_removed = mean(between_removed), within_kept = mean(within_kept),
                      n_left_smooth = sum(left_smooth_p < 0.05 / length(IND), na.rm = TRUE),
                      n_left_other  = sum(left_other_p  < 0.05 / length(IND), na.rm = TRUE),
                      .groups = "drop"), width = Inf)
    cat("\n(b) Final model under each adjustment:\n")
    print(step1b$model, width = Inf)
    cat("\n(c) Agreement of SH scores with lm_continuous:\n")
    print(step1b$sh_agreement)
  }

  # ---- Steps 2-4b: model selection (optional; primary dataset only) ---------
  step2 <- step3 <- step4b <- step4_raw <- step4_summ <- decision <- NULL
  FINAL_NAME   <- "A_equal_second_order"
  FINAL_SYNTAX <- MODEL_A_EQ
  if (!CFG$fixed_model && name == PRIMARY_DATASET) {
    cat("\n=== Step 2: candidate models ===\n")
    step2 <- compare_full(MODELS_STEP2, DATA); print(step2, width = Inf)
    cat("\n=== Step 3: ESEM diagnostic ===\n")
    step3 <- run_esem(DATA$orig)
    if (!is.null(step3)) { print(step3$loadings); cat("Heywood:", step3$heywood, "\n") }
    cat("\n=== Step 4: split-half comparison ===\n")
    step4_raw  <- compare_split_half(MODELS_STEP4, DATA$orig, CFG$n_splits, CFG$seed)
    step4_summ <- summarise_split_half(step4_raw)
    print(step4_summ$by_model, width = Inf); print(step4_summ$bic_wins)
    decision <- choose_model(step4_summ, bic_win = CFG$bic_win)
    print(decision$table, width = Inf)
    FINAL_NAME <- decision$chosen; FINAL_SYNTAX <- MODELS_STEP4[[FINAL_NAME]]
    if (FINAL_NAME == "A_original") {
      cat("\n=== Step 4b: equal second-order loadings ===\n")
      step4b <- test_equal_second_order(MODELS_STEP4$A_original, MODEL_A_EQ,
                                        DATA$orig, CFG$n_splits, CFG$seed)
      print(step4b$lrt); print(step4b$fit, width = Inf); print(step4b$convergence)
      if (!is.na(step4b$p) && step4b$p > CFG$eq_alpha) {
        FINAL_NAME <- "A_equal_second_order"; FINAL_SYNTAX <- MODEL_A_EQ
      }
    }
  } else {
    cat("\nModel selection skipped; using", FINAL_NAME, "\n")
  }

  # ---- Final model -------------------------------------------------------------
  fit_final <- fit_sem(FINAL_SYNTAX, DATA$orig)
  if (is.null(fit_final)) {
    warning("Run ", name, ": final model did not converge.")
    return(list(name = name, run = run, converged = FALSE))
  }
  print(summary(fit_final, standardized = TRUE, fit.measures = TRUE))

  # ---- Step 5: evaluation ------------------------------------------------------
  cat("\n=== Step 5a: residual correlations |r| > .10 ===\n")
  step5_misfit <- local_misfit(fit_final); print(step5_misfit)

  cat("\n=== Step 5b: composite reliability ===\n")
  step5_rel <- reliability(fit_final); print(step5_rel)

  cat("\n=== Step 5c: factor determinacy ===\n")
  step5_det <- factor_determinacy(fit_final); print(step5_det)

  # Step 5d: measurement invariance between groups, wherever the dataset holds
  # two groups with enough units: cropping system (annual vs perennial) and
  # region (US vs Mexico). Invariance is criterion 1 of the pooling decision.
  run_invariance <- function(group, label) {
    cat("\n=== Step 5d: measurement invariance by", label, "===\n")
    d <- DATA$orig %>%
      mutate(region = ifelse(grepl("^MX", sitecode), "Mexico", "US")) %>%
      group_by(.data[[group]]) %>% filter(n() >= CFG$inv_min_n) %>% ungroup()
    if (length(unique(d[[group]])) < 2) {
      cat("Fewer than two groups with >=", CFG$inv_min_n, "units: test skipped.\n")
      return(list(ok = c(skipped = FALSE)))
    }
    inv <- invariance_by_group(FINAL_SYNTAX, d, group = group)
    if (all(inv$ok)) {
      print(inv$sizes); print(inv$lrt)
      print(inv$fit %>% select(level, cfi.robust, rmsea.robust, srmr), width = Inf)
    } else {
      cat("Not all invariance models converged:\n"); print(inv$ok)
    }
    inv
  }
  step5_inv        <- run_invariance("system", "cropping system (annual vs perennial)")
  step5_inv_region <- run_invariance("region", "region (US vs Mexico)")

  cat("\n=== Step 5e: leave-site-out score stability (", CFG$K_scores, "folds) ===\n")
  step5_stab <- score_stability(FINAL_SYNTAX, DATA$orig, fit_final, CFG$K_scores, CFG$seed)
  print(step5_stab$folds_converged); print(step5_stab$summary)

  # ---- Step 6: export ----------------------------------------------------------
  scores <- export_scores(fit_final, DATA$orig)
  saveRDS(scores, files$scores)
  saveRDS(make_projection(dat_raw, fit_final, run$adjust),
          OUT(paste0(run$prefix, run$dataset, "_projection.rds")))
  saveRDS(list(name = FINAL_NAME, syntax = FINAL_SYNTAX, adjust = run$adjust,
               dataset = run$dataset, estimator = "MLR", missing = "fiml",
               cluster = "sitecode", indicators_standardized = TRUE),
          files$model)
  saveRDS(list(config = CFG, run = run, log_vars_tested = LOG_VARS,
               step1b_adjustment = step1b, step2 = step2,
               step3 = if (!is.null(step3)) step3[c("loadings", "heywood", "rmsea")],
               step4 = list(raw = step4_raw, summary = step4_summ, decision = decision),
               step4b = step4b, final_name = FINAL_NAME, final_fit = fit_final,
               misfit = step5_misfit, reliability = step5_rel, determinacy = step5_det,
               invariance = step5_inv[c("ok", "sizes", "lrt", "fit")],
               invariance_region = step5_inv_region[c("ok", "sizes", "lrt", "fit")],
               score_stability = step5_stab[c("summary", "folds_converged")],
               prepared_data = DATA$orig),
          files$report)
  cat("\nSaved:", unlist(files), sep = "\n  ")

  list(name = name, run = run, converged = TRUE, n = nrow(DATA$orig),
       sites = length(unique(DATA$orig$sitecode)), fit = fit_row(fit_final),
       determinacy = step5_det,
       stability_SH = step5_stab$summary$pearson[step5_stab$summary$latent == "SH"],
       invariance = step5_inv, invariance_region = step5_inv_region,
       step1b = step1b, scores = scores, files = files)
}


# ==============================================================================
# Run everything
# ==============================================================================
RESULTS <- imap(CFG$RUNS, function(run, name) run_sem(name, run))


# ==============================================================================
# Summary across runs
# ==============================================================================
cat("\n\n", strrep("=", 78), "\n SUMMARY ACROSS RUNS\n", strrep("=", 78), "\n", sep = "")

run_tab <- imap_dfr(RESULTS, function(r, name) {
  if (!isTRUE(r$converged)) return(tibble(run = name, converged = FALSE))
  tibble(run = name, dataset = r$run$dataset, adjust = r$run$adjust,
         units = r$n, sites = r$sites,
         cfi.robust = r$fit$cfi.robust, rmsea.robust = r$fit$rmsea.robust, srmr = r$fit$srmr,
         det_SH = r$determinacy[["SH"]], stability_SH = r$stability_SH)
})
cat("\nFit, determinacy and score stability:\n"); print(run_tab, width = Inf)

# Agreement of SH scores between runs, on the units they share
sh_of <- function(name) RESULTS[[name]]$scores %>% select(MGT_combo, SH_raw)
# every pair of datasets, plus each candidate vs its GAM-adjusted run
pairs <- c(combn(names(DATASETS), 2, simplify = FALSE),
           lapply(CANDIDATE_DATASETS, function(k) c(k, paste0("gam_", k))))
agree <- map_dfr(pairs, function(p) {
  if (!all(p %in% names(RESULTS)) || !all(sapply(RESULTS[p], function(r) isTRUE(r$converged))))
    return(tibble(comparison = paste(p, collapse = " vs ")))
  j <- inner_join(sh_of(p[1]), sh_of(p[2]), by = "MGT_combo", suffix = c("_1", "_2")) %>%
    filter(!is.na(SH_raw_1), !is.na(SH_raw_2))
  tibble(comparison = paste(p, collapse = " vs "), shared_units = nrow(j),
         pearson  = round(cor(j$SH_raw_1, j$SH_raw_2), 4),
         spearman = round(cor(j$SH_raw_1, j$SH_raw_2, method = "spearman"), 4))
})
cat("\nAgreement of SH scores between runs (shared units):\n"); print(agree)

# Invariance between cropping systems and between regions, in every dataset
# that holds two groups with enough units
inv_tab <- imap_dfr(RESULTS, function(r, name) {
  bind_rows(lapply(c(system = "invariance", region = "invariance_region"), function(el) {
    inv <- r[[el]]
    if (!isTRUE(r$converged) || is.null(inv$fit)) return(NULL)
    inv$fit %>% select(level, chisq.scaled, df, cfi.robust, rmsea.robust, srmr) %>%
      mutate(delta_cfi = cfi.robust - lag(cfi.robust), delta_rmsea = rmsea.robust - lag(rmsea.robust),
             delta_srmr = srmr - lag(srmr),
             groups = paste(names(inv$sizes), inv$sizes, sep = " = ", collapse = ", "))
  }), .id = "grouping") %>%
    { if (nrow(.)) mutate(., run = name, .before = 1) else . }
})
if (nrow(inv_tab)) { cat("\nInvariance between groups (criterion 1):\n"); print(inv_tab, width = Inf) }

# Climate/texture adjustment comparison, per dataset
adj_tab <- imap_dfr(RESULTS, function(r, name) {
  if (is.null(r$step1b)) return(NULL)
  r$step1b$per_indicator %>%
    group_by(method) %>%
    summarise(between_removed = mean(between_removed), within_kept = mean(within_kept),
              n_left_smooth = sum(left_smooth_p < 0.05 / length(IND), na.rm = TRUE),
              n_left_other  = sum(left_other_p  < 0.05 / length(IND), na.rm = TRUE),
              .groups = "drop") %>%
    left_join(r$step1b$sh_agreement %>% select(method, r_with_lm = pearson), by = "method") %>%
    mutate(dataset = name, .before = 1)
})
if (nrow(adj_tab)) { cat("\nClimate/texture adjustments by dataset:\n"); print(adj_tab, n = Inf, width = Inf) }

saveRDS(list(runs = run_tab, agreement = agree, invariance = inv_tab, adjustment = adj_tab,
             config = CFG), OUT("sem_summary_all_runs.rds"))
writexl::write_xlsx(list(sem_runs = run_tab, sh_agreement = agree,
                         invariance = if (nrow(inv_tab)) inv_tab else tibble(note = "not run"),
                         adjustment = if (nrow(adj_tab)) adj_tab else tibble(note = "not run")),
                    OUT("sem_summary_all_runs.xlsx"))
cat("\nSaved:", OUT("sem_summary_all_runs.rds"), "and .xlsx\n")
