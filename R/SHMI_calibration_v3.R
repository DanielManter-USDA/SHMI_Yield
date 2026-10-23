################################################################################
# SHMI_calibration_v3.R
#
# Final calibration of SHMI 1.2.0 weights against the SEM-derived soil health
# score (SH), using the weight structure supported by SHMI_weight_tests.R:
#
#   Seasons  two groups: cool (winter, spring) and warm (summer, fall).
#            Equal weights within each group; the warm share is estimated.
#   OrgInput amendments and animals weighted equally (0.5 / 0.5); the data
#            could not distinguish them.
#   Pillars  Cover, Diversity, InvDist, OrgInput weights estimated.
#            Diversity stays in the index structure a priori and is always
#            computed and reported. With div_floor = 0 its weight is
#            estimated freely (and may be zero); with div_floor = f > 0 its
#            pillar share is fixed at f and the other weights are estimated
#            given that constraint.
#   Scope    management units whose rotation includes at least one annual
#            crop (Inputs/in_scope_v120.rds).
#   Fixed    max_div = 10 (above the 99th percentile of effective species),
#            max_stir = 317.5 (annual STIR 80 = NRCS 345 reduced-till ceiling
#            at the top of the reduced-tillage TI classes), ti_rep = "max".
#   Selected hill in {0, 1, 2}, by nested leave-site-out CV with a
#            one-standard-error rule that keeps the a priori value (1)
#            unless another is clearly better.
#
# Why this is a regression: with the structure above,
#   SHMI = w_cov*[(1-ws)*Cool + ws*Warm] + w_div*Diversity + w_inv*InvDist + w_org*Org
# where Cool = (C_winter + C_spring)/2, Warm = (C_summer + C_fall)/2,
# Org = (O_amend + O_animal)/2. So SH ~ a + beta*SHMI is a linear regression
# of SH on (Cool, Warm, Diversity, InvDist, Org) with non-negative
# coefficients; every weight is a ratio of those coefficients. Non-negative
# least squares (NNLS) finds the optimum exactly.
#
# Validation (repeated nested leave-site-out CV):
#   M0  a priori weights (seasons, organic inputs, pillars all equal), hill 1
#   M1  weights estimated, hill fixed at 1
#   M2  weights estimated, hill selected inside each outer training set
#
# Outputs (prefix CFG$out_prefix):
#   results.rds            everything, for the manuscript
#   final_settings.rds     the settings to freeze in build_shmi()  (+ MD5)
#   oof_predictions.csv    out-of-fold predictions (use for honest R2)
#   optimum_shmi_scores.rds  SH scores joined with final SHMI indices
################################################################################

library(SHMI)
library(dplyr)
library(tidyr)
library(purrr)
library(tibble)
library(nnls)

stopifnot(packageVersion("SHMI") >= "1.2.0")

# ==============================================================================
# 0. Settings
# ==============================================================================
# Absolute paths, so the script does not depend on the working directory
BASE <- here::here("Outputs")   # all inputs to this script are outputs of earlier steps

CFG_default <- list(
  inputs         = file.path(BASE, "shmi_inputs_v120.rds"),      # from Prepare_SEM.R
  scores         = file.path(BASE, "sem_noMX_noGZ_scores.rds"),  # from SH_SEM_refinement_clean.R (overridden by the runner)
  in_scope       = file.path(BASE, "in_scope_v120.rds"),         # from Prepare_SEM.R
  target         = "SH",

  dist_meth      = "STIR",
  ti_rep         = "max",
  max_stir       = 317.5,
  max_div        = 8,       # gamma cap (SHMI 1.2.0 decision)
  animal_presence = "span",
  hill_grid      = c(0, 1, 2),     # used only for diversity = "gamma"

  # Test cases (development):
  #   diversity     : "gamma"    plant-day Shannon over the window (current SHMI)
  #                   "alpha"    within-year diversity (plant-day weighted mean)
  #                   "beta"     rotation diversity (gamma / alpha)
  #                   "combined" Alpha and Beta as two predictors; w_alpha estimated
  #                              (fixed at combined_w when a Diversity floor is set)
  #                   "combined_norm" as combined, but rotation diversity normalized to
  #                              the unit's own possible year-to-year change (Beta_norm)
  #                   "gamma_split" gamma's own within-year and rotation parts (log scale,
  #                              raw beta), each doubled so that w = 0.5 is exactly gamma:
  #                              Diversity = min(w * 2A + (1 - w) * 2B, 100)
  #   cover_seasons : "calendar" winter+spring vs summer+fall
  #                   "climate"  dormant vs active months (monthly mean temperature
  #                              below / at or above active_temp at the site)
  diversity      = "gamma",
  # Organic Inputs: "either" (SHMI 1.2.0) = share of years with an amendment or
  # animals; "weighted" = 0.5 x amendment share + 0.5 x animal share
  orginput_rule  = "either",
  combined_w     = 0.5,
  max_alpha      = 8,
  max_beta       = 4,
  cover_seasons  = "calendar",
  climate_file   = file.path(BASE, "site_climate_monthly.rds"),
  # Rule for an ACTIVE month under cover_seasons = "climate":
  #   "temperature"   mean temperature >= active_temp
  #   "temp_moisture" also not dry: precipitation (mm) >= 2 x mean temperature
  #                   (C), the Bagnouls-Gaussen dry-month criterion
  active_rule    = "temperature",
  # Sites whose months are never treated as dry (known irrigation); irrigation
  # is not recorded in the NAPESHM data, so empty by default
  irrigated_sites = character(0),
  # SH scores covering all units (incl. perennial-only), used only to evaluate
  # the calibrated index outside the calibration set
  explore_scores = file.path(BASE, "sem_all_scores.rds"),
  active_temp    = 5,
  hill_apriori   = 1,

  # Units to calibrate on. Default: the in-scope set. Sensitivity runs, e.g.
  #   function(d) filter(d, system == "annual")           (without fallow references)
  #   function(d) filter(d, MGT_combo %in% readRDS(CFG$in_scope), !grepl("^MX", sitecode))
  subset         = function(d) filter(d, MGT_combo %in% readRDS(file.path(BASE, "in_scope_v120.rds"))),

  # Diversity pillar share. 0 = estimate freely (first run). After the first
  # run, choose f inside the bootstrap 95% interval using a rule written down
  # beforehand, then rerun with that value.
  div_floor      = 0,
  floor_grid     = c(0, 0.05, 0.10, 0.15),   # Section 4b: CV cost of each floor
  # Warm-season share of the Cover weight. NA: estimated from the data.
  # A number in [0, 1]: fixed a priori (e.g. 0.75), only the rest estimated.
  warm_fix       = NA,
  warm_grid      = c(0.5, 0.6, 0.75, 0.9, 1),  # Section 4c: CV cost of each fixed warm share
  # Section 4d (combined Diversity types, with a floor): CV profile over the
  # spatial/temporal split, w_alpha = share of Diversity on within-year diversity
  w_alpha_grid   = c(0, 0.25, 0.5, 0.75, 1),
  # Section 4e: Organic Inputs folded into Diversity (Option A). Amendments
  # and animals each add up to x effective species, in proportion to the share
  # of years they occur (events per year do not matter):
  #   gamma* = gamma + x_amend * p_amend + x_animal * p_animal
  #   Diversity* = 100 * min(log(gamma*) / log(max_div), 1)
  # and the Organic Inputs pillar is dropped. org_x_grid gives the x values,
  # crossed for amendments and animals; NULL skips the profile.
  org_x_grid     = NULL,
  # Section 4f: the time-weighted version. A year with an organic amendment
  # counts like a plant species growing for D days that year (same for
  # animals): the amendment and animal pseudo-species enter the plant-day
  # pool, with D x (years with the input) days each, and
  #   Diversity* = 100 * min(H(plants + pseudo-species) / log(max_div), 1).
  # org_days_grid gives the D values, crossed for amendments and animals.
  org_days_grid  = NULL,
  # Section 4f also runs a "points" version: each year with the input adds a
  # fixed bonus on the 0-100 Diversity scale, whatever the plant diversity:
  #   Diversity* = min(100, 100 * log(D_q) / log(cap) + c_amend * p_amend + c_animal * p_animal)
  # org_points_grid gives the bonus values (points per input every year).
  org_points_grid = NULL,
  # Hill orders tested in Section 4f: 0 = richness, 1 = Shannon (plant-days),
  # 2 = inverse Simpson. The four-pillar reference is refitted for each order.
  org_hill_grid  = 1,
  # Section 4g: Diversity cap profile (four pillars). Diversity is rescored as
  # 100 * min(log(gamma) / log(cap), 1) for each cap; reports out-of-sample R2,
  # the share of units at 100 and whether any unit reaches 100. NULL skips it.
  maxdiv_grid    = NULL,
  # caps (effective species scoring 100) to test alongside k; each cap also gets
  # its own four-pillar reference. NULL = CFG$max_div only.
  org_maxdiv_grid = NULL,
  site_weighting = FALSE,   # TRUE: each site contributes equally to the fit

  K_inner        = 5,
  K_outer        = 5,
  R_outer        = getOption("shmi.R_outer", 20),   # Run_pipeline.R raises these
  B_boot         = getOption("shmi.B_boot", 1000),
  boot_reselect  = TRUE,    # re-select hill inside each bootstrap draw
  seed           = 2026,
  out_dir        = file.path(BASE, "Calibration"),
  out_prefix     = "calib_v3_"
)

# Scenario runs (Run_calibration_scenarios.R) set CFG_OVERRIDE before
# sourcing this script; a plain run uses the defaults above.
CFG <- utils::modifyList(CFG_default,
                         if (exists("CFG_OVERRIDE")) CFG_OVERRIDE else list())
dir.create(CFG$out_dir, showWarnings = FALSE, recursive = TRUE)
CFG$out_prefix <- file.path(CFG$out_dir, CFG$out_prefix)
cat("\n########## Calibration run:", basename(CFG$out_prefix), "##########\n")

# Diversity and Cover test cases other than the current definitions use a
# single hill value (hill only affects gamma)
if (CFG$diversity != "gamma") CFG$hill_grid <- CFG$hill_apriori

# Predictor columns. Cool/Warm are the two Cover groups: calendar
# (winter+spring / summer+fall) or climate (dormant / active months).
# Diversity: one column, or Alpha and Beta as two columns ("combined", free fit).
DIV_COMBINED <- CFG$diversity %in% c("combined", "combined_norm", "gamma_split")
DIV_FREE  <- if (DIV_COMBINED) c("DivA", "DivB") else "Diversity"
pred_cols <- function(floor) c("Cool", "Warm", "Org",
                               if (floor > 0) "Diversity" else DIV_FREE, "InvDist")
PRED <- pred_cols(CFG$div_floor)

# ==============================================================================
# 1. Data
# ==============================================================================
shmi_inputs <- readRDS(CFG$inputs)
if (!identical(attr(shmi_inputs, "rotation_window"), "calendar"))
  stop("Inputs were not prepared with the SHMI 1.2.0 calendar window.")

sh <- readRDS(CFG$scores) %>%
  transmute(MGT_combo, system, sitecode = as.character(sitecode), n_missing,
            y = .data[[CFG$target]]) %>%
  filter(!is.na(y))

# Nonlinear settings for a given hill value
struct_settings <- function(hill) {
  list(dist_meth = CFG$dist_meth, max_stir = CFG$max_stir, ti_rep = CFG$ti_rep,
       hill = hill, max_div = CFG$max_div, animal_presence = CFG$animal_presence,
       orginput_rule = CFG$orginput_rule)
}

# Cover groups from monthly components (climate seasons)
climate_cover <- function() {
  clim <- readRDS(CFG$climate_file)
  mo   <- compute_cover_monthly(shmi_inputs$crop, shmi_inputs$rot_bounds)
  site <- sub("^NAPESHM_([A-Za-z0-9]+)_.*$", "\\1", mo$MGT_combo)
  tav  <- as.matrix(clim[match(site, clim$sitecode), sprintf("tavg_%02d", 1:12)])
  prc  <- as.matrix(clim[match(site, clim$sitecode), sprintf("prec_%02d", 1:12)])
  if (anyNA(tav)) stop("Missing monthly temperature for: ",
                       paste(unique(site[rowSums(is.na(tav)) > 0]), collapse = ", "))
  active <- tav >= CFG$active_temp
  if (CFG$active_rule == "temp_moisture") {
    dry <- prc < 2 * tav                                   # Bagnouls-Gaussen dry month
    dry[site %in% CFG$irrigated_sites, ] <- FALSE          # irrigation removes water limits
    active <- active & !dry
  }
  P <- as.matrix(mo[, sprintf("plant_%02d", 1:12)]); D <- as.matrix(mo[, sprintf("days_%02d", 1:12)])
  a_days <- rowSums(D * active);  d_days <- rowSums(D * !active)
  warm <- ifelse(a_days > 0, 100 * rowSums(P * active)  / a_days, NA_real_)
  cool <- ifelse(d_days > 0, 100 * rowSums(P * !active) / d_days, NA_real_)
  # a site without dormant (or active) months: Cover is the other group alone
  tibble(MGT_combo = mo$MGT_combo,
         Cool = coalesce(cool, warm), Warm = coalesce(warm, cool),
         n_active_months = rowSums(active))
}

# Design matrix for one hill value: Cover groups, Org, Diversity column(s), InvDist
make_design <- function(hill) {
  comp <- compute_shmi_components(shmi_inputs, settings = struct_settings(hill))
  d <- comp %>%
    transmute(MGT_combo,
              Cool = (C_winter + C_spring) / 2,
              Warm = (C_summer + C_fall) / 2,
              Org  = if (CFG$orginput_rule == "either") O_any else (O_amend + O_animal) / 2,
              Diversity, InvDist)
  if (CFG$cover_seasons == "climate")
    d <- d %>% select(-Cool, -Warm) %>% inner_join(climate_cover(), by = "MGT_combo")
  if (!is.null(CFG$org_x_grid) || !is.null(CFG$org_days_grid) || !is.null(CFG$org_points_grid) ||
      !is.null(CFG$maxdiv_grid)) {
    d <- d %>%
      left_join(compute_diversity_components(shmi_inputs$crop, shmi_inputs$rot_bounds,
                                             max_gamma = CFG$max_div) %>%
                  select(MGT_combo, gamma_eff = gamma), by = "MGT_combo") %>%
      left_join(comp %>% transmute(MGT_combo, p_amend = O_amend / 100, p_animal = O_animal / 100),
                by = "MGT_combo") %>%
      mutate(gamma_eff = coalesce(gamma_eff, 1),
             p_amend = coalesce(p_amend, 0), p_animal = coalesce(p_animal, 0)) %>%
      left_join(shmi_inputs$rot_bounds %>%
                  transmute(MGT_combo, n_years = rot_end_yr - rot_start_yr + 1) %>% distinct(),
                by = "MGT_combo") %>%
      mutate(years_amend = round(p_amend * n_years), years_animal = round(p_animal * n_years))
  }
  if (CFG$diversity != "gamma") {
    dc <- compute_diversity_components(shmi_inputs$crop, shmi_inputs$rot_bounds,
                                       max_alpha = CFG$max_alpha, max_beta = CFG$max_beta,
                                       max_gamma = CFG$max_div)
    src <- switch(CFG$diversity,
                  combined_norm = c(DivA = "Alpha", DivB = "Beta_norm"),
                  gamma_split   = c(DivA = "Gamma_alpha_part", DivB = "Gamma_beta_part"),
                  c(DivA = "Alpha", DivB = "Beta"))
    mult <- if (CFG$diversity == "gamma_split") 2 else 1     # w = 0.5 reproduces gamma
    d <- d %>% select(-Diversity) %>%
      left_join(dc %>% transmute(MGT_combo, DivA = mult * .data[[src[["DivA"]]]],
                                 DivB = mult * .data[[src[["DivB"]]]]), by = "MGT_combo") %>%
      mutate(across(c(DivA, DivB), ~ coalesce(.x, 0)),
             Diversity = case_when(CFG$diversity == "alpha" ~ DivA,
                                   CFG$diversity == "beta"  ~ DivB,
                                   TRUE ~ pmin(CFG$combined_w * DivA + (1 - CFG$combined_w) * DivB, 100)))
  }
  d
}

base <- sh %>% CFG$subset() %>% arrange(MGT_combo)

X_list <- setNames(lapply(CFG$hill_grid, function(h) {
  d <- base %>% inner_join(make_design(h), by = "MGT_combo")
  if (anyNA(d[, unique(c(PRED, pred_cols(0), "Diversity"))])) stop("NA predictors for hill = ", h)
  d
}), paste0("hill_", CFG$hill_grid))

stopifnot(all(sapply(X_list, function(d) identical(d$MGT_combo, X_list[[1]]$MGT_combo))))
j0   <- match(CFG$hill_apriori, CFG$hill_grid)
y    <- X_list[[1]]$y
site <- X_list[[1]]$sitecode
n    <- length(y)
cat("Analysis set:", n, "units at", length(unique(site)), "sites\n")

cat("\nCorrelations among calibration predictors (hill =", CFG$hill_apriori, "):\n")
print(round(cor(X_list[[j0]][, unique(c(pred_cols(0), "Diversity"))]), 2))
cat("Test case: diversity =", CFG$diversity, "| cover seasons =", CFG$cover_seasons, "\n")


# ==============================================================================
# 2. Model fitting
# ==============================================================================
site_w <- function(clus) {
  if (!CFG$site_weighting) return(rep(1, length(clus)))
  as.numeric(1 / table(clus)[clus])
}

# NNLS with a free intercept (weighted centering is exact for LS with an
# intercept). With floor = f > 0 the Diversity pillar share is fixed at f:
# its coefficient must equal f/(1-f) times the sum of the other coefficients,
# so c*Diversity is added to every other predictor and the remaining weights
# are estimated given the floor. floor = 0 is the ordinary free fit.
#
# floor = 0: every predictor, including Diversity (or Alpha and Beta), is
# estimated freely. (Before this fix Diversity was left out at floor = 0,
# which forced its weight to zero.)
fit_nnls <- function(d, w, floor = CFG$div_floor, warm = CFG$warm_fix, drop = character(0)) {
  cols <- setdiff(pred_cols(floor), drop)
  free <- if (floor > 0) setdiff(cols, "Diversity") else cols
  Xf   <- as.matrix(d[, free])
  fixed_warm <- !is.null(warm) && !is.na(warm)
  if (fixed_warm) {
    # One Cover column with the seasonal split fixed: (1 - warm) * Cool + warm * Warm
    Xf <- cbind(Cover = (1 - warm) * d$Cool + warm * d$Warm,
                Xf[, setdiff(free, c("Cool", "Warm")), drop = FALSE])
  }
  cc <- if (floor > 0) floor / (1 - floor) else 0
  X  <- if (floor > 0) Xf + cc * d$Diversity else Xf
  yy <- d$y
  xm <- colSums(X * w) / sum(w); ym <- sum(yy * w) / sum(w)
  sw <- sqrt(w)
  f  <- nnls::nnls(sweep(X, 2, xm) * sw, (yy - ym) * sw)
  b_f <- setNames(f$x, colnames(Xf))
  a   <- ym - sum(b_f * xm)
  if (fixed_warm) {
    b_f <- c(Cool = (1 - warm) * b_f[["Cover"]], Warm = warm * b_f[["Cover"]],
             b_f[setdiff(names(b_f), "Cover")])
  }
  b <- if (floor > 0) c(b_f, Diversity = cc * sum(b_f))[cols] else b_f[cols]
  list(a = a, b = b)
}

# M0: a priori equal weights. Seasons equal -> Cool and Warm 1/2 each within
# Cover; pillars 1/4 each. Calibrated to SH by a simple linear fit.
DIR0 <- c(Cool = 1/8, Warm = 1/8, Org = 1/4, Diversity = 1/4, InvDist = 1/4)
fit_m0 <- function(d, w) {
  s <- as.vector(as.matrix(d[, names(DIR0)]) %*% DIR0)
  m <- lm(d$y ~ s, weights = w)
  list(a = unname(coef(m)[1]), b = unname(coef(m)[2]) * DIR0)
}

predict_lin <- function(fit, d) fit$a + as.vector(as.matrix(d[, names(fit$b)]) %*% fit$b)

# Coefficients -> build_shmi() settings
coef_to_settings <- function(b, hill) {
  div_b <- b[intersect(names(b), c("Diversity", "DivA", "DivB"))]
  beta <- c(cover = b[["Cool"]] + b[["Warm"]], diversity = sum(div_b),
            invdist = b[["InvDist"]], orginput = if ("Org" %in% names(b)) b[["Org"]] else 0)
  w_alpha <- if (all(c("DivA", "DivB") %in% names(b))) {
    if (sum(b[c("DivA", "DivB")]) > 0) b[["DivA"]] / sum(b[c("DivA", "DivB")]) else NA_real_
  } else if (DIV_COMBINED) CFG$combined_w else NA_real_
  if (sum(beta) <= 0) return(NULL)                     # SHMI carries no signal
  wp <- beta / sum(beta)
  ws <- if (beta[["cover"]] > 0) b[["Warm"]] / beta[["cover"]] else NA_real_
  ws_use <- if (is.na(ws)) 0.5 else ws                 # irrelevant when w_cover = 0
  settings <- c(struct_settings(hill),
                list(w_winter = (1 - ws_use) / 2, w_spring = (1 - ws_use) / 2,
                     w_summer = ws_use / 2,       w_fall   = ws_use / 2,
                     w_amend  = 0.5,              w_animal = 0.5,
                     w_cover  = wp[["cover"]],    w_diversity = wp[["diversity"]],
                     w_invdist = wp[["invdist"]], w_orginput  = wp[["orginput"]],
                     clip_to_rotation = TRUE))
  list(settings = settings, warm_share = ws, slope = sum(beta), w_alpha = w_alpha)
}

WEIGHT_NAMES <- c("warm_share", "w_cover", "w_diversity", "w_invdist", "w_orginput", "w_alpha")
weights_vec <- function(st) {
  if (is.null(st)) return(setNames(rep(NA_real_, length(WEIGHT_NAMES)), WEIGHT_NAMES))
  c(warm_share = st$warm_share, unlist(st$settings[WEIGHT_NAMES[2:5]]), w_alpha = st$w_alpha)
}


# ==============================================================================
# 3. Cross-validation machinery
# ==============================================================================
fold_by_cluster <- function(clus, K, seed) {
  u <- unique(clus)
  set.seed(seed)
  unname(setNames(sample(rep_len(seq_len(K), length(u))), u)[clus])
}

# Out-of-fold RMSE of NNLS for hill index j on rows idx
cv_hill <- function(j, idx, clus, K, seed) {
  fold <- fold_by_cluster(clus, K, seed)
  d    <- X_list[[j]][idx, ]
  pred <- rep(NA_real_, length(idx)); fold_rmse <- numeric(0)
  for (k in sort(unique(fold))) {
    tr <- fold != k
    f  <- fit_nnls(d[tr, ], site_w(clus[tr]))
    pred[!tr] <- predict_lin(f, d[!tr, ])
    fold_rmse <- c(fold_rmse, sqrt(mean((d$y[!tr] - pred[!tr])^2)))
  }
  list(rmse = sqrt(mean((d$y - pred)^2)), se = sd(fold_rmse) / sqrt(length(fold_rmse)))
}

# Hill selection, one-standard-error rule defaulting to the a priori value
select_hill <- function(idx, clus, K, seed) {
  tab <- map_dfr(seq_along(X_list), function(j) {
    cv <- cv_hill(j, idx, clus, K, seed)
    tibble(j = j, hill = CFG$hill_grid[j], rmse = cv$rmse, se = cv$se)
  })
  best  <- tab[which.min(tab$rmse), ]
  j_sel <- if (tab$rmse[j0] <= best$rmse + best$se) j0 else best$j
  list(j = j_sel, table = tab)
}


# ==============================================================================
# 4. Nested, repeated leave-site-out evaluation of M0 / M1 / M2
# ==============================================================================
cat("\n=== Nested CV:", CFG$R_outer, "repeats x", CFG$K_outer, "outer site folds ===\n")

outer <- lapply(seq_len(CFG$R_outer), function(r) {
  of <- fold_by_cluster(site, CFG$K_outer, CFG$seed + r)
  lapply(seq_len(CFG$K_outer), function(k) {
    tr <- which(of != k); te <- which(of == k)
    w  <- site_w(site[tr])

    f0  <- fit_m0(X_list[[j0]][tr, ], w)
    f1  <- fit_nnls(X_list[[j0]][tr, ], w)
    sel <- select_hill(tr, site[tr], CFG$K_inner, CFG$seed + 1000 * r + k)
    f2  <- fit_nnls(X_list[[sel$j]][tr, ], w)

    list(
      pred = tibble(rep_id = r, fold = k, row = te, y = y[te],
                    M0 = predict_lin(f0, X_list[[j0]][te, ]),
                    M1 = predict_lin(f1, X_list[[j0]][te, ]),
                    M2 = predict_lin(f2, X_list[[sel$j]][te, ])),
      sel  = tibble(rep_id = r, fold = k, hill = CFG$hill_grid[sel$j],
                    !!!as.list(weights_vec(coef_to_settings(f2$b, CFG$hill_grid[sel$j]))))
    )
  })
}) %>% unlist(recursive = FALSE)

oof     <- bind_rows(map(outer, "pred"))
sel_log <- bind_rows(map(outer, "sel"))

metrics <- oof %>%
  pivot_longer(c(M0, M1, M2), names_to = "model", values_to = "pred") %>%
  group_by(rep_id, model) %>%
  summarise(RMSE = sqrt(mean((y - pred)^2)),
            MAE  = mean(abs(y - pred)),
            R2   = 1 - sum((y - pred)^2) / sum((y - mean(y))^2),
            .groups = "drop")

cat("\nOut-of-sample performance across repeats - mean [2.5%, 97.5%]:\n")
print(metrics %>% group_by(model) %>%
        summarise(across(c(RMSE, MAE, R2),
                         ~ sprintf("%.3f [%.3f, %.3f]", mean(.x),
                                   quantile(.x, .025), quantile(.x, .975)))))

# Repeats share data, so report distributions and win rates, not t-tests
gains <- metrics %>% select(rep_id, model, RMSE) %>%
  pivot_wider(names_from = model, values_from = RMSE) %>%
  transmute(rep_id, weights_gain = M0 - M1, hill_gain = M1 - M2, total_gain = M0 - M2)
cat("\nRMSE reduction (positive = better):\n")
print(gains %>% summarise(across(-rep_id, list(mean = mean, win_rate = ~ mean(.x > 0)))),
      width = Inf)

cat("\nHill selected in outer training sets:\n")
print(count(sel_log, hill))

# Out-of-fold predictions averaged over repeats, and R2 within cropping system
oof_avg <- oof %>% group_by(row) %>%
  summarise(y = first(y), across(c(M0, M1, M2), mean), .groups = "drop") %>%
  mutate(MGT_combo = X_list[[1]]$MGT_combo[row],
         system    = X_list[[1]]$system[row],
         sitecode  = site[row])
cat("\nOut-of-fold R2 within cropping system:\n")
print(oof_avg %>% group_by(system) %>%
        summarise(n = n(),
                  across(c(M0, M1, M2), ~ 1 - sum((y - .x)^2) / sum((y - mean(y))^2)),
                  .groups = "drop"))


# ==============================================================================
# 4b. Cost of fixing the Diversity share (same outer folds as Section 4)
#     floor = 0 reproduces M1. Use this table, together with the bootstrap
#     interval for w_diversity from a div_floor = 0 run, to apply your
#     pre-written rule for choosing the floor.
# ==============================================================================
cat("\n=== Cross-validated cost of a fixed Diversity share (hill =", CFG$hill_apriori, ") ===\n")
floor_cost <- map_dfr(CFG$floor_grid, function(fl) map_dfr(seq_len(CFG$R_outer), function(r) {
  of <- fold_by_cluster(site, CFG$K_outer, CFG$seed + r)
  p  <- numeric(n)
  for (k in seq_len(CFG$K_outer)) {
    tr <- of != k
    f  <- fit_nnls(X_list[[j0]][tr, ], site_w(site[tr]), floor = fl)
    p[!tr] <- predict_lin(f, X_list[[j0]][!tr, ])
  }
  tibble(floor = fl, rep = r, RMSE = sqrt(mean((y - p)^2)),
         R2 = 1 - sum((y - p)^2) / sum((y - mean(y))^2))
}))
floor_table <- floor_cost %>%
  group_by(floor) %>%
  summarise(RMSE = mean(RMSE), R2 = mean(R2), .groups = "drop") %>%
  mutate(RMSE_increase = RMSE - RMSE[floor == 0]) %>%
  left_join(map_dfr(CFG$floor_grid, ~ as_tibble(as.list(weights_vec(coef_to_settings(
    fit_nnls(X_list[[j0]], site_w(site), floor = .x)$b, CFG$hill_apriori)))) %>%
      mutate(floor = .x)), by = "floor")
print(floor_table, width = Inf, digits = 3)

# ==============================================================================
# 4c. Cost of fixing the warm-season share of Cover (at the chosen floor)
#     Same folds as above; the other weights are re-estimated at each value.
# ==============================================================================
cat("\n=== Section 4c: cross-validated cost of a fixed warm-season share ===\n")
warm_cost <- map_dfr(c(NA, CFG$warm_grid), function(ws) map_dfr(seq_len(CFG$R_outer), function(r) {
  of <- fold_by_cluster(site, CFG$K_outer, CFG$seed + r)
  p  <- numeric(n)
  for (k in seq_len(CFG$K_outer)) {
    tr <- of != k
    f  <- fit_nnls(X_list[[j0]][tr, ], site_w(site[tr]), floor = CFG$div_floor, warm = ws)
    p[!tr] <- predict_lin(f, X_list[[j0]][!tr, ])
  }
  tibble(warm_fixed = ws, rep = r, RMSE = sqrt(mean((y - p)^2)),
         R2 = 1 - sum((y - p)^2) / sum((y - mean(y))^2))
}))
warm_table <- warm_cost %>%
  group_by(warm_fixed) %>%
  summarise(RMSE = mean(RMSE), R2 = mean(R2), .groups = "drop") %>%
  mutate(RMSE_vs_estimated = RMSE - RMSE[is.na(warm_fixed)]) %>%
  left_join(map_dfr(c(NA, CFG$warm_grid), ~ as_tibble(as.list(weights_vec(coef_to_settings(
    fit_nnls(X_list[[j0]], site_w(site), floor = CFG$div_floor, warm = .x)$b, CFG$hill_apriori)))) %>%
      mutate(warm_fixed = .x)), by = "warm_fixed") %>%
  arrange(!is.na(warm_fixed), warm_fixed)
cat("(warm_fixed = NA: warm share estimated from the data)\n")
print(warm_table, width = Inf, digits = 3)

# ==============================================================================
# 4d. Spatial vs temporal diversity: CV profile over the split w_alpha
#     (combined Diversity types with a floor only). Diversity is rebuilt as
#     w * Alpha + (1 - w) * Beta(_norm) for each w; same folds as above.
#     A flat profile means separating the two kinds of diversity neither
#     helps nor hurts SHMI.
# ==============================================================================
w_alpha_table <- NULL
if (DIV_COMBINED && CFG$div_floor > 0) {
  cat("\n=== Section 4d: profile over the spatial/temporal split of Diversity ===\n")
  with_w <- function(d, w) { d$Diversity <- pmin(w * d$DivA + (1 - w) * d$DivB, 100); d }
  w_alpha_table <- map_dfr(CFG$w_alpha_grid, function(wa) {
    Xw <- with_w(X_list[[j0]], wa)
    cvr <- map_dfr(seq_len(CFG$R_outer), function(r) {
      of <- fold_by_cluster(site, CFG$K_outer, CFG$seed + r)
      p  <- numeric(n)
      for (k in seq_len(CFG$K_outer)) {
        tr <- of != k
        f  <- fit_nnls(Xw[tr, ], site_w(site[tr]), floor = CFG$div_floor)
        p[!tr] <- predict_lin(f, Xw[!tr, ])
      }
      tibble(RMSE = sqrt(mean((y - p)^2)), R2 = 1 - sum((y - p)^2) / sum((y - mean(y))^2))
    })
    st <- coef_to_settings(fit_nnls(Xw, site_w(site), floor = CFG$div_floor)$b, CFG$hill_apriori)
    tibble(w_alpha = wa, RMSE = mean(cvr$RMSE), R2 = mean(cvr$R2),
           w_cover = st$settings$w_cover, w_invdist = st$settings$w_invdist,
           w_orginput = st$settings$w_orginput, warm_share = st$warm_share)
  })
  cat("(w_alpha = 0: rotation only; 1: within-year only)\n")
  print(w_alpha_table, width = Inf, digits = 3)
}

# ==============================================================================
# 4e. Organic Inputs folded into Diversity: CV profile over k (effective
#     species added by an amendment or animals in a year). Three pillars
#     (Cover, Diversity*, Inverse Disturbance) against the four-pillar model,
#     on the same folds, at the run's Diversity floor.
# ==============================================================================
org_k_table <- NULL
if (!is.null(CFG$org_x_grid) && CFG$diversity == "gamma") {
  cat("\n=== Section 4e: amendments and animals as effective species in Diversity ===\n")
  caps <- if (is.null(CFG$org_maxdiv_grid)) CFG$max_div else CFG$org_maxdiv_grid
  div_score <- function(d, xa, xn, cap)
    100 * pmin(log(d$gamma_eff + xa * d$p_amend + xn * d$p_animal) / log(cap), 1)
  cv_r2 <- function(Xk, drop) {
    r <- map_dfr(seq_len(CFG$R_outer), function(rr) {
      of <- fold_by_cluster(site, CFG$K_outer, CFG$seed + rr)
      p  <- numeric(n)
      for (kk in seq_len(CFG$K_outer)) {
        tr <- of != kk
        f  <- fit_nnls(Xk[tr, ], site_w(site[tr]), floor = CFG$div_floor, drop = drop)
        p[!tr] <- predict_lin(f, Xk[!tr, ])
      }
      tibble(RMSE = sqrt(mean((y - p)^2)), R2 = 1 - sum((y - p)^2) / sum((y - mean(y))^2))
    })
    c(RMSE = mean(r$RMSE), R2 = mean(r$R2))
  }
  weights_of <- function(Xk, drop) {
    st <- coef_to_settings(fit_nnls(Xk, site_w(site), floor = CFG$div_floor, drop = drop)$b,
                           CFG$hill_apriori)
    if (is.null(st)) return(c(w_cover = NA, w_diversity = NA, w_invdist = NA, w_orginput = NA, warm_share = NA))
    c(unlist(st$settings[c("w_cover", "w_diversity", "w_invdist", "w_orginput")]), warm_share = st$warm_share)
  }
  spread <- function(dv) c(div_max = max(dv), div_p99 = unname(quantile(dv, 0.99)),
                           share_at_cap = mean(dv >= 99.95), share_zero = mean(dv <= 0.05))
  X0 <- X_list[[j0]]
  grid <- tidyr::expand_grid(x_amend = CFG$org_x_grid, x_animal = CFG$org_x_grid)
  org_k_table <- map_dfr(caps, function(cap) {
    # four-pillar reference: gamma Diversity scored against this cap, Organic Inputs separate
    Xr <- X0; Xr$Diversity <- div_score(Xr, 0, 0, cap)
    ref <- tibble(max_div = cap, model = "four pillars (Organic Inputs separate)",
                  x_amend = NA_real_, x_animal = NA_real_,
                  !!!as.list(c(cv_r2(Xr, character(0)), weights_of(Xr, character(0)), spread(Xr$Diversity))))
    three <- purrr::pmap_dfr(grid, function(x_amend, x_animal) {
      Xk <- X0; Xk$Diversity <- div_score(Xk, x_amend, x_animal, cap)
      tibble(max_div = cap, model = "three pillars (inputs in Diversity)",
             x_amend = x_amend, x_animal = x_animal,
             !!!as.list(c(cv_r2(Xk, "Org"), weights_of(Xk, "Org"), spread(Xk$Diversity))))
    })
    bind_rows(ref, three) %>% mutate(R2_vs_four_pillars = R2 - R2[1],
                                     spans_0_100 = div_max >= 99.95)
  })
  print(org_k_table, n = Inf, width = Inf, digits = 3)
}

# ==============================================================================
# 4f. Time-weighted version: amendment and animal years as plant-day pseudo-
#     species (D days per year with the input). Three pillars against the
#     four-pillar reference for each cap, same folds.
# ==============================================================================
org_days_table <- NULL
if ((!is.null(CFG$org_days_grid) || !is.null(CFG$org_points_grid)) && CFG$diversity == "gamma") {
  cat("\n=== Section 4f: organic inputs in Diversity (days and points versions, by Hill order) ===\n")
  cap <- CFG$max_div
  X0 <- X_list[[j0]]
  spd <- compute_species_days(shmi_inputs$crop, shmi_inputs$rot_bounds)
  spd_list <- split(spd$days, spd$MGT_combo)
  plant <- lapply(X0$MGT_combo, function(u) if (is.null(spd_list[[u]])) numeric(0) else spd_list[[u]])
  hill_num <- function(d, q) {
    d <- d[d > 0]; if (!length(d)) return(1)
    p <- d / sum(d)
    if (q == 0) length(d) else if (q == 1) exp(-sum(p * log(p))) else 1 / sum(p^2)
  }
  score_of <- function(D) 100 * pmin(log(pmax(D, 1)) / log(cap), 1)
  D_days <- function(Da, Dn, q) vapply(seq_len(nrow(X0)), function(i)
    hill_num(c(plant[[i]], Da * X0$years_amend[i], Dn * X0$years_animal[i]), q), numeric(1))
  # check: Hill order 1 without inputs reproduces the Gamma score used elsewhere
  chk <- max(abs(score_of(D_days(0, 0, 1)) - X0$Diversity))
  cat("Check, plant-days only (q = 1) vs Gamma score: max |difference| =", signif(chk, 3), "\n")

  cv_r2 <- function(Xk, drop) {
    r <- map_dfr(seq_len(CFG$R_outer), function(rr) {
      of <- fold_by_cluster(site, CFG$K_outer, CFG$seed + rr)
      p  <- numeric(n)
      for (kk in seq_len(CFG$K_outer)) {
        tr <- of != kk
        f  <- fit_nnls(Xk[tr, ], site_w(site[tr]), floor = CFG$div_floor, drop = drop)
        p[!tr] <- predict_lin(f, Xk[!tr, ])
      }
      tibble(RMSE = sqrt(mean((y - p)^2)), R2 = 1 - sum((y - p)^2) / sum((y - mean(y))^2))
    })
    c(RMSE = mean(r$RMSE), R2 = mean(r$R2))
  }
  weights_of <- function(Xk, drop) {
    st <- coef_to_settings(fit_nnls(Xk, site_w(site), floor = CFG$div_floor, drop = drop)$b,
                           CFG$hill_apriori)
    if (is.null(st)) return(c(w_cover = NA, w_diversity = NA, w_invdist = NA, w_orginput = NA, warm_share = NA))
    c(unlist(st$settings[c("w_cover", "w_diversity", "w_invdist", "w_orginput")]), warm_share = st$warm_share)
  }
  spread <- function(dv) c(div_max = max(dv), div_p99 = unname(quantile(dv, 0.99)),
                           share_at_cap = mean(dv >= 99.95), share_zero = mean(dv <= 0.05))
  fit_row <- function(dv, drop, ...) {
    Xk <- X0; Xk$Diversity <- dv
    tibble(..., !!!as.list(c(cv_r2(Xk, drop), weights_of(Xk, drop), spread(dv))))
  }

  org_days_table <- map_dfr(CFG$org_hill_grid, function(q) {
    base <- score_of(D_days(0, 0, q))           # plants only, Hill order q
    ref <- fit_row(base, character(0), hill = q, mode = "reference",
                   model = "four pillars (Organic Inputs separate)",
                   x_amend = NA_real_, x_animal = NA_real_)
    days <- if (is.null(CFG$org_days_grid)) NULL else
      purrr::pmap_dfr(tidyr::expand_grid(x_amend = CFG$org_days_grid, x_animal = CFG$org_days_grid),
        function(x_amend, x_animal)
          fit_row(score_of(D_days(x_amend, x_animal, q)), "Org", hill = q, mode = "days",
                  model = "three pillars (inputs in Diversity)", x_amend = x_amend, x_animal = x_animal))
    pts <- if (is.null(CFG$org_points_grid)) NULL else
      purrr::pmap_dfr(tidyr::expand_grid(x_amend = CFG$org_points_grid, x_animal = CFG$org_points_grid),
        function(x_amend, x_animal)
          fit_row(pmin(base + x_amend * X0$p_amend + x_animal * X0$p_animal, 100), "Org",
                  hill = q, mode = "points", model = "three pillars (inputs in Diversity)",
                  x_amend = x_amend, x_animal = x_animal))
    bind_rows(ref, days, pts) %>%
      mutate(R2_vs_four_pillars = R2 - ref$R2, spans_0_100 = div_max >= 99.95, max_div = cap)
  })
  print(org_days_table, n = Inf, width = Inf, digits = 3)
}

# ==============================================================================
# 4g. Diversity cap profile (four pillars): the cap sets which effective number
#     of species scores 100. Same folds; the other weights re-estimated.
# ==============================================================================
maxdiv_table <- NULL
if (!is.null(CFG$maxdiv_grid) && CFG$diversity == "gamma") {
  cat("\n=== Section 4g: Diversity cap profile ===\n")
  X0 <- X_list[[j0]]
  maxdiv_table <- map_dfr(CFG$maxdiv_grid, function(cap) {
    Xc <- X0; Xc$Diversity <- 100 * pmin(log(Xc$gamma_eff) / log(cap), 1)
    cvr <- map_dfr(seq_len(CFG$R_outer), function(rr) {
      of <- fold_by_cluster(site, CFG$K_outer, CFG$seed + rr)
      p  <- numeric(n)
      for (kk in seq_len(CFG$K_outer)) {
        tr <- of != kk
        f  <- fit_nnls(Xc[tr, ], site_w(site[tr]), floor = CFG$div_floor)
        p[!tr] <- predict_lin(f, Xc[!tr, ])
      }
      tibble(RMSE = sqrt(mean((y - p)^2)), R2 = 1 - sum((y - p)^2) / sum((y - mean(y))^2))
    })
    st <- coef_to_settings(fit_nnls(Xc, site_w(site), floor = CFG$div_floor)$b, CFG$hill_apriori)
    tibble(max_div = cap, RMSE = mean(cvr$RMSE), R2 = mean(cvr$R2),
           w_cover = st$settings$w_cover, w_diversity = st$settings$w_diversity,
           w_invdist = st$settings$w_invdist, w_orginput = st$settings$w_orginput,
           div_mean = mean(Xc$Diversity), div_max = max(Xc$Diversity),
           share_at_100 = mean(Xc$Diversity >= 99.95), share_zero = mean(Xc$Diversity <= 0.05),
           spans_0_100 = max(Xc$Diversity) >= 99.95,
           cover_share_at_100 = mean(X0$Cool >= 99.95 & X0$Warm >= 99.95),
           invdist_share_at_100 = mean(X0$InvDist >= 99.95))
  })
  print(maxdiv_table, width = Inf, digits = 3)
}


# ==============================================================================
# 5. Final calibration on all data
# ==============================================================================
sel_final <- select_hill(seq_len(n), site, CFG$K_inner, CFG$seed)
jstar     <- sel_final$j
hill_star <- CFG$hill_grid[jstar]
fit_final <- fit_nnls(X_list[[jstar]], site_w(site))
final     <- coef_to_settings(fit_final$b, hill_star)
if (is.null(final)) stop("All pillar coefficients are zero in the final fit.")

cat("\n=== Final calibration (hill =", hill_star, ") ===\n")
cat("Hill profile (full-data inner CV):\n"); print(sel_final$table)
print(round(weights_vec(final), 4))
zero <- names(fit_final$b)[fit_final$b <= 1e-10]
if (length(zero)) cat("Predictors with zero coefficient:", paste(zero, collapse = ", "), "\n")
if (CFG$div_floor > 0) {
  cat("Diversity pillar share fixed at", CFG$div_floor, "(div_floor); other weights",
      "estimated given this constraint.\n")
} else if ("Diversity" %in% zero) {
  cat("Diversity weight is 0 in the free fit. Diversity is still computed and",
      "reported; see Section 4b and the bootstrap interval to set div_floor.\n")
}

# Round trip: build_shmi() with the final settings must reproduce the NNLS
# predictor. Only possible for definitions build_shmi() implements
# (calendar seasons and gamma Diversity); other test cases are skipped.
pkg <- NULL; err <- NA_real_
if (CFG$cover_seasons == "calendar" && CFG$diversity == "gamma") {
  pkg <- build_shmi(shmi_inputs, settings = final$settings, expert_mode = TRUE)
  shmi_pkg <- pkg$indicator_df$SHMI[match(X_list[[jstar]]$MGT_combo, pkg$indicator_df$MGT_combo)]
  err <- max(abs(fit_final$a + final$slope * shmi_pkg - predict_lin(fit_final, X_list[[jstar]])))
  cat("Package round-trip max |error|:", signif(err, 3), "\n")
  if (err > 1e-6) warning("build_shmi(final settings) does not reproduce the NNLS predictor.")
} else {
  cat("Round trip skipped: this test case (diversity =", CFG$diversity, ", cover =",
      CFG$cover_seasons, ") is not yet implemented in build_shmi().\n")
}


# ==============================================================================
# 6. Site-cluster bootstrap of the final weights
# ==============================================================================
cat("\n=== Site bootstrap: B =", CFG$B_boot,
    if (CFG$boot_reselect) "(hill re-selected in each draw)" else "", "===\n")
u_sites <- unique(site)
rows_of <- split(seq_len(n), site)

boot <- map_dfr(seq_len(CFG$B_boot), function(b) {
  set.seed(CFG$seed + 10000 + b)
  draw <- sample(u_sites, replace = TRUE)
  idx  <- unlist(rows_of[draw], use.names = FALSE)
  clus <- rep(paste0(draw, "_", seq_along(draw)), times = lengths(rows_of[draw]))
  j    <- if (CFG$boot_reselect) select_hill(idx, clus, CFG$K_inner, CFG$seed + b)$j else jstar
  f    <- fit_nnls(X_list[[j]][idx, ], site_w(clus))
  tibble(b = b, hill = CFG$hill_grid[j],
         !!!as.list(weights_vec(coef_to_settings(f$b, CFG$hill_grid[j]))))
})

weight_summary <- boot %>%
  pivot_longer(all_of(WEIGHT_NAMES), names_to = "weight", values_to = "value") %>%
  group_by(weight) %>%
  summarise(boot_median = median(value, na.rm = TRUE),
            lo95   = quantile(value, .025, na.rm = TRUE),
            hi95   = quantile(value, .975, na.rm = TRUE),
            p_zero = mean(value <= 1e-10, na.rm = TRUE),
            .groups = "drop") %>%
  mutate(estimate = weights_vec(final)[weight], .after = weight) %>%
  arrange(match(weight, WEIGHT_NAMES))
cat("\nFinal weights with site-bootstrap 95% intervals:\n")
print(weight_summary, digits = 3)
cat("\nHill selected across bootstrap draws:\n"); print(count(boot, hill))


# ==============================================================================
# 7. Freeze and save
# ==============================================================================
settings_file <- paste0(CFG$out_prefix, "final_settings.rds")
saveRDS(final$settings, settings_file)
md5 <- unname(tools::md5sum(settings_file))

if (!(CFG$cover_seasons == "calendar" && CFG$diversity == "gamma"))
  cat("\nNOTE: test case (diversity =", CFG$diversity, ", cover =", CFG$cover_seasons,
      "). The block below is indicative only: build_shmi() does not yet implement it.",
      "Warm share = share of Cover weight on", ifelse(CFG$cover_seasons == "climate",
      "active months;", "summer+fall;"), "w_alpha =", round(final$w_alpha, 3), "\n")
cat("\n=== Settings to freeze as official in build_shmi() ===\n")
s <- final$settings
cat(sprintf(paste0(
  "    # seasonal cover (two groups; warm share = %.4f)\n",
  "    w_winter = %.4f, w_spring = %.4f, w_summer = %.4f, w_fall = %.4f,\n",
  "    # diversity\n    hill = %d, max_div = %g,\n",
  "    # disturbance\n    dist_meth = \"EPA\", max_stir = %g, ti_rep = \"%s\",\n",
  "    # organic inputs\n    w_amend = %.4f, w_animal = %.4f, animal_presence = \"%s\",\n",
  "    # evaluation window\n    clip_to_rotation = TRUE,\n",
  "    # pillars\n    w_cover = %.4f, w_diversity = %.4f, w_invdist = %.4f, w_orginput = %.4f\n"),
  final$warm_share, s$w_winter, s$w_spring, s$w_summer, s$w_fall,
  as.integer(s$hill), s$max_div, s$max_stir, s$ti_rep,
  s$w_amend, s$w_animal, s$animal_presence,
  s$w_cover, s$w_diversity, s$w_invdist, s$w_orginput))
cat("\nNote: official scoring uses dist_meth = \"EPA\"; the calibration used STIR",
    "with max_stir =", s$max_stir, "\n")
cat("MD5 of", settings_file, ":", md5, "\n",
    "Record this checksum before analysing the independent yield data.\n")

write.csv(oof_avg, paste0(CFG$out_prefix, "oof_predictions.csv"), row.names = FALSE)

# Predictions of the final calibrated index for ALL units, including those
# outside the calibration set (e.g. Mexican and perennial units), with SH from
# the calibration target and from the all-units target. Used to explore how
# the calibrated weights transfer; correlations and within-site rankings are
# scale-free, so the two targets can be used side by side.
X_all <- make_design(hill_star)
pred_all <- X_all %>%
  transmute(MGT_combo, pred = predict_lin(fit_final, X_all),
            sitecode = sub("^NAPESHM_([A-Za-z0-9]+)_.*$", "\\1", MGT_combo),
            region = ifelse(grepl("^MX", sitecode), "Mexico", "US"),
            in_calibration = MGT_combo %in% X_list[[1]]$MGT_combo) %>%
  left_join(readRDS(CFG$scores) %>% transmute(MGT_combo, SH_target = .data[[CFG$target]]),
            by = "MGT_combo")
if (file.exists(CFG$explore_scores)) {
  pred_all <- pred_all %>%
    left_join(readRDS(CFG$explore_scores) %>%
                transmute(MGT_combo, system, SH_all = .data[[CFG$target]]), by = "MGT_combo")
}
write.csv(pred_all, paste0(CFG$out_prefix, "all_unit_predictions.csv"), row.names = FALSE)

if (!is.null(pkg)) {
  optimum <- readRDS(CFG$scores) %>% left_join(pkg$indicator_df, by = "MGT_combo")
  saveRDS(optimum, paste0(CFG$out_prefix, "optimum_shmi_scores.rds"))
}

saveRDS(list(config = CFG[setdiff(names(CFG), "subset")],
             final_settings = final$settings, final_md5 = md5,
             final_fit = fit_final, warm_share = final$warm_share, hill = hill_star,
             hill_profile = sel_final$table, round_trip_error = err,
             div_floor = CFG$div_floor, floor_table = floor_table,
             warm_fix = CFG$warm_fix, warm_table = warm_table, w_alpha_table = w_alpha_table,
             org_k_table = org_k_table, org_days_table = org_days_table,
             maxdiv_table = maxdiv_table,
             n_units = n, n_sites = length(unique(site)),
             metrics = metrics, gains = gains, selection_log = sel_log,
             bootstrap = boot, weight_summary = weight_summary,
             shmi_version = as.character(packageVersion("SHMI")),
             session = sessionInfo()),
        paste0(CFG$out_prefix, "results.rds"))

cat("\nSaved:", paste0(CFG$out_prefix, c("results.rds", "final_settings.rds",
                                        "oof_predictions.csv", "optimum_shmi_scores.rds")),
    sep = "\n  ")
