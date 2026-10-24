################################################################################
# SHMI_coarse_records.R
#
# Does SHMI keep its skill when the management record is only as detailed as
# remote-sensing and public data products (annual crop type, cover-crop presence,
# tillage class)? Real pixels cannot be used: NAPESHM and the long-term experiments
# are small plots, below a 30-m pixel. Instead the full records are degraded step
# by step (python/coarse_records/degrade_records.py) and rescored with the
# published weights, which are not refitted.
#
# Scenarios (cumulative):
#   D0  full record (the paper's scores)
#   D1  typical dates: every crop's dates set to the site's typical day of year
#   D2  + cover crop yes/no: one generic cover crop per gap between cash crops
#   D3  + tillage class: each crop interval's passes replaced by one pass at the
#         typical intensity of its class (no-till, reduced, conventional)
#   D4  + main crop only: mixtures, companion crops and intercrops dropped
#   D5  + no organic amendments or grazing (not visible from space)
#   D6  + last 2 years of record only (NAPESHM only)
# Noise (NAPESHM: CFG_C$n_noise seeds each; yield: CFG_C$yield_noise_seeds, 1 seed, 3 for J30):
#   J7, J14, J30  every crop and tillage date shifted at random by up to 7, 14, 30 days
#                 (each record on its own: replicate plots of a treatment no longer share dates)
#   S7, S14, S30  as J, but one shift per date and site: records on the same day move together
#   M10, M20      D3 with 10% or 20% of cover-crop gaps and tillage classes misclassified
#
# Measures
#   1. Agreement with the full-record score (NAPESHM units): RMS difference and mean
#      bias in SHMI points, median within-site Spearman rho, RMS difference by part.
#   2. Skill against soil health (NAPESHM, US units): out-of-fold R2, between- and
#      within-site r and the share of sites ranked correctly, on the calibration's
#      folds (10-fold leave-sites-out, repeated), with the PUBLISHED weights; only the
#      line from SHMI to lv_SH is fitted on the training sites, as in the calibration.
#      Secondary: the weights refitted on the coarse records (refit = TRUE).
#   3. Yield (D0-D5, and the noise scenarios): SHMI_yield_main.R rerun with each site scored from its degraded
#      workbook (SCORE_WB_DIR_OVERRIDE) and everything else (yields, N, climate) from
#      the full record: the within-site slope (Fig 1 model) and SHMI agreement.
#
# RUN: from the project (here::here()), after SHMI_calibration_plus.R (Part A) and
#      SHMI_yield_main.R (Part B; its site runs are reused). The coarse workbooks are written
#      with python (openpyxl) unless already in Outputs/Coarse/ (prebuilt: then python is not needed).
#      Settings: CFG_COARSE_OVERRIDE <- list(...) before sourcing (see CFG_C).
# Outputs (Outputs/Coarse/): coarse_agreement.csv, coarse_napeshm.csv,
#      coarse_napeshm_repeats.csv, coarse_yield.csv, Fig_coarse_records.png,
#      degrade_log.csv, the degraded workbooks, yield_<scenario>_agcros/ (main runs)
################################################################################
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(tibble); library(ggplot2); library(here); library(SHMI)
})
setwd(here::here())

CFG_C <- list(
  napeshm_wb   = here::here("Inputs", "NAPESHM_rebuilt.xlsx"),
  final        = here::here("Outputs", "Final_plus", "final_results.rds"),
  climate      = here::here("Outputs", "site_climate_monthly.rds"),
  treatments   = here::here("Inputs", "napeshm", "treatments.csv"),
  yield_wb_dir = here::here("Inputs", "yield_workbooks"),
  main_script  = here::here("R", "SHMI_yield_main.R"),
  script       = here::here("python", "coarse_records", "degrade_records.py"),
  out_dir      = here::here("Outputs", "Coarse"),
  python       = "",       # python executable; "" = python3, else python, on the PATH
  R            = NULL,     # CV repeats; NULL = the calibration's (100)
  n_noise      = 5,        # repeats (seeds) of each noise scenario (the prebuilt workbooks cover 5)
  refit        = TRUE,     # also refit the weights on the coarse records
  run_napeshm  = TRUE,
  napeshm_scenarios = NULL, # e.g. c("J7", "J14", "J30"): rerun only these (results merged with the saved ones)
  run_yield    = TRUE,
  yield_scenarios = NULL,  # e.g. c("D3", "D4", "D5"): run only these; the others are read from their saved runs
  # noise scenarios for yield (each seed is one run of SHMI_yield_main.R; prebuilt workbooks in
  # Outputs/Coarse/yield_workbooks/<scenario>_s0#/; finished seeds are read back, not rerun)
  yield_noise  = c("J7", "J14", "J30", "S7", "S14", "S30", "M10", "M20"),
  yield_noise_seeds = list(default = 1, J30 = 1:3),   # seeds per scenario (a vector applies to all)
  reuse_runs   = TRUE,     # read back noise seeds already run (FALSE: rerun them all, as run_all.Rmd does)
  rebuild      = FALSE     # rewrite degraded workbooks that already exist
)
if (exists("CFG_COARSE_OVERRIDE")) CFG_C <- utils::modifyList(CFG_C, CFG_COARSE_OVERRIDE)
dir.create(CFG_C$out_dir, recursive = TRUE, showWarnings = FALSE)
OUTC <- function(...) file.path(CFG_C$out_dir, ...)

SCEN <- tribble(
  ~scenario, ~ops,                                          ~label,                         ~yield,
  "D0",      "",                                            "Full record",                  TRUE,
  "D1",      "dates",                                       "Typical dates",                TRUE,
  "D2",      "dates,cover",                                 "+ cover crop yes/no",          TRUE,
  "D3",      "dates,cover,tillage",                         "+ tillage class",              TRUE,
  "D4",      "dates,cover,tillage,maincrop",                "+ main crop only",             TRUE,
  "D5",      "dates,cover,tillage,maincrop,noamend",        "+ no amendments or grazing",   TRUE,
  "D6",      "dates,cover,tillage,maincrop,noamend,last:2", "+ last 2 years only",          FALSE)
NOISE <- tribble(
  ~scenario, ~ops,                                ~label,
  "J7",      "jitter:7",                          "Dates shifted up to 7 days",
  "J14",     "jitter:14",                         "Dates shifted up to 14 days",
  "J30",     "jitter:30",                         "Dates shifted up to 30 days",
  "S7",      "sjitter:7",                         "Dates shifted up to 7 days, shared by date",
  "S14",     "sjitter:14",                        "Dates shifted up to 14 days, shared by date",
  "S30",     "sjitter:30",                        "Dates shifted up to 30 days, shared by date",
  "M10",     "dates,cover,tillage,misclass:0.1",  "D3, 10% of cover and tillage misclassified",
  "M20",     "dates,cover,tillage,misclass:0.2",  "D3, 20% of cover and tillage misclassified")

# ==============================================================================
# Degrading workbooks (python)
# ==============================================================================
PY <- local({
  if (nzchar(CFG_C$python)) CFG_C$python else {
    p <- unname(Sys.which(c("python3", "python", "py")))
    p <- p[nzchar(p)]
    ok <- vapply(p, function(x) identical(suppressWarnings(system2(shQuote(x), "--version", stdout = FALSE, stderr = FALSE)), 0L), logical(1))
    if (any(ok)) p[ok][1] else ""
  }
})
degrade <- function(inp, out, ops, tillage, seed = 1, site_regex = "NONE") {
  if (!nzchar(ops)) return(inp)
  if (file.exists(out) && !CFG_C$rebuild) return(out)
  if (!nzchar(PY)) stop("No python found to write ", basename(out), ". Unzip the prebuilt coarse workbooks into ",
                        "Outputs/Coarse/, or set CFG_COARSE_OVERRIDE$python", call. = FALSE)
  dir.create(dirname(out), recursive = TRUE, showWarnings = FALSE)
  args <- c(CFG_C$script, "--in", inp, "--out", out, "--ops", ops, "--tillage", tillage,
            "--seed", seed, "--site-regex", site_regex, "--log", OUTC("degrade_log.csv"))
  res <- suppressWarnings(system2(shQuote(PY), shQuote(args), stdout = TRUE, stderr = TRUE))
  if (!is.null(attr(res, "status")) && attr(res, "status") != 0)
    stop("degrade_records.py failed for ", basename(out), ":\n", paste(res, collapse = "\n"), call. = FALSE)
  message(basename(out), ": ", paste(trimws(res), collapse = " | "))
  out
}

# ==============================================================================
# Part A: NAPESHM, agreement and skill against soil health
# ==============================================================================
# shared helpers (as in SHMI_benchmark.R / SHMI_calibration_plus.R)
site_w <- function(s) 1 / as.numeric(table(s)[s])
wmean  <- function(x, w) sum(w * x) / sum(w)
r2     <- function(obs, pred) 1 - sum((obs - pred)^2) / sum((obs - mean(obs))^2)
safe_r <- function(a, b) if (length(a) > 2 && sd(a) > 0 && sd(b) > 0) cor(a, b) else NA_real_
line_fit <- function(s, y, w) {
  ms <- wmean(s, w); my <- wmean(y, w)
  b <- sum(w * (s - ms) * (y - my)) / sum(w * (s - ms)^2); c(a = my - b * ms, b = b)
}
q <- function(x, p) unname(stats::quantile(x, p, na.rm = TRUE))

if (CFG_C$run_napeshm) {
  res  <- readRDS(CFG_C$final)
  CFG  <- res$config                       # K, R, seed, fixed_div, graze_gap, amend_gap, div_window
  W    <- res$w_final                      # cover, diversity, disturbance, growing_share
  R_rep <- if (is.null(CFG_C$R)) CFG$R else CFG_C$R
  us0  <- tibble::as_tibble(res$us)
  OTHER <- c("NonGrow", "Grow", "InvDist")

  # units removed from all analyses (as in Prepare_SEM.R); the rest are matched to res$us below
  EXCLUDE_ALL <- c("NAPESHM_USMI04_306", "NAPESHM_USNC03_368", "NAPESHM_USNC03_369", "NAPESHM_USNC03_370",
                   "NAPESHM_USNC03_371", "NAPESHM_USNC03_372",
                   "NAPESHM_USAZ01_167", "NAPESHM_USAZ01_168", "NAPESHM_USAZ03_169", "NAPESHM_USAZ03_170",
                   "NAPESHM_USCO05_222", "NAPESHM_USCO05_223", "NAPESHM_USCO05_224",
                   "NAPESHM_USSD03_465", "NAPESHM_USSD03_466", "NAPESHM_USSD03_467",
                   "NAPESHM_USWY01_506", "NAPESHM_USWY01_507")
  irr <- read.csv(CFG_C$treatments, stringsAsFactors = FALSE) %>%
    transmute(MGT_combo = paste0("NAPESHM_", sitecode, "_", treatmentid),
              irrigated = !is.na(eu_irrigation) & !(trimws(tolower(eu_irrigation)) %in% c("none", "", "n_a")))
  clim <- readRDS(CFG_C$climate)
  options(SHMI.diversity_window = if (is.null(CFG$div_window)) 4 else CFG$div_window)
  org_parts <- getExportedValue("SHMI", "compute_orginput_components")

  # components exactly as in SHMI_calibration_plus.R (NAPESHM records carry STIR)
  components <- function(wb) {
    inp <- prepare_shmi_inputs(wb, exclude = EXCLUDE_ALL)
    climate <- inp$mgt %>%
      transmute(MGT_combo, sitecode = sub("^NAPESHM_(.*)_[^_]+$", "\\1", MGT_combo)) %>%
      left_join(clim, by = "sitecode") %>% left_join(irr, by = "MGT_combo") %>%
      mutate(irrigated = coalesce(irrigated, FALSE)) %>% select(-sitecode)
    oc <- org_parts(inp$rot_bounds, inp$amend, inp$animal) %>% select(MGT_combo, p_amend, p_animal)
    X <- shmi_components(inp, climate, dist_meth = "STIR")
    if (!"Diversity_plants" %in% names(X)) X$Diversity_plants <- X$Diversity
    X %>% transmute(MGT_combo, NonGrow = Cover_nongrowing, Grow = Cover_growing, Diversity_plants, Richness, InvDist) %>%
      left_join(oc, by = "MGT_combo") %>%
      mutate(across(c(p_amend, p_animal), ~ coalesce(.x, 0)),
             Diversity = Diversity_plants + CFG$graze_gap * p_animal * (100 - Diversity_plants))
  }
  # the published SHMI (weights fixed)
  cover_of <- function(D, w) (1 - w[["growing_share"]]) * D$NonGrow + w[["growing_share"]] * D$Grow
  shmi_w <- function(D, w) {
    base <- w[["cover"]] * cover_of(D, w) + w[["diversity"]] * D$Diversity + w[["disturbance"]] * D$InvDist
    base + CFG$amend_gap * D$p_amend * (100 - base)
  }
  # weights refitted (fit_plus() of the calibration, all-variation fit)
  fit_plus <- function(D, w) {
    k <- CFG$fixed_div / (1 - CFG$fixed_div); cp <- CFG$amend_gap * D$p_amend
    Xm <- (as.matrix(D[, OTHER]) + k * D$Diversity) * (1 - cp) + (1 + k) * 100 * cp
    xm <- colSums(Xm * w) / sum(w); ym <- wmean(D$SH, w)
    f <- nnls::nnls(sweep(Xm, 2, xm) * sqrt(w), (D$SH - ym) * sqrt(w))
    b <- setNames(f$x, OTHER); B <- (1 + k) * sum(b); cover <- b[["NonGrow"]] + b[["Grow"]]
    if (B <= 1e-12) return(NULL)
    c(cover = cover / B, diversity = CFG$fixed_div, disturbance = b[["InvDist"]] / B,
      growing_share = if (cover > 0) b[["Grow"]] / cover else 0.5)
  }

  rankable <- us0 %>% group_by(sitecode) %>% summarise(ok = n() >= 3 & sd(SH) > 0, .groups = "drop")
  N_RANKABLE <- sum(rankable$ok)
  perf <- function(obs, pred, site) {
    mo <- tapply(obs, site, mean); mp <- tapply(pred, site, mean)
    do <- obs - ave(obs, site); dp <- pred - ave(pred, site)
    rho <- vapply(split(seq_along(obs), site), function(i)
      if (length(i) < 3 || sd(obs[i]) == 0 || sd(pred[i]) == 0) NA_real_ else cor(obs[i], pred[i], method = "spearman"),
      numeric(1))
    rho <- rho[!is.na(rho)]
    c(R2 = r2(obs, pred), RMSE = sqrt(mean((obs - pred)^2)), between_r = safe_r(mo, mp),
      within_r = safe_r(do, dp), within_median_rho = if (length(rho)) median(rho) else NA_real_,
      sites_correct = sum(rho > 0) / N_RANKABLE)
  }
  eval_cv <- function(D) {
    sites <- unique(D$sitecode); y <- D$SH; out <- vector("list", R_rep)
    s_fixed <- shmi_w(D, W)
    for (r in seq_len(R_rep)) {
      set.seed(CFG$seed + r)                                       # the calibration's folds
      fold <- sample(rep_len(seq_len(CFG$K), length(sites)))[match(D$sitecode, sites)]
      pf <- pr <- rep(NA_real_, nrow(D))
      for (k in seq_len(CFG$K)) {
        tr <- fold != k; w <- site_w(D$sitecode[tr])
        l <- line_fit(s_fixed[tr], y[tr], w); pf[!tr] <- l[["a"]] + l[["b"]] * s_fixed[!tr]
        if (CFG_C$refit) {
          wt <- fit_plus(D[tr, ], w)
          if (!is.null(wt)) {
            s2 <- shmi_w(D, wt); l2 <- line_fit(s2[tr], y[tr], w); pr[!tr] <- l2[["a"]] + l2[["b"]] * s2[!tr]
          }
        }
      }
      out[[r]] <- bind_rows(tibble::as_tibble_row(c(repeat_id = r, weights = 1, perf(y, pf, D$sitecode))),
                            if (CFG_C$refit && !anyNA(pr))
                              tibble::as_tibble_row(c(repeat_id = r, weights = 2, perf(y, pr, D$sitecode))))
    }
    bind_rows(out) %>% mutate(weights = c("published", "refitted")[weights])
  }
  agreement <- function(D, D0) {
    s <- shmi_w(D, W); s0 <- shmi_w(D0, W)
    rho <- vapply(split(seq_along(s), D$sitecode), function(i)
      if (length(i) < 3 || sd(s[i]) == 0 || sd(s0[i]) == 0) NA_real_ else cor(s[i], s0[i], method = "spearman"), numeric(1))
    tibble(rms_diff = sqrt(mean((s - s0)^2)), bias = mean(s - s0), within_rho = median(rho, na.rm = TRUE),
           rms_cover = sqrt(mean((cover_of(D, W) - cover_of(D0, W))^2)),
           rms_diversity = sqrt(mean((D$Diversity - D0$Diversity)^2)),
           rms_invdist = sqrt(mean((D$InvDist - D0$InvDist)^2)),
           units_amended = sum(D$p_amend > 0), units_grazed = sum(D$p_animal > 0))
  }
  as_units <- function(X) us0 %>% select(MGT_combo, sitecode, system, SH) %>%
    left_join(X, by = "MGT_combo")

  # D0 from the same workbook and code, so every scenario is compared like with like
  X0 <- components(CFG_C$napeshm_wb)
  D0 <- as_units(X0)
  if (anyNA(D0$InvDist)) stop("Full-record components missing for some units")
  chk <- cor(shmi_w(D0, W), shmi_w(us0, W))
  message(sprintf("D0 check: SHMI from the workbook vs the calibration's components, r = %.4f", chk))
  if (chk < 0.99) warning("D0 differs from the calibration's components (r = ", round(chk, 3),
                          "): check that Inputs/NAPESHM_rebuilt.xlsx is the QC'd workbook used there")

  runs <- bind_rows(SCEN %>% mutate(seed = 1L),
                    NOISE %>% tidyr::crossing(seed = seq_len(CFG_C$n_noise)))
  if (!is.null(CFG_C$napeshm_scenarios)) runs <- runs %>% filter(scenario %in% CFG_C$napeshm_scenarios)
  rep_rows <- list(); agr_rows <- list()
  for (i in seq_len(nrow(runs))) {
    sc <- runs$scenario[i]; sd_ <- runs$seed[i]
    message(sprintf("[NAPESHM %d/%d] %s, seed %d", i, nrow(runs), sc, sd_))
    wb <- degrade(CFG_C$napeshm_wb, OUTC("napeshm", sprintf("NAPESHM_%s_s%02d.xlsx", sc, sd_)), runs$ops[i],
                  tillage = "stir", seed = sd_, site_regex = "^(NAPESHM_[^_]+)_")
    D <- if (sc == "D0") D0 else as_units(components(wb))
    if (anyNA(D$InvDist)) stop(sc, ": components missing for some units")
    agr_rows[[i]] <- agreement(D, D0) %>% mutate(scenario = sc, seed = sd_, .before = 1)
    rep_rows[[i]] <- eval_cv(D) %>% mutate(scenario = sc, seed = sd_, .before = 1)
  }
  reps <- bind_rows(rep_rows); agr <- bind_rows(agr_rows)
  if (!is.null(CFG_C$napeshm_scenarios)) {        # merge with the saved results of the other scenarios
    if (file.exists(OUTC("coarse_napeshm_repeats.csv")))
      reps <- bind_rows(read.csv(OUTC("coarse_napeshm_repeats.csv"), stringsAsFactors = FALSE) %>%
                          filter(!scenario %in% CFG_C$napeshm_scenarios), reps)
    old_agr <- if (file.exists(OUTC("coarse_agreement_seeds.csv")))
      read.csv(OUTC("coarse_agreement_seeds.csv"), stringsAsFactors = FALSE)
    else if (file.exists(OUTC("coarse_agreement.csv")))       # runs before the per-seed file existed
      read.csv(OUTC("coarse_agreement.csv"), stringsAsFactors = FALSE) %>% select(-any_of(c("label", "n_seeds"))) %>% mutate(seed = NA_integer_)
    if (!is.null(old_agr)) agr <- bind_rows(old_agr %>% filter(!scenario %in% CFG_C$napeshm_scenarios), agr)
  }
  write.csv(reps, OUTC("coarse_napeshm_repeats.csv"), row.names = FALSE)
  write.csv(agr, OUTC("coarse_agreement_seeds.csv"), row.names = FALSE)
  labs_all <- bind_rows(SCEN %>% select(scenario, label), NOISE %>% select(scenario, label))
  agr_sum <- agr %>% group_by(scenario) %>% summarise(across(-seed, ~ mean(.x, na.rm = TRUE)), n_seeds = n(), .groups = "drop") %>%
    left_join(labs_all, by = "scenario") %>% relocate(label, .after = scenario)
  write.csv(agr_sum %>% mutate(across(where(is.numeric), ~ signif(.x, 4))), OUTC("coarse_agreement.csv"), row.names = FALSE)
  nap_sum <- reps %>% group_by(scenario, weights) %>%
    summarise(across(c(R2, between_r, within_r, within_median_rho, sites_correct),
                     list(med = ~ median(.x, na.rm = TRUE), lo = ~ q(.x, 0.025), hi = ~ q(.x, 0.975))),
              n = n(), .groups = "drop") %>%
    left_join(labs_all, by = "scenario") %>% relocate(label, .after = scenario)
  write.csv(nap_sum %>% mutate(across(where(is.numeric), ~ round(.x, 4))), OUTC("coarse_napeshm.csv"), row.names = FALSE)
  print(nap_sum %>% filter(weights == "published") %>% select(scenario, label, R2_med, R2_lo, R2_hi, sites_correct_med), n = Inf)
}

# ==============================================================================
# Part B: yield, the within-site slope with coarse records
# ==============================================================================
if (CFG_C$run_yield) {
  leak <- intersect(c("SCORE_WB_DIR_OVERRIDE", "OUT_OVERRIDE", "SHMI_WEIGHTS_OVERRIDE"), ls(globalenv()))
  if (length(leak)) stop("Remove these from the global environment first: ", paste(leak, collapse = ", "))
  wbs <- list.files(CFG_C$yield_wb_dir, pattern = "^SHMI_.*\\.xlsx$", full.names = TRUE)
  ysc <- SCEN %>% filter(yield)
  y_rows <- list(); keep_corn <- list()
  # a scenario not in yield_scenarios is read back: from its saved units (coarse_units.rds), or, for runs
  # made before that file was written, from its results_main.csv (then without the agreement columns).
  # D0 is the base of the agreement, so it is rerun (fast: full-record site runs are cached) if not saved.
  read_saved <- function(sc, label) {
    d <- OUTC(paste0("yield_", sc)); if (!dir.exists(d)) d <- paste0(d, "_agcros"); rds <- file.path(d, "coarse_units.rds"); csv <- file.path(d, "results_main.csv")
    if (file.exists(rds)) return(readRDS(rds))
    if (sc == "D0" || !file.exists(csv)) return(NULL)
    r <- read.csv(csv, stringsAsFactors = FALSE); f1 <- r[r$figure == "Fig 1", ][1, ]
    xn <- r[grepl("^SHMI x N: change", r$estimate), ][1, ]
    list(row = tibble(scenario = sc, label = label, slope = f1$effect, lower = f1$lower, upper = f1$upper,
                      p = f1$p_val, SHMIxN = xn$effect, SHMIxN_p = xn$p_val, plot_years = NA_integer_), units = NULL)
  }
  for (i in seq_len(nrow(ysc))) {
    sc <- ysc$scenario[i]
    if (!is.null(CFG_C$yield_scenarios) && !(sc %in% CFG_C$yield_scenarios)) {
      sv <- read_saved(sc, ysc$label[i])
      if (!is.null(sv)) {
        message("Yield, ", sc, ": read from the saved run")
        y_rows[[i]] <- sv$row; if (!is.null(sv$units)) keep_corn[[sc]] <- sv$units
        next
      }
    }
    cat(sprintf("\n==================== Yield, %s: %s ====================\n", sc, ysc$label[i]))
    env <- new.env(parent = globalenv())
    env$OUT_OVERRIDE <- OUTC(paste0("yield_", sc))
    if (sc != "D0") {
      dir_s <- OUTC("yield_workbooks", sc)
      for (f in wbs) degrade(f, file.path(dir_s, basename(f)), ysc$ops[i], tillage = "epa", seed = 1, site_regex = "NONE")
      env$SCORE_WB_DIR_OVERRIDE <- dir_s
    }
    sys.source(CFG_C$main_script, envir = env)
    r <- read.csv(file.path(env$OUT, "results_main.csv"), stringsAsFactors = FALSE)
    f1 <- r[r$figure == "Fig 1", ][1, ]
    xn <- r[grepl("^SHMI x N: change", r$estimate), ][1, ]
    y_rows[[i]] <- tibble(scenario = sc, label = ysc$label[i], slope = f1$effect, lower = f1$lower, upper = f1$upper,
                          p = f1$p_val, SHMIxN = xn$effect, SHMIxN_p = xn$p_val, plot_years = nrow(env$corn))
    keep_corn[[sc]] <- env$corn %>% group_by(site, MGT_combo) %>% summarise(SHMI = mean(SHMI), .groups = "drop")
    saveRDS(list(row = y_rows[[i]], units = keep_corn[[sc]]), file.path(env$OUT, "coarse_units.rds"))
    rm(env); gc(verbose = FALSE)
  }
  yld <- bind_rows(y_rows)
  base <- if (is.null(keep_corn[["D0"]])) NULL else keep_corn[["D0"]] %>% rename(SHMI0 = SHMI)
  if (is.null(base)) stop("No D0 units for the agreement")
  yagr <- bind_rows(lapply(names(keep_corn), function(sc) {
    d <- keep_corn[[sc]] %>% inner_join(base, by = c("site", "MGT_combo"))
    rho <- d %>% group_by(site) %>% filter(n() >= 3, sd(SHMI) > 0, sd(SHMI0) > 0) %>%
      summarise(rho = cor(SHMI, SHMI0, method = "spearman"), .groups = "drop")
    tibble(scenario = sc, units = nrow(d), rms_diff = sqrt(mean((d$SHMI - d$SHMI0)^2)), bias = mean(d$SHMI - d$SHMI0),
           within_rho = median(rho$rho))
  }))
  yld <- yld %>% left_join(yagr, by = "scenario")
  write.csv(yld %>% mutate(across(where(is.numeric), ~ signif(.x, 4))), OUTC("coarse_yield.csv"), row.names = FALSE)
  print(yld)
}

# ==============================================================================
# Yield under noise (optional): one run of SHMI_yield_main.R per scenario and seed
# ==============================================================================
noise_row <- function(r, sc, seed) {
  f1 <- r[r$figure == "Fig 1", ][1, ]; xn <- r[grepl("^SHMI x N: change", r$estimate), ][1, ]
  tibble(scenario = sc, seed = seed, slope = f1$effect, lower = f1$lower, upper = f1$upper, p = f1$p_val,
         SHMIxN = xn$effect, SHMIxN_p = xn$p_val)
}
if (length(CFG_C$yield_noise)) {
  leak <- intersect(c("SCORE_WB_DIR_OVERRIDE", "OUT_OVERRIDE", "SHMI_WEIGHTS_OVERRIDE", "SHMI_FN_OVERRIDE"), ls(globalenv()))
  if (length(leak)) stop("Remove these from the global environment first: ", paste(leak, collapse = ", "))
  wbs <- list.files(CFG_C$yield_wb_dir, pattern = "^SHMI_.*\\.xlsx$", full.names = TRUE)
  n_rows <- list()
  seeds_for <- function(sc) {
    z <- CFG_C$yield_noise_seeds
    if (!is.list(z)) return(z)
    if (!is.null(z[[sc]])) z[[sc]] else if (!is.null(z$default)) z$default else 1
  }
  for (sc in CFG_C$yield_noise) for (seed in seeds_for(sc)) {
    tag <- sprintf("%s_s%02d", sc, seed)
    out <- OUTC(paste0("yield_", tag)); out_a <- paste0(out, "_agcros")
    done <- file.path(if (dir.exists(out_a)) out_a else out, "results_main.csv")
    if (CFG_C$reuse_runs && file.exists(done)) {                 # resume: a finished seed is not rerun
      message("Yield, ", tag, ": read from the saved run")
      n_rows[[tag]] <- noise_row(read.csv(done, stringsAsFactors = FALSE), sc, seed); next
    }
    cat(sprintf("\n==================== Yield, %s ====================\n", tag))
    dir_s <- OUTC("yield_workbooks", tag)
    ops <- NOISE$ops[NOISE$scenario == sc]
    for (f in wbs) degrade(f, file.path(dir_s, basename(f)), ops, tillage = "epa", seed = seed, site_regex = "NONE")
    env <- new.env(parent = globalenv())
    env$OUT_OVERRIDE <- out; env$SCORE_WB_DIR_OVERRIDE <- dir_s
    sys.source(CFG_C$main_script, envir = env)
    n_rows[[tag]] <- noise_row(read.csv(file.path(env$OUT, "results_main.csv"), stringsAsFactors = FALSE), sc, seed)
    rm(env); gc(verbose = FALSE)
  }
  new_runs <- bind_rows(n_rows)
  if (file.exists(OUTC("coarse_yield_noise_runs.csv"))) {      # keep the saved runs of other scenarios and seeds
    old_runs <- read.csv(OUTC("coarse_yield_noise_runs.csv"), stringsAsFactors = FALSE)
    new_runs <- bind_rows(old_runs %>% anti_join(new_runs, by = c("scenario", "seed")), new_runs)
  }
  write.csv(new_runs %>% arrange(scenario, seed) %>% mutate(across(where(is.numeric), ~ signif(.x, 4))),
            OUTC("coarse_yield_noise_runs.csv"), row.names = FALSE)
}

# ==============================================================================
# Figure: skill and yield slope by scenario
# ==============================================================================
pl <- list()
# yield under noise: mean over seeds of the slope and of its 95% CI limits
ynr <- if (file.exists(OUTC("coarse_yield_noise_runs.csv"))) read.csv(OUTC("coarse_yield_noise_runs.csv"), stringsAsFactors = FALSE) else NULL
if (!is.null(ynr) && nrow(ynr)) {
  yns <- ynr %>% group_by(scenario) %>%
    summarise(seeds = n(), slope_min = min(slope), slope_max = max(slope),
              across(c(slope, lower, upper, p, SHMIxN, SHMIxN_p), mean), .groups = "drop") %>%
    left_join(NOISE %>% select(scenario, label), by = "scenario")
  write.csv(yns %>% mutate(across(where(is.numeric), ~ signif(.x, 4))), OUTC("coarse_yield_noise.csv"), row.names = FALSE)
}
if (!exists("nap_sum") && file.exists(OUTC("coarse_napeshm.csv"))) nap_sum <- read.csv(OUTC("coarse_napeshm.csv"), stringsAsFactors = FALSE)
if (!exists("yld") && file.exists(OUTC("coarse_yield.csv"))) yld <- read.csv(OUTC("coarse_yield.csv"), stringsAsFactors = FALSE)
if (exists("nap_sum")) {
  ns <- nap_sum %>% filter(weights == "published")
  pl$a <- bind_rows(
    ns %>% transmute(scenario, label, panel = "a  Soil health: out-of-fold R²", est = R2_med, lo = R2_lo, hi = R2_hi),
    ns %>% transmute(scenario, label, panel = "b  Sites ranked correctly (%)", est = 100 * sites_correct_med,
                     lo = 100 * sites_correct_lo, hi = 100 * sites_correct_hi))
}
if (exists("yld")) pl$b <- yld %>% transmute(scenario, label, panel = "c  Corn yield per 10 SHMI (kg/ha)",
                                             est = slope, lo = lower, hi = upper)
if (exists("yns")) pl$c <- yns %>% transmute(scenario, label, panel = "c  Corn yield per 10 SHMI (kg/ha)",
                                             est = slope, lo = lower, hi = upper)
if (length(pl)) {
  fd <- bind_rows(pl)
  ord <- c(SCEN$scenario, NOISE$scenario)
  fd <- fd %>% mutate(scenario = factor(scenario, levels = rev(ord)),
                      lab = factor(paste0(scenario, "  ", label), levels = rev(paste0(c(SCEN$scenario, NOISE$scenario), "  ",
                                                                                  c(SCEN$label, NOISE$label)))),
                      kind = ifelse(scenario %in% NOISE$scenario, "Noise", ifelse(scenario == "D0", "Full record", "Coarser record")))
  ref <- fd %>% filter(scenario == "D0") %>% select(panel, ref = est)
  g <- ggplot(fd, aes(est, lab, colour = kind)) +
    geom_vline(data = ref, aes(xintercept = ref), linetype = "dashed", colour = "grey55") +
    geom_errorbar(aes(xmin = lo, xmax = hi), width = 0.25) + geom_point(size = 2.4) +
    facet_wrap(~ panel, scales = "free_x", nrow = 1) +
    scale_colour_manual(values = c("Full record" = "#1b7837", "Coarser record" = "#2166ac", "Noise" = "#b2182b"), name = NULL) +
    labs(x = NULL, y = NULL, title = "SHMI with coarser management records",
         caption = paste0("a, b: NAPESHM US units, published weights, the calibration's leave-sites-out folds; median and 95% range over repeats ",
                          "(noise: and seeds).\nb: D0-D6 have no range: a rising line from SHMI to soil health keeps the order of treatments within ",
                          "a site, so the share is the same in every repeat\n(58 rankable sites; binomial 95% CI at 81%: 69-89%). ",
                          "c: within-site slope of the Fig 1 model at 150 kg N/ha (95% CI); SHMI from the coarse record,\nyields and N from the ",
                          "full record (noise: mean over seeds). D6: NAPESHM only. J: each record shifted on its own; S: one shift per date and site. ",
                          "Dashed: full record.")) +
    theme_bw(base_size = 10) + theme(legend.position = "bottom", plot.caption = element_text(hjust = 0, size = 7.5))
  if (isTRUE(getOption("shmi.ms_figures", TRUE))) {     # manuscript: no title or caption (they go in the legend);
    g <- g + labs(title = NULL, subtitle = NULL, caption = NULL) +   # the panel heading (letter + measure, which
      theme(strip.background = element_blank(),                     # is also the x-axis meaning) stays, left-aligned
            strip.text = element_text(face = "bold", size = 11, hjust = 0))
  }
  ggsave(OUTC("Fig_coarse_records.png"), g, width = 11.5, height = 6.4, dpi = 300)
}
message("Coarse-record test done: ", CFG_C$out_dir)
