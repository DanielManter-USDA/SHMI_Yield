################################################################################
# B9_SHMI_coarse_records.R
#
# Does SHMI keep its skill when the management record is less detailed than the
# experiments' full records? Real remote-sensing pixels cannot be used: NAPESHM and
# the long-term experiments are small plots, below a 30-m pixel. Instead the full
# records are degraded (python/coarse_records/degrade_records.py) and rescored with
# the published weights, which are not refitted.
#
# a. Components: each piece of detail removed on its own from the full record, at
#    each relevant level
#      Crops            main crop only (mixtures, companion crops, intercrops dropped)
#      Cover crops      presence only (one generic rye cover per gap between cash
#                       crops); presence only with 10% or 20% of gaps misclassified
#      Dates            a typical date per crop and site (as public statistics give);
#                       every date shifted at random by up to 7, 14 or 30 days, either
#                       shared by everything recorded on the same day at a site (S, a
#                       farm calendar known approximately) or for each field and plot
#                       on its own (J, as detection from satellite images)
#      Tillage          class only (no-till, reduced, conventional: one pass per crop
#                       interval at the class's typical intensity); class with 10% or
#                       20% of intervals misclassified by one class
#      Manure, grazing  organic amendments and grazing removed
# b. Record types
#      Full record           the experiments' records (the paper's scores)
#      Typical farm record   crops and cover species kept; dates known to within 14
#                            days, shared across the farm (S14); tillage as class;
#                            manure and grazing recorded
#      Remote-sensing record main crop only; cover crops as presence; tillage class;
#                            dates detected for each field to within 14 days (J14);
#                            no manure or grazing (not visible)
#      Remote sensing + reported manure, grazing   the remote-sensing record with
#                            manure, compost and grazing added from the farmer
# Scenarios with a random element are repeated with several seeds (NAPESHM: n_noise;
# yield: yield_seeds).
#
# Measures
#   1. Agreement with the full-record score (NAPESHM units): RMS difference and mean
#      bias in SHMI points, median within-site Spearman rho, RMS difference by part.
#   2. Skill against soil health (NAPESHM, US units): out-of-fold R2, between- and
#      within-site r and the share of sites ranked correctly, on the calibration's
#      folds (10-fold leave-sites-out, repeated), with the PUBLISHED weights; only the
#      line from SHMI to lv_SH is fitted on the training sites, as in the calibration.
#      Secondary: the weights refitted on the coarse records (refit = TRUE).
#   3. Yield: B2_SHMI_yield_main.R rerun with each site scored from its degraded workbook
#      (SCORE_WB_DIR_OVERRIDE) and everything else (yields, N, climate) from the full
#      record: the within-site slope (Fig. 3a model) and SHMI agreement.
#
# RUN: from the project (here::here()), after A3a_SHMI_calibration_plus.R (Part A) and
#      B2_SHMI_yield_main.R (Part B; its site runs are reused). The degraded workbooks are
#      written with python (openpyxl) unless already in Outputs/Coarse/ (prebuilt).
#      Settings: CFG_COARSE_OVERRIDE <- list(...) before sourcing (see CFG_C).
# Outputs (Outputs/Coarse/): coarse_scenarios.csv, coarse_napeshm.csv,
#      coarse_napeshm_repeats.csv, coarse_agreement.csv, coarse_yield.csv,
#      coarse_yield_runs.csv, Fig_coarse_records.png, degrade_log.csv, the degraded
#      workbooks, yield_<scenario>[_s0#]_agcros/ (main runs)
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
  main_script  = here::here("R", "B2_SHMI_yield_main.R"),
  script       = here::here("python", "coarse_records", "degrade_records.py"),
  out_dir      = here::here("Outputs", "Coarse"),
  python       = "",       # python executable; "" = python3, else python, on the PATH
  R            = NULL,     # CV repeats; NULL = the calibration's (100)
  n_noise      = 5,        # NAPESHM seeds of each random scenario (the prebuilt workbooks cover 5)
  refit        = TRUE,     # also refit the weights on the degraded records
  run_napeshm  = TRUE,
  napeshm_scenarios = NULL, # e.g. c("J7", "J14"): rerun only these (results merged with the saved ones)
  run_yield    = TRUE,
  yield_scenarios = NULL,  # e.g. c("FARM", "RS"): run only these; the others are read from their saved runs
  yield_seeds  = list(default = 1, J30 = 1:3, FARM = 1:3, RSM = 1:3, RS = 1:3),   # yield seeds of each random scenario
  reuse_runs   = TRUE,     # read back yield runs already done (FALSE: rerun them all, as run_all.Rmd does)
  rebuild      = FALSE     # rewrite degraded workbooks that already exist
)
if (exists("CFG_COARSE_OVERRIDE")) CFG_C <- utils::modifyList(CFG_C, CFG_COARSE_OVERRIDE)
dir.create(CFG_C$out_dir, recursive = TRUE, showWarnings = FALSE)
OUTC <- function(...) file.path(CFG_C$out_dir, ...)

# one row per scenario: set (a components, b record types), component, level shown in
# the figure, operations for degrade_records.py, random (repeated with seeds), yield (rerun)
SCEN <- tribble(
  ~scenario, ~set, ~component,        ~level,                                   ~ops,                                       ~random, ~yield,
  "D0",      "b",  "Record type",     "Full record",                            "",                                         FALSE,   TRUE,
  "FARM",    "b",  "Record type",     "Typical farm record",                    "tillage,sjitter:14",                       TRUE,    TRUE,
  "RSM",     "b",  "Record type",     "Remote sensing + reported manure, grazing", "cover,maincrop,tillage,jitter:14",      TRUE,    TRUE,
  "RS",      "b",  "Record type",     "Remote-sensing record",                  "cover,maincrop,tillage,jitter:14,noamend", TRUE,    TRUE,
  "CROP",    "a",  "Crops",           "main crop only",                         "maincrop",                                 FALSE,   TRUE,
  "COVER",   "a",  "Cover crops",     "presence only",                          "cover",                                    FALSE,   TRUE,
  "COVER10", "a",  "Cover crops",     "presence, 10% misclassified",            "cover,misclass:0.1",                       TRUE,    TRUE,
  "COVER20", "a",  "Cover crops",     "presence, 20% misclassified",            "cover,misclass:0.2",                       TRUE,    TRUE,
  "D1",      "a",  "Dates",           "typical date per crop and site",         "dates",                                    FALSE,   TRUE,
  "S7",      "a",  "Dates",           "within 7 days, shared",                  "sjitter:7",                                TRUE,    TRUE,
  "S14",     "a",  "Dates",           "within 14 days, shared",                 "sjitter:14",                               TRUE,    TRUE,
  "S30",     "a",  "Dates",           "within 30 days, shared",                 "sjitter:30",                               TRUE,    TRUE,
  "J7",      "a",  "Dates",           "within 7 days, each field",              "jitter:7",                                 TRUE,    TRUE,
  "J14",     "a",  "Dates",           "within 14 days, each field",             "jitter:14",                                TRUE,    TRUE,
  "J30",     "a",  "Dates",           "within 30 days, each field",             "jitter:30",                                TRUE,    TRUE,
  "TILL",    "a",  "Tillage",         "class only",                             "tillage",                                  FALSE,   TRUE,
  "TILL10",  "a",  "Tillage",         "class, 10% misclassified",               "tillage,misclass:0.1",                     TRUE,    TRUE,
  "TILL20",  "a",  "Tillage",         "class, 20% misclassified",               "tillage,misclass:0.2",                     TRUE,    TRUE,
  "AMEND",   "a",  "Manure, grazing", "removed",                                "noamend",                                  FALSE,   TRUE) %>%
  mutate(label = ifelse(set == "b", level, paste0(component, ": ", level)))
write.csv(SCEN, OUTC("coarse_scenarios.csv"), row.names = FALSE)
labs_all <- SCEN %>% select(scenario, label)

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
# shared helpers (as in A4a_SHMI_benchmark.R / A3a_SHMI_calibration_plus.R)
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

  # units removed from all analyses (as in A1b_Prepare_SEM.R); the rest are matched to res$us below
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

  # components exactly as in A3a_SHMI_calibration_plus.R (NAPESHM records carry STIR)
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

  runs <- SCEN %>% select(scenario, ops, random) %>%
    tidyr::crossing(seed = seq_len(CFG_C$n_noise)) %>% filter(random | seed == 1)
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
                          filter(!scenario %in% CFG_C$napeshm_scenarios, scenario %in% SCEN$scenario), reps)
    if (file.exists(OUTC("coarse_agreement_seeds.csv")))
      agr <- bind_rows(read.csv(OUTC("coarse_agreement_seeds.csv"), stringsAsFactors = FALSE) %>%
                         filter(!scenario %in% CFG_C$napeshm_scenarios, scenario %in% SCEN$scenario), agr)
  }
  write.csv(reps, OUTC("coarse_napeshm_repeats.csv"), row.names = FALSE)
  write.csv(agr, OUTC("coarse_agreement_seeds.csv"), row.names = FALSE)
  agr_sum <- agr %>% group_by(scenario) %>% summarise(across(-seed, ~ mean(.x, na.rm = TRUE)), n_seeds = n(), .groups = "drop") %>%
    left_join(labs_all, by = "scenario") %>% relocate(label, .after = scenario)
  write.csv(agr_sum %>% mutate(across(where(is.numeric), ~ signif(.x, 4))), OUTC("coarse_agreement.csv"), row.names = FALSE)
  nap_sum <- reps %>% group_by(scenario, weights) %>%
    summarise(across(c(R2, between_r, within_r, within_median_rho, sites_correct),
                     list(med = ~ median(.x, na.rm = TRUE), lo = ~ q(.x, 0.025), hi = ~ q(.x, 0.975))),
              n = n(), .groups = "drop") %>%
    left_join(labs_all, by = "scenario") %>% relocate(label, .after = scenario) %>%
    arrange(match(scenario, SCEN$scenario))
  write.csv(nap_sum %>% mutate(across(where(is.numeric), ~ round(.x, 4))), OUTC("coarse_napeshm.csv"), row.names = FALSE)
  print(nap_sum %>% filter(weights == "published") %>% select(scenario, label, R2_med, R2_lo, R2_hi, sites_correct_med), n = Inf)
}

# ==============================================================================
# Part B: yield, the within-site slope with degraded records
# ==============================================================================
if (CFG_C$run_yield) {
  leak <- intersect(c("SCORE_WB_DIR_OVERRIDE", "OUT_OVERRIDE", "SHMI_WEIGHTS_OVERRIDE", "SHMI_FN_OVERRIDE"), ls(globalenv()))
  if (length(leak)) stop("Remove these from the global environment first: ", paste(leak, collapse = ", "))
  wbs <- list.files(CFG_C$yield_wb_dir, pattern = "^SHMI_.*\\.xlsx$", full.names = TRUE)
  seeds_for <- function(sc, random) {
    if (!random) return(1L)
    z <- CFG_C$yield_seeds
    if (!is.list(z)) return(z)
    if (!is.null(z[[sc]])) z[[sc]] else if (!is.null(z$default)) z$default else 1L
  }
  ysc <- SCEN %>% filter(yield) %>% rowwise() %>% mutate(seed = list(seeds_for(scenario, random))) %>%
    ungroup() %>% tidyr::unnest(seed) %>%
    mutate(tag = ifelse(random, sprintf("%s_s%02d", scenario, seed), scenario))
  run_row <- function(r, sc, seed, n_py = NA_integer_) {
    f1 <- r[r$figure == "Fig 1", ][1, ]; xn <- r[grepl("^SHMI x N: change", r$estimate), ][1, ]
    tibble(scenario = sc, seed = seed, slope = f1$effect, lower = f1$lower, upper = f1$upper, p = f1$p_val,
           SHMIxN = xn$effect, SHMIxN_p = xn$p_val, plot_years = n_py)
  }
  y_rows <- list(); keep_corn <- list()
  for (i in seq_len(nrow(ysc))) {
    sc <- ysc$scenario[i]; seed <- ysc$seed[i]; tag <- ysc$tag[i]
    out <- OUTC(paste0("yield_", tag)); out_dir <- if (dir.exists(paste0(out, "_agcros"))) paste0(out, "_agcros") else out
    rds <- file.path(out_dir, "coarse_units.rds"); csv <- file.path(out_dir, "results_main.csv")
    skip <- if (!is.null(CFG_C$yield_scenarios)) !(sc %in% CFG_C$yield_scenarios) else CFG_C$reuse_runs
    if (skip && file.exists(rds)) {                                  # resume: a finished run is read back
      sv <- readRDS(rds); y_rows[[tag]] <- sv$row; keep_corn[[tag]] <- sv$units
      message("Yield, ", tag, ": read from the saved run"); next
    }
    if (skip && sc != "D0" && file.exists(csv)) {                    # older runs: results only (no agreement)
      y_rows[[tag]] <- run_row(read.csv(csv, stringsAsFactors = FALSE), sc, seed)
      message("Yield, ", tag, ": read from the saved run (results only)"); next
    }
    cat(sprintf("\n==================== Yield, %s: %s ====================\n", tag, ysc$label[i]))
    env <- new.env(parent = globalenv())
    env$OUT_OVERRIDE <- out
    if (sc != "D0") {
      dir_s <- OUTC("yield_workbooks", tag)
      for (f in wbs) degrade(f, file.path(dir_s, basename(f)), ysc$ops[i], tillage = "epa", seed = seed, site_regex = "NONE")
      env$SCORE_WB_DIR_OVERRIDE <- dir_s
    }
    sys.source(CFG_C$main_script, envir = env)
    y_rows[[tag]] <- run_row(read.csv(file.path(env$OUT, "results_main.csv"), stringsAsFactors = FALSE), sc, seed, nrow(env$corn))
    keep_corn[[tag]] <- env$corn %>% group_by(site, MGT_combo) %>% summarise(SHMI = mean(SHMI), .groups = "drop")
    saveRDS(list(row = y_rows[[tag]], units = keep_corn[[tag]]), file.path(env$OUT, "coarse_units.rds"))
    rm(env); gc(verbose = FALSE)
  }
  yruns <- bind_rows(y_rows)
  base <- keep_corn[["D0"]]
  if (!is.null(base)) {
    base <- base %>% rename(SHMI0 = SHMI)
    yagr <- bind_rows(lapply(names(keep_corn), function(tg) {
      d <- keep_corn[[tg]] %>% inner_join(base, by = c("site", "MGT_combo"))
      rho <- d %>% group_by(site) %>% filter(n() >= 3, sd(SHMI) > 0, sd(SHMI0) > 0) %>%
        summarise(rho = cor(SHMI, SHMI0, method = "spearman"), .groups = "drop")
      tibble(tag = tg, units = nrow(d), rms_diff = sqrt(mean((d$SHMI - d$SHMI0)^2)), bias = mean(d$SHMI - d$SHMI0),
             within_rho = median(rho$rho))
    }))
    yruns <- yruns %>% mutate(tag = ifelse(SCEN$random[match(scenario, SCEN$scenario)], sprintf("%s_s%02d", scenario, seed), scenario)) %>%
      left_join(yagr, by = "tag") %>% select(-tag)
  }
  write.csv(yruns %>% mutate(across(where(is.numeric), ~ signif(.x, 4))), OUTC("coarse_yield_runs.csv"), row.names = FALSE)
  yld <- yruns %>% group_by(scenario) %>%
    summarise(seeds = n(), slope_min = min(slope), slope_max = max(slope),
              across(any_of(c("slope", "lower", "upper", "p", "SHMIxN", "SHMIxN_p", "rms_diff", "bias", "within_rho")),
                     ~ mean(.x, na.rm = TRUE)), .groups = "drop") %>%
    left_join(labs_all, by = "scenario") %>% relocate(label, .after = scenario) %>%
    mutate(pct_of_full = if (any(scenario == "D0")) 100 * slope / slope[scenario == "D0"][1] else NA_real_) %>%
    arrange(match(scenario, SCEN$scenario))
  write.csv(yld %>% mutate(across(where(is.numeric), ~ signif(.x, 4))), OUTC("coarse_yield.csv"), row.names = FALSE)
  print(yld %>% select(scenario, label, seeds, slope, lower, upper, pct_of_full), n = Inf)
}

# ==============================================================================
# Figure: a, components; b, record types (soil-health skill, ranking, yield slope)
# ==============================================================================
if (!exists("nap_sum") && file.exists(OUTC("coarse_napeshm.csv"))) nap_sum <- read.csv(OUTC("coarse_napeshm.csv"), stringsAsFactors = FALSE)
if (!exists("yld") && file.exists(OUTC("coarse_yield.csv"))) yld <- read.csv(OUTC("coarse_yield.csv"), stringsAsFactors = FALSE)
pl <- list()
if (exists("nap_sum")) {
  ns <- nap_sum %>% filter(weights == "published", scenario %in% SCEN$scenario)
  pl$r2 <- ns %>% transmute(scenario, measure = "Soil health: out-of-fold R²", est = R2_med, lo = R2_lo, hi = R2_hi)
  pl$sc <- ns %>% transmute(scenario, measure = "Sites ranked correctly (%)", est = 100 * sites_correct_med,
                            lo = 100 * sites_correct_lo, hi = 100 * sites_correct_hi)
}
if (exists("yld")) pl$y <- yld %>% filter(scenario %in% SCEN$scenario) %>%
  transmute(scenario, measure = "Corn yield per 10 SHMI (kg/ha)", est = slope, lo = lower, hi = upper)
if (length(pl)) {
  ms <- isTRUE(getOption("shmi.ms_figures", TRUE))
  fd <- bind_rows(pl) %>% left_join(SCEN %>% select(scenario, set, component, label), by = "scenario") %>%
    mutate(measure = factor(measure, c("Soil health: out-of-fold R²", "Sites ranked correctly (%)", "Corn yield per 10 SHMI (kg/ha)")),
           set = factor(set, c("a", "b"), if (ms) c("a", "b") else c("a  Components (each alone)", "b  Record types")),
           label = factor(label, rev(SCEN$label)))
  ref <- fd %>% filter(scenario == "D0") %>% select(measure, ref = est)
  g <- ggplot(fd, aes(est, label, colour = component)) +
    geom_vline(data = ref, aes(xintercept = ref), linetype = "dashed", colour = "grey55") +
    geom_errorbar(aes(xmin = lo, xmax = hi), width = 0.3) + geom_point(size = 2.2) +
    facet_grid(set ~ measure, scales = "free", space = "free_y") +
    scale_colour_manual(values = c("Record type" = "#1b7837", "Crops" = "#8c510a", "Cover crops" = "#35978f",
                                   "Dates" = "#2166ac", "Tillage" = "#b2182b", "Manure, grazing" = "#762a83"), name = NULL) +
    labs(x = NULL, y = NULL, title = "SHMI with less detailed management records",
         caption = paste0("Soil health: NAPESHM US units, published weights, the calibration's leave-site-out folds; median and 95% range over ",
                          "repeats and seeds. Sites ranked correctly: without a random element the share is the same in every repeat.\n",
                          "Yield: within-site slope of the Fig. 3a model at 150 kg N/ha (95% CI; mean over seeds); SHMI from the degraded record, ",
                          "yields and N from the full record. Dashed: full record.")) +
    theme_bw(base_size = 10) +
    theme(legend.position = "bottom", plot.caption = element_text(hjust = 0, size = 7.5),
          strip.text.y = element_text(angle = 0, face = "bold", hjust = 0))
  if (ms) g <- g + labs(title = NULL, subtitle = NULL, caption = NULL) +
    theme(strip.background = element_blank(), strip.text.x = element_text(face = "bold", size = 10, hjust = 0),
          strip.text.y = element_text(angle = 0, face = "bold", size = 13, hjust = 0))
  ggsave(OUTC("Fig_coarse_records.png"), g, width = 11.5, height = 7, dpi = 300)
}
message("Coarse-record test done: ", CFG_C$out_dir)
