################################################################################
# C1k_ceiling.R  (21 Nov 2026)
#
# How much of the explainable within-site variation in soil health does SHMI capture?
#
# The 0.35 within-site R2 for SHMI (Table 2) has no ceiling: part of the within-site
# variation in system-mean lv_SH is replicate (plot-to-plot) noise that no description
# of management could explain. The field replicates give that ceiling.
#
#   1. Every field replicate in soildata.csv is scored exactly as the system means
#      were (same code as C1i/A2c): the A1c outlier bounds, the A1c climate/texture
#      models (with the system's own clay+silt, so all replicates of a system get the
#      same adjustment), the A1c standardization, factor scores from the final US SEM
#      and the same 0-100 scale.
#   2. lmer(SH ~ 1 + (1 | site) + (1 | site:system)) on the replicate scores of the
#      354 systems used in the paper, REML. Within-site variation splits into
#      between-system variance (s2_sys) and replicate variance (s2_rep).
#   3. Ceilings, i.e. the most any description of management could explain:
#        system means (the scale of every SHMI analysis):
#          s2_sys / (s2_sys + s2_rep / n_h),  n_h = harmonic mean replicates per system
#        single replicates (shown for reference):  s2_sys / (s2_sys + s2_rep)
#   4. SHMI's within-site R2 on system means, in sample and at held-out sites
#      (leave-one-site-out, as C1j), and as a share of the system-mean ceiling.
#   5. Site bootstrap (B_BOOT resamples) for the ceiling and the shares.
#   SUBSET = "all" (default, 354 systems) or "annual" (annual cropland only;
#   writes to Outputs/Ceiling_annual). Set SUBSET before source().
#
# Inputs  Inputs/napeshm/soildata.csv (replicates), Outputs/semDF_prepped.rds (A1b),
#         Outputs/sem_noMX_projection.rds (A1c),
#         Outputs/SHAPE/SHAPE_by_unit_with_SHMI.csv (else data/...), as in C1j
# Writes  Outputs/Ceiling[_annual]/ceiling.csv and rep_scores.csv
# Runs in a few minutes (the bootstrap refits lmer B_BOOT times).
################################################################################
suppressPackageStartupMessages({ library(dplyr); library(lavaan); library(mgcv); library(lme4) })
setwd(here::here())
if (!exists("B_BOOT")) B_BOOT <- 1000
if (!exists("SUBSET")) SUBSET <- "all"          # "all" (354 systems) or "annual" (annual cropland only)
stopifnot(SUBSET %in% c("all", "annual"))
OUT <- if (SUBSET == "all") "Outputs/Ceiling" else "Outputs/Ceiling_annual"
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
set.seed(1)

F_REPS <- "Inputs/napeshm/soildata.csv"
F_SEM  <- "Outputs/semDF_prepped.rds"
F_PROJ <- "Outputs/sem_noMX_projection.rds"
F_UNIT <- if (file.exists("Outputs/SHAPE/SHAPE_by_unit_with_SHMI.csv")) "Outputs/SHAPE/SHAPE_by_unit_with_SHMI.csv" else "data/SHAPE_by_unit_with_SHMI.csv"
for (f in c(F_REPS, F_SEM, F_PROJ, F_UNIT)) if (!file.exists(f)) stop(f, " not found")

# ---- 1. replicate-level indicators and scores (as C1i sections 1-2) -------------------
make_combo <- function(sitecode, treatmentid) paste0("NAPESHM_", sitecode, "_", treatmentid)
reps <- read.csv(F_REPS, stringsAsFactors = FALSE) %>%
  mutate(MGT_combo = make_combo(sitecode, treatmentid)) %>%
  group_by(MGT_combo) %>% mutate(rep = row_number()) %>% ungroup() %>%
  transmute(MGT_combo, rep,
            Cmin = d_co2_c, Nmin = b_pmn_am_mgkg, ARYL = d_aryl, NAG = d_beta_minidase,
            SOC = b_soc, TN = b_total_n, POXC = a_act_c, ACE = a_ace,
            SSI = b_ssi, WAS = a_stabi, Ksat = shi_kfs_cmhr)

sem  <- readRDS(F_SEM) %>% mutate(sitecode = as.character(sitecode))
proj <- readRDS(F_PROJ)
cov_vars <- setdiff(proj$adj_vars, "sitecode")
reps <- reps %>% inner_join(sem %>% select(MGT_combo, sitecode, all_of(cov_vars)), by = "MGT_combo")

project_sh <- function(d, proj) {
  ind <- proj$indicators
  for (v in ind) { b <- proj$bounds[[v]]; x <- d[[v]]; x[!is.na(x) & (x < b[1] | x > b[2])] <- NA; d[[v]] <- x }
  for (v in ind) {
    keep <- complete.cases(d[, c(v, cov_vars)])
    out <- rep(NA_real_, nrow(d))
    if (any(keep)) {
      nd <- as.data.frame(d[keep, , drop = FALSE]); nd$site_f <- factor(nd$sitecode)
      tt <- mgcv::predict.gam(proj$models[[v]], newdata = nd, type = "terms")
      out[keep] <- d[[v]][keep] - rowSums(tt[, !grepl("site_f", colnames(tt)), drop = FALSE])
    }
    d[[v]] <- (out - proj$center[[v]]) / proj$sd[[v]]
  }
  has <- rowSums(!is.na(d[, ind])) > 0
  lv <- lavaan::lavPredict(proj$fit, newdata = as.data.frame(d[has, c(ind, "sitecode")]), type = "lv")
  sh <- rep(NA_real_, nrow(d)); sh[has] <- 100 * (lv[, "SH"] - proj$lv_min[["SH"]]) / proj$lv_range[["SH"]]
  d %>% select(MGT_combo, rep, sitecode) %>% mutate(n_missing = rowSums(is.na(d[, ind])), SH = sh)
}

# ---- 2. the 354 systems of the paper ---------------------------------------------------
units <- read.csv(F_UNIT, stringsAsFactors = FALSE) %>%
  filter(!is.na(SH), !is.na(SHMI), !is.na(SHAPE_mean)) %>%          # same filter as C1j
  mutate(sitecode = as.character(sitecode)) %>% select(MGT_combo, sitecode, SH, SHMI)
if (SUBSET == "annual") {
  annual <- sem %>% filter(system == "annual") %>% pull(MGT_combo)
  units  <- units %>% filter(MGT_combo %in% annual)
}
rep_sh <- project_sh(reps, proj) %>% filter(MGT_combo %in% units$MGT_combo, !is.na(SH))
write.csv(rep_sh, file.path(OUT, "rep_scores.csv"), row.names = FALSE)

chk <- rep_sh %>% group_by(MGT_combo) %>% summarise(SH_rep = mean(SH), n_rep = n(), .groups = "drop") %>%
  inner_join(units, by = "MGT_combo")
message(sprintf("C1k [%s]: %d replicates of %d systems at %d sites; replicates per system: %s",
                SUBSET, nrow(rep_sh), n_distinct(rep_sh$MGT_combo), n_distinct(rep_sh$sitecode),
                paste(names(table(chk$n_rep)), table(chk$n_rep), sep = ":", collapse = ", ")))
message(sprintf("C1k: replicate means vs system lv_SH: r = %.3f, mean |difference| = %.2f points",
                cor(chk$SH_rep, chk$SH), mean(abs(chk$SH_rep - chk$SH))))

# ---- 3. functions -----------------------------------------------------------------------
centre <- function(x, g) x - ave(x, g)
r2_within_in <- function(u) {                    # in sample, site-centred, one common slope
  y <- centre(u$SH, u$sitecode); x <- centre(u$SHMI, u$sitecode)
  1 - sum(resid(lm(y ~ 0 + x))^2) / sum(y^2)
}
r2_within_loso <- function(u) {                  # leave one site out, as C1j
  y <- centre(u$SH, u$sitecode); x <- centre(u$SHMI, u$sitecode); p <- numeric(length(y))
  for (s in unique(u$sitecode)) { m <- u$sitecode == s; b <- coef(lm(y[!m] ~ 0 + x[!m])); p[m] <- b * x[m] }
  1 - sum((y - p)^2) / sum(y^2)
}
ceiling_fit <- function(r) {
  r$sys <- interaction(r$sitecode, r$MGT_combo, drop = TRUE)
  m <- suppressMessages(lmer(SH ~ 1 + (1 | sitecode) + (1 | sys), data = r, REML = TRUE))
  vc <- as.data.frame(VarCorr(m))
  s2_sys <- vc$vcov[vc$grp == "sys"]; s2_rep <- vc$vcov[vc$grp == "Residual"]
  n_h <- 1 / mean(1 / as.numeric(table(droplevels(r$sys))))
  c(s2_sys = s2_sys, s2_rep = s2_rep, n_h = n_h,
    ceiling_means = s2_sys / (s2_sys + s2_rep / n_h),
    ceiling_single = s2_sys / (s2_sys + s2_rep))
}
estimate <- function(r, u) {
  cf <- ceiling_fit(r); r_in <- r2_within_in(u); r_lo <- r2_within_loso(u)
  c(cf, R2_in = r_in, R2_loso = r_lo,
    share_in = r_in / cf[["ceiling_means"]], share_loso = r_lo / cf[["ceiling_means"]])
}

# ---- 4. point estimate and site bootstrap ---------------------------------------------
u354 <- units %>% filter(MGT_combo %in% rep_sh$MGT_combo)
est <- estimate(rep_sh, u354)
sites <- unique(u354$sitecode)
boot <- t(replicate(B_BOOT, {
  pick <- sample(sites, length(sites), replace = TRUE)
  rb <- bind_rows(lapply(seq_along(pick), function(k) rep_sh %>% filter(sitecode == pick[k]) %>%
                           mutate(sitecode = paste0(sitecode, "_b", k), MGT_combo = paste0(MGT_combo, "_b", k))))
  ub <- bind_rows(lapply(seq_along(pick), function(k) u354 %>% filter(sitecode == pick[k]) %>%
                           mutate(sitecode = paste0(sitecode, "_b", k), MGT_combo = paste0(MGT_combo, "_b", k))))
  tryCatch(estimate(rb, ub), error = function(e) rep(NA_real_, length(est)))
}))
colnames(boot) <- names(est)
res <- data.frame(quantity = names(est), estimate = as.numeric(est),
                  lo = apply(boot, 2, quantile, 0.025, na.rm = TRUE),
                  hi = apply(boot, 2, quantile, 0.975, na.rm = TRUE),
                  subset = SUBSET, n_systems = nrow(u354), n_sites = length(sites), n_replicates = nrow(rep_sh),
                  B = sum(complete.cases(boot)), row.names = NULL)
write.csv(res, file.path(OUT, "ceiling.csv"), row.names = FALSE)
print(res, digits = 3)
message(sprintf("C1k [%s]: ceiling for system means %.2f (%.2f-%.2f); SHMI within-site R2 %.2f in sample (%.0f%% of ceiling), %.2f held out (%.0f%%)",
                SUBSET, est[["ceiling_means"]], res$lo[res$quantity == "ceiling_means"], res$hi[res$quantity == "ceiling_means"],
                est[["R2_in"]], 100 * est[["share_in"]], est[["R2_loso"]], 100 * est[["share_loso"]]))
