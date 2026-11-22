################################################################################
# C1i_pair_meta.R  (16 Nov 2026)
#
# What a conventional practice-by-practice synthesis concludes from the strict
# matched pairs (Supplementary Table 6), and whether SHMI explains what it leaves
# unexplained (Fig. 3a-c).
#
# Each strict pair is a "study": two systems at the same site that differ in one
# practice and match on every other recorded label. Effect sizes use the field
# replicates in soildata.csv (3 per treatment in the US):
#   a  SOC: log response ratio (escalc "ROM") of raw SOC. Climate is shared within a
#      pair and SOC has a true zero, so no adjustment is needed (adjusted SOC is a
#      residual with no true zero, so a ratio of it would be meaningless).
#   b  lv_SH: mean difference (escalc "MD"). Each replicate is scored exactly as the
#      treatment means were (A2c): outliers with the training bounds, the training
#      climate/texture models (with the treatment's own clay+silt, so every replicate
#      of a treatment gets the same adjustment), training standardization, factor
#      scores from the fitted US SEM, the training 0-100 scale. lv_SH is interval-
#      scale (zero fixed by convention), so differences, not ratios, are used.
# Models (metafor), separately for each practice:
#   rma.mv(yi, V, random = ~ 1 | site / pair), V built from the replicate data so
#   pairs that share a system are correlated. Pooled effect, 95% CI and prediction interval, Q, I2.
#   c  meta-regression of the lv_SH effects on the change in the targeted sub-index
#      and, separately, on the change in SHMI (per 10 points): likelihood-ratio test
#      (ML fits) and share of heterogeneity explained (REML: 1 - tau2_mod / tau2_0).
#      Also fitted across all practices (practice as a factor, one shared slope).
#
# Input   Inputs/napeshm/soildata.csv (replicates), Outputs/semDF_prepped.rds (A1b),
#         Outputs/sem_noMX_projection.rds (A1c), data/SHAPE_by_unit_with_SHMI.csv,
#         Outputs/Manuscript_work/SuppTable6_practice_pairs.csv (C1),
#         Inputs/napeshm/rr_lookup.csv and treatments.csv (section 7, SHI's pairs)
# Output  Outputs/Pair_meta/ rep_scores.csv, pair_effects.csv, pooled.csv,
#         metaregression.csv, pairs_label_vs_SHMI.csv, plotdata_pair_meta.rds, rr_pairs_classified.csv, rr_pooled.csv
#         Outputs/panels/ panel_meta_forest.rds, panel_meta_reg.rds
# Runs in about a minute.
################################################################################
suppressPackageStartupMessages({ library(dplyr); library(tidyr); library(lavaan); library(mgcv)
                                  library(metafor); library(ggplot2) })
setwd(here::here())
OUT <- "Outputs/Pair_meta"; PANEL_DIR <- "Outputs/panels"
dir.create(OUT, recursive = TRUE, showWarnings = FALSE); dir.create(PANEL_DIR, recursive = TRUE, showWarnings = FALSE)
F_REPS  <- "Inputs/napeshm/soildata.csv"
F_SEM   <- "Outputs/semDF_prepped.rds"
F_PROJ  <- "Outputs/sem_noMX_projection.rds"
F_UNITS <- "data/SHAPE_by_unit_with_SHMI.csv"
F_PAIRS <- "Outputs/Manuscript_work/SuppTable6_practice_pairs.csv"
F_RR    <- "Inputs/napeshm/rr_lookup.csv"          # the Soil Health Institute's pairs
F_TRT   <- "Inputs/napeshm/treatments.csv"
for (f in c(F_REPS, F_SEM, F_PROJ, F_UNITS, F_PAIRS)) if (!file.exists(f)) stop(f, " not found")

PRACTICES <- c("Cover crop vs none", "Organic vs synthetic inputs",
               "Reduced vs conventional tillage", "Rotation vs monoculture")
TARGET_LAB <- c("Cover crop vs none" = "Cover/Roots", "Organic vs synthetic inputs" = "Amendments",
                "Reduced vs conventional tillage" = "Inverse Disturbance", "Rotation vs monoculture" = "Diversity")

# ---- 1. replicate-level indicators (mapping as in A1b) --------------------------------
make_combo <- function(sitecode, treatmentid) paste0("NAPESHM_", sitecode, "_", treatmentid)
reps <- read.csv(F_REPS, stringsAsFactors = FALSE) %>%
  mutate(MGT_combo = make_combo(sitecode, treatmentid)) %>%
  group_by(MGT_combo) %>% mutate(rep = row_number()) %>% ungroup() %>%
  transmute(MGT_combo, rep,
            Cmin = d_co2_c, Nmin = b_pmn_am_mgkg, ARYL = d_aryl, NAG = d_beta_minidase,
            SOC = b_soc, TN = b_total_n, POXC = a_act_c, ACE = a_ace,
            SSI = b_ssi, WAS = a_stabi, Ksat = shi_kfs_cmhr)
SOC_RAW <- reps %>% select(MGT_combo, rep, SOC_raw = SOC)

# treatment-level covariates for the adjustment (site climate, the treatment's clay+silt)
sem <- readRDS(F_SEM) %>% mutate(sitecode = as.character(sitecode))
proj <- readRDS(F_PROJ)
cov_vars <- setdiff(proj$adj_vars, "sitecode")
reps <- reps %>% inner_join(sem %>% select(MGT_combo, sitecode, system, all_of(cov_vars)), by = "MGT_combo")

# ---- 2. score each replicate with the fitted model (as A2c) ----------------------------
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
rep_sh <- project_sh(reps, proj) %>% left_join(SOC_RAW, by = c("MGT_combo", "rep"))
write.csv(rep_sh, file.path(OUT, "rep_scores.csv"), row.names = FALSE)

# check: replicate means against the treatment-level lv_SH used everywhere else
units <- read.csv(F_UNITS, stringsAsFactors = FALSE) %>% select(MGT_combo, SH_unit = SH, SHMI_unit = SHMI)
chk <- rep_sh %>% group_by(MGT_combo) %>% summarise(SH_rep = mean(SH, na.rm = TRUE), .groups = "drop") %>%
  inner_join(units, by = "MGT_combo")
cat(sprintf("Replicate-mean vs treatment lv_SH: r = %.3f, mean |difference| = %.2f points (%d treatments)\n",
            cor(chk$SH_rep, chk$SH_unit, use = "complete.obs"), mean(abs(chk$SH_rep - chk$SH_unit), na.rm = TRUE), nrow(chk)))

# ---- 3. strict pairs and their effect sizes ---------------------------------------------
supp6 <- read.csv(F_PAIRS, check.names = FALSE, stringsAsFactors = FALSE)
# 16 Nov 2026: Supplementary Table 6 lists every with/without comparison; the meta-analysis
# uses the strict matches only (every other recorded label identical, N rate included)
pairs <- supp6 %>% filter(if ("Pair type" %in% names(.)) `Pair type` == "Strict match" else TRUE) %>%
  transmute(practice = Practice, site = Site, unit_A = `System with practice`, unit_B = `System without`,
            d_target = `Change in targeted part`, d_SHMI = `Change in SHMI`, d_SH_unit = `Change in soil health (lv_SH)`) %>%
  mutate(pair = paste0("P", row_number()))
stats <- rep_sh %>% group_by(MGT_combo) %>%
  summarise(sh_m = mean(SH, na.rm = TRUE), sh_sd = sd(SH, na.rm = TRUE), sh_n = sum(!is.na(SH)),
            soc_m = mean(SOC_raw, na.rm = TRUE), soc_sd = sd(SOC_raw, na.rm = TRUE), soc_n = sum(!is.na(SOC_raw)),
            .groups = "drop")
pe <- pairs %>% left_join(rename_with(stats, ~ paste0(.x, "_A")), by = c("unit_A" = "MGT_combo_A")) %>%
  left_join(rename_with(stats, ~ paste0(.x, "_B")), by = c("unit_B" = "MGT_combo_B"))
ok_sh  <- with(pe, sh_n_A >= 2 & sh_n_B >= 2 & sh_sd_A > 0 & sh_sd_B > 0)
ok_soc <- with(pe, soc_n_A >= 2 & soc_n_B >= 2 & soc_sd_A > 0 & soc_sd_B > 0 & soc_m_A > 0 & soc_m_B > 0)
ok_sh[is.na(ok_sh)] <- FALSE; ok_soc[is.na(ok_soc)] <- FALSE
cat(sprintf("Strict pairs: %d; usable for lv_SH %d, for SOC %d\n", nrow(pe), sum(ok_sh), sum(ok_soc)))

es_sh <- escalc("MD", m1i = sh_m_A, sd1i = sh_sd_A, n1i = sh_n_A, m2i = sh_m_B, sd2i = sh_sd_B, n2i = sh_n_B,
                data = pe[ok_sh, ])
es_soc <- escalc("ROM", m1i = soc_m_A, sd1i = soc_sd_A, n1i = soc_n_A, m2i = soc_m_B, sd2i = soc_sd_B, n2i = soc_n_B,
                 data = pe[ok_soc, ])
# Pairs that share a system are correlated. The covariance is built exactly from the
# replicate data: each contrast is (mean of A) - (mean of B) (or the log ratio), so two
# contrasts sharing system S have covariance +var(S) when S is on the same side in both
# and -var(S) when it is on opposite sides, where var(S) is the sampling variance of
# S's mean (MD: sd^2/n; ROM: sd^2/(n m^2), as in escalc). This matrix is always
# positive semi-definite (vcalc's approximation from sample sizes alone was not).
build_V <- function(e, pre, type) {
  g <- function(side) {
    sd <- e[[paste0(pre, "_sd_", side)]]; n <- e[[paste0(pre, "_n_", side)]]; m <- e[[paste0(pre, "_m_", side)]]
    if (type == "MD") sd^2 / n else sd^2 / (n * m^2)
  }
  vA <- g("A"); vB <- g("B"); k <- nrow(e)
  V <- diag(e$vi, k)
  if (k > 1) for (i in 1:(k - 1)) for (j in (i + 1):k) {
    if (e$site[i] != e$site[j]) next
    cv <- 0
    if (e$unit_A[i] == e$unit_A[j]) cv <- cv + vA[i]
    if (e$unit_B[i] == e$unit_B[j]) cv <- cv + vB[i]
    if (e$unit_A[i] == e$unit_B[j]) cv <- cv - vA[i]
    if (e$unit_B[i] == e$unit_A[j]) cv <- cv - vB[i]
    V[i, j] <- V[j, i] <- cv
  }
  V
}
write.csv(bind_rows(mutate(as.data.frame(es_sh), outcome = "lv_SH (mean difference)"),
                    mutate(as.data.frame(es_soc), outcome = "SOC (log response ratio)")) %>%
            select(outcome, practice, site, pair, unit_A, unit_B, yi, vi, d_target, d_SHMI, d_SH_unit),
          file.path(OUT, "pair_effects.csv"), row.names = FALSE)

# ---- 4. pooled effect per practice ------------------------------------------------------
# multilevel I2 (share of total variance due to heterogeneity). Pairs that share systems
# can be linearly dependent (e.g. 3 cover-crop systems x 2 controls at one site give 6
# contrasts from 5 means: 4 independent), so the exact sampling matrix V is singular;
# a generalized inverse is used here (the model itself inverts V + the random effects).
i2_ml <- function(m) {
  W <- MASS::ginv(m$V); X <- model.matrix(m)
  P <- W - W %*% X %*% MASS::ginv(t(X) %*% W %*% X) %*% t(X) %*% W
  df <- qr(m$V)$rank - m$p
  100 * sum(m$sigma2) / (sum(m$sigma2) + df / sum(diag(P)))
}
n_indep <- function(e, pre, type) qr(build_V(e, pre, type))$rank   # independent contrasts
fit_one <- function(es, pre, type, ...) {
  V <- build_V(es, pre, type)
  rma.mv(yi, V, random = ~ 1 | site / pair, data = es, ...)
}
pool <- function(es, pre, type, outcome, groups = PRACTICES) bind_rows(lapply(groups, function(p) {
  e <- es[es$practice == p, ]
  if (nrow(e) < 3) return(NULL)
  m <- tryCatch(fit_one(e, pre, type), error = function(err) { message(outcome, ", ", p, ": ", conditionMessage(err)); NULL })
  if (is.null(m)) return(NULL)
  pr <- predict(m)
  i2 <- tryCatch(i2_ml(m), error = function(err) NA_real_)
  tibble(outcome = outcome, practice = p, pairs = nrow(e), independent = n_indep(e, pre, type), sites = n_distinct(e$site),
         estimate = as.numeric(m$b), ci_lb = m$ci.lb, ci_ub = m$ci.ub, pval = m$pval,
         pi_lb = pr$pi.lb, pi_ub = pr$pi.ub, tau2 = sum(m$sigma2), Q = m$QE, Q_p = m$QEp, I2 = i2)
}))
pooled <- bind_rows(pool(es_soc, "soc", "ROM", "SOC (log response ratio)"),
                    pool(es_sh,  "sh",  "MD",  "lv_SH (mean difference)")) %>%
  mutate(pct = ifelse(grepl("SOC", outcome), 100 * (exp(estimate) - 1), NA_real_))   # LRR as % change
write.csv(mutate(pooled, across(where(is.double), ~ signif(.x, 4))), file.path(OUT, "pooled.csv"), row.names = FALSE)
cat("\n== Pooled effects per practice ==\n"); print(as.data.frame(pooled), digits = 3)

# ---- 5. meta-regression: what explains the heterogeneity in lv_SH effects? -------------
# moderators: the size of the change (in the practice's own sub-index, or in SHMI) and the
# system the practice joins (SHMI of the system without it: does a practice do more when
# added to a weaker system?). The baseline is management only, not the control's soil
# health, which would build in regression to the mean.
base_of <- function(e) units$SHMI_unit[match(e$unit_B, units$MGT_combo)]
es_sh  <- es_sh  %>% mutate(d_target10 = d_target / 10, d_SHMI10 = d_SHMI / 10, base_SHMI10 = base_of(es_sh) / 10)
es_soc <- es_soc %>% mutate(d_target10 = d_target / 10, d_SHMI10 = d_SHMI / 10, base_SHMI10 = base_of(es_soc) / 10)
metareg <- function(e, base, label, pre = "sh", type = "MD") {
  V <- build_V(e, pre, type)
  f <- function(mods, method) rma.mv(yi, V, mods = mods, random = ~ 1 | site / pair, data = e, method = method)
  m0r <- f(base, "REML"); m0m <- f(base, "ML")
  bind_rows(lapply(c(sub_index = "d_target10", SHMI = "d_SHMI10", base = "base_SHMI10"), function(x) {
    mods <- update(base, as.formula(paste("~ . +", x)))
    m1r <- f(mods, "REML"); m1m <- f(mods, "ML")
    lrt <- list(LRT = 2 * (as.numeric(logLik(m1m)) - as.numeric(logLik(m0m))))   # 1 d.f.
    lrt$pval <- pchisq(lrt$LRT, df = 1, lower.tail = FALSE)
    tibble(set = label, moderator = x, pairs = nrow(e), slope_per10 = coef(m1r)[[x]],
           slope_lb = m1r$ci.lb[names(coef(m1r)) == x], slope_ub = m1r$ci.ub[names(coef(m1r)) == x],
           LRT = lrt$LRT, LRT_p = lrt$pval, tau2_0 = sum(m0r$sigma2), tau2_mod = sum(m1r$sigma2),
           tau2_explained = max(0, 1 - sum(m1r$sigma2) / sum(m0r$sigma2)))
  })) %>% mutate(moderator = dplyr::recode(moderator, d_target10 = "Change in targeted sub-index",
                                             d_SHMI10 = "Change in SHMI", base_SHMI10 = "SHMI of the system without the practice"))
}
mr_one <- function(es, pre, type, outcome) bind_rows(
  bind_rows(lapply(PRACTICES, function(p) { e <- es[es$practice == p, ]
    if (nrow(e) < 5) return(NULL)
    tryCatch(metareg(e, ~ 1, p, pre, type), error = function(err) { message("meta-regression, ", p, ": ", conditionMessage(err)); NULL }) })),
  metareg(es, ~ 0 + practice, "All practices (shared slope)", pre, type)) %>% mutate(outcome = outcome, .before = 1)
mr <- bind_rows(mr_one(es_sh, "sh", "MD", "lv_SH (mean difference)"),
                mr_one(es_soc, "soc", "ROM", "SOC (log response ratio)"))
write.csv(mutate(mr, across(where(is.double), ~ signif(.x, 4))), file.path(OUT, "metaregression.csv"), row.names = FALSE)
cat("\n== Meta-regression of lv_SH effects ==\n"); print(as.data.frame(mr), digits = 3)

# ---- 5b. the analysis of Fig. 3d, on the strict pairs ---------------------------------
# Fig. 3d asks how much within-site variance in lv_SH a yes/no label, SHMI or both explain
# at held-out sites. For two systems at a site the within-site comparison is their
# difference, so on the strict pairs the same question is: how much of the pair differences
# in lv_SH is explained by a constant practice effect (the label), by a slope on the change
# in SHMI, or by both? R2 is uncentred (a difference of zero is the null) and out of fold,
# each site held out in turn; unweighted least squares as in Fig. 3d. Shared and unique
# parts as in C1g: shared = label + SHMI - both, SHMI only = both - label, label only =
# both - SHMI. Per practice, and across practices (label = one effect per practice,
# SHMI = one slope).
r2u <- function(y, yhat) 1 - sum((y - yhat)^2) / sum(y^2)
fit_pred <- function(tr, te, model) {
  if (model == "label") {
    m <- tapply(tr$yi, tr$practice, mean); out <- unname(m[te$practice])
    out[is.na(out)] <- mean(tr$yi); out
  } else if (model == "SHMI") {
    sum(tr$d_SHMI * tr$yi) / sum(tr$d_SHMI^2) * te$d_SHMI
  } else {
    tr$pf <- factor(tr$practice); te$pf <- factor(te$practice, levels(tr$pf))
    f <- if (nlevels(tr$pf) > 1) lm(yi ~ 0 + pf + d_SHMI, data = tr) else lm(yi ~ d_SHMI, data = tr)
    out <- predict(f, newdata = te); out[is.na(out)] <- 0; out
  }
}
MODS <- c("label", "SHMI", "both")
pair_r2 <- function(e, label) {
  e <- as.data.frame(e)
  ins <- sapply(MODS, function(m) r2u(e$yi, fit_pred(e, e, m)))
  pred <- matrix(NA_real_, nrow(e), 3, dimnames = list(NULL, MODS))
  for (s in unique(e$site)) { te <- e$site == s
    for (m in MODS) pred[te, m] <- fit_pred(e[!te, ], e[te, ], m) }
  oof <- apply(pred, 2, function(p) r2u(e$yi, p))
  tibble(set = label, pairs = nrow(e), sites = n_distinct(e$site),
         R2_label_in = ins[["label"]], R2_SHMI_in = ins[["SHMI"]], R2_both_in = ins[["both"]],
         R2_label_oof = oof[["label"]], R2_SHMI_oof = oof[["SHMI"]], R2_both_oof = oof[["both"]]) %>%
    mutate(shared_oof = R2_label_oof + R2_SHMI_oof - R2_both_oof,
           SHMI_only_oof = R2_both_oof - R2_label_oof, label_only_oof = R2_both_oof - R2_SHMI_oof)
}
pr2 <- bind_rows(
  bind_rows(lapply(PRACTICES, function(p) { e <- es_sh[es_sh$practice == p, ]
    if (n_distinct(e$site) < 3) NULL else pair_r2(e, p) })),
  pair_r2(es_sh, "All practices"))
write.csv(mutate(pr2, across(where(is.double), ~ round(.x, 3))), file.path(OUT, "pairs_label_vs_SHMI.csv"), row.names = FALSE)
cat("\n== Strict pairs: label vs SHMI (as Fig. 3d) ==\n"); print(as.data.frame(pr2), digits = 2)

# ---- 6. plot data and panels --------------------------------------------------------------
pd <- list(pooled = pooled, metareg = mr, pair_r2 = pr2,
           pairs = bind_rows(mutate(as.data.frame(es_soc), outcome = "SOC (log response ratio)"),
                             mutate(as.data.frame(es_sh), outcome = "lv_SH (mean difference)")))
saveRDS(pd, file.path(OUT, "plotdata_pair_meta.rds"))

lab <- pooled %>% mutate(row = factor(practice, rev(PRACTICES)),
                         txt = sprintf("k = %d, %d sites, I² = %.0f%%", pairs, sites, I2),
                         panel = ifelse(grepl("SOC", outcome), "SOC, log response ratio", "Soil health (lv_SH), mean difference"))
pts <- pd$pairs %>% mutate(row = factor(practice, rev(PRACTICES)),
                           panel = ifelse(grepl("SOC", outcome), "SOC, log response ratio", "Soil health (lv_SH), mean difference"))
p_forest <- ggplot(lab, aes(y = row)) +
  geom_vline(xintercept = 0, linetype = 2, colour = "grey50") +
  geom_point(data = pts, aes(x = yi), colour = "grey60", alpha = 0.5, size = 1,
             position = position_jitter(height = 0.15, seed = 1)) +
  geom_errorbarh(aes(xmin = pi_lb, xmax = pi_ub), height = 0, linewidth = 0.4, colour = "grey30") +
  geom_errorbarh(aes(xmin = ci_lb, xmax = ci_ub), height = 0, linewidth = 1.4) +
  geom_point(aes(x = estimate), shape = 21, fill = "white", size = 2.6, stroke = 1) +
  geom_text(aes(x = Inf, label = txt), hjust = 1.02, vjust = -0.9, size = 2.4, colour = "grey30") +
  facet_wrap(~ panel, scales = "free_x") +
  labs(x = "Effect of the practice (pooled; thick, 95% CI; thin, 95% prediction interval)", y = NULL) +
  theme_bw(base_size = 9) + theme(panel.grid.minor = element_blank())
saveRDS(p_forest, file.path(PANEL_DIR, "panel_meta_forest.rds"))

# panel c: share of between-pair heterogeneity explained, on the same "% variance" axis as
# panel d; each bar is labelled with its slope per 10 points and the LRT P value, so a
# small share with a clear slope reads correctly (larger changes help, but dose explains
# little of the spread within pairs that hold the rest of the system fixed)
mrp <- mr %>% filter(outcome == "lv_SH (mean difference)", moderator != "SHMI of the system without the practice") %>%
  mutate(set = factor(set, rev(c(PRACTICES, "All practices (shared slope)"))),
                     moderator = factor(moderator, c("Change in targeted sub-index", "Change in SHMI")),
                     lab = sprintf("%+.1f per 10, P %s", slope_per10,
                                   ifelse(LRT_p < 0.001, "< 0.001", ifelse(LRT_p < 0.01, sprintf("= %.3f", LRT_p), sprintf("= %.2f", LRT_p)))))
p_reg <- ggplot(mrp, aes(x = 100 * tau2_explained, y = set, fill = moderator)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.75) +
  geom_text(aes(label = lab, x = 100 * tau2_explained), position = position_dodge(width = 0.8),
            hjust = -0.05, size = 2.3, colour = "grey25") +
  scale_fill_manual(values = c("Change in targeted sub-index" = "grey65", "Change in SHMI" = "#2f7fb5"), name = NULL) +
  scale_x_continuous(limits = c(0, 50), expand = expansion(mult = c(0, 0.02))) +
  labs(x = "Between-pair heterogeneity in soil-health effects explained (%)", y = NULL) +
  theme_bw(base_size = 9) + theme(panel.grid.minor = element_blank(), legend.position = "top")
saveRDS(p_reg, file.path(PANEL_DIR, "panel_meta_reg.rds"))
# ---- 7. the Soil Health Institute's pairs (rr_lookup.csv) against ours ---------------------
# rr_lookup lists SHI's control/treatment pairs by comparison type (no effect sizes). Each
# pair is classed against our pairing: strict (one of our strict pairs, either order),
# loose (same site, our with/without groups differ for that practice, other labels free),
# variant (both systems in the same group, e.g. two cover-crop mixtures), or not coded
# (comparison type we do not code, or a label missing). SOC log response ratios and lv_SH
# mean differences are then pooled per comparison type exactly as above, for all of SHI's
# US pairs and for the subset that are our strict pairs.
if (file.exists(F_RR)) {
  LABEL_OF <- function(ct) dplyr::case_when(grepl("cover", ct, ignore.case = TRUE) ~ "eu_cover_crops",
                                           grepl("till", ct, ignore.case = TRUE) ~ "eu_tillage",
                                           grepl("rotat|divers", ct, ignore.case = TRUE) ~ "eu_crop_rot",
                                           grepl("organic|manure|amend|compost|nutri", ct, ignore.case = TRUE) ~ "eu_nutr_type",
                                           TRUE ~ NA_character_)
  group_of <- list(
    eu_cover_crops = function(x) ifelse(x %in% c("grasses", "mix", "legume"), "A", ifelse(x == "none", "B", NA)),
    eu_nutr_type   = function(x) ifelse(x %in% c("organic", "synthetic_and_organic"), "A", ifelse(x == "synthetic", "B", NA)),
    eu_tillage     = function(x) ifelse(x %in% c("none", "minimum"), "A", ifelse(x %in% c("intense", "moderate"), "B", NA)),
    eu_crop_rot    = function(x) ifelse(x %in% c("2", "3", "greater_than_3"), "A", ifelse(x == "1", "B", NA)))
  trt <- read.csv(F_TRT, stringsAsFactors = FALSE) %>%
    mutate(MGT_combo = make_combo(sitecode, treatmentid),
           across(c(eu_tillage, eu_cover_crops, eu_crop_rot, eu_nutr_type), ~ tolower(trimws(as.character(.x))))) %>%
    select(MGT_combo, eu_tillage, eu_cover_crops, eu_crop_rot, eu_nutr_type)
  rr <- read.csv(F_RR, stringsAsFactors = FALSE) %>%
    transmute(practice = comparisontype, site = sitecode,
              unit_A = make_combo(sitecode, treatmentid_treatment), unit_B = make_combo(sitecode, treatmentid_control),
              label = LABEL_OF(comparisontype), pair = paste0("R", row_number()))
  cat("\nSHI comparison types (pairs):\n"); print(table(rr$practice, useNA = "ifany"))
  strict_keys <- with(pairs, c(paste(unit_A, unit_B), paste(unit_B, unit_A)))
  lab_A <- trt[match(rr$unit_A, trt$MGT_combo), ]; lab_B <- trt[match(rr$unit_B, trt$MGT_combo), ]
  gA <- gB <- rep(NA_character_, nrow(rr))
  for (i in seq_len(nrow(rr))) if (!is.na(rr$label[i])) {
    gA[i] <- group_of[[rr$label[i]]](lab_A[[rr$label[i]]][i]); gB[i] <- group_of[[rr$label[i]]](lab_B[[rr$label[i]]][i])
  }
  rr <- rr %>% mutate(region = ifelse(grepl("^MX", site), "MX", "US"),
                      in_ours = unit_A %in% units$MGT_combo & unit_B %in% units$MGT_combo,
                      class = case_when(paste(unit_A, unit_B) %in% strict_keys ~ "strict",
                                        is.na(label) | is.na(gA) | is.na(gB) ~ "not coded",
                                        gA != gB ~ "loose",
                                        TRUE ~ "variant"),
                      direction = case_when(class %in% c("strict", "loose") & gA == "A" ~ "with vs without",
                                            class %in% c("strict", "loose") ~ "reversed (control has the practice)",
                                            TRUE ~ NA_character_))
  write.csv(rr, file.path(OUT, "rr_pairs_classified.csv"), row.names = FALSE)
  cat("\nSHI pairs whose control has the practice (re-oriented to with vs without below):",
      sum(rr$direction == "reversed (control has the practice)", na.rm = TRUE), "\n")
  # orient strict and loose pairs as ours (with the practice vs without) so effects compare
  flip <- rr$direction %in% "reversed (control has the practice)"
  tmp <- rr$unit_A[flip]; rr$unit_A[flip] <- rr$unit_B[flip]; rr$unit_B[flip] <- tmp
  cat("\n== SHI pairs against our pairing (US) ==\n")
  print(as.data.frame(rr %>% filter(region == "US") %>% count(practice, class) %>%
                        pivot_wider(names_from = class, values_from = n, values_fill = 0)))
  cat("\nOur strict pairs found among SHI's:", sum(paste(pairs$unit_A, pairs$unit_B) %in% c(paste(rr$unit_A, rr$unit_B), paste(rr$unit_B, rr$unit_A))),
      "of", nrow(pairs), "\n")

  # effect sizes for SHI's US pairs (strict and loose: with vs without; variants: SHI's treatment vs control)
  rpe <- rr %>% filter(region == "US") %>%
    left_join(rename_with(stats, ~ paste0(.x, "_A")), by = c("unit_A" = "MGT_combo_A")) %>%
    left_join(rename_with(stats, ~ paste0(.x, "_B")), by = c("unit_B" = "MGT_combo_B"))
  okS <- with(rpe, soc_n_A >= 2 & soc_n_B >= 2 & soc_sd_A > 0 & soc_sd_B > 0 & soc_m_A > 0 & soc_m_B > 0); okS[is.na(okS)] <- FALSE
  okH <- with(rpe, sh_n_A >= 2 & sh_n_B >= 2 & sh_sd_A > 0 & sh_sd_B > 0); okH[is.na(okH)] <- FALSE
  rr_soc <- escalc("ROM", m1i = soc_m_A, sd1i = soc_sd_A, n1i = soc_n_A, m2i = soc_m_B, sd2i = soc_sd_B, n2i = soc_n_B, data = rpe[okS, ])
  rr_sh  <- escalc("MD",  m1i = sh_m_A,  sd1i = sh_sd_A,  n1i = sh_n_A,  m2i = sh_m_B,  sd2i = sh_sd_B,  n2i = sh_n_B,  data = rpe[okH, ])
  # pool only SHI's practice contrasts (types we code); its "crop count" and "cropping
  # system" lists are mostly every system against every other at a site (not coded here),
  # whose contrasts are largely linearly dependent and not practice comparisons
  rr_soc <- filter(rr_soc, class != "not coded"); rr_sh <- filter(rr_sh, class != "not coded")
  types <- sort(unique(rr_sh$practice))
  rr_pool <- bind_rows(
    mutate(pool(rr_soc, "soc", "ROM", "SOC (log response ratio)", types), pairs_used = "all SHI pairs"),
    mutate(pool(filter(rr_soc, class == "strict"), "soc", "ROM", "SOC (log response ratio)", types), pairs_used = "SHI pairs that are our strict pairs"),
    mutate(pool(rr_sh, "sh", "MD", "lv_SH (mean difference)", types), pairs_used = "all SHI pairs"),
    mutate(pool(filter(rr_sh, class == "strict"), "sh", "MD", "lv_SH (mean difference)", types), pairs_used = "SHI pairs that are our strict pairs")) %>%
    mutate(pct = ifelse(grepl("SOC", outcome), 100 * (exp(estimate) - 1), NA_real_)) %>%
    relocate(pairs_used, .after = practice)
  write.csv(mutate(rr_pool, across(where(is.double), ~ signif(.x, 4))), file.path(OUT, "rr_pooled.csv"), row.names = FALSE)
  cat("\n== SHI's pairs, pooled per comparison type (compare with our strict-pair results above) ==\n")
  print(as.data.frame(rr_pool %>% select(outcome, practice, pairs_used, pairs, independent, sites, estimate, ci_lb, ci_ub, pct, I2)), digits = 3)
}

# ---- 7b. scores against the meta-analysis (16 Nov 2026) --------------------------------------
# For the systems in each practice's strict matches, and in its practice matches (the systems of
# Fig. 2d), lv_SH is regressed within site (site means removed) on the yes/no label, on SHMI, and
# on the practice's own part of SHMI (Cover/Roots for cover crops, the amendment bonus for organic
# inputs, Inverse Disturbance for reduced tillage, Diversity for rotation). Reported: in-sample
# and leave-one-site-out within-site R2, the slopes, and the change each model implies for the
# strict pairs (slope x mean change in the pairs; for the label, its within-site difference),
# next to the pooled meta-analysis effect. Published SHMI weights; system means, not replicates.
OWN <- c("Cover crop vs none" = "Cover", "Organic vs synthetic inputs" = "Amend_bonus",
         "Reduced vs conventional tillage" = "InvDist", "Rotation vs monoculture" = "Diversity")
usys <- read.csv(F_UNITS, stringsAsFactors = FALSE)
within_fit <- function(z, x) {
  y <- z$SH - ave(z$SH, z$sitecode); xx <- z[[x]] - ave(z[[x]], z$sitecode)
  b <- sum(xx * y) / sum(xx^2)
  oof <- vapply(unique(z$sitecode), function(st) {           # leave one site out
    tr <- z$sitecode != st; bt <- sum(xx[tr] * y[tr]) / sum(xx[tr]^2)
    sum((y[!tr] - bt * xx[!tr])^2) }, 0)
  c(b = b, R2_in = 1 - sum((y - b * xx)^2) / sum(y^2), R2_oof = 1 - sum(oof) / sum(y^2))
}
sc_rows <- list()
for (p in names(OWN)) for (set in c("Strict match", "Practice match")) {
  sp <- supp6 %>% filter(Practice == p, set == "Practice match" | `Pair type` == "Strict match")
  A <- unique(sp$`System with practice`); B <- unique(sp$`System without`)
  z <- usys %>% filter(MGT_combo %in% c(A, B)) %>% mutate(label = as.numeric(MGT_combo %in% A))
  f <- lapply(c(label = "label", SHMI = "SHMI", part = OWN[[p]]), function(x) within_fit(z, x))
  m <- pooled %>% filter(outcome == "lv_SH (mean difference)", practice == p)
  sc_rows[[length(sc_rows) + 1]] <- tibble(
    practice = p, set, systems = nrow(z), sites = n_distinct(z$sitecode), own_part = OWN[[p]],
    R2_label_in = f$label[["R2_in"]], R2_SHMI_in = f$SHMI[["R2_in"]], R2_part_in = f$part[["R2_in"]],
    R2_label_oof = f$label[["R2_oof"]], R2_SHMI_oof = f$SHMI[["R2_oof"]], R2_part_oof = f$part[["R2_oof"]],
    slope10_SHMI = 10 * f$SHMI[["b"]], slope10_part = 10 * f$part[["b"]],
    mean_dSHMI_pairs = mean(sp$`Change in SHMI`), mean_dpart_pairs = mean(sp$`Change in targeted part`),
    implied_label = f$label[["b"]], implied_SHMI = f$SHMI[["b"]] * mean(sp$`Change in SHMI`),
    implied_part = f$part[["b"]] * mean(sp$`Change in targeted part`),
    meta_estimate = if (set == "Strict match") m$estimate else NA_real_,
    meta_ci_lb = if (set == "Strict match") m$ci_lb else NA_real_,
    meta_ci_ub = if (set == "Strict match") m$ci_ub else NA_real_)
}
scores_vs_meta <- bind_rows(sc_rows)
write.csv(mutate(scores_vs_meta, across(where(is.double), ~ round(.x, 3))), file.path(OUT, "scores_vs_meta.csv"), row.names = FALSE)
cat("\n== Label, SHMI and the practice's own part within site, against the meta-analysis ==\n")
print(as.data.frame(scores_vs_meta %>% select(practice, set, systems, R2_label_oof, R2_SHMI_oof, R2_part_oof,
                                              implied_label, implied_SHMI, implied_part, meta_estimate)), digits = 2)

# ---- 8. Supplementary Table 6 with the pairs used in Fig. 3a,b (16 Nov 2026) ---------------
used_key <- function(ok) paste(pe$unit_A[ok], pe$unit_B[ok])
k6 <- paste(supp6$`System with practice`, supp6$`System without`)
supp6$`Used in Fig. 2a,b (SOC / lv_SH)` <- ifelse(
  supp6$`Pair type` %in% "Strict match",
  paste(ifelse(k6 %in% used_key(ok_soc), "yes", "no"), "/", ifelse(k6 %in% used_key(ok_sh), "yes", "no")), "")
supp6 <- relocate(supp6, `Used in Fig. 2a,b (SOC / lv_SH)`, .after = `Pair type`)
write.csv(supp6, file.path(OUT, "SuppTable6_practice_pairs.csv"), row.names = FALSE, na = "")
cat("Supplementary Table 6 (with the pairs used in Fig. 3a,b) written to", OUT, "\n")

cat("\nSaved:", OUT, "and the two panels in", PANEL_DIR, "\n")
