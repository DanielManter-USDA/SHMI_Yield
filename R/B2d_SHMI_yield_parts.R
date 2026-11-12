################################################################################
# B2d_SHMI_yield_parts.R
#
# Does the whole score predict corn yield better than its parts? The yield version
# of Fig. 2a. B2's within-site model is refitted with SHMI replaced by:
#   each part alone (Cover/Roots, Diversity, Inverse Disturbance, the share of years
#   amended), equal weights, and, if B2c has been run, SHMI from the current season only.
# The four parts are also fitted together and reported as the weights yield would choose
# (in-sample, so not compared with the fixed scores).
# Same plot-years, same terms: site x year, N per site, other factors, X_w, X_w x N_c,
# X_w x age_c; random treatment/plot; AR(1). X_w = X minus its site mean.
# Compared by:
#   AIC (maximum likelihood; lower is better; delta from SHMI), with its range when each
#     experiment is left out in turn (set B2D_NO_LOO <- TRUE to skip; that part is slow)
#   the standardized within-site slope: kg/ha per within-site SD of the score, at
#     150 kg N/ha and the average age (for the free model: the combined change in
#     yield from a one-SD rise in its fitted combination)
#   experiments with a positive slope in their own model (of 10)
# B2D_REDRAW <- TRUE before sourcing redraws the figure from the saved yield_parts.csv (seconds).
#
# Input   Outputs/yield_main_default_total_agcros/B2_data.rds (B2)
#         Outputs/yield_record_length/rolling_shmi.rds (B2c; optional)
# Outputs Outputs/yield_parts/ yield_parts.csv, yield_parts_sites.csv, yield_parts_loo.csv,
#         yield_weights_free.csv (the four parts fitted together, in-sample: the weights yield
#         would choose, not a competitor), Fig_yield_parts.png/.pdf
#
# Rotation check (B2D_ROTATION <- TRUE before sourcing): every model also gets the
# preceding legume (legume_class from B2: soybean, alfalfa, ... or none) as a factor, so the
# direct rotation effect on yield (beyond the legume N credit already in N) is not left to
# the scores. Outputs go to Outputs/yield_parts_rotation/. If Diversity's advantage over
# SHMI in the main run is the rotation effect, it shrinks here.
################################################################################
suppressPackageStartupMessages({ library(dplyr); library(tidyr); library(purrr); library(nlme); library(ggplot2) })
setwd(here::here())
IN   <- "Outputs/yield_main_default_total_agcros/B2_data.rds"
ROLL <- "Outputs/yield_record_length/rolling_shmi.rds"
ROTATION <- exists("B2D_ROTATION") && isTRUE(B2D_ROTATION)
OUT  <- if (ROTATION) "Outputs/yield_parts_rotation" else "Outputs/yield_parts"
if (!file.exists(IN)) stop(IN, " not found: run B2_SHMI_yield_main.R first")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
ctrl <- lmeControl(opt = "optim", maxIter = 200, msMaxIter = 200)
W <- c(Cover = 0.687, Diversity = 0.150, InvDist = 0.163); AMEND_GAP <- 0.60
# 11 Nov 2026: the ~170 model fits run in parallel. B2D_WORKERS sets the number of R
# processes (default: cores - 1, at most 6; 1 = one after another, as before). The fits and
# results are the same either way.
WORKERS <- if (exists("B2D_WORKERS")) B2D_WORKERS else max(1, min(6, parallel::detectCores() - 1))
par_map <- function(X, FUN) {
  if (WORKERS <= 1 || length(X) < 2) return(lapply(X, FUN))
  cl <- parallel::makeCluster(min(WORKERS, length(X)))
  on.exit(parallel::stopCluster(cl))
  parallel::clusterEvalQ(cl, suppressPackageStartupMessages({ library(dplyr); library(purrr); library(nlme); library(tibble) }))
  e <- environment(FUN)
  parallel::clusterExport(cl, ls(e, all.names = TRUE), envir = e)
  parallel::parLapply(cl, X, FUN)
}

REDRAW <- exists("B2D_REDRAW") && isTRUE(B2D_REDRAW)   # TRUE: only redraw the figure from yield_parts.csv
if (!REDRAW) {
D <- readRDS(IN); corn <- D$corn
if (ROTATION) {
  if (!"legume_class" %in% names(corn)) stop("B2's corn data has no legume_class column: rerun B2_SHMI_yield_main.R")
  corn <- corn %>% mutate(prev_crop = factor(ifelse(is.na(legume_class), "none", as.character(legume_class))))
  message("Rotation check: preceding legume added to every model; ", paste(capture.output(print(table(corn$prev_crop))), collapse = " "))
}
xtra <- function(dd) if (ROTATION) fac_term(dd, "prev_crop") else character()

# ---- the parts: column names as build_shmi() returns them (the first match is used) --
pick <- function(cands, what) {
  hit <- cands[cands %in% names(corn)][1]
  if (is.na(hit)) stop("No column for ", what, " in B2's corn data. Columns are:\n  ",
                       paste(names(corn), collapse = ", "), call. = FALSE)
  hit
}
col_cover <- pick(c("Cover", "Cover_Roots", "CoverRoots", "cover"), "Cover/Roots")
col_div   <- pick(c("Diversity", "diversity"), "Diversity")
col_inv   <- pick(c("InvDist", "Inverse_Disturbance", "invdist"), "Inverse Disturbance")
col_amend <- pick(c("p_amend", "amend_share", "OrgInput_amend"), "share of years amended")
message("Parts: ", paste(c(col_cover, col_div, col_inv, col_amend), collapse = ", "))
corn <- corn %>% mutate(P_cover = .data[[col_cover]], P_div = .data[[col_div]], P_inv = .data[[col_inv]],
                        P_amend = 100 * .data[[col_amend]],
                        EQ_base = (P_cover + P_div + P_inv) / 3,
                        P_equal = EQ_base + AMEND_GAP * .data[[col_amend]] * (100 - EQ_base),
                        # 10 Nov 2026: the variants of Supplementary Table 1, with the published
                        # weights (the remaining parts' weights rescaled to sum to 1; the amendment
                        # term applied as in SHMI). SHMI_base: SHMI without the amendment term.
                        P_base   = W[["Cover"]] * P_cover + W[["Diversity"]] * P_div + W[["InvDist"]] * P_inv,
                        B_nocov  = (W[["Diversity"]] * P_div + W[["InvDist"]] * P_inv) / (W[["Diversity"]] + W[["InvDist"]]),
                        B_nodiv  = (W[["Cover"]] * P_cover + W[["InvDist"]] * P_inv) / (W[["Cover"]] + W[["InvDist"]]),
                        B_noinv  = (W[["Cover"]] * P_cover + W[["Diversity"]] * P_div) / (W[["Cover"]] + W[["Diversity"]]),
                        P_nocov  = B_nocov + AMEND_GAP * .data[[col_amend]] * (100 - B_nocov),
                        P_nodiv  = B_nodiv + AMEND_GAP * .data[[col_amend]] * (100 - B_nodiv),
                        P_noinv  = B_noinv + AMEND_GAP * .data[[col_amend]] * (100 - B_noinv))
# check: SHMI rebuilt from its parts with these formulas (should be ~1; if not, the variants
# above are approximations of the package's scoring and the log says so)
chk <- cor(corn$SHMI, corn$P_base + AMEND_GAP * corn[[col_amend]] * (100 - corn$P_base), use = "complete.obs")
message(sprintf("SHMI rebuilt from its parts: r = %.4f", chk))

# current-season SHMI (B2c, K = 1), if available
has_k1 <- file.exists(ROLL)
if (has_k1) {
  k1 <- readRDS(ROLL) %>% filter(K == 1) %>% select(MGT_combo, year, P_season = SHMI)
  corn <- corn %>% left_join(k1, by = c("MGT_combo", "year"))
}
CAND <- c("SHMI" = "SHMI", "Cover/Roots alone" = "P_cover", "Diversity alone" = "P_div",
          "Inverse Disturbance alone" = "P_inv", "Share of years amended alone" = "P_amend",
          "SHMI, equal weights" = "P_equal", "SHMI_base (no amendment term)" = "P_base",
          "SHMI without Cover/Roots" = "P_nocov", "SHMI without Diversity" = "P_nodiv",
          "SHMI without Inverse Disturbance" = "P_noinv",
          if (has_k1) c("SHMI, current season only" = "P_season"))
d <- corn %>% filter(if_all(all_of(unname(CAND)), ~ !is.na(.x))) %>% droplevels()
message(nrow(d), " of ", nrow(corn), " corn plot-years have every score")
for (v in unname(CAND)) d <- d %>% group_by(site) %>% mutate("{v}_w" := .data[[v]] - mean(.data[[v]])) %>% ungroup()

# ---- model (B2's terms; aliased columns dropped as in B2b/B2c) ------------------------
varies <- function(d, v) { d <- droplevels(d); any(tapply(as.numeric(d[[v]]), d$site_year, function(x) diff(range(x))) > 1e-8, na.rm = TRUE) }
n_terms <- function(d) { v <- grep("^N[12]_", names(d), value = TRUE); v[vapply(v, function(k) varies(d, k), logical(1))] }
fac_term <- function(d, v) if (n_distinct(d[[v]]) > 1 && varies(d, v)) v else character()
nm <- function(x) make.names(gsub("[^A-Za-z0-9_.]", "_", x), unique = TRUE)
fit <- function(d, rhs, method = "ML") {
  X <- model.matrix(reformulate(rhs), d)
  qx <- qr(X, tol = 1e-7); keep <- sort(qx$pivot[seq_len(qx$rank)])
  X <- X[, keep, drop = FALSE]; colnames(X) <- nm(colnames(X))
  dd <- data.frame(yield_kgha = d$yield_kgha, trt = d$trt, plot_id = d$plot_id, year = d$year, X, check.names = FALSE)
  lme(reformulate(c("0", colnames(X)), "yield_kgha"), random = ~ 1 | trt/plot_id,
      correlation = corCAR1(form = ~ year | trt/plot_id), data = dd, method = method, control = ctrl)
}
base_rhs <- c("site_year", n_terms(d), fac_term(d, "other"), xtra(d))
inter <- function(x) c(x, if (varies(d, "N_c")) paste0(x, ":N_c"), paste0(x, ":age_c"))

# one candidate on data dd: AIC (ML) and the yield change per within-site SD (single scores)
# 10 Nov 2026: also likelihood-ratio tests (ML) of nested models against SHMI: does the
# candidate add to SHMI (SHMI vs SHMI + candidate), and does SHMI add to the candidate
# (candidate vs SHMI + candidate)? Each score enters with its N and age interactions.
NESTED <- !exists("B2D_NO_NESTED")
terms_of <- function(xw, dd) c(xw, if (varies(dd, "N_c")) paste0(xw, ":N_c"), paste0(xw, ":age_c"))
lrt <- function(m0, m1) {
  a <- anova(m0, m1); tibble(chisq = a$L.Ratio[2], df = a$df[2] - a$df[1], p = a$`p-value`[2])
}
one <- function(label, x, dd = d, nested = FALSE) {
  xw <- paste0(x, "_w")
  rhs <- c("site_year", n_terms(dd), fac_term(dd, "other"), xtra(dd))
  m  <- fit(dd, c(rhs, terms_of(xw, dd)), "ML")
  b  <- fixef(m); V <- vcov(m); cn <- nm(xw); sdx <- sd(dd[[xw]])
  est <- b[[cn]] * sdx; se <- sqrt(V[cn, cn]) * sdx
  out <- tibble(score = label, AIC = AIC(m), slope_per_sd = est, se = se, lower = est - 1.96 * se, upper = est + 1.96 * se,
                p_slope = 2 * pnorm(-abs(est / se)))
  if (nested && x != "SHMI") {
    m_ref  <- fit(dd, c(rhs, terms_of("SHMI_w", dd)), "ML")
    m_both <- tryCatch(fit(dd, c(rhs, terms_of("SHMI_w", dd), terms_of(xw, dd)), "ML"), error = function(e) NULL)
    if (!is.null(m_both) && length(fixef(m_both)) > length(fixef(m_ref))) {
      a1 <- lrt(m_ref, m_both); a2 <- lrt(m, m_both)
      out <- out %>% mutate(p_adds_to_SHMI = a1$p, df_adds = a1$df, p_SHMI_adds = a2$p, df_SHMI_adds = a2$df)
    } else out <- out %>% mutate(p_adds_to_SHMI = NA_real_, p_SHMI_adds = NA_real_)   # aliased with SHMI
  }
  out
}
message("Fitting ", length(CAND), " scores (", WORKERS, " R processes)")
res <- bind_rows(par_map(seq_along(CAND), function(i) one(names(CAND)[i], unname(CAND)[i], nested = NESTED))) %>%
  mutate(dAIC_vs_SHMI = AIC - AIC[score == "SHMI"],
         pct_of_SHMI = 100 * slope_per_sd / slope_per_sd[score == "SHMI"])

# weights yield would choose: all four parts entered together (in-sample; reported, not plotted)
xw4 <- c("P_cover_w", "P_div_w", "P_inv_w", "P_amend_w")
m_free <- fit(d, c(base_rhs, unlist(lapply(xw4, inter))), "REML")
tt <- summary(m_free)$tTable[nm(xw4), ]
free <- tibble(part = c("Cover/Roots", "Diversity", "Inverse Disturbance", "Share of years amended (per % point)"),
               per_point = tt[, "Value"], se = tt[, "Std.Error"], p = tt[, "p-value"],
               per_sd = tt[, "Value"] * sapply(xw4, function(x) sd(d[[x]])))
three <- free$per_point[1:3]
free$implied_weight <- c(three / sum(three), NA)
free$shmi_weight <- c(W[["Cover"]], W[["Diversity"]], W[["InvDist"]], NA)
write.csv(mutate(free, across(where(is.double), ~ signif(.x, 4))), file.path(OUT, "yield_weights_free.csv"), row.names = FALSE)
cat("\nParts entered together (in-sample), and the weights they imply:\n"); print(free, width = Inf)

# each experiment's own model: is the slope positive?
sites <- bind_rows(par_map(seq_along(CAND), function(i) { v <- unname(CAND)[i]; lab <- names(CAND)[i]
 split(d, d$site) %>% imap_dfr(function(ds, s) {
  ds <- droplevels(ds); x <- paste0(v, "_w")
  if (sd(ds[[x]]) < 1e-8) return(NULL)
  m <- tryCatch(fit(ds, c("site_year", n_terms(ds), fac_term(ds, "other"), xtra(ds), x), "REML"),
                error = function(e) { message("  ", s, " / ", lab, ": ", conditionMessage(e)); NULL })
  if (is.null(m)) return(NULL)
  tt <- summary(m)$tTable[nm(x), ]
  tibble(score = lab, site = s, slope_per_10 = 10 * tt[["Value"]], p = tt[["p-value"]])
 })
}))
pos <- sites %>% group_by(score) %>% summarise(experiments = n(), positive = sum(slope_per_10 > 0),
                                               positive_p05 = sum(slope_per_10 > 0 & p < 0.05), .groups = "drop")

# leave one experiment out: how much each comparison depends on a single experiment
# (every candidate refitted without it; 10 x candidates ML fits, the slow part)
DO_LOO <- !exists("B2D_NO_LOO")
if (DO_LOO) {
  message("Leaving out each experiment in turn (", WORKERS, " R processes)")
  loo <- bind_rows(par_map(sort(unique(as.character(d$site))), function(s) {
    dd <- droplevels(filter(d, site != s))
    bind_rows(unname(imap(CAND, function(v, lab) one(lab, v, dd)))) %>%
      mutate(left_out = s, dAIC_vs_SHMI = AIC - AIC[score == "SHMI"])
  }))
  write.csv(mutate(loo, across(where(is.double), ~ signif(.x, 4))), file.path(OUT, "yield_parts_loo.csv"), row.names = FALSE)
  rng <- loo %>% group_by(score) %>%
    summarise(dAIC_loo_min = min(dAIC_vs_SHMI), dAIC_loo_max = max(dAIC_vs_SHMI),
              slope_loo_min = min(slope_per_sd), slope_loo_max = max(slope_per_sd), .groups = "drop")
  res <- res %>% left_join(rng, by = "score")
}
res <- res %>% left_join(pos, by = "score")
cat("\nFor the Extended Data Fig. 4 legend, AIC minus SHMI's (leaving out each experiment in turn):\n")
for (i in seq_len(nrow(res))) cat(sprintf("  %-30s %+6.1f%s\n", res$score[i], res$dAIC_vs_SHMI[i],
    if (DO_LOO) sprintf(" (%+.0f to %+.0f)", res$dAIC_loo_min[i], res$dAIC_loo_max[i]) else ""))
write.csv(mutate(res, across(where(is.double), ~ signif(.x, 4))), file.path(OUT, "yield_parts.csv"), row.names = FALSE)
write.csv(mutate(sites, across(where(is.double), ~ signif(.x, 4))), file.path(OUT, "yield_parts_sites.csv"), row.names = FALSE)
cat("\nCorn yield within sites: SHMI against its parts (same", nrow(d), "plot-years)\n"); print(res, width = Inf)

} else {                                   # redraw only, from the saved results
  res <- read.csv(file.path(OUT, "yield_parts.csv"), check.names = FALSE, stringsAsFactors = FALSE)
  DO_LOO <- "dAIC_loo_min" %in% names(res)
  message("Figure redrawn from ", file.path(OUT, "yield_parts.csv"))
}

# ---- figure (Extended Data Fig. 4b, assembled by B7) ---------------------------------------------------
# two panels; the AIC comparison (yield_parts.csv) is reported in the legend
# bars: range when each experiment is left out in turn (if run), else the slope's 95% CI
# 6 Nov 2026: the figure drops the current-season score (B2D_SHOW_SEASON) and the count of
# experiments with a positive slope (B2D_SHOW_COUNT; it raised the question of significance and
# Extended Data Fig. 5a now shows each experiment's slope); yield_parts.csv keeps every row.
if (!exists("B2D_SHOW_SEASON")) B2D_SHOW_SEASON <- FALSE
if (!exists("B2D_SHOW_COUNT"))  B2D_SHOW_COUNT  <- FALSE
if (!exists("B2D_SHOW_EQUAL"))  B2D_SHOW_EQUAL  <- FALSE   # equal weights: reported in the legend / Supp. Data Table 1
res_fig <- res %>% filter(score != "All parts entered freely", B2D_SHOW_SEASON | score != "SHMI, current season only",
                          B2D_SHOW_EQUAL | score != "SHMI, equal weights",
                          !score %in% c("SHMI_base (no amendment term)", "SHMI without Cover/Roots",
                                        "SHMI without Diversity", "SHMI without Inverse Disturbance"))
res <- res_fig
lev <- rev(res$score)
lo_s <- if (DO_LOO) res$slope_loo_min else res$lower; hi_s <- if (DO_LOO) res$slope_loo_max else res$upper
P1 <- "Experiments with a positive slope"
P3 <- "Yield change per within-site SD of the score (kg/ha)"
fd <- bind_rows(
  res %>% transmute(score, panel = P1, est = positive, lo = NA, hi = NA, txt = sprintf("%d of %d", positive, experiments)),
  res %>% transmute(score, panel = P3, est = slope_per_sd, lo = lo_s, hi = hi_s)) %>%
  filter(B2D_SHOW_COUNT | panel != P1) %>%
  mutate(panel = factor(panel, c(P1, P3)), score = factor(score, lev),
         grp = case_when(score == "SHMI" ~ "SHMI", grepl("alone", score) ~ "Single part", TRUE ~ "Other"))
g <- ggplot(fd, aes(est, score, colour = grp)) +
  geom_vline(data = filter(fd, score == "SHMI"), aes(xintercept = est), colour = "#1b6ca8", linewidth = 0.3) +
  geom_errorbar(aes(xmin = lo, xmax = hi), width = 0, orientation = "y", na.rm = TRUE) +
  geom_point(size = 2.4) +
  geom_text(aes(label = txt), hjust = -0.35, size = 3, na.rm = TRUE, show.legend = FALSE) +
  facet_wrap(~ panel, nrow = 1, scales = "free_x") +
  scale_x_continuous(expand = expansion(mult = c(0.08, 0.25))) +
  scale_colour_manual(values = c(SHMI = "#1b6ca8", `Single part` = "#8a8a8a", Other = "#e6a141"),
                      breaks = c("SHMI", "Single part", "Other"), name = NULL) +
  labs(x = NULL, y = NULL, caption = if (DO_LOO) "Bars: range with each experiment left out in turn" else "Bars: 95% CI") +
  theme_bw(base_size = 10) +
  theme(legend.position = "bottom", strip.background = element_blank(), strip.text = element_text(face = "bold"),
        panel.spacing.x = unit(1.5, "lines"))                     # keeps axis labels of adjacent panels apart
W_FIG <- if (B2D_SHOW_COUNT) 8 else 5
ggsave(file.path(OUT, "Fig_yield_parts.png"), g, width = W_FIG, height = 3.3, dpi = 300)
ggsave(file.path(OUT, "Fig_yield_parts.pdf"), g, width = W_FIG, height = 3.3)
message("Yield parts comparison done: ", OUT)
