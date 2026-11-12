################################################################################
# B2b_SHMI_soybean.R
#
# Does SHMI also predict soybean yield? The corn analysis (B2_SHMI_yield_main.R)
# repeated for soybean at the experiments that grew it, with the same scores, plots,
# treatments and model, except that soybean receives no N, so the model has no N
# terms and no SHMI x N interaction:
#   yield ~ site x year + other + N level (NEMLTCRS, a treatment factor applied to
#           the other crops) + SHMI_w + SHMI_w:age_c + (1 | treatment/plot), AR(1)
# SHMI_w = SHMI minus the site mean over the units that grew soybean.
#
# Why it matters: if the corn benefit came mainly through N supply (largest at low
# fertilizer N), a legume that fixes its own N should gain less per SHMI point,
# relative to its yield, than corn does. Both crops are compared as a percentage of
# their mean yield, at the same experiments, with the corn model refitted on those
# experiments (its slope at 150 kg N/ha, as in Fig. 3).
#
# Input   Outputs/yield_main_default_total_agcros/B2_data.rds (written by B2)
# Outputs Outputs/yield_soybean/
#   soy_results.csv      pooled slopes, soybean and corn at the same experiments
#                        (kg/ha and % of mean yield per 10 SHMI, 95% CI, P)
#   soy_site_slopes.csv  each experiment's own slope, both crops; pooled by random
#                        effects (Knapp-Hartung) with I2
#   soy_data_summary.csv experiments, systems, plot-years, SHMI range for soybean
#   soy_treatment_yields.csv soybean yield by treatment (used by B8b)
#   Fig_soybean_a.png    panel a alone, at the corn panel's size (B7 uses it as Fig. 2b)
#   Fig_soybean.png/.pdf a: soybean yield against SHMI within sites (as Fig. 3a);
#                        b: per-experiment and pooled slopes, % of mean yield per 10 SHMI,
#                        soybean and corn
################################################################################
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(purrr); library(nlme); library(metafor); library(ggplot2)
})
setwd(here::here())
IN  <- "Outputs/yield_main_default_total_agcros/B2_data.rds"
OUT <- "Outputs/yield_soybean"
if (exists("SOY_IN_OVERRIDE")) IN <- SOY_IN_OVERRIDE
if (!file.exists(IN)) stop(IN, " not found: run B2_SHMI_yield_main.R first")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

RESID_CUT   <- 4      # as in B2
MIN_LEVELS  <- 2      # an experiment needs at least 2 SHMI levels among its soybean units ...
MIN_PY      <- 30     # ... at least this many soybean plot-years ...
MIN_RANGE   <- 2      # ... and an SHMI range of at least this many points among them (PAUP's soybean
                      #     units span 85.3-85.6: no SHMI contrast to estimate a slope from)
N_REF       <- 150
ctrl <- lmeControl(opt = "optim", maxIter = 200, msMaxIter = 200)

D <- readRDS(IN)
runs <- D$runs; site_start <- D$site_start; corn <- D$corn

# ---- soybean rows, built as B2 builds the corn rows ----------------------------------
# treatment code -> parts, per experiment (the regexes of B2)
TRT <- list(NELITCSE  = list(re = "^(CH|DK|MP|NT|RT|SS)(CC|CSB|SBC|SB)$", parts = c("tillage", "rotation")),
            NEMLTCRS  = list(re = "^ROT([0-9])([0-9])$",               parts = c("rotation", "N_code")),
            OHTVD_Ho  = list(re = "(R[0-9]+)(T[0-9]+)(P[0-9]+)",        parts = c("rotation", "tillage", "phase")),
            OHTVD_Wo  = list(re = "(R[0-9]+)(T[0-9]+)(P[0-9]+)",        parts = c("rotation", "tillage", "phase")),
            WICST_ARS = list(re = "(CS[0-9]+)(P[0-9]+)",                parts = c("rotation", "phase")),
            NEMEIRR   = list(re = "(NT|DK)([0-9]+)",                    parts = c("tillage", "residue")),
            NECCIRR   = list(re = "(Full|Dfct)([XCM])([0-9]+)([NA])",   parts = c("irrigation", "amelioration", "N_rate", "residue")))
crop_rows <- function(r, re, parts, crop = "^soybean") {
  r$inputs$yield %>%
    filter(grepl(crop, tolower(trimws(CD_name))), !is.na(yield_kg_ha)) %>%
    group_by(MGT_combo, year) %>% summarise(yield_kgha = mean(yield_kg_ha), .groups = "drop") %>%
    mutate(year = as.integer(year)) %>%
    inner_join(r$shmi, by = "MGT_combo", relationship = "many-to-many") %>%
    filter(year >= period_start, year <= period_end) %>%
    tidyr::extract(MGT_trt, into = parts, regex = re, remove = FALSE) %>%
    mutate(MGT_field = as.character(MGT_field), MGT_trt_base = sub("P[0-9]+$", "", MGT_trt))
}
agcros_desc <- function(loc) {
  u <- read.csv(file.path(D$AGCROS_DIR, paste0("SHMI_", loc, "_units.csv")), stringsAsFactors = FALSE)
  u %>% transmute(MGT_combo, treatment_descriptor,
                  system = paste(treatment_descriptor, rotation_descriptor, tillage_descriptor, n_treatment_descriptor,
                                 ifelse(is.na(cover_crop) | cover_crop == "", "no cover", cover_crop),
                                 paste("residue removal", residue_removal),
                                 ifelse(organic_management %in% "Yes", "organic", "conventional"), sep = " | "),
                  residue = residue_removal, organic = organic_management %in% "Yes")
}
soy <- imap_dfr(runs, function(r, s) {
  if (s %in% names(TRT)) {
    d <- crop_rows(r, TRT[[s]]$re, TRT[[s]]$parts)
    if (!nrow(d)) return(NULL)
    d <- switch(s,
      NELITCSE  = d %>% filter(!is.na(tillage)) %>%
        mutate(MGT_trt_base = paste0(tillage, ifelse(rotation %in% c("CSB", "SBC"), "CS", rotation)), other = "none"),
      NEMLTCRS  = d %>% filter(year <= 2006) %>% mutate(other = "none"),
      # Ohio establishment years excluded as in B2 (de Camargo Santos et al. 2025)
      OHTVD_Wo  = d %>% filter(year >= 1964) %>% mutate(other = "none"),
      OHTVD_Ho  = d %>% filter(year >= 1965, year != 1981) %>% mutate(other = "none"),   # 1981 flooded, replanted (field notes)
      WICST_ARS = d %>% mutate(other = ifelse(rotation %in% c("CS3", "CS5"), "WICST_organic", "none")),
      d %>% mutate(other = "none"))
  } else {                                                     # AgCROS sites: systems from the descriptors
    d <- crop_rows(r, "(.+)", "trt_code")
    if (!nrow(d)) return(NULL)
    d <- d %>% left_join(agcros_desc(s), by = "MGT_combo") %>%
      mutate(MGT_trt_base = system,
             other = case_when(organic & residue %in% c("Partial", "Yes") ~ paste0(s, "_organic_res", residue),
                               organic ~ paste0(s, "_organic"),
                               residue %in% c("Partial", "Yes") ~ paste0(s, "_res", residue),
                               TRUE ~ "none"))
    if (s %in% D$AGCROS_CONVENTIONAL_ONLY) d <- filter(d, !organic)
  }
  d %>% mutate(site = s)
})
if (!nrow(soy)) stop("No soybean yields in the site runs")

soy <- soy %>%
  mutate(nlev = ifelse(site == "NEMLTCRS", paste0("NEMLTCRS_N", N_code), "none"),   # N level of the other crops
         other = relevel(factor(other), ref = "none"), nlev = relevel(factor(nlev), ref = "none"),
         plot_id = factor(paste(site, MGT_field, sep = "_")),
         site_year = factor(paste(site, year, sep = "_")),
         trt = factor(paste(site, MGT_trt_base)),
         age = unname(year - site_start[as.character(site)])) %>%
  group_by(site) %>%
  mutate(SHMI_w = SHMI - mean(SHMI[!duplicated(MGT_combo)])) %>% ungroup()

summ <- soy %>% group_by(site) %>%
  summarise(systems = n_distinct(MGT_trt_base), shmi_levels = n_distinct(round(SHMI)), plot_years = n(),
            years = paste(range(year), collapse = "-"), SHMI_min = min(SHMI), SHMI_max = max(SHMI),
            mean_yield = mean(yield_kgha), .groups = "drop") %>%
  mutate(used = shmi_levels >= MIN_LEVELS & plot_years >= MIN_PY & SHMI_max - SHMI_min >= MIN_RANGE)
print(summ, width = Inf)
soy <- soy %>% filter(site %in% summ$site[summ$used]) %>% droplevels()
if (any(duplicated(soy[c("plot_id", "year")]))) soy$plot_id <- factor(soy$MGT_combo)

# ---- models ---------------------------------------------------------------------------
# (droplevels: a split keeps every site-year level, whose empty groups would give NA)
varies <- function(d, v) { d <- droplevels(d); any(tapply(as.numeric(d[[v]]), d$site_year, function(x) diff(range(x))) > 1e-8, na.rm = TRUE) }
fac_term <- function(d, v) if (n_distinct(d[[v]]) > 1 && varies(d, v)) v else character()
# The fixed effects are expanded to a model matrix and any column that is aliased with
# the others (e.g. an "other" or N level that coincides with a site-year or with SHMI
# differences in the soybean data) is dropped, with a message, so lme does not fail with
# "Singularity in backsolve". SHMI terms are never dropped: if one is aliased, the
# experiment has no SHMI contrast left and the fit stops.
# Coefficient names: ":" becomes "_" (SHMI_w:age_c -> SHMI_w_age_c).
nm <- function(x) make.names(gsub("[^A-Za-z0-9_.]", "_", x), unique = TRUE)
fit <- function(d, rhs) {
  d <- droplevels(d)
  X <- model.matrix(reformulate(rhs), d)
  qx <- qr(X, tol = 1e-7)
  keep <- sort(qx$pivot[seq_len(qx$rank)])
  drop <- setdiff(colnames(X), colnames(X)[keep])
  shmi_cols <- grep("^SHMI_w", colnames(X), value = TRUE)
  if (any(shmi_cols %in% drop)) stop("SHMI term aliased with the other fixed effects: ", paste(intersect(shmi_cols, drop), collapse = ", "))
  if (length(drop)) message("  dropped aliased columns: ", paste(drop, collapse = ", "))
  X <- X[, keep, drop = FALSE]
  colnames(X) <- nm(colnames(X))
  dd <- data.frame(yield_kgha = d$yield_kgha, trt = d$trt, plot_id = d$plot_id, year = d$year, X, check.names = FALSE)
  f <- reformulate(c("0", colnames(X)), "yield_kgha")
  lme(f, random = ~ 1 | trt/plot_id, correlation = corCAR1(form = ~ year | trt/plot_id), data = dd,
      method = "REML", control = ctrl)
}
soy_rhs <- function(d, extra = character()) c("site_year", fac_term(d, "other"), fac_term(d, "nlev"), "SHMI_w", extra)

# residual screening, per experiment, as in B2
soy <- split(soy, soy$site) %>% map_dfr(function(ds) {
  ds <- arrange(ds, plot_id, year)
  m <- tryCatch(fit(ds, soy_rhs(ds)), error = function(e) NULL)
  ds$resid_n <- if (is.null(m)) 0 else as.numeric(residuals(m, type = "normalized"))
  ds
})
message(sum(abs(soy$resid_n) > RESID_CUT), " soybean plot-years set aside (|normalized residual| > ", RESID_CUT, ")")
soy <- soy %>% filter(abs(resid_n) <= RESID_CUT) %>% select(-resid_n) %>% droplevels()
soy$age_c <- soy$age - mean(soy$age)

# soybean yield by treatment (for the treatment table, B8b): plain mean and mean relative to
# the site-year mean of all soybean treatments (removes year)
soy_trt <- soy %>% group_by(site, year) %>% mutate(dev = yield_kgha - mean(yield_kgha)) %>% ungroup() %>%
  mutate(trt = as.character(trt)) %>% group_by(site, trt) %>%
  summarise(soy_MGT_trt = paste(sort(unique(as.character(MGT_trt))), collapse = ", "),
            soy_SHMI = round(mean(SHMI[!duplicated(MGT_combo)]), 1), soy_years = n_distinct(year),
            soy_plot_years = n(), soy_yield_mean = round(mean(yield_kgha)), soy_yield_sd = round(sd(yield_kgha)),
            soy_yield_vs_siteyr = round(mean(dev)), .groups = "drop")
write.csv(soy_trt, file.path(OUT, "soy_treatment_yields.csv"), row.names = FALSE)

per10 <- function(m, term = "SHMI_w") {
  tt <- summary(m)$tTable
  b <- 10 * tt[term, "Value"]; se <- 10 * tt[term, "Std.Error"]; df <- tt[term, "DF"]
  tibble(effect = b, se = se, df = df, lower = b - qt(0.975, df) * se, upper = b + qt(0.975, df) * se,
         p = tt[term, "p-value"])
}
m_soy <- fit(soy, soy_rhs(soy, "SHMI_w:age_c"))
soy_row <- per10(m_soy) %>% mutate(crop = "Soybean", sites = n_distinct(soy$site), plot_years = nrow(soy),
                                   mean_yield = mean(soy$yield_kgha))

# corn at the same experiments: the B2 model (N terms, SHMI x N, SHMI x age), slope at
# 150 kg N/ha and the average age of these data
cs <- corn %>% filter(site %in% unique(soy$site)) %>% droplevels() %>% mutate(age_c = age - mean(age))
n_terms <- function(d) { v <- grep("^N[12]_", names(d), value = TRUE); v[vapply(v, function(k) varies(d, k), logical(1))] }
corn_rhs <- c("site_year", n_terms(cs), fac_term(cs, "other"), "SHMI_w",
              if (varies(cs, "N_c")) "SHMI_w:N_c", "SHMI_w:age_c")
m_corn <- fit(cs, corn_rhs)
corn_row <- per10(m_corn) %>% mutate(crop = "Corn", sites = n_distinct(cs$site), plot_years = nrow(cs),
                                     mean_yield = mean(cs$yield_kgha))

# the soybean slope is reported at the average age of the soybean data (age_c = 0);
# the age interaction is given too
age_row <- { tt <- summary(m_soy)$tTable["SHMI_w_age_c", ]
  tibble(crop = "Soybean", term = "SHMI x age (per 10 SHMI per year)", effect = 10 * tt[["Value"]],
         se = 10 * tt[["Std.Error"]], p = tt[["p-value"]]) }

res <- bind_rows(soy_row, corn_row) %>%
  mutate(pct = 100 * effect / mean_yield, pct_lower = 100 * lower / mean_yield, pct_upper = 100 * upper / mean_yield) %>%
  select(crop, sites, plot_years, mean_yield, effect, lower, upper, p, pct, pct_lower, pct_upper, se, df)

# ---- each experiment's own slope, both crops -------------------------------------------
site_slope <- function(d, rhs_fun) {
  m <- tryCatch(fit(d, rhs_fun(d)), error = function(e) { message("  ", d$site[1], ": site model failed: ", conditionMessage(e)); NULL })
  if (is.null(m)) return(NULL)
  per10(m) %>% mutate(mean_yield = mean(d$yield_kgha), plot_years = nrow(d))
}
soy_sites <- split(soy, soy$site) %>% imap_dfr(function(d, s) {
  r <- site_slope(d, soy_rhs); if (is.null(r)) NULL else mutate(r, site = s, crop = "Soybean") })
corn_sites <- split(cs, cs$site) %>% imap_dfr(function(d, s) {
  rhs <- function(d) c("site_year", n_terms(d), fac_term(d, "other"), "SHMI_w")   # slope at the site's own N
  r <- site_slope(d, rhs); if (is.null(r)) NULL else mutate(r, site = s, crop = "Corn")
})
sites <- bind_rows(soy_sites, corn_sites) %>%
  mutate(pct = 100 * effect / mean_yield, pct_se = 100 * se / mean_yield,
         pct_lower = 100 * lower / mean_yield, pct_upper = 100 * upper / mean_yield)
pool <- sites %>% group_by(crop) %>% group_modify(function(g, k) {
  if (nrow(g) < 2) return(tibble())
  m <- rma(yi = pct, sei = pct_se, data = g, method = "REML", test = "knha")
  tibble(pooled_pct = as.numeric(m$b), pooled_lower = m$ci.lb, pooled_upper = m$ci.ub, pooled_p = m$pval,
         I2 = m$I2, k = m$k)
}) %>% ungroup()

write.csv(mutate(res, across(where(is.double), ~ signif(.x, 4))), file.path(OUT, "soy_results.csv"), row.names = FALSE)
write.csv(mutate(sites, across(where(is.double), ~ signif(.x, 4))), file.path(OUT, "soy_site_slopes.csv"), row.names = FALSE)
write.csv(mutate(pool, across(where(is.double), ~ signif(.x, 4))), file.path(OUT, "soy_site_pooled.csv"), row.names = FALSE)
write.csv(mutate(summ, across(where(is.double), ~ round(.x, 1))), file.path(OUT, "soy_data_summary.csv"), row.names = FALSE)
cat("\nPer 10 SHMI points (same experiments):\n"); print(res, width = Inf)
cat("\nSoybean SHMI x age:\n"); print(age_row, width = Inf)
cat("\nPooled per-experiment slopes (% of mean yield per 10 SHMI):\n"); print(pool, width = Inf)

# ---- figure: % of mean yield per 10 SHMI, per experiment and pooled -------------------
lev <- rev(c(sort(unique(sites$site)), "Pooled (one model)"))
fd <- bind_rows(sites %>% select(site, crop, pct, pct_lower, pct_upper),
                res %>% transmute(site = "Pooled (one model)", crop, pct, pct_lower, pct_upper)) %>%
  mutate(site = factor(site, lev), crop = factor(crop, c("Corn", "Soybean")))
g <- ggplot(fd, aes(pct, site, colour = crop)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey55") +
  geom_errorbar(aes(xmin = pct_lower, xmax = pct_upper), width = 0.25, position = position_dodge(width = 0.55)) +
  geom_point(aes(shape = site == "Pooled (one model)"), size = 2.6, position = position_dodge(width = 0.55)) +
  scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 18), guide = "none") +
  scale_colour_manual(values = c(Corn = "#b8860b", Soybean = "#2e7d32"), name = NULL) +
  labs(x = "Yield change per 10 SHMI (% of the crop's mean yield)", y = NULL) +
  theme_bw(base_size = 10) + theme(legend.position = "bottom", panel.grid.minor = element_blank())
# ---- figure a: soybean yield against SHMI within sites (as Fig. 3a for corn) ----------
# points: each experiment's soybean plot-years at one SHMI level (whole points), yields
# with the model's other terms removed (partial residuals at the average age) plus the
# mean soybean yield; line: the pooled within-site slope with its 95% band, pivoting on
# the mean SHMI (site x year sets each experiment's level, so only the slope is estimated)
site_cols <- c(NELITCSE = "#0072B2", NEMLTCRS = "#D55E00", OHTVD_Ho = "#56B4E9", OHTVD_Wo = "#CC79A7",
               WICST_ARS = "#009E73", NEMEIRR = "#E69F00", NECCIRR = "#8C510A", MNMOFS = "#999999",
               INWLACRE = "#000000", INWLTPAC = "#882255", IAAM7071 = "#F0E442", PAUP = "#6A3D9A")   # as in B2
tt_s <- summary(m_soy)$tTable["SHMI_w", ]
b_s <- tt_s[["Value"]]; se_s <- tt_s[["Std.Error"]]; df_s <- tt_s[["DF"]]
piv <- mean(soy$SHMI); y0 <- mean(soy$yield_kgha)
soy$y_adj <- as.numeric(resid(m_soy, level = 0)) + b_s * (soy$SHMI - piv) + y0
# Points as in Fig. 2a (10 Nov 2026): one per management system at each site, its mean SHMI
# and mean adjusted yield over its plot-years; treatments with the same SHMI record (differing
# only in factors SHMI does not score) merged. FIG3A_POINTS = "level" restores whole-point bins.
if (!exists("FIG3A_POINTS")) FIG3A_POINTS <- "system"
pts <- if (FIG3A_POINTS == "system") {
  soy %>% group_by(site, MGT_trt_base) %>% mutate(SHMI_level = round(mean(SHMI), 1)) %>%
    group_by(site, SHMI_level) %>%
    summarise(SHMI = mean(SHMI), y = mean(y_adj), se = sd(y_adj) / sqrt(n()), .groups = "drop")
} else {
  soy %>% group_by(site, SHMI_level = round(SHMI)) %>%
    summarise(SHMI = mean(SHMI), y = mean(y_adj), se = sd(y_adj) / sqrt(n()), .groups = "drop")
}
ln <- tibble(SHMI = seq(min(soy$SHMI), max(soy$SHMI), length.out = 200)) %>%
  mutate(fit = y0 + b_s * (SHMI - piv), half = qt(0.975, df_s) * se_s * abs(SHMI - piv),
         lower = fit - half, upper = fit + half)
lab_a <- sprintf("Yield change per 10 SHMI: %+.0f kg/ha (95%% CI %.0f\u2013%.0f)",
                 10 * b_s, 10 * (b_s - qt(0.975, df_s) * se_s), 10 * (b_s + qt(0.975, df_s) * se_s))
ga <- ggplot() +
  geom_ribbon(data = ln, aes(SHMI, ymin = lower, ymax = upper), fill = "#2e7d32", alpha = 0.15) +
  geom_line(data = ln, aes(SHMI, fit), colour = "#2e7d32", linewidth = 1.1) +
  geom_errorbar(data = pts, aes(SHMI, ymin = y - se, ymax = y + se), width = 0, colour = "grey35") +
  geom_point(data = pts, aes(SHMI, y, fill = site), shape = 21, size = 2.6, colour = "black") +
  scale_fill_manual(values = site_cols, name = NULL) +
  # corner labels as in B2's corn panel (slope top-left, P bottom-right, size 2.8)
  annotate("text", x = -Inf, y = Inf, hjust = -0.1, vjust = 1.1, size = 2.8, label = lab_a) +
  annotate("text", x = Inf, y = -Inf, hjust = 1.1, vjust = -1, size = 2.8,
           label = ifelse(tt_s[["p-value"]] < 0.001, "p < 0.001", sprintf("p = %.3f", tt_s[["p-value"]]))) +
  labs(x = "Soil Health Management Index (SHMI)",
       y = "Soybean yield adjusted for site, year\nand other treatments (kg ha\u207B\u00B9)") +
  theme_bw(base_size = 10)
write.csv(pts, file.path(OUT, "Fig_soybean_points.csv"), row.names = FALSE)
# plot data for R/figures/F3_yield.R (Fig. 3b), 11 Nov 2026
saveRDS(list(points = pts, line = ln, points_mode = FIG3A_POINTS, site_cols = site_cols,
             slope = 10 * b_s, lower = 10 * (b_s - qt(0.975, df_s) * se_s), upper = 10 * (b_s + qt(0.975, df_s) * se_s),
             p = tt_s[["p-value"]], n_sites = n_distinct(soy$site)),
        file.path(OUT, "plotdata_soybean.rds"))
# panel a on its own at the corn panel's size and fonts (B2: 3.6 x 2.6 in panel, 6.8 x 4.6 in,
# base size 10), so B7 can set it beside corn as Fig. 2b without cropping or shrinking
ga_alone <- ga + if (requireNamespace("ggh4x", quietly = TRUE))
  ggh4x::force_panelsizes(rows = unit(2.6, "in"), cols = unit(3.6, "in"))
ggsave(file.path(OUT, "Fig_soybean_a.png"), ga_alone, width = 6.8, height = 4.6, dpi = 300, bg = "white")

# b: the comparison with corn (% of mean yield), panel g above
gb <- g + theme(legend.position = "bottom")
if (requireNamespace("patchwork", quietly = TRUE)) {
  fig <- patchwork::wrap_plots(ga, gb, ncol = 1, heights = c(1, 0.9)) + patchwork::plot_annotation(tag_levels = "a") &
    theme(plot.tag = element_text(face = "bold", size = 14))
  h <- 4.2 + 1.6 + 0.42 * length(lev)
} else { fig <- ga; h <- 4.2 }
ggsave(file.path(OUT, "Fig_soybean.png"), fig, width = 7, height = h, dpi = 300)
ggsave(file.path(OUT, "Fig_soybean.pdf"), fig, width = 7, height = h)
message("Soybean test done: ", OUT)
