################################################################################
# SHMI_yield_stability.R
# Is corn yield more stable where SHMI is higher? Research sites only.
#
# Uses the data and figure helpers of SHMI_yield_main.R (its sections 1-4), so the N
# accounting (LEGUME_SCENARIO, MANURE_SCENARIO) is whatever the main script is set to.
# Minnesota farms are not used: 3-5 corn years per farm are too few for stability.
#
# Every test extends the Fig 1 model (within sites, same year):
#   yield ~ site x year + N per site + other treatment factors
#           + SHMI_w + SHMI_w x N + SHMI_w x age + (1 | treatment/plot), AR(1)
# so N rate, system age and the non-SHMI treatments are accounted for throughout, and
# every slope is the within-site SHMI slope at 150 kg N/ha and the average system age.
#
#   S1  Sensitivity to the year (Finlay-Wilkinson): + SHMI_w x county yield deviation.
#       The year's conditions = the county's corn yield that year minus the site's mean
#       county yield (external to the plots, so no plot is regressed on itself). A
#       negative interaction = the SHMI benefit is larger in poor years (buffering);
#       it also means a smaller benefit in good years. Fig S1: the SHMI slope against
#       the county yield deviation, as Figs 2-3.
#   S2  Downside risk: for each plot with at least MIN_YEARS corn years, the 10th
#       percentile of its yield deviations from the site-year mean of all plots (how
#       far it falls in its worst years), and the same around the plot's own mean
#       deviation (the low tail without the plot's average advantage). Plot-level
#       model: q10 ~ site + N curve + other + SHMI_w + (1 | treatment). Fig S2: as
#       Fig 1, points = treatments with the same SHMI at a site.
#   S3  Drought: + SHMI_w x SPEI (continuous), and SHMI_w x drought year
#       (SPEI <= SPEI_DROUGHT) at sites with at least 3 drought years, for two SPEI-3
#       windows (May-July and June-August) and two site sets (all sites; rainfed sites
#       without the irrigated NEMEIRR and NECCIRR): four analyses. SPEI from
#       Outputs/yield/drought_years.csv (filled by Get_SPEI.R from SPEIbase). A
#       negative SPEI interaction = more SHMI benefit in drier years. Fig S3: the SHMI
#       slope against SPEI for each analysis, plus a summary panel.
#   S4  Variance: the Fig 1 model with residual SD = sigma x exp(delta x SHMI_w / 10); the
#       ratio exp(delta) is the change in year-to-year scatter per 10 SHMI (< 1 =
#       steadier), after the mean effects (unlike CV or SD, it does not mix the mean
#       yield into the variability). Results table only.
#
# Output: Outputs/yield_stability_<legume>_<manure> (+ "_agcros" with INCLUDE_AGCROS):
#   stability_results.csv, FigS1_sensitivity.png, FigS2_downside.png,
#   FigS3_drought.png (May-July, all sites), FigS3_drought_<mjj|jja>_<all|rainfed>.png,
#   FigS3_drought_summary.png/.csv, FigS2_points.csv
################################################################################

# ---- data and figure helpers from SHMI_yield_main.R (sections 1-4) -----------------
# Always reloaded, so the sites and settings (INCLUDE_AGCROS, AGCROS_SITES, N scenarios)
# match SHMI_yield_main.R rather than whatever `corn` is left in memory. Fast: the site
# runs come from the cache (Outputs/yield/site_runs).
local({
  txt <- readLines(here::here("R", "SHMI_yield_main.R"), encoding = "UTF-8")
  stop_at <- grep("^# 6\\. FIG 1-3", txt)[1] - 2   # sections 1-4
  eval(parse(text = txt[1:stop_at], encoding = "UTF-8"), envir = globalenv())
})
message("Stability: ", n_distinct(corn$site), " sites: ", paste(sort(unique(corn$site)), collapse = ", "))
OUTS <- paste0("Outputs/yield_stability_", LEGUME_SCENARIO, "_", MANURE_SCENARIO)
if (isTRUE(INCLUDE_AGCROS)) OUTS <- paste0(OUTS, "_agcros")
dir.create(OUTS, recursive = TRUE, showWarnings = FALSE)
message("Stability: N accounting legume = ", LEGUME_SCENARIO, ", manure = ", MANURE_SCENARIO, "; writing to ", OUTS)

MIN_YEARS    <- 8        # S2: corn years a plot needs for a 10th percentile
SPEI_DROUGHT <- -1.0     # S3: drought year = SPEI at or below this (moderate drought or worse)
set.seed(2026)

# ---- the Fig 1 model, extended ------------------------------------------------------
SHMI_TERMS <- c("SHMI_w", "SHMI_w:N_c", "SHMI_w:age_c")
fit_fig1 <- function(d, extra = NULL, weights = NULL, ar1 = TRUE) {
  d <- droplevels(d)
  f <- reformulate(c("site_year", n_terms(d), o_term(d), SHMI_TERMS, mei_terms(d), extra), "yield_kgha")
  lme(f, random = ~ 1 | trt/plot_id, correlation = if (ar1) corCAR1(form = ~ year | trt/plot_id),
      weights = weights, data = d, method = "REML", control = ctrl)
}
# SHMI slope per 10 SHMI at x of the moderator `int` (N_REF, average age)
slope_x <- function(m, int, x) {
  b <- fixef(m); V <- vcov(m)
  tibble(x = x, b = 10 * (b[["SHMI_w"]] + x * b[[int]]),
         se = 10 * sqrt(V["SHMI_w", "SHMI_w"] + x^2 * V[int, int] + 2 * x * V["SHMI_w", int]))
}
row_int <- function(m, int, estimate, figure, n) {
  tt <- summary(m)$tTable
  coef_row(figure, estimate, tt[int, "Value"], tt[int, "Std.Error"], tt[int, "DF"], tt[int, "p-value"], n)
}
row_at <- function(m, int, x, estimate, figure, n) {
  s <- slope_x(m, int, x); df <- summary(m)$tTable["SHMI_w", "DF"]
  coef_row(figure, estimate, s$b / 10, s$se / 10, df, 2 * pt(-abs(s$b / s$se), df), n)
}
# the slope against a moderator, with its 95% band (as Figs 2-3)
plot_mod <- function(m, int, xs, marks, xlab, title, cap, col, rug, xscale = 1, p_lab) {
  df <- summary(m)$tTable["SHMI_w", "DF"]
  sl <- slope_x(m, int, xs / xscale) %>% mutate(x = x * xscale,
          lower = b - qt(0.975, df) * se, upper = b + qt(0.975, df) * se)
  mk <- slope_x(m, int, marks / xscale) %>% mutate(x = x * xscale, lab = sprintf("%+.0f", b))
  ggplot(sl, aes(x, b)) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey60") +
    geom_ribbon(aes(ymin = lower, ymax = upper), fill = col, alpha = 0.18) +
    geom_line(colour = col, linewidth = 1.2) +
    geom_point(data = mk, size = 2.6, colour = col) +
    geom_text(data = mk, aes(label = lab), vjust = -1.1, size = 3) +
    geom_rug(data = tibble(x = rug), aes(x = x), inherit.aes = FALSE, sides = "b", alpha = 0.25, colour = "grey30") +
    annotate("text", x = -Inf, y = Inf, hjust = -0.05, vjust = 1.4, size = 2.8,
             label = paste0("Yield change per 10 SHMI, 95% band; ", lab_p(summary(m)$tTable[int, "p-value"], p_lab))) +
    labs(x = xlab, y = "Corn yield change per 10 SHMI (kg ha⁻¹)", title = title, caption = cap) + theme_fig
}
res <- list()

################################################################################
# S1. SENSITIVITY TO THE YEAR (Finlay-Wilkinson, county yield as the year's conditions)
################################################################################
d1 <- corn %>% filter(!is.na(county_kgha)) %>%
  group_by(site) %>% mutate(env_c = (county_kgha - mean(county_kgha)) / 1000) %>% ungroup()  # t/ha
m_fw <- fit_fig1(d1, "SHMI_w:env_c")
q_env <- quantile(d1$env_c, c(0.1, 0.5, 0.9))
res$fw <- bind_rows(
  row_int(m_fw, "SHMI_w:env_c", "SHMI x year conditions: change in slope per 1,000 kg/ha county yield", "S1", n_distinct(d1$site)),
  map_dfr(seq_along(q_env), function(k) row_at(m_fw, "SHMI_w:env_c", q_env[[k]],
    sprintf("Slope in a %s year (county yield %+.0f kg/ha vs the site's mean)", c("poor", "typical", "good")[k], 1000 * q_env[[k]]),
    "S1", n_distinct(d1$site))))
# sensitivity: the site-year mean of all plots as the year's conditions (includes the plot itself)
d1b <- corn %>% group_by(site_year) %>% mutate(sy_mean = mean(yield_kgha)) %>%
  group_by(site) %>% mutate(env2_c = (sy_mean - mean(sy_mean)) / 1000) %>% ungroup()
m_fw2 <- fit_fig1(d1b, "SHMI_w:env2_c")
res$fw2 <- row_int(m_fw2, "SHMI_w:env2_c", "Sensitivity: SHMI x site-year mean yield, per 1,000 kg/ha", "S1", n_distinct(d1b$site))

figS1 <- plot_mod(m_fw, "SHMI_w:env_c", seq(min(d1$env_c), max(d1$env_c), length.out = 100) * 1000,
                  round(1000 * q_env, -1), "County corn yield that year minus the site's mean (kg ha⁻¹)",
                  "SHMI benefit and the year's conditions (research sites)",
                  paste0("Within-site SHMI slope (Fig 1 model + SHMI × county yield deviation), at ", N_REF,
                         " kg N/ha and\naverage system age; 95% band. Points: 10th, 50th and 90th percentile years. ",
                         "A falling line = the SHMI\nbenefit is larger in poor years. Rug: site-years in the data."),
                  "#08519c", 1000 * unique(select(d1, site_year, env_c))$env_c, xscale = 1000, p_lab = "SHMI × year")
if (!exists("MS")) MS <- if (isTRUE(getOption("shmi.ms_figures", TRUE))) ggplot2::labs(title = NULL, subtitle = NULL, caption = NULL) else NULL
ggsave(file.path(OUTS, "FigS1_sensitivity.png"), figS1 + MS + PANEL, width = 6.8, height = 4.4, dpi = 300)

################################################################################
# S2. DOWNSIDE RISK (10th percentile of each plot's deviations from the site-year mean)
################################################################################
u <- corn %>% group_by(site_year) %>% mutate(dev = yield_kgha - mean(yield_kgha)) %>% ungroup() %>%
  group_by(site, trt, MGT_trt_base, plot_id) %>%
  summarise(n = n(), q10 = quantile(dev, 0.10), q10o = quantile(dev - mean(dev), 0.10), mean_dev = mean(dev),
            SHMI = first(SHMI), SHMI_w = first(SHMI_w), other = first(as.character(other)), N_m = mean(N100),
            .groups = "drop") %>%
  filter(n >= MIN_YEARS) %>%
  mutate(other = relevel(factor(other), ref = "none"))
cat("S2:", nrow(u), "plots with at least", MIN_YEARS, "corn years at", n_distinct(u$site), "sites\n")
for (s_ in unique(as.character(u$site))) {   # as the main model: linear N per site, quadratic at the N-rate sites
  if (diff(range(u$N_m[u$site == s_])) < 1e-8) next
  u[[paste0("U1_", s_)]] <- ifelse(u$site == s_, u$N_m, 0)
  if (s_ %in% N_QUAD_SITES) u[[paste0("U2_", s_)]] <- u[[paste0("U1_", s_)]]^2
}
fit_q <- function(resp) {
  rhs <- c("site", grep("^U[12]_", names(u), value = TRUE), if (nlevels(droplevels(u$other)) > 1) "other", "SHMI_w")
  lmer(reformulate(c(rhs, "(1 | trt)"), resp), data = droplevels(u), REML = TRUE)
}
mq  <- fit_q("q10"); mqo <- fit_q("q10o")
q_row <- function(m, estimate) {
  cf <- summary(m)$coefficients["SHMI_w", ]
  coef_row("S2", estimate, cf[["Estimate"]], cf[["Std. Error"]], cf[["df"]], cf[["Pr(>|t|)"]], n_distinct(u$site))
}
res$q10 <- bind_rows(q_row(mq,  "Worst years (10th percentile of deviations from the site-year mean), kg/ha"),
                     q_row(mqo, "Worst years around the plot's own mean deviation (low tail without the average advantage), kg/ha"))

# figure as Fig 1: partial residuals, points = treatments with the same SHMI at a site
pr_points <- function(m, resp, measure) {
  b <- fixef(m); X <- getME(m, "X")
  keep <- setdiff(colnames(X), c("(Intercept)", "SHMI_w"))
  Xc <- sweep(X[, keep, drop = FALSE], 2, colMeans(X[, keep, drop = FALSE]))
  dd <- droplevels(u) %>% mutate(y = .data[[resp]] - drop(Xc %*% b[keep]))
  pts <- dd %>% group_by(site, SHMI_level = round(SHMI)) %>%
    summarise(SHMI = mean(SHMI), y_m = mean(y), se = sd(y) / sqrt(n()), n_plots = n(), .groups = "drop") %>%
    rename(y = y_m) %>% mutate(measure = measure)
  cf <- summary(m)$coefficients["SHMI_w", ]; df <- cf[["df"]]
  x0 <- mean(dd$SHMI); y0 <- mean(dd$y)
  ln <- tibble(SHMI = seq(min(dd$SHMI), max(dd$SHMI), length.out = 100)) %>%
    mutate(fit = y0 + cf[["Estimate"]] * (SHMI - x0), se = cf[["Std. Error"]] * abs(SHMI - x0),
           lower = fit - qt(0.975, df) * se, upper = fit + qt(0.975, df) * se, measure = measure)
  lab <- tibble(measure = measure,
                txt = paste0(lab_slope(10 * cf[["Estimate"]], 10 * (cf[["Estimate"]] - qt(0.975, df) * cf[["Std. Error"]]),
                                       10 * (cf[["Estimate"]] + qt(0.975, df) * cf[["Std. Error"]]), "Change per 10 SHMI"),
                             "; ", lab_p(cf[["Pr(>|t|)"]])))
  list(pts = pts, ln = ln, lab = lab)
}
M1 <- "Worst years vs the site-year mean"; M2 <- "Worst years vs the plot's own average"
pp <- list(pr_points(mq, "q10", M1), pr_points(mqo, "q10o", M2))
pts <- bind_rows(map(pp, "pts")) %>% mutate(measure = factor(measure, levels = c(M1, M2)))
lns <- bind_rows(map(pp, "ln"))  %>% mutate(measure = factor(measure, levels = c(M1, M2)))
labs_ <- bind_rows(map(pp, "lab")) %>% mutate(measure = factor(measure, levels = c(M1, M2)))
write.csv(pts, file.path(OUTS, "FigS2_points.csv"), row.names = FALSE)
figS2 <- ggplot() +
  geom_hline(yintercept = 0, colour = "grey80") +
  geom_ribbon(data = lns, aes(SHMI, ymin = lower, ymax = upper), fill = "darkgreen", alpha = 0.15) +
  geom_line(data = lns, aes(SHMI, fit), colour = "darkgreen", linewidth = 1.1) +
  geom_errorbar(data = pts, aes(SHMI, ymin = y - se, ymax = y + se), width = 0, colour = "grey35") +
  geom_point(data = pts, aes(SHMI, y, fill = site), shape = 21, size = 2.4, colour = "black") +
  geom_text(data = labs_, aes(x = -Inf, y = Inf, label = txt), hjust = -0.04, vjust = 1.4, size = 2.5) +
  facet_wrap(~ measure, nrow = 1) +
  scale_fill_manual(values = site_cols, name = NULL) +
  labs(x = XLAB, y = "Worst-year deviation (kg ha⁻¹)",
       title = "SHMI and yield in the worst years (research sites)",
       caption = paste0("Each plot with at least ", MIN_YEARS, " corn years: the 10th percentile of its yield minus the mean of all plots at the site that\n",
                        "year (left), and the same around the plot's own mean deviation (right). Model: site + N per site + other treatment\n",
                        "factors + SHMI + (1 | treatment). Points: treatments with the same SHMI at a site (± SE over plots); higher = smaller losses.")) +
  theme_fig + theme(strip.background = element_blank(), strip.text = element_text(face = "bold"))
ggsave(file.path(OUTS, "FigS2_downside.png"),
       figS2 + MS + force_panelsizes(rows = unit(2.6, "in"), cols = unit(3.2, "in")), width = 9, height = 4.6, dpi = 300)

################################################################################
# S3. DROUGHT: two SPEI windows x all sites / rainfed sites only
################################################################################
# Windows (SPEIbase SPEI-3, filled by Get_SPEI.R):
#   mjj  May-July   (SPEI-3 for July; to silking)
#   jja  June-August (SPEI-3 for August; silking and early grain fill)
# Site sets: all sites, and rainfed sites only (NEMEIRR and NECCIRR are irrigated,
# where rainfall-based SPEI says little about crop water stress; the AgCROS sites
# INWLACRE, INWLTPAC, IAAM7071 and PAUP are rainfed).
# For each of the four analyses: the continuous SHMI x SPEI interaction (slopes at
# SPEI -1.5, 0, +1.5), and SHMI x drought year (SPEI <= SPEI_DROUGHT) at sites with at
# least 3 drought years. FigS3_drought.png keeps the original analysis (May-July,
# all sites); FigS3_drought_<window>_<set>.png for each analysis and
# FigS3_drought_summary.png comparing them.
SPEI_WINDOWS <- c(mjj = "May-July", jja = "June-August")
IRRIGATED    <- c("NEMEIRR", "NECCIRR")
SITE_SETS    <- list(all = NULL, rainfed = IRRIGATED)        # sites to leave out
dfile <- file.path(CACHE, "drought_years.csv")
if (!file.exists(dfile)) {
  write.csv(corn %>% distinct(site, year) %>% arrange(site, year) %>%
              mutate(spei_mjj = NA_real_, spei_jja = NA_real_), dfile, row.names = FALSE)
  message("S3 skipped: run Get_SPEI.R to fill ", dfile, " (SPEI per site-year), then rerun.")
} else {
  sp_all <- read.csv(dfile)
  no_spei <- setdiff(unique(corn$site), unique(sp_all$site[!is.na(sp_all$spei_jja)]))
  if (length(no_spei))
    warning("No SPEI for ", paste(no_spei, collapse = ", "), " in ", dfile,
            ": these sites drop out of S3. Rerun Get_SPEI.R with INCLUDE_AGCROS <- TRUE.", call. = FALSE)
  q_sp <- c(-1.5, 0, 1.5)
  s3_rows <- list(); s3_slopes <- list()
  for (w in names(SPEI_WINDOWS)) {
    col <- paste0("spei_", w)
    if (!col %in% names(sp_all) || all(is.na(sp_all[[col]]))) {
      message("S3: no ", col, " in ", dfile, " - rerun Get_SPEI.R"); next
    }
    for (set in names(SITE_SETS)) {
      lab <- sprintf("[%s SPEI, %s]", SPEI_WINDOWS[[w]], if (set == "all") "all sites" else "rainfed sites")
      d3 <- corn %>% filter(!site %in% SITE_SETS[[set]]) %>%
        inner_join(sp_all %>% filter(!is.na(.data[[col]])) %>% select(site, year, spei = all_of(col)),
                   by = c("site", "year")) %>%
        mutate(drought = as.integer(spei <= SPEI_DROUGHT))
      cat("\nS3", lab, "\n")
      print(d3 %>% distinct(site, year, drought) %>% group_by(site) %>% summarise(years = n(), drought_years = sum(drought)))
      m_sp <- tryCatch(fit_fig1(d3, "SHMI_w:spei"), error = function(e) { message("S3 ", lab, " failed: ", conditionMessage(e)); NULL })
      if (is.null(m_sp)) next
      ns <- n_distinct(d3$site)
      s3_rows[[length(s3_rows) + 1]] <- bind_rows(
        row_int(m_sp, "SHMI_w:spei", paste(lab, "SHMI x SPEI: change in slope per SPEI unit (negative = more benefit when drier)"), "S3", ns),
        map_dfr(q_sp, function(x) row_at(m_sp, "SHMI_w:spei", x, sprintf("%s Slope at SPEI %+.1f", lab, x), "S3", ns)))
      s3_slopes[[length(s3_slopes) + 1]] <- slope_x(m_sp, "SHMI_w:spei", q_sp) %>%
        mutate(window = SPEI_WINDOWS[[w]], set = if (set == "all") "All sites" else "Rainfed sites",
               df = summary(m_sp)$tTable["SHMI_w", "DF"], p_int = summary(m_sp)$tTable["SHMI_w:spei", "p-value"],
               p_main = summary(m_sp)$tTable["SHMI_w", "p-value"])   # SHMI slope at SPEI 0
      keep <- d3 %>% distinct(site, year, drought) %>% group_by(site) %>% summarise(k = sum(drought)) %>% filter(k >= 3) %>% pull(site)
      if (length(keep)) {
        m_dr <- tryCatch(fit_fig1(filter(d3, site %in% keep), "SHMI_w:drought"), error = function(e) NULL)
        if (!is.null(m_dr)) s3_rows[[length(s3_rows) + 1]] <- bind_rows(
          row_int(m_dr, "SHMI_w:drought", sprintf("%s SHMI x drought year (SPEI <= %.1f): drought minus other years", lab, SPEI_DROUGHT), "S3", length(keep)),
          row_at(m_dr, "SHMI_w:drought", 0, paste(lab, "Slope in other years"), "S3", length(keep)),
          row_at(m_dr, "SHMI_w:drought", 1, paste(lab, "Slope in drought years"), "S3", length(keep)))
      }
      fig <- plot_mod(m_sp, "SHMI_w:spei", seq(min(d3$spei), max(d3$spei), length.out = 100), q_sp,
                      paste(SPEI_WINDOWS[[w]], "SPEI (negative = drier than usual)"),
                      sprintf("SHMI benefit and drought (%s)", if (set == "all") "all research sites" else "rainfed research sites"),
                      paste0("Within-site SHMI slope (Fig 1 model + SHMI × ", SPEI_WINDOWS[[w]], " SPEI), at ", N_REF,
                             " kg N/ha and\naverage system age; 95% band. A line rising to the left = more SHMI benefit in dry years.\n",
                             "SPEIbase, 3-month SPEI for ", if (w == "mjj") "July" else "August", ". Rug: site-years."),
                      "#b35806", unique(select(d3, site_year, spei))$spei, p_lab = "SHMI × SPEI")
      ggsave(file.path(OUTS, sprintf("FigS3_drought_%s_%s.png", w, set)), fig + MS + PANEL, width = 6.8, height = 4.4, dpi = 300)
      if (w == "mjj" && set == "all") ggsave(file.path(OUTS, "FigS3_drought.png"), fig + MS + PANEL, width = 6.8, height = 4.4, dpi = 300)
    }
  }
  res$spei <- bind_rows(s3_rows)
  if (length(s3_slopes)) {
    # Panel headers: the SHMI slope at SPEI 0 (P for the SHMI term; all four are well above 0).
    # Reference: the overall slope of the Fig 1 model (green line and band), so a reader sees
    # at once that every estimate sits on it; the SHMI x SPEI P values go in the legend
    # (smallest printed in the console).
    m0 <- fit_fig1(corn)
    s0 <- slope_x(m0, "SHMI_w:N_c", 0); df0 <- summary(m0)$tTable["SHMI_w", "DF"]
    ref <- tibble(b = s0$b, lower = s0$b - qt(0.975, df0) * s0$se, upper = s0$b + qt(0.975, df0) * s0$se)
    ss <- bind_rows(s3_slopes) %>%
      mutate(lower = b - qt(0.975, df) * se, upper = b + qt(0.975, df) * se,
             cond = factor(sprintf("SPEI %+.1f", x), sprintf("SPEI %+.1f", q_sp)),
             analysis = factor(paste0(window, "\n", set, "\n", lab_p(p_main, "SHMI"))))
    message("S3: SHMI x SPEI interaction P values: ", paste(sprintf("%.3f", unique(ss$p_int)), collapse = ", "),
            " (smallest ", sprintf("%.3f", min(ss$p_int)), ")")
    figS3s <- ggplot(ss, aes(cond, b, colour = cond)) +
      annotate("rect", xmin = -Inf, xmax = Inf, ymin = ref$lower, ymax = ref$upper, fill = "darkgreen", alpha = 0.10) +
      geom_hline(yintercept = ref$b, colour = "darkgreen", linewidth = 0.6) +
      geom_hline(yintercept = 0, colour = "grey25", linewidth = 0.5) +
      geom_pointrange(aes(ymin = lower, ymax = upper), size = 0.4) +
      facet_wrap(~ analysis, nrow = 1) +
      scale_colour_manual(values = c("#b35806", "grey40", "#2166ac"), guide = "none") +
      labs(x = NULL, y = "Corn yield change per 10 SHMI (kg ha⁻¹)", title = "SHMI benefit in dry, average and wet growing seasons",
           caption = paste0("Within-site SHMI slope at ", N_REF, " kg N/ha and average system age, 95% CI. SPEI −1.5 = severe drought, +1.5 = very wet.\n",
                            "Green line and band: the overall slope (Fig 1 model). Headers: P for the SHMI slope at SPEI 0. ",
                            "Rainfed sites exclude the irrigated NEMEIRR and NECCIRR.")) + theme_fig +
      theme(strip.background = element_blank(), strip.text = element_text(face = "bold"))
    ggsave(file.path(OUTS, "FigS3_drought_summary.png"), figS3s + MS, width = 9, height = 4.4, dpi = 300)
    write.csv(ss %>% select(window, set, cond, b, se, lower, upper, p_main, p_int), file.path(OUTS, "FigS3_drought_summary.csv"), row.names = FALSE)
  }
}

################################################################################
# S4. VARIANCE: does year-to-year scatter change with SHMI?
################################################################################
# SHMI_w is scaled to tens of points (SHMI10) so the variance parameter stays small: on
# the raw scale exp(delta x SHMI_w) over +/- 25 points overflowed the optimizer. If the
# model with AR(1) still fails, it is refitted without AR(1) (the variance question is
# about the residual scatter, which AR(1) only reweights). p: likelihood-ratio test
# against the same model with constant variance.
corn4 <- corn %>% mutate(SHMI10 = SHMI_w / 10)
s4 <- NULL
for (ar in c(TRUE, FALSE)) {
  m_var <- tryCatch(fit_fig1(corn4, weights = varExp(form = ~ SHMI10), ar1 = ar),
                    error = function(e) { message("S4 (AR1 = ", ar, ") failed: ", conditionMessage(e)); NULL })
  if (is.null(m_var)) next
  m_con <- tryCatch(fit_fig1(corn4, ar1 = ar), error = function(e) NULL)
  ci <- tryCatch(intervals(m_var, which = "var-cov")$varStruct, error = function(e) NULL)
  dl <- coef(m_var$modelStruct$varStruct, unconstrained = FALSE)[["expon"]]
  lr <- if (!is.null(m_con)) anova(m_con, m_var) else NULL
  s4 <- tibble(figure = "S4", estimate = paste0("Residual SD ratio per 10 SHMI (< 1 = steadier yields)",
                                              if (!ar) ", model without AR(1)" else ""),
               effect = exp(dl), se = NA_real_, df = NA_real_,
               lower = if (!is.null(ci)) exp(ci[1, "lower"]) else NA_real_,
               upper = if (!is.null(ci)) exp(ci[1, "upper"]) else NA_real_,
               p_val = if (!is.null(lr)) lr[2, "p-value"] else NA_real_, n = n_distinct(corn$site))
  cat("S4 (AR1 =", ar, "): SD ratio per 10 SHMI", round(exp(dl), 3), "\n")
  break
}
res$var <- s4

################################################################################
# RESULTS
################################################################################
out <- bind_rows(res) %>% mutate(across(c(effect, se, lower, upper), ~ round(.x, 2)), p_val = signif(p_val, 3),
                                 legume = LEGUME_SCENARIO, manure = MANURE_SCENARIO)
print(out, n = Inf, width = Inf)
write.csv(out, file.path(OUTS, "stability_results.csv"), row.names = FALSE)
# units: S1-S3 kg/ha per 10 SHMI (interactions: per 10 SHMI per unit of the moderator);
# S2 kg/ha of the 10th-percentile deviation per 10 SHMI; S4 a ratio
message("Done: stability results and figures in ", OUTS)
