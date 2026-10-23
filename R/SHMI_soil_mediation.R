################################################################################
# SHMI_soil_mediation.R
#
# Is part of the within-site SHMI -> corn yield slope carried by measured soil health?
#
#   a path   mediator ~ SHMI          (treatments within a site; site fixed effects)
#   b path   yield ~ ... + SHMI + mediator   (the Fig. 3a model plus the mediator)
#   c        yield ~ ... + SHMI           (same model and same plot-years, without it)
#   indirect = a x b; direct = c'; share mediated = a x b / (a x b + c')
# CIs for the indirect effect: Monte Carlo method (20,000 draws of a and b from their
# estimates and standard errors; Preacher & Selig 2012). Results are per 10 SHMI points
# at N_REF (150 kg N/ha) and the average system age, as in the manuscript.
#
# Soil data
#  * NAPESHM sampled three of the yield experiments: Hoytville (USOH02) and Wooster
#    (USOH03), 9 tillage x rotation treatments each, and WICST (USWI01), 5 cropping
#    systems with corn. These give lv_SH, SHAPE and the individual indicators (0-15 cm,
#    one sampling ~2019) for 23 treatments. napeshm_mediators_crosswalk.csv holds the
#    values; fill its 'unit' column with the yield data's MGT_trt_base.
#  * AgCROS: SOC by plot, layer and year at the other sites (soil_c.csv; optional).
# Mediators: lv_SH (primary: the construct SHMI was built to predict), SHAPE mean,
# SOC 0-15 cm (pooled NAPESHM + AgCROS, the only one available at all sites),
# POXC and aggregate stability; SOC 0-30 cm and stock as AgCROS sensitivities.
#
# Interpretation: SHMI was assigned by the experiments, but soil properties were not,
# so the mediation is descriptive ("consistent with"), not a causal decomposition.
# Higher yields also return more residue and can raise SOC (reverse path). lv_SH at
# the three NAPESHM sites was part of SHMI's calibration, so its a path is partly
# in-sample; the b path (soil -> yield) is not. Sites are few -> exploratory.
#
# INPUT (b) Outputs/yield/soil_c.csv, one row per sample or plot x date x layer:
#   site            site code as in SHMI_yield_main.R (e.g. NEMEIRR)
#   unit            plot (MGT_field) or treatment (MGT_trt_base) as in the SHMI workbook;
#                   set UNIT_LEVEL below to match
#   year            sampling year
#   top_cm, bottom_cm   layer depths (cm)
#   soc_gkg         soil ORGANIC C, g/kg (use total C only where soils are not calcareous)
#   bd_gcm3         bulk density (optional; for the stock sensitivity)
# build_soil_c_from_agcros() below helps build this file from AgCROS/GRACEnet
# downloads; check its column matching on the first site.
#
# RUN: after SHMI_yield_main.R in the same R session (uses corn, n_terms, o_term,
#      mei_terms, ctrl, N_REF, OUT); sources it if corn is missing.
# OUTPUT (OUT): mediation_results.csv, mediation_units.csv, Fig_mediation.png
################################################################################
if (!exists("corn") || !exists("m_full")) source(here::here("R", "SHMI_yield_main.R"))
suppressPackageStartupMessages({ library(dplyr); library(tidyr); library(ggplot2); library(nlme); library(readxl) })

NAPESHM_FILE <- "Outputs/yield/napeshm_mediators_crosswalk.csv"  # NAPESHM soil data at OHTVD_Ho, OHTVD_Wo, WICST_ARS
SOIL_FILE   <- "Outputs/yield/soil_c.csv"                      # AgCROS SOC by plot and layer: NEMLTCRS, NELITCSE, NECCIRR, NEMEIRR
UNIT_LEVEL  <- "plot"        # AgCROS rows: "plot" = MGT_field ; "treatment" = MGT_trt_base
DEPTH_MAIN  <- 15            # cm: primary layer 0-15 (SHAPE convention; surface responds to tillage)
DEPTH_SENS  <- 30            # cm: sensitivity 0-30 (concentration, and stock if bd is available)
MIN_COVER   <- 0.8           # share of 0-DEPTH that sampled layers must cover
MIN_UNITS   <- 3             # units with SOC per site
YEAR_WINDOW <- NA            # e.g. 5: use only corn years within 5 years of a SOC sampling
NMC         <- 20000
set.seed(7)

# ============================================================================
# Helper: tidy AgCROS / GRACEnet soil tables into soil_c.csv
# Exports differ by site and template version, so columns are matched by pattern and
# printed; check them. 'crosswalk' maps the AgCROS experimental unit ID to the SHMI
# workbook unit: data.frame(site, eu_id, unit).
# ============================================================================
build_soil_c_from_agcros <- function(file, sheet, site, crosswalk, out = SOIL_FILE, append = TRUE) {
  x <- if (grepl("\\.csv$", file)) read.csv(file, check.names = FALSE) else read_excel(file, sheet = sheet)
  pick <- function(pat, required = TRUE, exclude = "^$") {
    k <- grep(pat, names(x), ignore.case = TRUE, value = TRUE)
    k <- k[!grepl(exclude, k, ignore.case = TRUE)]
    if (!length(k)) { if (required) stop("no column matching ", pat, " in: ", paste(names(x), collapse = " | ")) else return(NA) }
    k[1]
  }
  cols <- c(eu = pick("exp.*unit|^eu|plot.*id"), date = pick("date|year"),
            top = pick("upper|top"), bot = pick("lower|bottom"),
            oc = pick("organic.*c|^soc", FALSE, exclude = "nitrogen|inorg"),
            tc = pick("total.*c(arbon)?", FALSE, exclude = "nitrogen|inorg|\\bN\\b"),
            bd = pick("bulk.*dens", FALSE))
  cat("Columns used for", site, ":\n"); print(cols)
  carbon <- if (!is.na(cols[["oc"]])) cols[["oc"]] else { warning(site, ": no organic C column; using total C"); cols[["tc"]] }
  yr <- x[[cols[["date"]]]]
  yr <- if (inherits(yr, c("Date", "POSIXt"))) as.integer(format(yr, "%Y")) else as.integer(substr(as.character(yr), 1, 4))
  d <- data.frame(site = site, eu_id = as.character(x[[cols[["eu"]]]]), year = yr,
                  top_cm = as.numeric(x[[cols[["top"]]]]), bottom_cm = as.numeric(x[[cols[["bot"]]]]),
                  soc_gkg = as.numeric(x[[carbon]]),
                  bd_gcm3 = if (!is.na(cols[["bd"]])) as.numeric(x[[cols[["bd"]]]]) else NA_real_)
  # units: g C/kg expected; % values (< 10 throughout) are converted
  if (max(d$soc_gkg, na.rm = TRUE) < 10) { message(site, ": SOC looks like %, converting to g/kg"); d$soc_gkg <- 10 * d$soc_gkg }
  d <- inner_join(d, crosswalk %>% mutate(eu_id = as.character(eu_id)), by = c("site", "eu_id")) %>%
    select(site, unit, year, top_cm, bottom_cm, soc_gkg, bd_gcm3) %>% filter(!is.na(soc_gkg))
  if (append && file.exists(out)) d <- bind_rows(read.csv(out) %>% filter(site != !!site), d)
  write.csv(d, out, row.names = FALSE); invisible(d)
}

# ============================================================================
# 1. unit-level mediator tables
#    (a) NAPESHM: one row per treatment (napeshm_mediators_crosswalk.csv; fill 'unit'
#        with the MGT_trt_base used in the yield data, e.g.
#        unique(corn$MGT_trt_base[corn$site == "OHTVD_Ho"]); rows with blank unit are dropped)
#    (b) AgCROS: SOC by layer (soil_c.csv), depth-weighted to 0-D cm, as the deviation
#        from the site x sampling-year mean, averaged over sampling years
#    Every mediator enters the models as its deviation from the site mean (M).
# ============================================================================
napeshm <- if (file.exists(NAPESHM_FILE)) read.csv(NAPESHM_FILE) %>%
  mutate(unit = trimws(sub("\\[.*\\]", "", as.character(unit)))) %>% filter(!is.na(unit), unit != "") else NULL
if (!is.null(napeshm)) message("NAPESHM mediators: ", nrow(napeshm), " treatments at ", paste(unique(napeshm$site), collapse = ", "))
soil <- if (file.exists(SOIL_FILE)) read.csv(SOIL_FILE) else NULL

layer_mean <- function(s, D) {
  # layers lying wholly within 0-D cm (a 0-30 cm composite is not used for 0-15 cm)
  s %>% filter(bottom_cm <= D) %>% mutate(w = bottom_cm - pmax(top_cm, 0)) %>% filter(w > 0) %>%
    group_by(site, unit, year, top_cm, bottom_cm, w) %>%
    summarise(soc_gkg = mean(soc_gkg), bd_gcm3 = mean(bd_gcm3), .groups = "drop") %>%   # replicate cores
    group_by(site, unit, year) %>%
    summarise(cover = sum(w) / D, soc = sum(w * soc_gkg) / sum(w),
              stock = if (all(!is.na(bd_gcm3))) sum(w * soc_gkg * bd_gcm3) / 10 else NA_real_,   # Mg C/ha
              .groups = "drop") %>% filter(cover >= MIN_COVER)
}
agcros_units <- function(D) {
  if (is.null(soil)) return(NULL)
  layer_mean(soil, D) %>% group_by(site, year) %>%
    mutate(soc_dev = soc - mean(soc), stock_dev = stock - mean(stock)) %>%
    group_by(site, unit) %>%
    summarise(SOC_gkg = mean(soc_dev), STOCK = mean(stock_dev), .groups = "drop") %>%
    mutate(unit = as.character(unit), level = UNIT_LEVEL, source = "AgCROS")
}
# mediator table for one variable: site, unit, level, value
med_table <- function(var, source = c("napeshm", "agcros"), D = DEPTH_MAIN) {
  out <- list()
  if ("napeshm" %in% source && !is.null(napeshm) && var %in% names(napeshm))
    out$n <- napeshm %>% transmute(site, unit = as.character(unit), level = "treatment", value = .data[[var]], source = "NAPESHM")
  if ("agcros" %in% source) { a <- agcros_units(D)
    if (!is.null(a) && var %in% names(a)) out$a <- a %>% transmute(site, unit, level, value = .data[[var]], source) }
  bind_rows(out)
}

prep <- function(u, drop = NULL) {
  corn_u <- bind_rows(corn %>% mutate(unit = as.character(MGT_field), level = "plot"),
                      corn %>% mutate(unit = as.character(MGT_trt_base), level = "treatment")) %>%
    filter(!site %in% drop)
  key <- function(x) sub("^0+(?=.)", "", trimws(as.character(x)), perl = TRUE)
  d <- corn_u %>% mutate(unit = key(unit)) %>%
    inner_join(u %>% mutate(unit = key(unit)), by = c("site", "level", "unit")) %>% filter(!is.na(value))
  if (!is.na(YEAR_WINDOW) && !is.null(soil)) {
    sy <- soil %>% distinct(site, year) %>% group_by(site) %>% summarise(sy = list(year), .groups = "drop")
    d <- d %>% left_join(sy, by = "site") %>%
      filter(source != "AgCROS" | mapply(function(y, s) is.null(s) || min(abs(y - s)) <= YEAR_WINDOW, year, sy)) %>% select(-sy)
  }
  ok <- d %>% group_by(site) %>% summarise(n_u = n_distinct(unit), n_s = n_distinct(round(SHMI)), .groups = "drop") %>%
    filter(n_u >= MIN_UNITS, n_s >= 2)
  missing <- setdiff(unique(u$site), ok$site)
  if (length(missing)) message("  not used (too few matched units): ", paste(missing, collapse = ", "))
  d %>% filter(site %in% ok$site) %>% droplevels() %>% group_by(site) %>%
    mutate(SHMI_w = SHMI - mean(SHMI[!duplicated(MGT_combo)]),
           M = value - mean(value[!duplicated(unit)])) %>% ungroup() %>%
    arrange(plot_id, year)
}

# ============================================================================
# 2. models
# ============================================================================
yield_fit <- function(d, extra = NULL) {
  nv <- any(tapply(d$N_c, d$site_year, function(x) diff(range(x))) > 1e-8)
  rhs <- c("site_year", n_terms(d), o_term(d), "SHMI_w", mei_terms(d), if (nv) "SHMI_w:N_c",
           if (n_distinct(d$age_c) > 1) "SHMI_w:age_c", extra)
  lme(reformulate(rhs, "yield_kgha"), random = ~ 1 | trt/plot_id,
      correlation = corCAR1(form = ~ year | trt/plot_id), data = d, method = "REML", control = ctrl)
}
coef_se <- function(m, k) { tt <- summary(m)$tTable; c(b = tt[k, "Value"], se = tt[k, "Std.Error"], p = tt[k, "p-value"]) }

mediate <- function(u, label, drop = NULL) {
  if (is.null(u) || !nrow(u)) { message("No data for ", label); return(NULL) }
  message(label)
  d <- prep(u, drop)
  if (!nrow(d)) { message("No site has enough units with ", label); return(NULL) }
  units <- d %>% group_by(site, unit) %>%
    summarise(SHMI = mean(SHMI), SHMI_w = mean(SHMI_w), M = mean(M), N = 100 * mean(N100), .groups = "drop") %>%
    group_by(site) %>% mutate(N_w = N - mean(N)) %>% ungroup()
  # a path: SOC on SHMI within sites (+ N where it varies)
  fa <- reformulate(c("SHMI_w", if (sd(units$N_w) > 1e-8) "N_w", if (n_distinct(units$site) > 1) "site"), "M")
  ma <- lm(fa, data = units)
  a <- summary(ma)$coefficients["SHMI_w", ]
  # c and c' (+ b) on the same plot-years
  mc  <- yield_fit(d)
  mcp <- yield_fit(d, extra = "M")
  c_  <- coef_se(mc, "SHMI_w"); cp <- coef_se(mcp, "SHMI_w"); b_ <- coef_se(mcp, "M")
  # Monte Carlo indirect effect (per 10 SHMI)
  ab <- 10 * rnorm(NMC, a[["Estimate"]], a[["Std. Error"]]) * rnorm(NMC, b_[["b"]], b_[["se"]])
  cpd <- 10 * rnorm(NMC, cp[["b"]], cp[["se"]])
  share <- ab / (ab + cpd)
  res <- tibble(
    mediator = label, sites = paste(sort(unique(d$site)), collapse = ", "),
    n_sites = n_distinct(d$site), n_units = nrow(units), plot_years = nrow(d),
    r_SHMI_M = cor(units$SHMI_w, units$M),
    a_per10 = 10 * a[["Estimate"]], a_lo = 10 * (a[["Estimate"]] - qt(.975, ma$df.residual) * a[["Std. Error"]]),
    a_hi = 10 * (a[["Estimate"]] + qt(.975, ma$df.residual) * a[["Std. Error"]]), a_p = a[["Pr(>|t|)"]],
    b_per_unit = b_[["b"]], b_se = b_[["se"]], b_p = b_[["p"]],
    c_total_per10 = 10 * c_[["b"]], c_se = 10 * c_[["se"]],
    c_direct_per10 = 10 * cp[["b"]], c_direct_se = 10 * cp[["se"]],
    indirect_per10 = 10 * a[["Estimate"]] * b_[["b"]],
    indirect_lo = unname(quantile(ab, .025)), indirect_hi = unname(quantile(ab, .975)),
    share_mediated = (10 * a[["Estimate"]] * b_[["b"]]) / (10 * a[["Estimate"]] * b_[["b"]] + 10 * cp[["b"]]),
    share_lo = unname(quantile(share, .025)), share_hi = unname(quantile(share, .975)))
  attr(res, "units") <- units %>% mutate(mediator = label)
  res
}

MEDIATORS <- list(
  list(var = "lv_SH",      src = "napeshm", label = "Soil health, lv_SH (0-100; NAPESHM)"),
  list(var = "SHAPE_mean", src = "napeshm", label = "SHAPE mean score (0-1; NAPESHM)"),
  list(var = "SOC_gkg",    src = c("napeshm", "agcros"), D = DEPTH_MAIN, label = sprintf("SOC 0-%d cm (g/kg; NAPESHM + AgCROS)", DEPTH_MAIN)),
  list(var = "POXC_mgkg",  src = "napeshm", label = "Active carbon, POXC (mg/kg; NAPESHM)"),
  list(var = "WAS",        src = "napeshm", label = "Wet aggregate stability (NAPESHM)"),
  list(var = "SOC_gkg",    src = "agcros", D = DEPTH_SENS, label = sprintf("SOC 0-%d cm (g/kg; AgCROS)", DEPTH_SENS)),
  list(var = "STOCK",      src = "agcros", D = DEPTH_SENS, label = sprintf("SOC stock 0-%d cm (Mg C/ha; AgCROS)", DEPTH_SENS)))
tabs <- lapply(MEDIATORS, function(m) med_table(m$var, m$src, if (is.null(m$D)) DEPTH_MAIN else m$D))
runs <- Map(function(m, u) tryCatch(mediate(u, m$label), error = function(e) { message("  failed: ", conditionMessage(e)); NULL }),
            MEDIATORS, tabs)
keep <- !vapply(runs, is.null, logical(1))
runs <- runs[keep]; tabs <- tabs[keep]
stopifnot(length(runs) > 0)
res <- bind_rows(runs)

# leave one site out, for every mediator with at least 3 sites
for (k in seq_along(runs)) {
  if (runs[[k]]$n_sites < 3) next
  for (s in strsplit(runs[[k]]$sites, ", ")[[1]])
    res <- bind_rows(res, tryCatch(mediate(tabs[[k]], paste0(runs[[k]]$mediator, ", without ", s), drop = s),
                                   error = function(e) NULL))
}
print(as.data.frame(res %>% mutate(across(where(is.numeric), ~ signif(.x, 3)))), row.names = FALSE)
write.csv(res, file.path(OUT, "mediation_results.csv"), row.names = FALSE)
units_all <- bind_rows(lapply(runs, attr, "units"))
write.csv(units_all, file.path(OUT, "mediation_units.csv"), row.names = FALSE)

# ============================================================================
# 3. figure: a path (units) and the decomposition (primary mediator)
# ============================================================================
u1 <- attr(runs[[1]], "units"); r1 <- runs[[1]]
pa <- ggplot(u1, aes(SHMI_w, M, colour = site)) + geom_point(size = 2) +
  geom_smooth(aes(group = 1), method = "lm", formula = y ~ x, colour = "black", linewidth = 0.7) +
  labs(x = "SHMI minus site mean", y = paste(r1$mediator, "minus site mean"),
       title = sprintf("a  Mediator change per 10 SHMI: %.2f (P = %.2g)", r1$a_per10, r1$a_p)) +
  theme_bw(base_size = 10)
dec <- tibble(part = factor(c("Total (c)", "Direct (c')", "Indirect via soil (a x b)"),
                            c("Total (c)", "Direct (c')", "Indirect via soil (a x b)")),
              est = c(r1$c_total_per10, r1$c_direct_per10, r1$indirect_per10),
              lo = c(r1$c_total_per10 - 1.96 * r1$c_se, r1$c_direct_per10 - 1.96 * r1$c_direct_se, r1$indirect_lo),
              hi = c(r1$c_total_per10 + 1.96 * r1$c_se, r1$c_direct_per10 + 1.96 * r1$c_direct_se, r1$indirect_hi))
pb <- ggplot(dec, aes(part, est)) + geom_hline(yintercept = 0, linetype = 2, colour = "grey60") +
  geom_pointrange(aes(ymin = lo, ymax = hi), colour = "#08519c") +
  labs(x = NULL, y = "Corn yield per 10 SHMI (kg/ha)",
       title = sprintf("b  Share carried by the mediator: %.0f%% (%.0f to %.0f%%)", 100 * r1$share_mediated,
                       100 * r1$share_lo, 100 * r1$share_hi)) +
  theme_bw(base_size = 10)
if (requireNamespace("patchwork", quietly = TRUE)) {
  fig <- patchwork::wrap_plots(pa, pb, widths = c(1.2, 1)) +
    patchwork::plot_annotation(caption = sprintf("%s; %d units at %s; %d plot-years. Slopes at %d kg N/ha and average age.",
                                                 r1$mediator, r1$n_units, r1$sites, r1$plot_years, N_REF))
  ggsave(file.path(OUT, "Fig_mediation.png"), fig, width = 9.5, height = 4, dpi = 300)
}
cat("\nSaved mediation_results.csv, mediation_units.csv and Fig_mediation.png in", OUT, "\n")
