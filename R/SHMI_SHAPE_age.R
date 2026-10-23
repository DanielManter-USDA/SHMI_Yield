################################################################################
# SHMI_SHAPE_age.R
#
# 1. Does SHMI track SHAPE? Correlations of SHMI (and SHMI_base) with the SHAPE
#    scores (mean of SOC, AC, Resp, KSSL_WAS, and each one), for all units, site
#    means and within sites (each unit minus its site mean), as in Fig 5. SH (the
#    measured soil-health score SHMI was calibrated on) is shown alongside for
#    reference. SHAPE scores are put on 0-100 (x 100) so slopes read as SHAPE
#    points per SHMI point.
#
# 2. Does treatment age change how well SHMI tracks SH and SHAPE?
#    Age = sampling year - treatment start year. Start years from treatments.csv,
#    missing ones filled from sites_AI_startdate.csv; any still missing take the
#    site median (as in SHMI_age_disturbance.R) and are flagged.
#    For y = SH and SHAPE_mean (0-100), mixed model with SHMI split into its site
#    mean (between) and the unit's deviation from it (within):
#      y ~ SHMI_between + SHMI_within * age_c + (1 | site)
#    SHMI_within:age_c = how much the within-site slope changes per year of age.
#    Sensitivity: random within-site slope per site, (1 + SHMI_within | site).
#    Also: each site's own within-site slope against the site's age (Fig B).
#
# Inputs
#   Outputs/SHAPE/SHAPE_by_unit_with_SHMI.csv   (Merge_SHAPE_SHMI.R)
#   Inputs/CSV Files/treatments.csv             (treatment_start_year)
#   Inputs/CSV Files/sites_AI_startdate.csv     (updated start dates; site- or
#                                                treatment-level, see read_starts)
# Outputs: Outputs/SHAPE_age/
#   shape_correlations.csv, age_models.csv, age_slopes.csv, site_slopes.csv,
#   unit_ages.csv, FigA_SHMI_vs_SHAPE.png, FigB_slope_by_age.png
################################################################################
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(purrr); library(ggplot2); library(here); library(lme4)
})
if (!requireNamespace("patchwork", quietly = TRUE)) stop("Install patchwork: install.packages('patchwork')")
HAS_LMERTEST <- requireNamespace("lmerTest", quietly = TRUE)

in_csv <- function(f) {                      # Inputs/ or inputs/
  p <- here::here("Inputs", "CSV Files", f)
  if (!file.exists(p)) p <- here::here("inputs", "CSV Files", f)
  p
}
CFG <- list(
  shape_file  = here::here("Outputs", "SHAPE", "SHAPE_by_unit_with_SHMI.csv"),
  treatments  = in_csv("treatments.csv"),
  start_file  = in_csv("sites_AI_startdate.csv"),
  prefer_new  = TRUE,       # TRUE: the new file's start year replaces treatments.csv where both exist
  sample_year = 2019,       # NAPESHM sampling year, as in SHMI_age_disturbance.R
  min_units   = 3,          # Fig B: sites with at least this many units ...
  min_range   = 5,          # ... and an SHMI range of at least this many points
  ages        = NULL,       # ages at which slopes are reported; NULL = every 5 years over the data
  age_max_sens = 60,        # sensitivity: refit without sites older than this
  out_dir     = here::here("Outputs", "SHAPE_age")
)
if (exists("CFG_OVERRIDE")) CFG <- utils::modifyList(CFG, CFG_OVERRIDE)
dir.create(CFG$out_dir, recursive = TRUE, showWarnings = FALSE)
OUTF <- function(x) file.path(CFG$out_dir, x)

# ==============================================================================
# Data
# ==============================================================================
d0 <- read.csv(CFG$shape_file, stringsAsFactors = FALSE)
SHAPE_IND <- grep("^SHAPE_(SOC|AC|ACE|Resp|KSSL_WAS)$", names(d0), value = TRUE)
d0 <- d0 %>% mutate(across(all_of(c("SHAPE_mean", SHAPE_IND)), ~ 100 * .x))     # 0-1 -> 0-100

year_of <- function(x) {                     # a year or a date in any common format -> year
  x <- as.character(x)
  y <- suppressWarnings(as.numeric(regmatches(x, regexpr("(18|19|20)[0-9]{2}", x))))
  out <- rep(NA_real_, length(x)); out[grepl("(18|19|20)[0-9]{2}", x)] <- y
  ifelse(out > 1900 & out <= CFG$sample_year, out, NA_real_)
}

trt <- read.csv(CFG$treatments, stringsAsFactors = FALSE) %>%
  transmute(MGT_combo = paste0("NAPESHM_", sitecode, "_", treatmentid), sitecode,
            start_trt = year_of(treatment_start_year))

# sites_AI_startdate.csv: keyed by MGT_combo, by sitecode + treatmentid, or by
# sitecode alone (one start for the site). Column names are lowercased and made
# unique (Start_Year + start_year -> start_year, start_year_1); metadata columns
# (confidence, note, source) are skipped and the most complete start column is used.
read_starts <- function(f) {
  if (!file.exists(f)) { warning("No start-date file at ", f, "; using treatments.csv only"); return(NULL) }
  a <- read.csv(f, stringsAsFactors = FALSE, check.names = FALSE)
  orig <- names(a)
  names(a) <- make.unique(tolower(trimws(names(a))), sep = "_")
  cand <- grep("start", names(a), value = TRUE)
  cand <- cand[!grepl("confidence|note|source", cand)]
  if (!length(cand)) stop("No start-year column in ", f, ": ", paste(orig, collapse = ", "))
  yrs <- sapply(cand, function(cn) sum(!is.na(year_of(a[[cn]]))))
  message("Start-date file, valid years per column: ",
          paste(sprintf("%s = %d", orig[match(cand, names(a))], yrs), collapse = "; "))
  sc <- cand[which.max(yrs)]
  message("Using '", orig[match(sc, names(a))], "'")
  a$start_new <- year_of(a[[sc]])
  if ("mgt_combo" %in% names(a)) {
    a %>% transmute(MGT_combo = mgt_combo, start_new) %>% mutate(level = "unit")
  } else if (all(c("sitecode", "treatmentid") %in% names(a))) {
    a %>% transmute(MGT_combo = paste0("NAPESHM_", sitecode, "_", treatmentid), start_new) %>% mutate(level = "unit")
  } else if ("sitecode" %in% names(a)) {
    a %>% group_by(sitecode) %>% summarise(start_new = first(na.omit(start_new)), .groups = "drop") %>% mutate(level = "site")
  } else stop("Start-date file needs MGT_combo, sitecode + treatmentid, or sitecode: ", paste(orig, collapse = ", "))
}
new <- read_starts(CFG$start_file)
ages <- trt
if (!is.null(new)) {
  ages <- if (new$level[1] == "unit") left_join(ages, select(new, MGT_combo, start_new), by = "MGT_combo")
          else left_join(ages, select(new, sitecode, start_new), by = "sitecode")
} else ages$start_new <- NA_real_
dis <- ages %>% filter(!is.na(start_trt), !is.na(start_new), start_trt != start_new)
if (nrow(dis)) {
  message(nrow(dis), " units have different start years in the two files (",
          if (CFG$prefer_new) "using the new file" else "keeping treatments.csv; set prefer_new = TRUE to switch", ")")
  write.csv(dis, OUTF("start_year_disagreements.csv"), row.names = FALSE)
}
ages <- ages %>%
  mutate(start = if (CFG$prefer_new) coalesce(start_new, start_trt) else coalesce(start_trt, start_new),
         source = case_when(!is.na(start_trt) & (!CFG$prefer_new | is.na(start_new)) ~ "treatments.csv",
                            !is.na(start_new) ~ "sites_AI_startdate.csv", TRUE ~ NA_character_))

d <- d0 %>% left_join(ages %>% select(MGT_combo, start, source), by = "MGT_combo") %>%
  filter(!is.na(SHMI), !is.na(SH), !is.na(SHAPE_mean)) %>%
  mutate(age = CFG$sample_year - start) %>%
  group_by(sitecode) %>%
  mutate(source = ifelse(is.na(age) & any(!is.na(age)), "site median (imputed)", source),
         age = ifelse(is.na(age), median(age, na.rm = TRUE), age)) %>%
  ungroup()
if (any(is.na(d$age))) {
  warning(sum(is.na(d$age)), " units at ", n_distinct(d$sitecode[is.na(d$age)]),
          " sites still have no start year and are dropped from part 2: ",
          paste(unique(d$sitecode[is.na(d$age)]), collapse = ", "))
}
write.csv(d %>% select(MGT_combo, sitecode, system, start, age, source), OUTF("unit_ages.csv"), row.names = FALSE)
cat("Units:", nrow(d), "| sites:", n_distinct(d$sitecode), "\nAge sources:\n"); print(table(d$source, useNA = "ifany"))
cat("Age (years): median", median(d$age, na.rm = TRUE), "| range", paste(range(d$age, na.rm = TRUE), collapse = "-"), "\n")

# ==============================================================================
# 1. SHMI vs SHAPE (and SH): all units, site means, within sites
# ==============================================================================
# Part 1 uses every unit with SHMI, SH and SHAPE, whether or not it has a start year;
# only part 2 (da below) drops units without an age
stopifnot(nrow(d) == sum(!is.na(d0$SHMI) & !is.na(d0$SH) & !is.na(d0$SHAPE_mean)))
cat("Part 1 uses all", nrow(d), "units (", sum(is.na(d$age)), "without a start year are kept here, dropped in part 2)\n")
YS <- c("SHAPE_mean", SHAPE_IND, "SH")
XS <- c("SHMI", "SHMI_base")
lev <- function(df, level) switch(level,
  all   = df,
  sites = df %>% group_by(sitecode) %>% summarise(across(all_of(c(XS, YS)), mean), .groups = "drop"),
  within = df %>% group_by(sitecode) %>% filter(n() > 1) %>% mutate(across(all_of(c(XS, YS)), ~ .x - mean(.x))) %>% ungroup())
cor_row <- function(df, x, y) {
  ok <- complete.cases(df[c(x, y)]); a <- df[[x]][ok]; b <- df[[y]][ok]
  ct <- cor.test(a, b)
  tibble(n = sum(ok), r = unname(ct$estimate), lo = ct$conf.int[1], hi = ct$conf.int[2], p = ct$p.value,
         rho = suppressWarnings(cor(a, b, method = "spearman")), slope = unname(coef(lm(b ~ a))[2]))
}
cors <- expand_grid(level = c("all", "sites", "within"), x = XS, y = YS) %>%
  mutate(res = pmap(list(level, x, y), function(level, x, y) cor_row(lev(d, level), x, y))) %>%
  unnest(res)
write.csv(cors %>% mutate(across(where(is.numeric), ~ signif(.x, 3))), OUTF("shape_correlations.csv"), row.names = FALSE)
cat("\n=== 1. Correlations with SHMI (r; slope = points per SHMI point) ===\n")
print(cors %>% filter(x == "SHMI") %>% select(level, y, n, r, lo, hi, slope) %>%
        mutate(across(where(is.double), ~ round(.x, 2))) %>%
        pivot_wider(id_cols = y, names_from = level, values_from = c(r, slope)), width = Inf)
cat("Is SHMI closer to SH or to SHAPE? (all units: SHAPE_mean vs SH, r =",
    round(cor(d$SHAPE_mean, d$SH), 2), "; within sites r =",
    round(with(lev(d, "within"), cor(SHAPE_mean, SH)), 2), ")\n")

# Fig A: rows SHAPE_mean and SH, columns all / site means / within (style of Fig 5)
panel <- function(level, y, ylab) {
  df <- lev(d, level); cr <- cors %>% filter(level == !!level, x == "SHMI", y == !!y)
  xl <- switch(level, all = "SHMI", sites = "Site mean SHMI", within = "SHMI minus site mean")
  yl <- switch(level, all = ylab, sites = paste("Site mean", ylab), within = paste(ylab, "minus site mean"))
  ttl <- switch(level, all = sprintf("all %d units", nrow(df)), sites = sprintf("%d site means", nrow(df)), within = "within sites")
  ggplot(df, aes(SHMI, .data[[y]])) +
    geom_point(colour = "steelblue", alpha = 0.6, size = 1) +
    geom_smooth(method = "lm", formula = y ~ x, se = FALSE, colour = "#C8414B", linewidth = 1) +
    annotate("text", x = -Inf, y = Inf, hjust = -0.05, vjust = 1.2, size = 3.3,
             label = sprintf("r = %.2f\nR² = %.2f\nslope = %.2f", cr$r, cr$r^2, cr$slope)) +
    labs(x = xl, y = yl, title = sprintf("%s vs SHMI: %s", ylab, ttl)) +
    theme_classic(base_size = 10) + theme(plot.title = element_text(size = 10))
}
figA <- patchwork::wrap_plots(c(map(c("all", "sites", "within"), panel, y = "SHAPE_mean", ylab = "SHAPE"),
                                map(c("all", "sites", "within"), panel, y = "SH", ylab = "Soil health (SH)")), ncol = 3) +
  patchwork::plot_annotation(caption = paste0("SHAPE: mean of the SOC, active C, respiration and aggregate stability scores (0-100). ",
                                              "SH: measured soil-health score SHMI was calibrated on.\nWithin sites: each unit minus its site mean."))
ggsave(OUTF("FigA_SHMI_vs_SHAPE.png"), figA, width = 13, height = 8, dpi = 300)

# ==============================================================================
# 2. Age: does the within-site slope change with treatment age?
# ==============================================================================
da <- d %>% filter(!is.na(age)) %>%
  group_by(sitecode) %>% mutate(SHMI_between = mean(SHMI), SHMI_within = SHMI - SHMI_between) %>% ungroup() %>%
  mutate(age_c = age - mean(age), SHMI_between_c = SHMI_between - mean(SHMI_between))
AGE_MEAN <- mean(da$age)
AGES <- if (is.null(CFG$ages)) seq(5, 5 * floor(max(da$age) / 5), by = 5) else CFG$ages
fmt_p <- function(p) ifelse(p < 0.001, "p < 0.001", sprintf("p = %.3f", p))
age_varies <- da %>% group_by(sitecode) %>% summarise(v = n_distinct(age) > 1) %>% pull(v)
cat(sprintf("\n=== 2. Age (mean %.1f years; age varies within %d of %d sites) ===\n", AGE_MEAN, sum(age_varies), length(age_varies)))

fit_lmm <- function(f, data, REML = TRUE) {
  if (HAS_LMERTEST) lmerTest::lmer(f, data = data, REML = REML) else lmer(f, data = data, REML = REML)
}
coef_tab <- function(m) {
  cf <- summary(m)$coefficients
  p <- if ("Pr(>|t|)" %in% colnames(cf)) cf[, "Pr(>|t|)"] else 2 * pnorm(-abs(cf[, "t value"]))
  tibble(term = rownames(cf), estimate = cf[, "Estimate"], se = cf[, "Std. Error"], p = p)
}
slopes_at <- function(m, y, model) {             # within-site slope (per 10 SHMI) at chosen ages
  b <- fixef(m); V <- as.matrix(vcov(m)); k <- c("SHMI_within", "SHMI_within:age_c")
  map_dfr(AGES, function(a) {
    g <- c(1, a - AGE_MEAN); est <- sum(g * b[k]); se <- sqrt(drop(t(g) %*% V[k, k] %*% g))
    tibble(y, model, age = a, slope10 = 10 * est, lo = 10 * (est - 1.96 * se), hi = 10 * (est + 1.96 * se))
  })
}
models <- list(); slopes <- list()
for (y in c("SH", "SHAPE_mean")) {
  f0 <- reformulate(c("SHMI_between_c", "SHMI_within", "(1 | sitecode)"), y)
  f1 <- reformulate(c("SHMI_between_c", "SHMI_within * age_c", "(1 | sitecode)"), y)
  f2 <- reformulate(c("SHMI_between_c", "SHMI_within * age_c", "(1 + SHMI_within | sitecode)"), y)
  m0 <- fit_lmm(f0, da); m1 <- fit_lmm(f1, da)
  m2 <- tryCatch(fit_lmm(f2, da), error = function(e) { message("  ", y, " random slope: ", conditionMessage(e)); NULL })
  ds <- filter(da, age <= CFG$age_max_sens)                 # sensitivity: without the oldest sites
  m1s <- fit_lmm(f1, ds)
  m2s <- tryCatch(fit_lmm(f2, ds), error = function(e) NULL)
  sens_lab <- sprintf("age <= %d", CFG$age_max_sens)
  # likelihood-ratio test of the two age terms (ML fits)
  l0 <- logLik(fit_lmm(f0, da, REML = FALSE)); l1 <- logLik(fit_lmm(f1, da, REML = FALSE))
  lrt_p <- pchisq(2 * (as.numeric(l1) - as.numeric(l0)), df = attr(l1, "df") - attr(l0, "df"), lower.tail = FALSE)
  models[[y]] <- bind_rows(
    coef_tab(m1) %>% mutate(y, model = "random intercept"),
    if (!is.null(m2)) coef_tab(m2) %>% mutate(y, model = "random slope by site"),
    coef_tab(m1s) %>% mutate(y, model = paste("random intercept,", sens_lab)),
    if (!is.null(m2s)) coef_tab(m2s) %>% mutate(y, model = paste("random slope by site,", sens_lab))) %>%
    mutate(LRT_age_p = lrt_p, singular_random_slope = if (!is.null(m2)) isSingular(m2) else NA)
  slopes[[y]] <- bind_rows(slopes_at(m1, y, "random intercept"),
                           if (!is.null(m2)) slopes_at(m2, y, "random slope by site"))
  cat(sprintf("%s: %d units at %d sites older than %d years left out in the sensitivity fits\n", y,
              nrow(da) - nrow(ds), n_distinct(da$sitecode) - n_distinct(ds$sitecode), CFG$age_max_sens))
  tt <- coef_tab(m1) %>% filter(term %in% c("SHMI_within", "SHMI_within:age_c", "age_c", "SHMI_between_c"))
  cat(sprintf("\n%s: within slope %.2f per 10 SHMI at %.0f years; slope change per 10 years of age %+.2f (p = %.3f); LRT for age terms p = %.3f\n",
              y, 10 * tt$estimate[tt$term == "SHMI_within"], AGE_MEAN,
              100 * tt$estimate[tt$term == "SHMI_within:age_c"], tt$p[tt$term == "SHMI_within:age_c"], lrt_p))
}
models <- bind_rows(models); slopes <- bind_rows(slopes)
write.csv(models %>% mutate(across(where(is.numeric), ~ signif(.x, 4))), OUTF("age_models.csv"), row.names = FALSE)
write.csv(slopes %>% mutate(across(where(is.numeric), ~ round(.x, 2))), OUTF("age_slopes.csv"), row.names = FALSE)
cat("\nWithin-site slope (points per 10 SHMI) by age:\n")
print(slopes %>% mutate(across(where(is.double), ~ round(.x, 2))), n = Inf)

# each site's own within-site slope against its age (site-level view of the same question)
site_sl <- da %>% group_by(sitecode) %>%
  filter(n() >= CFG$min_units, diff(range(SHMI)) >= CFG$min_range) %>%
  group_modify(~ map_dfr(c("SH", "SHAPE_mean"), function(y) {
    m <- lm(reformulate("SHMI", y), data = .x)
    tibble(y, slope10 = 10 * coef(m)[["SHMI"]], se10 = 10 * summary(m)$coefficients["SHMI", "Std. Error"],
           n = nrow(.x), age = median(.x$age))
  })) %>% ungroup() %>% filter(is.finite(se10), se10 > 0)
write.csv(site_sl %>% mutate(across(where(is.double), ~ round(.x, 3))), OUTF("site_slopes.csv"), row.names = FALSE)
meta <- site_sl %>% group_by(y) %>%
  group_modify(~ { m <- lm(slope10 ~ age, data = .x, weights = n); s <- summary(m)$coefficients
                   tibble(sites = nrow(.x), per10yr = 10 * s["age", 1], p = s["age", 4]) }) %>% ungroup()
# weighted by units per site: inverse-variance weights are dominated by sites whose 3-4 units
# happen to sit almost on a line (SE near 0), which made the trend look far stronger than it is
cat("\nSite slopes vs site age (weighted by units per site):\n"); print(meta)
cat("\nInteraction (change in within-site slope per 10 SHMI per 10 years), all models:\n")
print(models %>% filter(term == "SHMI_within:age_c") %>% transmute(y, model, per10yr = 100 * estimate, se = 100 * se, p) %>%
        mutate(across(where(is.double), ~ round(.x, 3))))

ylab_of <- c(SH = "Soil health (SH)", SHAPE_mean = "SHAPE")
figB <- ggplot(site_sl, aes(age, slope10)) +
  geom_hline(yintercept = 0, colour = "grey70") +
  geom_errorbar(aes(ymin = slope10 - se10, ymax = slope10 + se10), width = 0, colour = "grey60") +
  geom_point(aes(size = n), shape = 21, fill = "steelblue", colour = "black", alpha = 0.8) +
  geom_line(data = slopes %>% filter(model == "random slope by site"), aes(age, slope10), colour = "#C8414B", linewidth = 1) +
  geom_ribbon(data = slopes %>% filter(model == "random slope by site"), aes(age, ymin = lo, ymax = hi),
              inherit.aes = FALSE, fill = "#C8414B", alpha = 0.12) +
  geom_text(data = meta %>%
              left_join(models %>% filter(model == "random slope by site", term == "SHMI_within:age_c") %>%
                          select(y, est = estimate, p_rs = p), by = "y") %>%
              left_join(models %>% filter(model == "random intercept", term == "SHMI_within:age_c") %>%
                          select(y, est_ri = estimate, p_ri = p), by = "y"),
            aes(x = -Inf, y = Inf, label = sprintf(paste0("Change in slope per 10 years: %+.2f (%s; random intercept only %+.2f, %s)\n",
                                                          "Site slopes vs age: %+.2f per 10 years (%s, %d sites)"),
                                                   100 * est, fmt_p(p_rs), 100 * est_ri, fmt_p(p_ri), per10yr, fmt_p(p), sites)),
            hjust = -0.03, vjust = 1.3, size = 3.1, inherit.aes = FALSE) +
  facet_wrap(~ y, scales = "free_y", labeller = labeller(y = ylab_of)) +
  scale_size_continuous(range = c(1.5, 5), name = "Units at site") +
  labs(x = "Treatment age at sampling (years)", y = "Within-site slope (points per 10 SHMI)",
       caption = paste0("Points: each site's own slope (± SE; sites with ≥ ", CFG$min_units, " units and an SHMI range ≥ ",
                        CFG$min_range, "). Line and band: within-site slope from the mixed model\n",
                        "y ~ SHMI between + SHMI within × age + (1 + SHMI within | site), with 95% CI. Site-slope trend weighted by units. Age = ", CFG$sample_year,
                        " minus treatment start year.")) +
  theme_bw(base_size = 11) + theme(panel.grid.minor = element_blank())
ggsave(OUTF("FigB_slope_by_age.png"), figB, width = 11, height = 5, dpi = 300)
message("Saved results and figures to ", CFG$out_dir)
