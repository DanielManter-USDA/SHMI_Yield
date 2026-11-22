################################################################################
# C2_export_tables.R
#
# Every computed table in the paper, from the outputs of the steps that compute it, in one
# workbook (one sheet per table) plus a CSV per table (17 Nov 2026 layout):
#   Table 2  SHMI against its parts       C1j Parts/SHMI_vs_parts.csv + B2d yield_parts/yield_parts.csv
#   Table 3  practice effects, strict pairs  C1i Pair_meta/pooled.csv + C1f Practice_additivity/data_budget.csv
#   Table 1 (file ED_Table1.csv) scores by system type          data/SHAPE_by_unit_with_SHMI.csv
#   (C1's example-systems Table1.csv left the paper on 17 Nov 2026 and is not exported.)
#   Supplementary Data Table 1   the ten experiments            B2 site_SHMI_ranges.csv
#   Supplementary Data Table 2   coarse-record scenarios        B9 Coarse/coarse_*.csv (definitions and results)
# (Table 4, "Scoring a field with SHMI", is text and has no computed values.)
# Retired 17 Nov 2026: SuppTable3 (B5 mediation) and SuppTable4 (A4a benchmark; now Table 2).
# A table whose source is missing (its step was not run) is skipped with a message.
#
# Output  Outputs/Manuscript_work/Tables.xlsx and Outputs/Manuscript_work/<table>.csv
################################################################################
suppressPackageStartupMessages({ library(dplyr); library(tidyr) })
OUT  <- here::here("Outputs", "Manuscript_work")
MAIN <- here::here("Outputs", "yield_main_default_total_agcros")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
TABLES <- list()
add <- function(name, src, build) {
  if (!file.exists(src)) { message("  ", name, ": skipped (", src, " not found)"); return(invisible()) }
  t <- build(src)
  TABLES[[name]] <<- t
  write.csv(t, file.path(OUT, paste0(name, ".csv")), row.names = FALSE, fileEncoding = "UTF-8")
  message("  ", name, ": ", nrow(t), " rows")
}
rng <- function(x, digits = 0) {
  f <- function(v) formatC(v, format = "f", digits = digits)
  if (diff(range(x)) == 0) f(median(x)) else sprintf("%s (%s–%s)", f(median(x)), f(min(x)), f(max(x)))
}


# ---- shared formatting (17 Nov 2026) -----------------------------------------------------
MINUS <- "\u2212"
num <- function(x, d = 2, plus = FALSE) {
  v <- formatC(round(x, d) + 0, format = "f", digits = d, flag = if (plus) "+" else "")
  gsub("-", MINUS, v)
}
ci_to <- function(lo, hi, d) paste0(num(lo, d, TRUE), " to ", num(hi, d, TRUE))

# Table 2: SHMI against its parts (C1j for soil health; B2d for corn yield)
if (!exists("R2_DIGITS"))   R2_DIGITS   <- 3       # held-out R2, its difference from SHMI and the CI
if (!exists("CEILING_COL")) CEILING_COL <- FALSE   # TRUE adds "% of replicate ceiling" (needs C1k)
add("Table2", here::here("Outputs", "Parts", "SHMI_vs_parts.csv"), function(f) {
  p <- read.csv(f, stringsAsFactors = FALSE)
  LEV <- c("SHMI", "Cover/Roots", "Diversity", "Inverse Disturbance", "Amendments")
  LAB <- c("SHMI", "Cover/Roots alone", "Diversity alone", "Inverse Disturbance alone", "Amendments alone")
  soil <- function(o) {
    q <- p[p$outcome == o, ][match(LEV, p$score[p$outcome == o]), ]
    data.frame(check.names = FALSE,
      a = ifelse(q$score == "SHMI", "0", num(q$dAICc, 0, TRUE)),
      # held-out R2 and its difference from SHMI are both computed from the unrounded values in
      # C1j and shown to 3 decimals, so the printed difference matches the printed R2 values
      b = ifelse(q$score == "SHMI", num(q$R2_site_out, R2_DIGITS),
                 sprintf("%s (%s, %s)", num(q$R2_site_out, R2_DIGITS), num(q$d_R2, R2_DIGITS, TRUE),
                         ci_to(q$d_lo, q$d_hi, R2_DIGITS))),
      r2 = q$R2_site_out)
  }
  t <- data.frame(Score = LAB, check.names = FALSE)
  s1 <- soil("lv_SH"); s2 <- soil("SHAPE")
  t[["lv_SH: ΔAICc"]] <- s1$a
  t[["lv_SH: R² at held-out sites (difference from SHMI, 95% CI)"]] <- s1$b
  # optional: held-out R2 as a share of the replicate ceiling (C1k; lv_SH only)
  cf <- here::here("Outputs", "Ceiling", "ceiling.csv")
  if (isTRUE(CEILING_COL) && file.exists(cf)) {
    cl <- read.csv(cf, stringsAsFactors = FALSE)
    ceil <- cl$estimate[cl$quantity == "ceiling_means"]
    t[["lv_SH: % of replicate ceiling"]] <- num(100 * s1$r2 / ceil, 0)
  }
  t[["SHAPE: ΔAICc"]] <- s2$a
  t[["SHAPE: R² at held-out sites (difference from SHMI, 95% CI)"]] <- s2$b
  yf <- here::here("Outputs", "yield_parts", "yield_parts.csv")
  if (file.exists(yf)) {
    yp <- read.csv(yf, stringsAsFactors = FALSE)
    YK <- c("SHMI", "Cover/Roots alone", "Diversity alone", "Inverse Disturbance alone", "Share of years amended alone")
    y <- yp[match(YK, yp$score), ]
    if (!"pct_of_SHMI" %in% names(y)) y$pct_of_SHMI <- 100 * y$slope_per_sd / y$slope_per_sd[1]
    t[["Corn yield: ΔAIC"]] <- ifelse(YK == "SHMI", "0", num(y$dAIC_vs_SHMI, 0, TRUE))
    t[["Corn yield: gain, kg ha-1 per within-site SD (95% CI)"]] <-
      sprintf("%s (%s to %s)", num(y$slope_per_sd, 0), num(y$lower, 0), num(y$upper, 0))
    t[["Corn yield: % of SHMI"]] <- num(y$pct_of_SHMI, 0)
    if ("experiments" %in% names(y)) t[["Corn yield: experiments in which the score varied"]] <- y$experiments
    if (all(c("dAIC_loo_min", "dAIC_loo_max") %in% names(y)))
      t[["Corn yield: ΔAIC, one experiment left out (range)"]] <-
        ifelse(YK == "SHMI", "", paste0(num(y$dAIC_loo_min, 0, TRUE), " to ", num(y$dAIC_loo_max, 0, TRUE)))
  } else message("  Table2: corn yield columns skipped (", yf, " not found; run B2d)")
  t
})

# Table 3: practice effects from strict matched pairs (C1i pooled meta-analyses; C1f data budget)
add("Table3", here::here("Outputs", "Pair_meta", "pooled.csv"), function(f) {
  pl <- read.csv(f, stringsAsFactors = FALSE)
  PR <- c("Cover crop vs none" = "Cover crop", "Organic vs synthetic inputs" = "Organic inputs",
          "Reduced vs conventional tillage" = "Reduced tillage", "Rotation vs monoculture" = "Rotation")
  soc <- pl[grepl("^SOC", pl$outcome), ]; sh <- pl[grepl("^lv_SH", pl$outcome), ]
  soc <- soc[match(names(PR), soc$practice), ]; sh <- sh[match(names(PR), sh$practice), ]
  pct <- function(x) 100 * (exp(x) - 1)
  t <- data.frame(check.names = FALSE, Practice = unname(PR))
  bf <- here::here("Outputs", "Practice_additivity", "data_budget.csv")
  if (file.exists(bf)) {
    bd <- read.csv(bf, stringsAsFactors = FALSE)
    bd <- bd[bd$pairing == "strict", ][match(unname(PR), bd$practice[bd$pairing == "strict"]), ]
    t[["Systems usable in strict pairs (% of all)"]] <- sprintf("%d (%.0f%%)", bd$systems, bd$pct_systems)
  }
  t[["Pairs (sites)"]] <- sprintf("%d (%d)", sh$pairs, sh$sites)
  t[["SOC change, % (95% CI)"]] <- sprintf("%s (%s)", num(pct(soc$estimate), 1, TRUE), ci_to(pct(soc$ci_lb), pct(soc$ci_ub), 1))
  t[["lv_SH change, points (95% CI)"]] <- sprintf("%s (%s)", num(sh$estimate, 1, TRUE), ci_to(sh$ci_lb, sh$ci_ub, 1))
  t[["I², lv_SH (%)"]] <- sprintf("%.0f", sh$I2)
  t[["lv_SH effect expected in a new pair (95% PI)"]] <- ci_to(sh$pi_lb, sh$pi_ub, 1)
  t
})

# Table 1 (written as ED_Table1.csv): median (lowest-highest) by system type
add("ED_Table1", here::here("data", "SHAPE_by_unit_with_SHMI.csv"), function(f) {
  d <- read.csv(f)
  lab <- c(annual = "Annual", perennial = "Perennial", reference = "Reference", fallow = "Fallow")
  vars <- c(Cover = "Cover/Roots", Diversity = "Diversity (with livestock)", InvDist = "Inverse Disturbance",
            SHMI_base = "SHMI_base", Amend_bonus = "Amendment bonus", SHMI = "SHMI")
  out <- lapply(names(lab), function(s) {
    g <- d[d$system == s, ]
    setNames(data.frame(vapply(names(vars), function(v) rng(g[[v]]), "")), sprintf("%s (%d)", lab[[s]], nrow(g)))
  })
  cbind(Score = unname(vars), do.call(cbind, out))
})

# Supplementary Data Table 1: the ten experiments (B2, section 7b)
add("SuppTable1", file.path(MAIN, "site_SHMI_ranges.csv"), function(f) read.csv(f, check.names = FALSE))

# Supplementary Data Table 2: every coarse-record scenario of B9 with its results (the four
# record types are also drawn in Extended Data Fig. 5; the components removed one at a time
# appear only here). Soil health: NAPESHM, out-of-fold, published weights, median (95% range);
# yield: within-site slope per 10 SHMI (95% CI); agreement with the full-record score.
add("SuppTable2", here::here("Outputs", "Coarse", "coarse_scenarios.csv"), function(f) {
  CO <- function(x) here::here("Outputs", "Coarse", x)
  sc <- read.csv(f, stringsAsFactors = FALSE) %>% filter(set %in% c("a", "b"))
  mr <- function(m, lo, hi, d = 2) ifelse(is.na(m), "", ifelse(abs(hi - lo) < 1e-9, formatC(m, format = "f", digits = d),
                  sprintf("%s (%s–%s)", formatC(m, format = "f", digits = d), formatC(lo, format = "f", digits = d), formatC(hi, format = "f", digits = d))))
  t <- sc %>% transmute(scenario, `Set` = ifelse(set == "b", "Record type", "One detail removed"),
                        `Detail` = ifelse(set == "b", "", component), `Record` = level)
  if (file.exists(CO("coarse_agreement.csv"))) {
    ag <- read.csv(CO("coarse_agreement.csv"), stringsAsFactors = FALSE)
    t <- t %>% left_join(ag %>% transmute(scenario, `RMS difference from full score (SHMI points)` = round(rms_diff, 1),
                                          `Mean bias (SHMI points)` = round(bias, 1),
                                          `Within-site rank agreement (Spearman)` = round(within_rho, 2)), by = "scenario")
  }
  if (file.exists(CO("coarse_napeshm.csv"))) {
    ns <- read.csv(CO("coarse_napeshm.csv"), stringsAsFactors = FALSE) %>% filter(weights == "published")
    t <- t %>% left_join(ns %>% transmute(scenario,
           `Within-site r` = mr(within_r_med, within_r_lo, within_r_hi),
           `Sites ranked in the right direction (%)` = mr(100 * sites_correct_med, 100 * sites_correct_lo, 100 * sites_correct_hi, 0),
           `Between-site r` = mr(between_r_med, between_r_lo, between_r_hi),
           `Overall R2` = mr(R2_med, R2_lo, R2_hi)), by = "scenario")
  }
  if (file.exists(CO("coarse_yield.csv"))) {
    y <- read.csv(CO("coarse_yield.csv"), stringsAsFactors = FALSE)
    full <- y$slope[y$scenario == "D0"][1]
    t <- t %>% left_join(y %>% transmute(scenario,
           `Corn yield per 10 SHMI, kg ha-1 (95% CI)` = sprintf("%.0f (%.0f to %.0f)", slope, lower, upper),
           `Yield slope, % of full record` = if (length(full) && !is.na(full)) round(100 * slope / full) else NA), by = "scenario")
  }
  t %>% arrange(match(Set, c("Record type", "One detail removed"))) %>% rename(Scenario = scenario)
})

if (length(TABLES)) {
  if (requireNamespace("writexl", quietly = TRUE)) {
    writexl::write_xlsx(TABLES, file.path(OUT, "Tables.xlsx"))
    message("Saved Tables.xlsx (", length(TABLES), " sheets) and one CSV per table in ", OUT)
  } else message("writexl not installed: CSVs only, in ", OUT)
}
