################################################################################
# C2_export_tables.R
#
# Every computed table in the paper, from the outputs of the steps that compute it, in one
# workbook (one sheet per table) plus a CSV per table (17 Nov 2026 layout):
#   Table 1 (file ED_Table1.csv)  scores by system type     data/SHAPE_by_unit_with_SHMI.csv
#   Table 2  SHMI against its parts   C1j Parts/SHMI_vs_parts.csv + B2d yield_parts/yield_parts.csv
#   SuppTable1.csv  the ten experiments (Supplementary Table 2 in the paper)   B2 site_SHMI_ranges.csv
# The other Supplementary Tables are written by C1m (1, 4, 5, 6) and B8/B8b (3); C3 combines
# them into Outputs/Manuscript/numbered/Supplementary_Tables.xlsx.
# Retired 21 Nov 2026: Table 3 (practice effects from strict pairs) and the coarse-record table.
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

if (length(TABLES)) {
  if (requireNamespace("writexl", quietly = TRUE)) {
    writexl::write_xlsx(TABLES, file.path(OUT, "Tables.xlsx"))
    message("Saved Tables.xlsx (", length(TABLES), " sheets) and one CSV per table in ", OUT)
  } else message("writexl not installed: CSVs only, in ", OUT)
}
