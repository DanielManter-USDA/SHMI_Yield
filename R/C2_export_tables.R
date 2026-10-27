################################################################################
# C2_export_tables.R
#
# Every table in the paper, from the outputs of the steps that compute it, in one
# workbook (one sheet per table) plus a CSV per table:
#   Table 1                      C1_Fig1_and_Table1.R          Outputs/Manuscript/Table1.csv
#   Extended Data Table 1        scores by system type          data/SHAPE_by_unit_with_SHMI.csv
#   Supplementary Data Table 1   the ten experiments            B2 site_SHMI_ranges.csv
#   Supplementary Data Table 2   coarse-record scenarios        B9 Coarse/coarse_scenarios.csv
#   Supplementary Data Table 3   soil properties as carriers    B5 mediation_results.csv
#   Supplementary Data Table 4   SHMI against its parts          A4a Benchmark/benchmark_bootstrap.csv
# (Table 2, "Scoring a field with SHMI", is text and has no computed values.)
# A table whose source is missing (its step was not run) is skipped with a message.
#
# Output  Outputs/Manuscript/Tables.xlsx and Outputs/Manuscript/<table>.csv
################################################################################
suppressPackageStartupMessages({ library(dplyr); library(tidyr) })
OUT  <- here::here("Outputs", "Manuscript")
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

# Table 1 (written by C1)
add("Table1", file.path(OUT, "Table1.csv"), function(f) read.csv(f, check.names = FALSE, encoding = "UTF-8"))

# Extended Data Table 1: median (lowest-highest) by system type
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

# Supplementary Data Table 2: record types and scenarios of the coarse-record test (B9)
add("SuppTable2", here::here("Outputs", "Coarse", "coarse_scenarios.csv"), function(f) read.csv(f, check.names = FALSE))

# Supplementary Data Table 3: soil properties as carriers of the SHMI-yield slope (B5)
add("SuppTable3", file.path(MAIN, "mediation_results.csv"), function(f) {
  m <- read.csv(f)
  ci <- function(e, lo, hi, d = 1) sprintf("%s (%s to %s)", formatC(e, format = "f", digits = d),
                                           formatC(lo, format = "f", digits = d), formatC(hi, format = "f", digits = d))
  main <- m %>% filter(!grepl(", without", mediator))
  loo  <- m %>% filter(grepl(", without", mediator)) %>%
    mutate(base = sub(", without.*$", "", mediator)) %>%
    group_by(base) %>% summarise(loo = sprintf("%.0f to %.0f", 100 * min(share_mediated), 100 * max(share_mediated)), .groups = "drop")
  main %>% left_join(loo, by = c(mediator = "base")) %>%
    transmute(`Soil measure` = mediator,
              `Experiments (systems; plot-years)` = sprintf("%d (%d; %s)", n_sites, n_units, format(plot_years, big.mark = ",")),
              `r with SHMI` = round(r_SHMI_M, 2),
              `Change per 10 SHMI (95% CI)` = ci(a_per10, a_lo, a_hi, 2),
              `Yield per unit, kg ha-1 (P)` = sprintf("%.1f (P = %.2f)", b_per_unit, b_p),
              `Total SHMI slope` = round(c_total_per10),
              `Indirect effect (95% CI)` = ci(indirect_per10, indirect_lo, indirect_hi, 0),
              `Share carried, % (95% CI)` = ci(100 * share_mediated, 100 * share_lo, 100 * share_hi, 0),
              `Share, leaving out one site, %` = loo)
})

# Supplementary Data Table 4: SHMI against its variants, parts and unconstrained models (A4a)
add("SuppTable4", here::here("Outputs", "Benchmark", "benchmark_bootstrap.csv"), function(f) {
  b <- read.csv(f)
  lab <- c("SHMI" = "SHMI (weights refitted in each fold)",
           "SHMI, published weights (fixed)" = "SHMI, published weights held fixed",
           "SHMI_base (no amendment term)" = "SHMI_base (no amendment term)",
           "SHMI without Cover/Roots" = "SHMI without Cover/Roots",
           "SHMI without Diversity" = "SHMI without Diversity",
           "SHMI without InvDist" = "SHMI without Inverse Disturbance",
           "SHMI without amendments" = "SHMI without amendment term (refitted)",
           "SHMI, equal weights" = "SHMI, equal weights",
           "Cover/Roots" = "Cover/Roots alone", "Diversity" = "Diversity alone",
           "InvDist" = "Inverse Disturbance alone", "Amendments" = "Share of years amended alone",
           "OLS on all inputs" = "Linear model (OLS)", "Random forest on all inputs" = "Random forest")
  metric <- c(R2 = "Overall R²", between_r = "r between site means", within_r = "r within sites")
  b %>% filter(metric %in% names(.env$metric), model %in% names(lab)) %>%
    mutate(cell_value = sprintf("%.2f", value_avg_oof),
           cell_diff = ifelse(model == "SHMI", "—",
                              sprintf("%+.3f (%+.3f to %+.3f)%s", diff_from_SHMI, diff_lo, diff_hi,
                                      ifelse(diff_lo > 0 | diff_hi < 0, "*", ""))),
           model = factor(lab[model], unname(lab)), metric = .env$metric[metric]) %>%
    select(model, metric, Value = cell_value, `Difference from SHMI (95% CI)` = cell_diff) %>%
    pivot_wider(names_from = metric, values_from = c(Value, `Difference from SHMI (95% CI)`), names_glue = "{metric}: {.value}") %>%
    arrange(model) %>% rename(Score = model) %>% mutate(Score = as.character(Score))
})

if (length(TABLES)) {
  if (requireNamespace("writexl", quietly = TRUE)) {
    writexl::write_xlsx(TABLES, file.path(OUT, "Tables.xlsx"))
    message("Saved Tables.xlsx (", length(TABLES), " sheets) and one CSV per table in ", OUT)
  } else message("writexl not installed: CSVs only, in ", OUT)
}
