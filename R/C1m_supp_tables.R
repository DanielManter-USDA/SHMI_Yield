################################################################################
# C1m_supp_tables.R  (21 Nov 2026)
#
# The Supplementary Tables not written by another step, as CSV files that C3 combines
# into Outputs/Manuscript/numbered/Supplementary_Tables.xlsx:
#   SuppTable_sites.csv        Supplementary Table 1  the 74 NAPESHM sites
#   SuppTable_sem.csv          Supplementary Table 4  the structural equation model (A1c)
#   SuppTable_commonality.csv  Supplementary Table 5  commonality analysis of lv_SH
#   SuppTable_N_scenarios.csv  Supplementary Table 6  N-accounting alternatives (B6)
# (Supplementary Table 2, the ten experiments, is written by C2 from B2; Supplementary
# Table 3, every treatment, by B8/B8b.)
#
# Inputs  Outputs/SHAPE/SHAPE_by_unit_with_SHMI.csv (else data/), the 354 systems with
#           lv_SH (SH), SHMI, its parts, MAT and MAP
#         Outputs/semDF_prepped.rds (A1b: state, treatment start year, system type)
#         Inputs/napeshm/treatments.csv (NAPESHM treatment labels, for the contrasts)
#         Outputs/sem_noMX_report.rds or sem_noMX_projection.rds (A1c: the final US model)
#         Outputs/yield_main_default_total_agcros/N_scenarios_all_sites.csv (B6)
# A table whose input is missing is skipped with a message.
################################################################################
suppressPackageStartupMessages({ library(dplyr) })
setwd(here::here())
WORK <- "Outputs/Manuscript_work"; dir.create(WORK, recursive = TRUE, showWarnings = FALSE)
EN <- "–"
rng <- function(x) sprintf("%.0f%s%.0f", min(x), EN, max(x))
save_csv <- function(d, f) { write.csv(d, file.path(WORK, f), row.names = FALSE, na = "", fileEncoding = "UTF-8"); message("C1m: wrote ", f) }

F_UNIT <- if (file.exists("Outputs/SHAPE/SHAPE_by_unit_with_SHMI.csv")) "Outputs/SHAPE/SHAPE_by_unit_with_SHMI.csv" else "data/SHAPE_by_unit_with_SHMI.csv"
units <- read.csv(F_UNIT, stringsAsFactors = FALSE) %>%
  filter(!is.na(SH), !is.na(SHMI), !is.na(SHAPE_mean)) %>%          # the 354 systems (same filter as C1j)
  mutate(sitecode = as.character(sitecode))

# ---- Supplementary Table 1: the NAPESHM sites --------------------------------------------
if (all(file.exists(c("Outputs/semDF_prepped.rds", "Inputs/napeshm/treatments.csv")))) {
  sem <- readRDS("Outputs/semDF_prepped.rds") %>%
    transmute(MGT_combo, state_prov, start = suppressWarnings(as.integer(treatment_start_year)))
  trt <- read.csv("Inputs/napeshm/treatments.csv", stringsAsFactors = FALSE) %>%
    transmute(MGT_combo = paste0("NAPESHM_", sitecode, "_", treatmentid),
              till = tolower(trimws(eu_tillage)), cover = tolower(trimws(eu_cover_crops)),
              rot = tolower(trimws(eu_crop_rot)),
              organic = tolower(trimws(eu_nutr_type)) %in% c("organic", "synthetic_and_organic"))
  u <- units %>% left_join(sem, by = "MGT_combo") %>% left_join(trt, by = "MGT_combo") %>%
    mutate(state = state.abb[match(state_prov, state.name)])
  # a contrast: annual systems at the site differ in the NAPESHM label (organic: organic or
  # synthetic-and-organic nutrient source vs not); perennial vs annual: both kinds present
  contrasts <- function(g) {
    a <- g[g$system == "annual", ]
    k <- c(tillage = n_distinct(a$till) > 1, `cover crop` = n_distinct(a$cover) > 1,
           rotation = n_distinct(a$rot) > 1, `organic amendment` = n_distinct(a$organic) > 1,
           `perennial vs annual` = any(g$system == "annual") && any(g$system %in% c("perennial", "reference")))
    k
  }
  rows <- lapply(split(u, u$sitecode), function(g) {
    k <- contrasts(g); first <- min(g$start, na.rm = TRUE)
    data.frame(check.names = FALSE, Site = g$sitecode[1], State = g$state[1], Systems = nrow(g),
               Annual = sum(g$system == "annual"), Perennial = sum(g$system == "perennial"),
               Reference = sum(g$system == "reference"), Fallow = sum(g$system == "fallow"),
               `First treatment year` = as.character(first), `Years under treatment in 2019` = 2019 - first,
               `MAT (°C)` = round(g$MATC[1], 1), `MAP (mm)` = round(g$MAPmm[1]),
               `Within-site contrasts` = if (any(k)) paste(names(k)[k], collapse = "; ") else "—",
               `SHMI range` = rng(g$SHMI), `lv_SH range` = rng(g$SH), k = I(list(k)))
  })
  st1 <- bind_rows(rows)
  kk <- colSums(do.call(rbind, st1$k)); st1$k <- NULL
  all_row <- data.frame(check.names = FALSE, Site = "All sites", State = sprintf("%d states", n_distinct(u$state)),
                        Systems = nrow(u), Annual = sum(u$system == "annual"), Perennial = sum(u$system == "perennial"),
                        Reference = sum(u$system == "reference"), Fallow = sum(u$system == "fallow"),
                        `First treatment year` = paste0(min(as.integer(st1$`First treatment year`)), EN, max(as.integer(st1$`First treatment year`))),
                        `Years under treatment in 2019` = paste("median", format(median(st1$`Years under treatment in 2019`))),
                        `MAT (°C)` = NA, `MAP (mm)` = NA,
                        `Within-site contrasts` = paste0(paste(names(kk), kk, collapse = "; "), " sites"),
                        `SHMI range` = rng(u$SHMI), `lv_SH range` = rng(u$SH))
  st1$`Years under treatment in 2019` <- as.character(st1$`Years under treatment in 2019`)
  save_csv(bind_rows(st1, all_row), "SuppTable_sites.csv")
} else message("C1m: Supplementary Table 1 skipped (needs Outputs/semDF_prepped.rds and Inputs/napeshm/treatments.csv)")

# ---- Supplementary Table 4: the structural equation model --------------------------------
# The final fit is read from A1c's report if present, otherwise from A1c's projection (the same
# fitted model, also read by C1k and A8); reliability and determinacy are then recomputed exactly
# as A1c computes them.
REP  <- "Outputs/sem_noMX_report.rds"
PROJ <- "Outputs/sem_noMX_projection.rds"
determinacy <- function(fit) {   # as factor_determinacy() in A1c
  L     <- lavaan::lavInspect(fit, "est")$lambda
  Phi   <- lavaan::lavInspect(fit, "cov.lv")[colnames(L), colnames(L)]
  Sigma <- lavaan::lavInspect(fit, "cov.ov")[rownames(L), rownames(L)]
  C     <- Phi %*% t(L)
  round(sqrt(diag(C %*% solve(Sigma, t(C))) / diag(Phi)), 3)
}
rep <- if (file.exists(REP)) readRDS(REP) else if (file.exists(PROJ)) {
  f <- readRDS(PROJ)$fit
  list(final_fit = f, determinacy = determinacy(f),
       reliability = if (requireNamespace("semTools", quietly = TRUE))
         tryCatch(semTools::compRelSEM(f, higher = "SH"), error = function(e) NULL))
}
if (!is.null(rep) && requireNamespace("lavaan", quietly = TRUE)) {
  fit <- rep$final_fit
  LAT <- c(SH = "Soil health (lv_SH)", SH.Bio = "Biological activity", SH.CN = "Carbon and nitrogen pools",
           SH.Phys = "Physical condition")
  nm <- function(x) ifelse(x %in% names(LAT), LAT[x], x)
  pe <- lavaan::parameterEstimates(fit, standardized = TRUE)
  fmt_p <- function(p) ifelse(is.na(p), "", ifelse(p < 0.001, "<0.001", sprintf("%.3f", p)))
  par_row <- function(d, type, note = "") data.frame(check.names = FALSE, Parameter = type, From = nm(d$lhs), To = nm(d$rhs),
    Estimate = sprintf("%.3f", d$est), SE = ifelse(is.na(d$se) | d$se == 0, "", sprintf("%.3f", d$se)),
    z = ifelse(is.na(d$z), "", sprintf("%.2f", d$z)), P = fmt_p(d$pvalue), Standardized = sprintf("%.2f", d$std.all),
    Note = ifelse(nzchar(d$label), paste0("constrained equal (label ", d$label, ")"), note))
  so  <- pe %>% filter(op == "=~", lhs == "SH")
  fo  <- pe %>% filter(op == "=~", lhs != "SH")
  rc  <- pe %>% filter(op == "~~", lhs != rhs)
  r2  <- lavaan::inspect(fit, "r2"); r2 <- r2[intersect(names(r2), unique(fo$rhs))]
  fm  <- lavaan::fitMeasures(fit, c("chisq.scaled", "df", "cfi.robust", "tli.robust", "rmsea.robust",
                                     "rmsea.ci.lower.robust", "rmsea.ci.upper.robust", "srmr"))
  val <- function(type, from, v, note = "") data.frame(check.names = FALSE, Parameter = type, From = from, To = "",
    Estimate = v, SE = "", z = "", P = "", Standardized = "", Note = note)
  rel <- unlist(rep$reliability); det <- unlist(rep$determinacy)
  st4 <- bind_rows(
    par_row(so, "Second-order loading"), par_row(fo, "First-order loading"),
    par_row(rc, "Residual covariance", "standardized value is the residual correlation"),
    data.frame(check.names = FALSE, Parameter = "Indicator R²", From = names(r2), To = "", Estimate = "", SE = "", z = "", P = "",
               Standardized = sprintf("%.2f", r2), Note = "share of the indicator's variance explained by its factor"),
    val("Model fit", "Scaled χ² (df)", sprintf("%.1f (%d)", fm[["chisq.scaled"]], as.integer(fm[["df"]])),
        "354 systems; MLR, FIML, standard errors clustered by site"),
    val("Model fit", "Robust CFI", sprintf("%.3f", fm[["cfi.robust"]])),
    val("Model fit", "Robust TLI", sprintf("%.3f", fm[["tli.robust"]])),
    val("Model fit", "Robust RMSEA (90% CI)", sprintf("%.3f (%.3f%s%.3f)", fm[["rmsea.robust"]], fm[["rmsea.ci.lower.robust"]], EN, fm[["rmsea.ci.upper.robust"]])),
    val("Model fit", "SRMR", sprintf("%.3f", fm[["srmr"]])),
    val("Model fit", "Log10-transformed indicators (skewness > 1), same model",
        "robust CFI 0.924, TLI 0.890, RMSEA 0.135, SRMR 0.048", "Supplementary Methods"),
    if (length(rel)) val("Composite reliability", nm(names(rel)), sprintf("%.2f", rel)),
    if (length(det)) val("Factor-score determinacy", nm(names(det)), sprintf("%.2f", det),
                         "correlation of factor scores with the factor"))
  save_csv(st4, "SuppTable_sem.csv")
} else message("C1m: Supplementary Table 4 skipped (needs ", REP, " or ", PROJ, ", and lavaan)")

# ---- Supplementary Table 5: commonality analysis ------------------------------------------
# OLS of lv_SH on the three sub-indices and amendment frequency, across all systems and within
# sites (everything centred on its site mean); Mobius inversion of the R2 of all 15 subsets.
PRED <- c(`Cover/Roots` = "Cover", Diversity = "Diversity", `Inverse Disturbance` = "InvDist", Amendments = "p_amend")
commonality <- function(y, X) {
  v <- names(X)
  r2 <- function(S) if (!length(S)) 0 else summary(lm(y ~ ., data = X[, S, drop = FALSE]))$r.squared
  subsets <- unlist(lapply(seq_along(v), function(k) combn(v, k, simplify = FALSE)), recursive = FALSE)
  R2 <- setNames(sapply(subsets, r2), sapply(subsets, paste, collapse = "|"))
  getR2 <- function(S) if (!length(S)) 0 else R2[[paste(v[v %in% S], collapse = "|")]]
  total <- getR2(v)
  parts <- sapply(subsets, function(S) {
    s <- 0
    for (j in 0:length(S)) for (Tt in if (j == 0) list(character()) else combn(S, j, simplify = FALSE))
      s <- s + (-1)^(length(S) - length(Tt) + 1) * getR2(setdiff(v, Tt))
    s
  })
  lab <- sapply(subsets, function(S) if (length(S) == 1) paste("Unique to", S) else paste("Shared:", paste(S, collapse = " & ")))
  list(parts = setNames(parts, lab), total = total)
}
if (all(PRED %in% names(units))) {
  X <- setNames(units[, PRED], names(PRED))
  centre <- function(x) x - ave(x, units$sitecode)
  a <- commonality(units$SH, X)
  w <- commonality(centre(units$SH), as.data.frame(lapply(X, centre), check.names = FALSE))
  st5 <- data.frame(check.names = FALSE, Component = c(names(a$parts), "Total (full-model R²)"),
                    `All systems: R² part` = round(c(a$parts, a$total), 3),
                    `All systems: % of explained` = round(100 * c(a$parts, a$total) / a$total, 1),
                    `Within sites: R² part` = round(c(w$parts, w$total), 3),
                    `Within sites: % of explained` = round(100 * c(w$parts, w$total) / w$total, 1))
  save_csv(st5, "SuppTable_commonality.csv")
  message(sprintf("C1m: Cover/Roots & Diversity shared %.0f%% of the explained variation; Diversity unique %.0f%%",
                  100 * a$parts[["Shared: Cover/Roots & Diversity"]] / a$total, 100 * a$parts[["Unique to Diversity"]] / a$total))
} else message("C1m: Supplementary Table 5 skipped (", F_UNIT, " lacks ", paste(setdiff(PRED, names(units)), collapse = ", "), ")")

# ---- Supplementary Table 6: N-accounting alternatives (B6) --------------------------------
NS <- "Outputs/yield_main_default_total_agcros/N_scenarios_all_sites.csv"
if (file.exists(NS)) {
  ns <- read.csv(NS, stringsAsFactors = FALSE) %>% arrange(dAIC_best)
  st6 <- data.frame(check.names = FALSE, `Legume credit` = ns$legume,
                    `Manure-N availability` = ifelse(ns$manure == "total", "total N", ns$manure),
                    `Mean N, kg ha⁻¹` = ns$N_mean, `ΔAIC from best` = ns$dAIC_best,
                    `SHMI slope, kg ha⁻¹ per 10 points` = ns$SHMI_slope, SE = ns$SHMI_se,
                    `SHMI × N, kg ha⁻¹ per 10 points per 100 kg N ha⁻¹` = ns$SHMIxN)
  save_csv(st6, "SuppTable_N_scenarios.csv")
} else message("C1m: Supplementary Table 6 skipped (", NS, " not found; run B6)")
