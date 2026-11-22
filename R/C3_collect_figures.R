################################################################################
# C3_collect_figures.R  (21 Nov 2026: submission layout)
#
# Copies every display item from the folder of the step that makes it into
# Outputs/Manuscript/numbered, under its number in the paper, with a manifest
# (manifest_numbered.csv: file, source, present). Only that folder is shared in git.
#   Figs 1-4     R/figures/F1-F4 (Outputs/figures)
#   Tables 1-2   C2 (Outputs/Manuscript_work: ED_Table1.csv, Table2.csv)
#   Supplementary Tables 1-6, combined into Supplementary_Tables.xlsx (one sheet per table,
#   a Contents sheet with each title and its notes):
#     1 NAPESHM sites (C1m), 2 the ten experiments (C2 from B2), 3 every treatment (B8, B8b),
#     4 structural equation model (C1m from A1c), 5 commonality analysis (C1m),
#     6 N-accounting alternatives (C1m from B6)
# Anything still written straight to Outputs/Manuscript is moved to Outputs/Manuscript_work first.
################################################################################
suppressPackageStartupMessages(library(dplyr))
OUT  <- here::here("Outputs", "Manuscript")        # holds only numbered/ (shared in git)
WORK <- here::here("Outputs", "Manuscript_work")   # where C1, C1m and C2 write their working copies
O    <- function(...) here::here("Outputs", ...)
MAIN <- "yield_main_default_total_agcros"
FIGS <- "manuscript_figures_default_total_agcros"
dir.create(OUT, recursive = TRUE, showWarnings = FALSE); dir.create(WORK, recursive = TRUE, showWarnings = FALSE)
# Scripts not yet pointed at Outputs/Manuscript_work (or older runs) may have left files
# directly in Outputs/Manuscript: move them to the working folder so only numbered/ remains.
loose <- setdiff(list.files(OUT, full.names = TRUE), file.path(OUT, "numbered"))
if (length(loose)) {
  file.copy(loose, WORK, overwrite = TRUE, recursive = TRUE, copy.date = TRUE)
  unlink(loose, recursive = TRUE)
  message("Moved ", length(loose), " file(s) from Outputs/Manuscript to Outputs/Manuscript_work")
}

# ---- the paper's numbering (submission, 21 Nov 2026) -------------------------------------------
NUM <- file.path(OUT, "numbered"); dir.create(NUM, showWarnings = FALSE)
unlink(list.files(NUM, full.names = TRUE))   # start clean, so items renumbered or dropped leave no stale copies
NUMBERED <- data.frame(stringsAsFactors = FALSE, rbind(
  c("Fig1.png",        O("figures", "F1_soil_health_model.png")),  # a SEM | b path diagram
  c("Fig2.png",        O("figures", "F2_system_score.png")),       # SHMI vs soil health: all, between, within
  c("Fig3.png",        O("figures", "F3_yield.png")),              # corn and soybean | each experiment
  c("Fig4.png",        O("figures", "F4_yield_moderators.png")),   # N | age | good and bad years
  c("Table1.csv",      file.path(WORK, "ED_Table1.csv")),          # scores by system type (C2)
  c("Table2.csv",      file.path(WORK, "Table2.csv")),             # SHMI against its parts (C2)
  # Supplementary Tables 1-6, staged here, then combined into one workbook below
  c("SuppTable1.csv",  file.path(WORK, "SuppTable_sites.csv")),            # NAPESHM sites (C1m)
  c("SuppTable2.csv",  file.path(WORK, "SuppTable1.csv")),                 # the ten experiments (C2 from B2)
  c("SuppTable3.xlsx", O("SHMI_study_treatments_with_yield.xlsx")),        # every treatment (B8, B8b)
  c("SuppTable4.csv",  file.path(WORK, "SuppTable_sem.csv")),              # SEM (C1m from A1c)
  c("SuppTable5.csv",  file.path(WORK, "SuppTable_commonality.csv")),      # commonality (C1m)
  c("SuppTable6.csv",  file.path(WORK, "SuppTable_N_scenarios.csv"))))     # N accounting (C1m from B6)
names(NUMBERED) <- c("file", "source")
ok <- file.exists(NUMBERED$source)
file.copy(NUMBERED$source[ok], file.path(NUM, NUMBERED$file[ok]), overwrite = TRUE, copy.date = TRUE)
# CSVs are UTF-8 (en dashes, superscripts); Excel on Windows reads a CSV without a byte-order
# mark as Windows-1252 and shows "–" as "â€“". Add the UTF-8 byte-order mark so it opens cleanly.
BOM <- as.raw(c(0xef, 0xbb, 0xbf))
for (f in file.path(NUM, NUMBERED$file[ok & grepl("\\.csv$", NUMBERED$file)])) {
  x <- readBin(f, "raw", file.size(f))
  if (length(x) < 3 || !identical(x[1:3], BOM)) writeBin(c(BOM, x), f)
}
NUMBERED$present <- file.exists(file.path(NUM, NUMBERED$file))
# ---- Supplementary Tables 1-6 as one workbook ------------------------------------------------
# Nature takes Supplementary tables only as data files. The manuscript lists the five titles;
# the tables go in numbered/Supplementary_Tables.xlsx, one sheet per table, with a Contents
# sheet giving each title and its notes. Supplementary Table 1 also gets the corn-yield
# comparison from B2d (yield_parts.csv): AIC of the yield model with each score in place of
# SHMI minus SHMI's, its range with each experiment left out, and the P of the score's slope.
SUPP <- list(
  `Table 1` = c(
    "Supplementary Table 1: The 74 NAPESHM sites used for SHMI and lv_SH.",
    "Systems: NAPESHM treatments kept for SHMI and lv_SH (354 at 74 US sites), grouped as in the main-text Table 1. First treatment year: earliest treatment start among a site's systems; years under treatment: 2019 (NAPESHM sampling) minus that year. MAT, MAP: site mean annual temperature and precipitation from the NAPESHM site metadata.",
    "Within-site contrasts: annual systems at the site differ in the NAPESHM treatment label for tillage, cover crops, rotation length or organic nutrient source; perennial vs annual: the site has both annual and perennial or reference systems. SHMI and lv_SH ranges: lowest-highest system value at the site."),
  `Table 2` = c(
    "Supplementary Table 2: The ten long-term experiments used for the yield analysis.",
    "Design and records: treatment structure and data sources. Years, treatments, plots and corn plot-years are those analysed (treatments with corn; rotation phases combined). N: fertilizer plus manure (total N) plus legume credit.",
    "Share of SHMI spread: the share of each experiment's within-site SHMI variance from each part of SHMI (weight x Cover/Roots, Diversity and Inverse Disturbance, and the amendment bonus; cov(part, SHMI)/var(SHMI)); the four add to 100%, and a negative share means that part moved against SHMI. Experiments are ordered by their slope at 150 kg N ha⁻¹, as in Fig. 3b. Slope, experiment's own N and age: each experiment's own within-site model, averaged over its N rates and years. Slope at 150 kg N ha⁻¹ and the average age: the same model standardized with the pooled SHMI × N and SHMI × age coefficients, as plotted in Fig. 3b."),
  `Table 3` = c(
    "Supplementary Table 3: Every treatment in the ten experiments.",
    "Treatment design, SHMI and its sub-indices, and mean corn (grain dry matter) and soybean yields over the analysis period. Treatments that differ only in rotation phase are listed together. Treatments without corn (corn plot-years = 0) are not in the corn-yield analysis. Other factors: treatment factors SHMI does not score, entered as a covariate in the yield models. Full management records and the code that builds them are in the repository (Code availability)."),
  `Table 4` = c(
    "Supplementary Table 4: Structural equation model for soil health (lv_SH).",
    "Estimates from the final model fitted to the 354 US NAPESHM systems (robust maximum likelihood, full-information maximum likelihood for missing values, standard errors clustered by site). Latent variances were fixed at 1; the two second-order loadings labelled a were constrained equal (likelihood-ratio test, P = 0.95). Standardized: completely standardized loadings; for residual covariances, the residual correlation. Indicators were adjusted for climate and texture and standardized before fitting (Methods)."),
  `Table 5` = c(
    "Supplementary Table 5: Commonality analysis of soil health (lv_SH) on the SHMI components.",
    "Ordinary least-squares regression of lv_SH on Cover/Roots, Diversity, Inverse Disturbance and amendment frequency (share of years amended), across all 354 systems and within sites (lv_SH and every predictor centred on its site mean). Each part is a unique or shared contribution from Möbius inversion of the R² values of all 15 subset regressions; % of explained is the part divided by the full-model R². Negative parts indicate suppression between predictors."),
  `Table 6` = c(
    "Supplementary Table 6: Nitrogen-accounting alternatives for the corn-yield model.",
    "Each row is the within-site corn-yield model refitted by maximum likelihood on the same 8,674 plot-years with N supply computed under one combination of legume credits and manure-N availability (Supplementary Methods). Mean N: mean N supply per corn plot-year. ΔAIC: difference from the best-fitting combination (default legume credits, manure N counted in full). SHMI slope at 150 kg N ha⁻¹ and the mean age. Models were fitted by maximum likelihood so that AICs are comparable; the SHMI slope for the selected accounting (204) therefore differs slightly from the REML estimate in the main text (208 kg ha⁻¹ per 10 points)."))
# Readable column names for the tables B2 writes in raw form (Supplementary Tables 2 and 6)
ci_txt <- function(x, lo, hi, d = 0) {
  r <- function(v) round(v, d) + 0   # + 0 turns -0 into 0
  x <- r(x); lo <- r(lo); hi <- r(hi)
  ifelse(is.na(x), "", sprintf(paste0("%.", d, "f (%.", d, "f to %.", d, "f)"), x, lo, hi))
}
tidy_supp <- function(k, tab) {
  if (k == 2 && all(c("site", "share_cover", "slope_adj") %in% names(tab))) {   # the ten experiments (Supplementary Table 2)
    pct <- function(x) sprintf("%.0f", 100 * x)
    tab <- tab[order(tab$slope_adj, tab$site), ]     # 17 Nov 2026: Fig. 3b order (slope)
    tab <- data.frame(check.names = FALSE,
      Experiment = tab$site, Location = tab$location, Years = paste0(tab$first_year, "\u2013", tab$last_year),
      Treatments = tab$treatments, Plots = tab$plots, `Corn plot-years` = tab$corn_plot_years,
      `N, kg ha-1 (range)` = paste0(tab$N_min, "\u2013", tab$N_max),
      `SHMI range (mean)` = sprintf("%.0f\u2013%.0f (%.0f)", tab$SHMI_min, tab$SHMI_max, tab$SHMI_mean),
      `Share of SHMI spread, %: Cover / Diversity / Disturbance / Amendment` =
        paste(pct(tab$share_cover), pct(tab$share_div), pct(tab$share_dist), pct(tab$share_amend), sep = " / "),
      `Slope per 10 SHMI, kg ha-1 (95% CI), experiment's own N and age` = ci_txt(tab$slope, tab$slope_lower, tab$slope_upper),
      `Mean N, kg ha-1` = tab$N_mean, `Mean system age, years` = tab$age_mean,
      `Slope at 150 kg N ha-1 and the average age (95% CI); Fig. 3b` = ci_txt(tab$slope_adj, tab$slope_adj_lower, tab$slope_adj_upper))
  }
  tab
}
if (requireNamespace("writexl", quietly = TRUE) && requireNamespace("readxl", quietly = TRUE)) {
  sheets <- list(); contents <- character()
  for (k in seq_along(SUPP)) {
    nm <- names(SUPP)[k]; stem <- paste0("SuppTable", k)
    f_csv <- file.path(NUM, paste0(stem, ".csv")); f_x <- file.path(NUM, paste0(stem, ".xlsx"))
    contents <- c(contents, SUPP[[k]], "")
    if (file.exists(f_csv)) {
      tab <- read.csv(f_csv, check.names = FALSE, stringsAsFactors = FALSE, fileEncoding = "UTF-8-BOM")
      sheets[[nm]] <- tidy_supp(k, tab); unlink(f_csv)
    } else if (file.exists(f_x)) {                     # a workbook: its first sheet is the table, the rest follow
      sh <- readxl::excel_sheets(f_x)
      for (i in seq_along(sh)) {
        lab <- if (i == 1) nm else substr(paste(nm, sh[i]), 1, 31)
        sheets[[lab]] <- as.data.frame(readxl::read_excel(f_x, sheet = sh[i]))
      }
      unlink(f_x)
    } else { contents <- c(contents[seq_len(length(contents) - 1)], "  (not yet produced)", ""); message("  missing: ", stem) }
  }
  wb <- c(list(Contents = data.frame(`Supplementary Tables` = contents, check.names = FALSE)), sheets)
  writexl::write_xlsx(wb, file.path(NUM, "Supplementary_Tables.xlsx"))
  NUMBERED$file[grepl("^SuppTable", NUMBERED$file)] <- paste0("Supplementary_Tables.xlsx [", names(SUPP), "]")
  NUMBERED$present[grepl("^Supplementary_Tables", NUMBERED$file)] <- names(SUPP) %in% names(sheets)
} else message("writexl or readxl not installed: Supplementary tables left as separate files")
write.csv(NUMBERED, file.path(NUM, "manifest_numbered.csv"), row.names = FALSE)
message("Numbered as in the paper (Outputs/Manuscript/numbered): ", sum(NUMBERED$present), " of ", nrow(NUMBERED))
if (any(!NUMBERED$present)) message("  missing: ", paste(NUMBERED$file[!NUMBERED$present], collapse = "; "))
