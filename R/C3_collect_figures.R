################################################################################
# C3_collect_figures.R
#
# Copies every display item from the folder of the step that makes it into
# Outputs/Manuscript/numbered, under its number in the paper (manuscript v5,
# 6 Nov 2026: items renumbered in order of first citation; 7 Nov 2026: Fig. 3 corn and
# soybean, Fig. 4 SHMI x N, x age and the experiments, stability to Extended Data Fig. 5,
# soybean Supplementary Fig. removed), with a manifest
# (manifest_numbered.csv: file, source, present). Only that folder is shared in git.
# Fig. 1, Fig. 2, Table 1, Extended Data Figs 2 and 6 (written as ED_Fig2.png and
# ED_Fig5.png), Extended Data Table 1 and the Supplementary Data Tables (written by C2
# under their earlier numbers) are made by their own scripts in Outputs/Manuscript_work
# (not shared); this script only copies them. Anything still written straight to
# Outputs/Manuscript is moved to Outputs/Manuscript_work first.
#
# Extended Data Fig. 1 (the soil health model diagram) is drawn by hand from the SEM
# loadings exported by A8 (Outputs/sem_noMX_standardized.csv), so the loadings file is
# copied in its place.
################################################################################
suppressPackageStartupMessages(library(dplyr))
OUT  <- here::here("Outputs", "Manuscript")        # holds only numbered/ (shared in git)
WORK <- here::here("Outputs", "Manuscript_work")   # where C1, C1b, C1c, A4b and C2 write their copies
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

# ---- the paper's current numbering (manuscript v5, 6 Nov 2026) --------------------------------
# Display items were renumbered in order of first citation; the files above keep their
# script names. This copies each one to Outputs/Manuscript/numbered/ under the number it
# has in the paper, so what goes to the journal can be read off the file names.
NUM <- file.path(OUT, "numbered"); dir.create(NUM, showWarnings = FALSE)
unlink(list.files(NUM, full.names = TRUE))   # start clean, so items renumbered or dropped leave no stale copies
NUMBERED <- data.frame(stringsAsFactors = FALSE, rbind(
  # 7 Nov 2026: three main figures (Nature 8-page budget). Fig. 1 merges old Figs 1b-d and 2a
  # (C1); Fig. 2 = corn and soybean (B7 writes it as Fig3.png); Fig. 3 = N, age and the
  # experiments (B7 writes it as Fig4.png). The path diagram and the validation scatterplots
  # move to Extended Data, which is renumbered by first citation.
  # 10 Nov 2026 layout: Fig. 1 SHMI and soil health; Fig. 2 label vs SHMI, data budget, parts
  # (both from C1); Fig. 3 corn, soybean and good/bad years (B7's Fig3.png); Fig. 4 N, age and
  # experiments (B7's Fig4.png). Extended Data: 1 SEM (drawn by hand; loadings), 2 practices do
  # not add up (C1), 3 SHMI within label categories (C1g), 4 SHMI and SHAPE (C1h), 5 coarse records.
  c("Fig1.png",                  file.path(WORK, "Fig1.png")),
  c("Fig2.png",                  file.path(WORK, "Fig2.png")),
  c("Fig3.png",                  O("figures", "F3_yield.png")),                       # corn | soybean | years (R/figures/F3_yield.R)
  c("Fig4.png",                  O("figures", "F4_yield_moderators.png")),            # SHMI x N | x age | experiments (R/figures/F4_yield_moderators.R)
  c("ED_Fig1_loadings.csv",      O("sem_noMX_standardized.csv")),                    # SEM; the diagram is drawn by hand
  c("ED_Fig2.png",               file.path(WORK, "ED_Fig_additivity.png")),          # practices do not add up (C1)
  c("ED_Fig3.png",               O("Practice_additivity", "Fig_label_within.png")),  # SHMI within label categories (C1g)
  c("ED_Fig4.png",               O("Practice_additivity", "Fig_SHMI_vs_SHAPE.png")), # SHMI and SHAPE (C1h)
  c("ED_Fig5.png",               file.path(WORK, "ED_Fig5.png")),                      # coarse records (C1c)
  c("Table1.csv",                file.path(WORK, "Table1.csv")),
  c("ED_Table1.csv",             file.path(WORK, "ED_Table1.csv")),
  # Supplementary Tables 1-5, renumbered by first citation (soil properties as carriers, B5,
  # cut to a sentence); copied here as staging files, then combined into one workbook below.
  c("SuppTable1.csv",            file.path(WORK, "SuppTable4.csv")),                   # SHMI vs variants (A4a) + yield (B2d)
  c("SuppTable2.csv",            file.path(WORK, "SuppTable1.csv")),                   # the ten experiments (B2)
  c("SuppTable3.csv",            O(MAIN, "SHMI_x_other_by_site.csv")),                # slope by level of unscored factors (B2)
  c("SuppTable4.xlsx",           O("SHMI_study_treatments_with_yield.xlsx")),         # every treatment (B8, B8b)
  c("SuppTable5.csv",            file.path(WORK, "SuppTable2.csv")),                   # coarse-record scenarios (B9)
  c("SuppTable6.csv",            file.path(WORK, "SuppTable6_practice_pairs.csv"))))   # matched practice pairs (C1)
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
# Extended Data Fig. 4 (9 Nov 2026): a = practice labels vs SHMI for soil health (C1e);
# b = practice pairs for corn yield (B2h). Stacked here when both exist; else a alone.
PAIRS_PNG <- O("yield_paired_budget", "Fig_practice_pairs.png")
ED4_A     <- O("Practice_additivity", "Fig_labels_vs_SHMI.png")
if (!exists("ED_PAIRS_PANEL")) ED_PAIRS_PANEL <- FALSE   # 10 Nov 2026: practice-pairs panel cut from the paper
if (ED_PAIRS_PANEL && file.exists(PAIRS_PNG) && file.exists(ED4_A) && requireNamespace("magick", quietly = TRUE)) {
  lab_ed <- function(img, letter) {
    w <- magick::image_info(img)$width
    img <- magick::image_append(c(magick::image_blank(w, 70, "white"), img), stack = TRUE)
    magick::image_annotate(img, letter, size = 80, weight = 700, gravity = "northwest", location = "+15+0")
  }
  a  <- magick::image_read(ED4_A); b <- magick::image_read(PAIRS_PNG)
  b  <- magick::image_scale(b, as.character(magick::image_info(a)$width))
  ed4 <- magick::image_append(c(lab_ed(a, "a"), lab_ed(b, "b")), stack = TRUE)
  magick::image_write(magick::image_background(ed4, "white"), file.path(NUM, "ED_Fig3.png"), format = "png", density = 300)
  file.copy(O("yield_paired_budget", "practice_pairs_summary.csv"), file.path(NUM, "ED_Fig3b_practice_pairs.csv"), overwrite = TRUE)
} else if (ED_PAIRS_PANEL && !file.exists(PAIRS_PNG)) message("Extended Data Fig. 4 without panel b: ", PAIRS_PNG, " not found (run B2h_practice_pairs.R)")
NUMBERED$present <- file.exists(file.path(NUM, NUMBERED$file))
# ---- Supplementary Tables 1-5 as one workbook (8 Nov 2026) ----------------------------------
# Nature takes Supplementary tables only as data files. The manuscript lists the five titles;
# the tables go in numbered/Supplementary_Tables.xlsx, one sheet per table, with a Contents
# sheet giving each title and its notes. Supplementary Table 1 also gets the corn-yield
# comparison from B2d (yield_parts.csv): AIC of the yield model with each score in place of
# SHMI minus SHMI's, its range with each experiment left out, and the P of the score's slope.
SUPP <- list(
  `Table 1` = c(
    "Supplementary Table 1: SHMI against its variants and its parts, for soil health and corn yield.",
    "Soil health: out-of-fold performance at US sites (354 systems, 74 sites) on the calibration's folds (10-fold leave-sites-out cross-validation repeated 100 times), from the repeat-averaged out-of-fold predictions (overall R2; Pearson r of site means; Pearson r of deviations from site means). SHMI's row gives its own value; the other rows give the candidate minus SHMI, with 95% intervals from a site-level bootstrap (2,000 resamples); negative means worse than SHMI; * interval excludes zero.",
    "Weights were refitted on the training sites in every fold for SHMI and its variants, except the published-weights version, which refits only the line mapping the score to soil health. Equal weights: one third each to Cover/Roots, Diversity and Inverse Disturbance, equal season shares, amendment term unchanged. Cover/Roots alone: season split fitted on the training sites. Every score was mapped to lv_SH by a site-weighted line fitted on the training sites.",
    "Corn yield: the yield model of Fig. 2a on the same data (ten experiments) with each score in place of SHMI, using the published weights (variants without a part rescale the remaining weights to sum to 1). Slope: the yield change per within-site standard deviation of the score, so scores on different scales compare; % of SHMI's: that slope relative to SHMI's. P, score adds to SHMI: likelihood-ratio test of SHMI against SHMI plus the score (maximum likelihood); P, SHMI adds to score: the score against SHMI plus the score. A small P in the first means the score carries yield information SHMI misses; in the second, that SHMI carries information the score misses."),
  `Table 2` = c(
    "Supplementary Table 2: The ten experiments and what drove their SHMI range.",
    "Corn years, treatments and plots are those analysed (treatments with corn; rotation phases combined). N: fertilizer plus manure (total N) plus legume credit.",
    "Share of SHMI spread: the share of each experiment's within-site SHMI variance from each part of SHMI (weight x Cover, Diversity and Inverse Disturbance, and the amendment bonus; cov(part, SHMI)/var(SHMI)); the four add to 100%, and a negative share means that part moved against SHMI. Experiments are ordered by their Diversity share, as in Fig. 3c. Slope: each experiment's own within-site model, averaged over its N rates and years."),
  `Table 3` = c(
    "Supplementary Table 3: The SHMI slope by level of the treatment factors SHMI does not score.",
    "For each experiment with residue removal, irrigation, organic management or a planter treatment: the slope per 10 SHMI points in each level, at the experiment's own N and age and moved to 150 kg N/ha and the average age, with the SHMI x factor interaction F test."),
  `Table 4` = c(
    "Supplementary Table 4: Every treatment in the ten experiments.",
    "Description, SHMI and its parts (points from each sub-index = weight x sub-index), N, and mean corn and soybean yields, with notes on data corrections (IAAM7071 cover crops and harvest dates; INWLACRE 2004 planting dates). Further sheets of the treatment workbook follow as 'Table 4 ...'."),
  `Table 5` = c(
    "Supplementary Table 5: Every coarse-record scenario and its results.",
    "For each scenario, what was degraded, soil-health skill on the calibration's folds (correlation within sites, share of sites ranked in the right direction, correlation between site means, out-of-fold R2) and the within-site corn yield slope, as in Extended Data Fig. 5."),
  `Table 6` = c(
    "Supplementary Table 6: Every matched practice pair in the soil-health data (Fig. 1a).",
    "Two annual systems at the same NAPESHM site that differ in one practice (cover crop, organic inputs, reduced tillage or rotation) and share every other recorded label. Shared columns give the background the pair holds constant; it differs from pair to pair, so each pair measures the practice within one system. Changes: the system with the practice minus the one without (points)."))
# Readable column names for the tables B2 writes in raw form (Supplementary Tables 2 and 3)
ci_txt <- function(x, lo, hi, d = 0) {
  r <- function(v) round(v, d) + 0   # + 0 turns -0 into 0
  x <- r(x); lo <- r(lo); hi <- r(hi)
  ifelse(is.na(x), "", sprintf(paste0("%.", d, "f (%.", d, "f to %.", d, "f)"), x, lo, hi))
}
tidy_supp <- function(k, tab) {
  if (k == 2 && all(c("site", "share_cover", "slope_adj") %in% names(tab))) {
    pct <- function(x) sprintf("%.0f", 100 * x)
    tab <- tab[order(tab$share_div, tab$site), ]     # 10 Nov 2026: Fig. 3c order (Diversity share); groups dropped
    tab <- data.frame(check.names = FALSE,
      Experiment = tab$site, Location = tab$location, Years = paste0(tab$first_year, "\u2013", tab$last_year),
      Treatments = tab$treatments, Plots = tab$plots, `Corn plot-years` = tab$corn_plot_years,
      `N, kg ha-1 (range)` = paste0(tab$N_min, "\u2013", tab$N_max),
      `SHMI range (mean)` = sprintf("%.0f\u2013%.0f (%.0f)", tab$SHMI_min, tab$SHMI_max, tab$SHMI_mean),
      `Share of SHMI spread, %: Cover / Diversity / Disturbance / Amendment` =
        paste(pct(tab$share_cover), pct(tab$share_div), pct(tab$share_dist), pct(tab$share_amend), sep = " / "),
      `Slope per 10 SHMI, kg ha-1 (95% CI), experiment's own N and age` = ci_txt(tab$slope, tab$slope_lower, tab$slope_upper),
      `Mean N, kg ha-1` = tab$N_mean, `Mean system age, years` = tab$age_mean,
      `Slope at 150 kg N ha-1 and the average age (95% CI); Fig. 3c` = ci_txt(tab$slope_adj, tab$slope_adj_lower, tab$slope_adj_upper))
  }
  if (k == 3 && all(c("factor_level", "F_int", "p_int") %in% names(tab))) {
    tab <- data.frame(check.names = FALSE,
      Experiment = tab$site, Level = tab$factor_level, Treatments = tab$n_trt, `Plot-years` = tab$plot_years,
      `SHMI range` = sub("-", "\u2013", tab$SHMI_range),
      `Slope per 10 SHMI, kg ha-1 (95% CI), experiment's own N and age` = ci_txt(tab$slope, tab$lower, tab$upper),
      `Slope at 150 kg N ha-1 and the average age` = ifelse(is.na(tab$slope_adj), "", sprintf("%.0f", tab$slope_adj)),
      `SHMI x factor F` = ifelse(is.na(tab$F_int), "", sprintf("%.2f", tab$F_int)),
      `df` = ifelse(is.na(tab$df1), "", paste0(tab$df1, ", ", round(tab$df2))),
      P = ifelse(is.na(tab$p_int), "", sprintf("%.2f", tab$p_int)),
      Note = ifelse(is.na(tab$note), "", tab$note))
  }
  tab
}
if (requireNamespace("writexl", quietly = TRUE) && requireNamespace("readxl", quietly = TRUE)) {
  yp_f <- O("yield_parts", "yield_parts.csv")
  sheets <- list(); contents <- character()
  for (k in seq_along(SUPP)) {
    nm <- names(SUPP)[k]; stem <- paste0("SuppTable", k)
    f_csv <- file.path(NUM, paste0(stem, ".csv")); f_x <- file.path(NUM, paste0(stem, ".xlsx"))
    contents <- c(contents, SUPP[[k]], "")
    if (file.exists(f_csv)) {
      tab <- read.csv(f_csv, check.names = FALSE, stringsAsFactors = FALSE, fileEncoding = "UTF-8-BOM")
      if (k == 1) {
        # 10 Nov 2026: unconstrained models dropped; one column per soil-health measure (SHMI's
        # own value in its row, differences from SHMI below); corn yield as the slope per
        # within-site SD, its share of SHMI's, and nested likelihood-ratio tests against SHMI.
        sc <- names(tab)[1]
        tab <- tab[!grepl("Linear model|OLS|Random forest", tab[[sc]]), ]
        ref <- grepl("^SHMI \\(weights refitted", tab[[sc]])
        for (m in c("Overall R²", "r between site means", "r within sites")) {
          v <- paste0(m, ": Value"); dcol <- paste0(m, ": Difference from SHMI (95% CI)")
          if (all(c(v, dcol) %in% names(tab))) {
            tab[[dcol]][ref] <- paste0(tab[[v]][ref], " (reference)")
            tab[[v]] <- NULL
            names(tab)[names(tab) == dcol] <- paste0("Soil health: ", m, ", difference from SHMI (95% CI)")
          }
        }
        if (file.exists(yp_f)) {
          yp <- read.csv(yp_f, stringsAsFactors = FALSE)
          if (!"p_adds_to_SHMI" %in% names(yp)) yp$p_adds_to_SHMI <- yp$p_SHMI_adds <- NA_real_
          if (!"pct_of_SHMI" %in% names(yp)) yp$pct_of_SHMI <- 100 * yp$slope_per_sd / yp$slope_per_sd[yp$score == "SHMI"]
          fp <- function(p) ifelse(is.na(p), "", format.pval(p, digits = 2, eps = 0.001))
          yp <- yp %>% transmute(ykey = score,
            `Corn yield: slope per within-site SD of the score, kg ha-1 (95% CI)` =
              sprintf("%.0f (%.0f to %.0f)", slope_per_sd, lower, upper),
            `Corn yield: % of SHMI's slope` = ifelse(score == "SHMI", "100 (reference)", sprintf("%.0f", pct_of_SHMI)),
            `Corn yield: P, score adds to SHMI` = ifelse(score == "SHMI", "", fp(p_adds_to_SHMI)),
            `Corn yield: P, SHMI adds to score` = ifelse(score == "SHMI", "", fp(p_SHMI_adds)))
          # the yield model uses the published weights; the refitted-without-amendment variant
          # scores yield as SHMI_base
          key <- tab[[sc]]
          key[grepl("^SHMI \\(weights refitted|published weights", key)] <- "SHMI"
          key[grepl("without amendment term", key)] <- "SHMI_base (no amendment term)"
          tab <- tab %>% mutate(ykey = key) %>% left_join(yp, by = "ykey") %>% select(-ykey)
          pub <- grepl("published weights", tab[[sc]])
          for (cn in grep("^Corn yield", names(tab), value = TRUE)) tab[[cn]][pub] <- ifelse(cn == names(tab)[grep("^Corn yield", names(tab))[1]], "as SHMI (yield uses the published weights)", "")
        }
      }
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
