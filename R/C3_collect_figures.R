################################################################################
# C3_collect_figures.R
#
# 17 Nov 2026: Figs 1-5, Tables 1-3 (Table 4 is text), Supplementary Tables 1-3; no Supplementary Figs.
# 16 Nov 2026: Figs 1-5 and Supplementary Figs 1-3 and Tables 1-6, numbered by first citation (see NUMBERED below).
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
  # 12 Nov 2026 (Nature Communications): six main figures, one per Results section, drawn by
  # R/figures/F1-F6 from saved plot data and panels; the former Extended Data Figs 1-5 are now
  # panels of Figs 1-3 and 6. Extended Data Table 1 keeps its file name until the
  # Supplementary Information is renumbered.
  # 16 Nov 2026: five main figures; the soil-health model (F1) and the principles by outcome
  # (F6) move to the Supplementary Information. Supplementary Figs and Tables are numbered by
  # first citation in the main text (tables: old 5 -> 3, old 6 -> 5, old 3 -> 6).
  # 17 Nov 2026 (streamlined layout): Fig. 1 a SEM, b path diagram; Fig. 2 SHMI vs lv_SH (all,
  # between, within); Fig. 3 corn and soybean, each experiment; Fig. 4 N, age, good and bad years;
  # Fig. 5 simplified records. Table 2 SHMI against its parts (C1j + B2d), Table 3 practice
  # effects from strict pairs (C1i + C1f). No Supplementary Figures. Supplementary Tables:
  # 1 the ten experiments, 2 every treatment, 3 coarse-record scenarios.
  c("Fig1.png",                  O("figures", "F1_soil_health_model.png")),  # SEM | path diagram
  c("Fig2.png",                  O("figures", "F2_system_score.png")),       # SHMI vs soil health: all, between, within
  c("Fig3.png",                  O("figures", "F3_yield.png")),              # corn and soybean | each experiment
  c("Fig4.png",                  O("figures", "F4_yield_moderators.png")),   # N | age | good and bad years
  c("Fig5.png",                  O("figures", "F5_records.png")),            # simplified records
  # Table 1 = scores by system type (C2 still writes it as ED_Table1.csv); C1's example-systems
  # Table1.csv left the paper on 17 Nov 2026 and is no longer copied.
  c("Table1.csv",                file.path(WORK, "ED_Table1.csv")),          # scores by system type
  c("Table2.csv",                file.path(WORK, "Table2.csv")),             # SHMI against its parts
  c("Table3.csv",                file.path(WORK, "Table3.csv")),             # practice effects, strict pairs
  # Supplementary Tables 1-3, staged here, then combined into one workbook below.
  c("SuppTable1.csv",            file.path(WORK, "SuppTable1.csv")),                   # the ten experiments (B2)
  c("SuppTable2.xlsx",           O("SHMI_study_treatments_with_yield.xlsx")),         # every treatment (B8, B8b)
  c("SuppTable3.csv",            file.path(WORK, "SuppTable2.csv"))))                  # coarse-record scenarios (B9)
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
    "Supplementary Table 1: The ten experiments and what drove their SHMI range.",
    "Corn years, treatments and plots are those analysed (treatments with corn; rotation phases combined). N: fertilizer plus manure (total N) plus legume credit.",
    "Share of SHMI spread: the share of each experiment's within-site SHMI variance from each part of SHMI (weight x Cover, Diversity and Inverse Disturbance, and the amendment bonus; cov(part, SHMI)/var(SHMI)); the four add to 100%, and a negative share means that part moved against SHMI. Experiments are ordered by their slope at 150 kg N ha-1, as in Fig. 3b. Slope: each experiment's own within-site model, averaged over its N rates and years."),
  `Table 2` = c(
    "Supplementary Table 2: Every treatment in the ten experiments.",
    "Description, SHMI and its parts (points from each sub-index = weight x sub-index), N, and mean corn and soybean yields, with notes on data corrections (IAAM7071 cover crops and harvest dates; INWLACRE 2004 planting dates). Further sheets of the treatment workbook follow as 'Table 2 ...'."),
  `Table 3` = c(
    "Supplementary Table 3: Every coarse-record scenario and its results.",
    "For each scenario, what was degraded, soil-health skill on the calibration's folds (correlation within sites, share of sites ranked in the right direction, correlation between site means, out-of-fold R2) and the within-site corn yield slope, as in Fig. 5."))
# Readable column names for the tables B2 writes in raw form (Supplementary Tables 2 and 6)
ci_txt <- function(x, lo, hi, d = 0) {
  r <- function(v) round(v, d) + 0   # + 0 turns -0 into 0
  x <- r(x); lo <- r(lo); hi <- r(hi)
  ifelse(is.na(x), "", sprintf(paste0("%.", d, "f (%.", d, "f to %.", d, "f)"), x, lo, hi))
}
tidy_supp <- function(k, tab) {
  if (k == 1 && all(c("site", "share_cover", "slope_adj") %in% names(tab))) {   # 17 Nov 2026: was Table 2
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
