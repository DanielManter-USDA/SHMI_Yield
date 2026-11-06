################################################################################
# C3_collect_figures.R
#
# Copies every display item from the folder of the step that makes it into
# Outputs/Manuscript/numbered, under its number in the paper (manuscript v5,
# 6 Nov 2026: items renumbered in order of first citation), with a manifest
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
NUMBERED <- data.frame(stringsAsFactors = FALSE, rbind(
  c("Fig1.png",                  file.path(WORK, "Fig1.png")),
  c("Fig2.png",                  file.path(WORK, "Fig2.png")),
  c("Fig3.png",                  O(FIGS, "Fig3.png")),
  c("Fig4.png",                  O(FIGS, "Fig4.png")),
  c("ED_Fig1_loadings.csv",      O("sem_noMX_standardized.csv")),                    # loadings; the diagram is drawn by hand
  c("ED_Fig2.png",               file.path(WORK, "ED_Fig2.png")),                      # calibrated weights
  c("ED_Fig3.png",               O("Practice_additivity", "Fig_labels_vs_SHMI.png")), # checklist vs SHMI (C1e)
  c("ED_Fig4.png",               O("SHAPE_age", "FigA_SHMI_vs_SHAPE.png")),           # SHMI and SHAPE
  c("ED_Fig5.png",               O(FIGS, "ED_Fig4_sites.png")),                       # experiments, own slopes, parts
  c("ED_Fig6.png",               file.path(WORK, "ED_Fig5.png")),                      # coarse records (C1c)
  c("SuppFig1.png",              O("yield_soybean", "Fig_soybean.png")),              # soybean
  c("SuppFig2.png",              O("yield_paired_budget", "Fig_yield_budget.png")),   # matched pairs, yield (B2f)
  c("SuppFig3.png",              O("Calib_scenarios", "SuppFig2_weights_by_comparison.png")), # weights by comparison
  c("SuppFig4.png",              O(MAIN, "Fig_effect_sizes.png")),                    # effect sizes
  c("Table1.csv",                file.path(WORK, "Table1.csv")),
  c("ED_Table1.csv",             file.path(WORK, "ED_Table1.csv")),
  c("SuppTable1.csv",            file.path(WORK, "SuppTable4.csv")),                   # SHMI vs variants (A4a)
  c("SuppTable2.csv",            file.path(WORK, "SuppTable3.csv")),                   # soil properties as carriers (B5)
  c("SuppTable3.csv",            file.path(WORK, "SuppTable1.csv")),                   # the ten experiments (B2)
  c("SuppTable4.csv",            file.path(WORK, "SuppTable2.csv"))))                  # coarse-record scenarios (B9)
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
write.csv(NUMBERED, file.path(NUM, "manifest_numbered.csv"), row.names = FALSE)
message("Numbered as in the paper (Outputs/Manuscript/numbered): ", sum(NUMBERED$present), " of ", nrow(NUMBERED))
if (any(!NUMBERED$present)) message("  missing: ", paste(NUMBERED$file[!NUMBERED$present], collapse = "; "))
