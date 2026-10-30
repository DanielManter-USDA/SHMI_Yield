################################################################################
# C3_collect_figures.R
#
# Copies every figure in the paper from the folder of the step that makes it into
# Outputs/Manuscript, under the paper's numbering, and writes a manifest of display
# items (Outputs/Manuscript/manifest.csv: item, file, source script, present or not).
# Fig. 1 (C1), Fig. 2 (A4b) and Extended Data Figs. 2 and 5 (C1b, C1c) are written to Outputs/Manuscript directly; tables
# are written there by C2.
#
# Extended Data Fig. 1 (the soil health model diagram) is drawn by hand from the
# SEM loadings exported by A8 (Outputs/sem_noMX_standardized.csv); it has no plotting
# script, so the manifest lists the loadings file as its source.
################################################################################
OUT  <- here::here("Outputs", "Manuscript")
O    <- function(...) here::here("Outputs", ...)
MAIN <- "yield_main_default_total_agcros"
FIGS <- "manuscript_figures_default_total_agcros"
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

ITEMS <- data.frame(stringsAsFactors = FALSE, rbind(
  #  item                          file in Outputs/Manuscript    source file                                  script
  c("Fig. 1",                      "Fig1.png",                   file.path(OUT, "Fig1.png"),                  "C1_Fig1_and_Table1.R"),
  c("Fig. 2",                      "Fig2.png",                   file.path(OUT, "Fig2.png"),                  "A4b_Fig2_soil_health.R"),
  c("Fig. 3",                      "Fig3.png",                   O(FIGS, "Fig3.png"),                         "B7_Assemble_manuscript_figures.R"),
  c("Fig. 4",                      "Fig4.png",                   O(FIGS, "Fig4.png"),                         "B7_Assemble_manuscript_figures.R"),
  c("Extended Data Fig. 1",        "ED_Fig1_loadings.csv",       O("sem_noMX_standardized.csv"),              "A8_Export_appendix_tables.R (diagram drawn by hand)"),
  c("Extended Data Fig. 2",        "ED_Fig2.png",                file.path(OUT, "ED_Fig2.png"),               "C1b_ED_Fig2.R (from A3a's final_results.rds)"),
  c("Extended Data Fig. 3",        "ED_Fig3.png",                O("SHAPE_age", "FigA_SHMI_vs_SHAPE.png"),    "A6c_SHMI_SHAPE_age.R"),
  c("Extended Data Fig. 4",        "ED_Fig4.png",                O(FIGS, "ED_Fig4_sites.png"),                "B7_Assemble_manuscript_figures.R"),
  c("Extended Data Fig. 5",        "ED_Fig5.png",                file.path(OUT, "ED_Fig5.png"),               "C1c_ED_Fig5.R (from B9's saved results)"),
  c("Extended Data Fig. 6",        "ED_Fig6.png",                O("Practice_additivity", "Fig_labels_vs_SHMI.png"), "C1e_labels_vs_SHMI.R"),
  c("Supplementary Data Fig. 1",   "SuppFig1.png",               O("yield_soybean", "Fig_soybean.png"),       "B2b_SHMI_soybean.R"),
  c("Supplementary Data Fig. 2",   "SuppFig2.png",               O("Calib_scenarios", "SuppFig2_weights_by_comparison.png"), "A5b_SuppFig2_weights_by_comparison.R"),
  c("Supplementary Data Fig. 3",   "SuppFig3.png",               O(MAIN, "Fig_effect_sizes.png"),             "B2_SHMI_yield_main.R"),
  c("Table 1",                     "Table1.csv",                 file.path(OUT, "Table1.csv"),                "C1_Fig1_and_Table1.R"),
  c("Extended Data Table 1",       "ED_Table1.csv",              file.path(OUT, "ED_Table1.csv"),             "C2_export_tables.R"),
  c("Supplementary Data Table 1",  "SuppTable1.csv",             file.path(OUT, "SuppTable1.csv"),            "C2_export_tables.R (from B2)"),
  c("Supplementary Data Table 2",  "SuppTable2.csv",             file.path(OUT, "SuppTable2.csv"),            "C2_export_tables.R (from B9)"),
  c("Supplementary Data Table 3",  "SuppTable3.csv",             file.path(OUT, "SuppTable3.csv"),            "C2_export_tables.R (from B5)"),
  c("Supplementary Data Table 4",  "SuppTable4.csv",             file.path(OUT, "SuppTable4.csv"),            "C2_export_tables.R (from A4a)")))
names(ITEMS) <- c("item", "file", "source", "script")

for (i in seq_len(nrow(ITEMS))) {
  src <- ITEMS$source[i]; dst <- file.path(OUT, ITEMS$file[i])
  if (file.exists(src) && normalizePath(src, mustWork = FALSE) != normalizePath(dst, mustWork = FALSE))
    file.copy(src, dst, overwrite = TRUE, copy.date = TRUE)
}
ITEMS$present  <- file.exists(file.path(OUT, ITEMS$file))
ITEMS$modified <- ifelse(ITEMS$present, format(file.mtime(file.path(OUT, ITEMS$file)), "%Y-%m-%d %H:%M"), NA)
ITEMS$source   <- sub(paste0("^", here::here(), "/?"), "", ITEMS$source)
write.csv(ITEMS, file.path(OUT, "manifest.csv"), row.names = FALSE)
message("Display items in Outputs/Manuscript: ", sum(ITEMS$present), " of ", nrow(ITEMS))
if (any(!ITEMS$present))
  message("  missing (their step was not run): ", paste(ITEMS$item[!ITEMS$present], collapse = "; "))
