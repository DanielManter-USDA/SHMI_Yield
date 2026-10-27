################################################################################
# run_all.R: reruns every analysis in Moore & Manter, in order, and nothing else.
#
#   source("run_all.R")                      # the default run (yield analyses)
#
# or set RUN_SETTINGS first to change any setting below for one run, e.g.
#   RUN_SETTINGS <- list(run_calibration = TRUE, run_n_scenarios = FALSE)
#   source("run_all.R")
#
# The console shows one line per step (start, done or FAILED, minutes). Each step's
# full output goes to its own log in Outputs/logs/<run time>/, with a table of the
# steps in steps.csv. To see the results afterwards, render report.Rmd, which only
# reads Outputs/ and takes seconds:
#   rmarkdown::render("report.Rmd", output_dir = "Outputs")
#
# Inputs that cannot be shipped with the code are listed in Inputs/README.md; the
# check below reports any that are missing for the parts selected.
################################################################################

SETTINGS <- list(
  # Optional rebuilds: show how the shipped workbooks were made from the raw records, and
  # check that they are reproduced. The analyses use the shipped workbooks either way.
  # Part A (optional)
  run_napeshm_rebuild = FALSE,  # A0  NAPESHM workbook from the public NAPESHM tables in Inputs/napeshm (python);
                                #     writes Outputs/NAPESHM_rebuilt.xlsx (before QC) to compare with the
                                #     distributed copy, Inputs/NAPESHM_rebuilt.xlsx (after QC)
  # Part A
  run_calibration     = FALSE,  # A1-A8 soil health model, SHMI calibration, benchmark, SHAPE (hours; frozen for the paper)

  # Part B (optional)
  run_conversion      = FALSE,  # B0  AgCROS workbooks and soil C from the raw AgCROS tables (python). The tables
                                #     are not redistributed here; if you have them, set agcros_tables to where you
                                #     keep them. The paper used the snapshot exported 2026-04-21; a newer export
                                #     can change the workbooks and results.
                                #     Overwrites the files in Inputs/yield_workbooks (git shows any change)
  agcros_tables       = "Inputs/agcros/extracted_tables",   # raw AgCROS tables (run_conversion only)

  # Part B
  run_yield           = TRUE,   # B1-B5, B7-B8 yield analyses, figures and treatment table
  run_weights         = FALSE,  # B5b yield slope with weights fitted between or within sites (reruns B2 three times)
  run_n_scenarios     = TRUE,   # B6  N-accounting sensitivity: 20 legume x manure scenarios (30-60 min)
  run_coarse          = FALSE,  # B9  coarse-record test (hours)

  # Configuration
  in_session          = FALSE,  # TRUE: run each step in this R session instead of a fresh one
                                #   (use it if a step fails with "callr subprocess failed", e.g. on Windows)
  python              = ""      # python executable; "" = RETICULATE_PYTHON (RStudio's Python option), else the PATH
)
# Random steps use fixed seeds in the scripts (cross-validation folds 2026 + repeat,
# bootstraps, random forest seed 1 with one thread, coarse-record seeds), so a rerun
# with the same R and package versions reproduces the paper's numbers.
if (exists("RUN_SETTINGS")) {
  bad <- setdiff(names(RUN_SETTINGS), names(SETTINGS))
  if (length(bad)) stop("Unknown setting(s) in RUN_SETTINGS: ", paste(bad, collapse = ", "))
  SETTINGS <- utils::modifyList(SETTINGS, RUN_SETTINGS)
}
p <- SETTINGS

# ---- setup ---------------------------------------------------------------------
ROOT <- here::here()
setwd(ROOT)
source(file.path(ROOT, "R", "run_helpers.R"))
options(shmi.in_session = isTRUE(p$in_session))
start_log(ROOT)
saveRDS(p, file.path(LOG_DIR, "settings.rds"))
writeLines(capture.output(sessionInfo()), file.path(LOG_DIR, "sessionInfo.txt"))
message("SHMI_Yield run ", basename(LOG_DIR), " | logs: ", LOG_DIR)
message("Settings: ", paste(names(p)[vapply(p, isTRUE, logical(1))], collapse = ", "))

seed_outputs()     # shipped copies of derived inputs (data/), used when the local ones are missing
PY <- find_python(p$python)
PY_WHY <- attr(PY, "why"); PY <- as.character(PY)
py_need <- python_needed(p)
if (is.na(PY) && length(py_need))
  stop("Python is needed to build: ", paste(py_need, collapse = "; "),
       ". Install python with pandas, openpyxl and numpy, set RUN_SETTINGS$python, or put those files in place (see Inputs/README.md).",
       "\n  Pythons tried:\n    ", paste(PY_WHY, collapse = "\n    "))
message("Python: ", if (is.na(PY)) paste0("not available (R-only run): ", paste(PY_WHY, collapse = "; ")) else PY,
        " | SHMI package: ", as.character(utils::packageVersion("SHMI")))
need <- required_inputs(p)
write.csv(need, file.path(LOG_DIR, "inputs_check.csv"), row.names = FALSE)
if (any(!need$present & need$required)) {
  show_inputs(need[need$required, ])
  stop("Required inputs are missing (above). See Inputs/README.md.")
}

# ---- Part A: SHMI calibration (frozen for the paper; off by default) -------------
if (isTRUE(p$run_napeshm_rebuild)) {
  # A0. NAPESHM management records -> SHMI workbook (rules: python/napeshm_rebuild/conversion_rules.md)
  run_py(PY, c("python/napeshm_rebuild/convert.py", "Inputs/napeshm", "Outputs/napeshm_converted.json"),
         "A0_napeshm_convert")
  run_py(PY, c("python/napeshm_rebuild/write_wb.py", "Inputs/napeshm/NAPESHM_test.xlsx",
               "Outputs/napeshm_converted.json", "Outputs/NAPESHM_rebuilt.xlsx"), "A0_napeshm_write")
  message("  Compare Outputs/NAPESHM_rebuilt.xlsx with Inputs/NAPESHM_rebuilt.xlsx (the workbook used, after QC).")
}

if (isTRUE(p$run_calibration)) {
  # A1. soils, SEM data, soil health (lv_SH), climate normals (Extended Data Fig. 1)
  if (file.exists(Sys.getenv("GSSURGO_GPKG", "Inputs/gSSURGO_CONUS.gpkg"))) {
    run_r("A1a_Get_site_ssurgo.R")
  } else message("  gSSURGO not found: using Outputs/site_soils_ssurgo.csv (see Inputs/README.md)")
  run_r("A1b_Prepare_SEM.R")
  run_r("A1c_soil_health_SEM.R")
  if (!file.exists("Outputs/site_climate_monthly.rds")) run_r("A1d_Get_site_climate.R")
  # A2. US model applied to the Mexican sites (a check, not reported), the four-sub-index
  # calibration kept for comparison, summary tables
  for (s in c("A2c_Project_SH.R", "A2d_SHMI_calibration_final.R", "A2e_Summarize_results.R"))
    run_r(s)
  # A3. final calibration (Extended Data Fig. 2; Extended Data Table 1)
  run_r("A3a_SHMI_calibration_plus.R")
  run_r("A3b_calibration_report_figures.R")
  # A4. benchmark (Supplementary Data Table 4) and Fig. 2
  run_r("A4a_SHMI_benchmark.R")
  run_r("A4b_Fig2_soil_health.R")
  # A5. weights by comparison (Supplementary Data Fig. 1)
  run_r("A5a_Run_within_calibration_scenarios.R")
  run_r("A5b_SuppFig1_weights_by_comparison.R")
  # A6. SHAPE (Extended Data Fig. 3); the SHAPE app runs outside R
  run_r("A6a_Prepare_SHAPE.R")
  dir.create("Outputs/SHAPE", showWarnings = FALSE, recursive = TRUE)
  if (file.exists("Inputs/SHAPE/batch_Combined_scores.csv")) {
    file.copy("Inputs/SHAPE/batch_Combined_scores.csv", "Outputs/SHAPE/batch_Combined_scores.csv", overwrite = TRUE)
    run_r("A6b_Merge_SHAPE_SHMI.R")
    run_r("A6c_SHMI_SHAPE_age.R")
  } else message("  SHAPE: run Outputs/SHAPE/SHAPE_multi_indicator_input.csv through the SHAPE app, put ",
                 "batch_Combined_scores.csv in Inputs/SHAPE/, then rerun.")
  # A7-A8. climate source check (Methods); SEM loadings and adjustment slopes (Extended Data Fig. 1)
  run_r("A7_Check_climate_source.R")
  run_r("A8_Export_appendix_tables.R")
}

# ---- Part B: corn yield --------------------------------------------------------
# B0. AgCROS workbooks and soil C, converted from the raw AgCROS tables only with
# run_conversion = TRUE (both are shipped); NAPESHM crosswalk for the mediation test
TEMPLATE <- system.file("extdata", "SHMI_template.xlsx", package = "SHMI")
if (!nzchar(TEMPLATE)) stop("SHMI_template.xlsx not found in the installed SHMI package (extdata)")
WB <- "Inputs/yield_workbooks"
dir.create(WB, showWarnings = FALSE, recursive = TRUE)
conv <- function(sites, extra = character(), label)
  run_py(PY, c("python/agcros2shmi/agcros2shmi.py", p$agcros_tables, sites,
               "--template", TEMPLATE, "--out", WB, extra), label)
if (conv_needed(p)) {
  conv(c("NELITCSE", "NEMLTCRS", "NEMEIRR", "NECCIRR"),
       c("--supplement", "python/agcros2shmi/supplements/NECCIRR_amendments.csv",
         "--corrections", "python/agcros2shmi/corrections/NEMLTCRS_corrections.csv"), "B0_convert_nebraska")
  conv(c("INWLACRE", "IAAM7071", "PAUP"), label = "B0_convert_added")
}
dir.create("Outputs/yield", showWarnings = FALSE, recursive = TRUE)
if (conv_needed(p) || !file.exists("Outputs/yield/soil_c.csv"))
  run_py(PY, c("python/build_soil_c.py", p$agcros_tables, "Outputs/yield"), "B0_soil_c")
run_r("B0_build_napeshm_crosswalk.R")

if (isTRUE(p$run_yield)) {
  run_r("B1_Get_SPEI.R")                    # growing-season SPEI for every site-year
  run_r("B2_SHMI_yield_main.R")             # Fig. 3, Extended Data Fig. 4, Supp. Data Table 1, Supp. Data Fig. 2
  run_r("B3_SHMI_yield_stability.R")        # Fig. 4
  run_r("B4_SHMI_N_equivalent.R")           # fertilizer-N equivalent (reruns B2 in its own session)
  run_r("B5_SHMI_soil_mediation.R")         # Supp. Data Table 3 (reruns B2 in its own session)
}
if (isTRUE(p$run_weights)) run_r("B5b_SHMI_yield_weights_sensitivity.R")   # Methods: weights by comparison
if (isTRUE(p$run_n_scenarios)) run_r("B6_SHMI_N_scenarios.R")               # Methods: N accounting
if (isTRUE(p$run_yield)) {
  run_r("B7_Assemble_manuscript_figures.R") # Figs. 3-4, Extended Data Fig. 4
  if (!is.na(PY)) {
    run_py(PY, c("python/build_study_treatments.py", "--wb-dir", "Inputs/yield_workbooks",
                 "--units-dir", "Inputs/yield_workbooks", "--old", "data/SHMI_study_treatments_descriptions.xlsx",
                 "--out", "Outputs/SHMI_study_treatments.xlsx"), "B8_study_treatments")
  } else message("  No python: treatment table (B8) not rebuilt")
}
if (isTRUE(p$run_coarse)) {
  # B9. coarse records (Extended Data Fig. 5; Supp. Data Table 2); python writes the coarse
  # workbooks only when they are not already in Outputs/Coarse/
  run_r("B9_SHMI_coarse_records.R", vars = list(CFG_COARSE_OVERRIDE = list(python = if (is.na(PY)) "" else PY,
                                                                        reuse_runs = FALSE)))
}

# ---- Part C: display items into Outputs/Manuscript (seconds) ----------------------
run_r("C1_Fig1_and_Table1.R")       # Fig. 1, Table 1
run_r("C1b_ED_Fig2.R")              # Extended Data Fig. 2, from the saved calibration (Outputs/Final_plus)
# When Part A was not run, its drawing steps (seconds each) are rerun from the saved
# Part A outputs, so every display item is current; a step whose input is missing is skipped
if (!isTRUE(p$run_calibration)) {
  REDRAW <- c("A4b_Fig2_soil_health.R"                = "Outputs/Final_plus/us_oof_predictions.csv",          # Fig. 2
              "A6c_SHMI_SHAPE_age.R"                  = "Outputs/SHAPE/SHAPE_by_unit_with_SHMI.csv",          # Extended Data Fig. 3
              "A5b_SuppFig1_weights_by_comparison.R"  = "Outputs/Calib_scenarios/scenario_summary.csv",       # Supplementary Data Fig. 1
              "A8_Export_appendix_tables.R"           = "Outputs/sem_noMX_projection.rds")                    # Extended Data Fig. 1 loadings
  if (!file.exists("Outputs/Benchmark/benchmark_summary.csv")) REDRAW <- REDRAW[names(REDRAW) != "A4b_Fig2_soil_health.R"]
  for (s in names(REDRAW)) {
    if (file.exists(REDRAW[[s]])) run_r(s) else message("  ", s, " skipped: ", REDRAW[[s]], " not found (run Part A)")
  }
}
run_r("C2_export_tables.R")         # Extended Data Table 1, Supplementary Data Tables 1-4, Tables.xlsx
run_r("C3_collect_figures.R")       # the other figures, under the paper's numbering; manifest.csv

# ---- checks ----------------------------------------------------------------------
# a rerun of Part A must reproduce the weights frozen in the SHMI package
fr <- "Outputs/Final_plus/final_results.rds"
if (isTRUE(p$run_calibration) && file.exists(fr)) {
  wf <- readRDS(fr)$w_final
  ref <- c(cover = SHMI::shmi_weights()[["Cover"]], diversity = SHMI::shmi_weights()[["Diversity"]],
           disturbance = SHMI::shmi_weights()[["InvDist"]], growing_share = SHMI::shmi_parameters()$growing_share)
  got <- unlist(wf[names(ref)])
  message("Calibrated weights ", paste(sprintf("%s %.3f", names(got), got), collapse = ", "),
          if (isTRUE(all.equal(unname(round(got, 3)), unname(ref), tolerance = 1e-3))) " match" else " DIFFER FROM",
          " the SHMI package (", paste(sprintf("%.3f", ref), collapse = ", "), ")")
}
mf <- "Outputs/Manuscript/manifest.csv"
if (file.exists(mf)) { m <- utils::read.csv(mf); message("Display items in Outputs/Manuscript: ", sum(m$present), " of ", nrow(m)) }

# ---- summary ---------------------------------------------------------------------
steps <- utils::read.csv(file.path(LOG_DIR, "steps.csv"))
message("\nAll steps finished (", round(sum(steps$minutes), 1), " min). Logs: ", LOG_DIR)
print(steps[, c("step", "ok", "minutes")], row.names = FALSE)
message("Results: rmarkdown::render(\"report.Rmd\", output_dir = \"Outputs\")")
