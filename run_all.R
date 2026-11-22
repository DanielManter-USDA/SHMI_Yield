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
  run_calibration     = FALSE,  # A1-A8 soil health model, SHMI calibration, benchmark, fixed-value grid, SHAPE (hours; frozen for the paper)

  # Part B (optional)
  run_conversion      = FALSE,  # B0  AgCROS workbooks from the raw AgCROS tables (python). The tables
                                #     are not redistributed here; if you have them, set agcros_tables to where you
                                #     keep them. The paper used the snapshot exported 2026-04-21; a newer export
                                #     can change the workbooks and results.
                                #     Overwrites the files in Inputs/yield_workbooks (git shows any change)
  agcros_tables       = "Inputs/agcros",   # raw AgCROS tables (run_conversion only)

  # Part B
  run_yield           = TRUE,   # B1-B4, B8 yield analyses, figures and treatment table
  run_n_scenarios     = TRUE,   # B6  N-accounting sensitivity: 20 legume x manure scenarios (30-60 min)
  run_record_length   = TRUE,   # B2c rolling (current-season) scores, used by B2d
  run_yield_parts     = TRUE,   # B2d yield: SHMI against its parts, Table 2 (20-40 min)

  # Configuration
  in_session          = TRUE,  # TRUE: run each step in this R session instead of a fresh one
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
  # A2. the US soil health model applied to the Mexican sites (read by A3a; not in the paper)
  run_r("A2c_Project_SH.R")
  # A3. final calibration (weights, growing-season share; out-of-fold predictions)
  run_r("A3a_SHMI_calibration_plus.R")
  # A4. benchmark: out-of-fold R2 with the weights refitted in each fold (Fig. 1b, Results)
  run_r("A4a_SHMI_benchmark.R")
  # A5. fixed-value grid: Diversity weight, livestock and amendment shares (Methods, Supplementary Methods)
  run_r("A5a_Run_within_calibration_scenarios.R")
  # A6. SHAPE (Table 2); the SHAPE app runs outside R
  run_r("A6a_Prepare_SHAPE.R")
  dir.create("Outputs/SHAPE", showWarnings = FALSE, recursive = TRUE)
  if (file.exists("Inputs/SHAPE/batch_Combined_scores.csv")) {
    file.copy("Inputs/SHAPE/batch_Combined_scores.csv", "Outputs/SHAPE/batch_Combined_scores.csv", overwrite = TRUE)
    run_r("A6b_Merge_SHAPE_SHMI.R")
  } else message("  SHAPE: run Outputs/SHAPE/SHAPE_multi_indicator_input.csv through the SHAPE app, put ",
                 "batch_Combined_scores.csv in Inputs/SHAPE/, then rerun.")
  # A7-A8. climate source check (Supplementary Methods); SEM loadings for Fig. 1a
  run_r("A7_Check_climate_source.R")
  run_r("A8_Export_appendix_tables.R")
}

# ---- Part B: corn yield --------------------------------------------------------
# B0. AgCROS workbooks, converted from the raw AgCROS tables only with
# run_conversion = TRUE (the workbooks are shipped)
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
  conv(c("INWLACRE", "IAAM7071", "PAUP"),
       c("--supplement", "python/agcros2shmi/supplements/IAAM7071_covers.csv"), "B0_convert_added")   # IAAM7071 rye and clover (not in AgCROS)
}
dir.create("Outputs/yield", showWarnings = FALSE, recursive = TRUE)

if (isTRUE(p$run_yield)) {
  run_r("B1_Get_SPEI.R")                    # growing-season SPEI for every site-year
  run_r("B2_SHMI_yield_main.R")             # Figs 3 and 4, Supplementary Table 2
  run_r("B2b_SHMI_soybean.R")               # soybean (Fig. 3a)
  if (isTRUE(p$run_record_length)) run_r("B2c_SHMI_record_length.R")   # rolling (current-season) scores for B2d
  if (isTRUE(p$run_yield_parts)) run_r("B2d_SHMI_yield_parts.R")       # Table 2, corn yield (after B2c)
  run_r("B2i_SHMI_variance_shares.R")       # share of within-site SHMI variance from each part (Fig. 3b, Supplementary Table 2)
  run_r("B3_SHMI_yield_stability.R")        # Fig. 4c
  run_r("B4_SHMI_N_equivalent.R")           # fertilizer-N equivalent (reads B2_model.rds)
}
if (isTRUE(p$run_n_scenarios)) run_r("B6_SHMI_N_scenarios.R")               # Supplementary Table 6
if (isTRUE(p$run_yield)) {
  if (!is.na(PY)) {
    run_py(PY, c("python/build_study_treatments.py", "--wb-dir", "Inputs/yield_workbooks",
                 "--units-dir", "Inputs/yield_workbooks",
                 "--old", "data/SHMI_study_treatments_descriptions.xlsx",
                 "--out", "Outputs/SHMI_study_treatments.xlsx"), "B8_study_treatments")
  } else message("  No python: treatment table (B8) not rebuilt")
  run_r("B8b_study_treatments_with_yield.R")   # Supplementary Table 3: every treatment, with SHMI and mean yields
}

# ---- Part C: soil-health analyses, figures and tables (minutes) -----------------
run_r("C1_Fig1_path.R")             # Fig. 1b, the path diagram
run_r("C1f_paired_data_budget.R")   # share of systems usable by strict single-practice pairs (Discussion)
run_r("C1h_SHMI_vs_soil_health.R")  # Fig. 2; SHAPE
run_r("C1j_SHMI_vs_parts.R")        # Table 2, soil-health columns
run_r("C1k_ceiling.R")              # ceiling on explainable within-site variation (Results)
run_r("C1l_rf_within.R")            # within-site random forests (Results, Supplementary Methods)
# When Part A was not run, A8 (seconds) is rerun from the saved Part A outputs for Fig. 1a
if (!isTRUE(p$run_calibration) && file.exists("Outputs/sem_noMX_projection.rds")) run_r("A8_Export_appendix_tables.R")
run_r("C1m_supp_tables.R")          # Supplementary Tables 1, 4, 5 and 6
source(here::here("R", "figures", "run_figures.R"))   # Figs 1-4 from the saved plot data and panels
run_r("C2_export_tables.R")         # Tables 1-2, Supplementary Table 2
run_r("C3_collect_figures.R")       # every display item under the paper's numbering; manifest_numbered.csv
run_r("C4_bundle_for_review.R")     # one zip of every small results file, for checking the text

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
mf <- "Outputs/Manuscript/numbered/manifest_numbered.csv"
if (file.exists(mf)) { m <- utils::read.csv(mf); message("Display items in Outputs/Manuscript/numbered: ", sum(m$present), " of ", nrow(m)) }

# ---- summary ---------------------------------------------------------------------
steps <- utils::read.csv(file.path(LOG_DIR, "steps.csv"))
message("\nAll steps finished (", round(sum(steps$minutes), 1), " min). Logs: ", LOG_DIR)
print(steps[, c("step", "ok", "minutes")], row.names = FALSE)
message("Results: rmarkdown::render(\"report.Rmd\", output_dir = \"Outputs\")")
