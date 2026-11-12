################################################################################
# _targets.R  (11 Nov 2026)   The SHMI_Yield pipeline
#
#   targets::tar_make()          run every step that is out of date (and only those)
#   targets::tar_visnetwork()    the pipeline as a graph: what depends on what, what is stale
#   targets::tar_outdated()      the steps tar_make() would run now
#   targets::tar_read(B2_main)   the files a step wrote (or any target's value)
#   targets::tar_invalidate(B3)  force a step to rerun next time
#
# Once per project: targets::tar_config_set(reporter_make = "timestamp")
#   (plain timestamped lines instead of the live progress bar, which the scripts' own
#   progress messages can clear when they run inside the pipeline's R process)
#
# Each step runs one script from R/ (its log goes to Outputs/logs/<run>/) and returns the
# files later steps read. A step reruns when its script, its input files or an upstream
# step's files change. Changing a figure script in R/figures/ redraws only that figure.
#
# Part A (the calibration) is frozen: its outputs are tracked as input files. B9 (the
# coarse-record test, hours) reruns only when forced: Sys.setenv(SHMI_FORCE = "B9").
# Retired (not run): B7 (Figs 3 and 4 now come from R/figures), B0 (crosswalk for B5 only), B2e, B2g, B2h, B5, C1b, A5b and A6c redraws.
# Each step lists only the upstream steps whose files its script actually reads.
################################################################################
library(targets)
source("R/targets_helpers.R")
# run each script inside the pipeline's own R process (tar_make() already runs apart from
# RStudio), so B2 does not need a second R session; FALSE runs each step in its own process
options(shmi.in_session = TRUE)
tar_option_set(error = "stop")

MAIN <- "Outputs/yield_main_default_total_agcros"
STAB <- "Outputs/yield_stability_default_total_agcros"
PA   <- "Outputs/Practice_additivity"
WORK <- "Outputs/Manuscript_work"
sc   <- function(f) file.path("R", f)

list(
  # ---- the SHMI package: reinstalling it reruns everything that scores records ----------
  # (tracks the installed package's DESCRIPTION file, which changes with every install, so
  # tar_outdated() stays accurate; B2's site cache is keyed on the install time in the same way)
  tar_target(shmi_version, system.file("DESCRIPTION", package = "SHMI", mustWork = TRUE), format = "file"),

  # ---- inputs ------------------------------------------------------------------------------
  tar_target(in_workbooks, list.files("Inputs/yield_workbooks", full.names = TRUE), format = "file"),
  tar_target(in_napeshm, c("data/SHAPE_by_unit_with_SHMI.csv", "Inputs/napeshm/treatments.csv"), format = "file"),
  tar_target(in_spei, "Inputs/climate/spei03.nc", format = "file"),
  tar_target(in_cache, present("data/cache/nass_site_county_corn_yield.csv", "data/cache/site_soils_ssurgo.csv",
                               "data/napeshm_yield_site_units.csv", "data/soil_c.csv"), format = "file"),
  tar_target(in_partA, present("Outputs/Final_plus/us_oof_predictions.csv", "Outputs/Final_plus/final_results.rds",
                               "Outputs/Benchmark/benchmark_summary.csv", "Outputs/Benchmark/benchmark_bootstrap.csv",
                               "Outputs/sem_noMX_standardized.csv", "Outputs/site_climate_monthly.rds",
                               "Outputs/SHAPE/SHAPE_by_unit_with_SHMI.csv"), format = "file"),
  tar_target(in_helpers, sc(c("SHMI_legume_credit.R", "SHMI_manure_n.R", "run_helpers.R")), format = "file"),

  # ---- scripts (an edited script reruns its step) -------------------------------------------
  tar_target(s_B1,  sc("B1_Get_SPEI.R"), format = "file"),
  tar_target(s_B2,  sc("B2_SHMI_yield_main.R"), format = "file"),
  tar_target(s_B2b, sc("B2b_SHMI_soybean.R"), format = "file"),
  tar_target(s_B2c, sc("B2c_SHMI_record_length.R"), format = "file"),
  tar_target(s_B2d, sc("B2d_SHMI_yield_parts.R"), format = "file"),
  tar_target(s_B2f, sc("B2f_SHMI_yield_paired_budget.R"), format = "file"),
  tar_target(s_B3,  sc("B3_SHMI_yield_stability.R"), format = "file"),
  tar_target(s_B4,  sc("B4_SHMI_N_equivalent.R"), format = "file"),
  tar_target(s_B5b, sc("B5b_SHMI_yield_weights_sensitivity.R"), format = "file"),
  tar_target(s_B6,  sc("B6_SHMI_N_scenarios.R"), format = "file"),
  tar_target(s_B8,  "python/build_study_treatments.py", format = "file"),
  tar_target(s_B8b, sc("B8b_study_treatments_with_yield.R"), format = "file"),
  tar_target(s_B9,  c(sc("B9_SHMI_coarse_records.R"), "python/coarse_records/degrade_records.py"), format = "file"),
  tar_target(s_A4b, sc("A4b_Fig2_soil_health.R"), format = "file"),
  tar_target(s_C1d, sc("C1d_practice_additivity.R"), format = "file"),
  tar_target(s_C1e, sc("C1e_labels_vs_SHMI.R"), format = "file"),
  tar_target(s_C1f, sc("C1f_paired_data_budget.R"), format = "file"),
  tar_target(s_C1g, sc("C1g_label_vs_SHMI.R"), format = "file"),
  tar_target(s_C1h, sc("C1h_SHMI_vs_soil_health.R"), format = "file"),
  tar_target(s_C1,  sc("C1_Fig1_and_Table1.R"), format = "file"),
  tar_target(s_C1c, sc("C1c_ED_Fig5.R"), format = "file"),
  tar_target(s_C2,  sc("C2_export_tables.R"), format = "file"),
  tar_target(s_C3,  sc("C3_collect_figures.R"), format = "file"),
  tar_target(s_fig, sc(c("figures/fig_common.R", "figures/F3_yield.R", "figures/F4_yield_moderators.R")), format = "file"),

  # ---- yield analyses (Part B) -------------------------------------------------------------
  tar_target(B1_spei, step("B1_Get_SPEI.R", c("data/cache/site_spei.csv", "Outputs/yield/drought_years.csv"),
                           s_B1, in_workbooks, in_spei), format = "file"),
  tar_target(B2_main, step("B2_SHMI_yield_main.R",
                           file.path(MAIN, c("B2_data.rds", "B2_model.rds", "plotdata_yield_main.rds", "results_main.csv",
                                             "site_SHMI_ranges.csv", "SHMI_x_other_by_site.csv")),
                           s_B2, in_helpers, in_workbooks, in_cache, shmi_version), format = "file"),
  tar_target(B2b_soybean, step("B2b_SHMI_soybean.R",
                               file.path("Outputs/yield_soybean", c("plotdata_soybean.rds", "soy_results.csv", "soy_treatment_yields.csv")),
                               s_B2b, B2_main), format = "file"),
  tar_target(B2c_record, step("B2c_SHMI_record_length.R",
                              "Outputs/yield_record_length/rolling_shmi.rds",
                              s_B2c, B2_main), format = "file"),
  tar_target(B2d_parts, step("B2d_SHMI_yield_parts.R", file.path("Outputs/yield_parts", c("yield_parts.csv", "yield_parts_loo.csv")),
                             s_B2d, B2_main, B2c_record), format = "file"),
  tar_target(B2f_budget, step("B2f_SHMI_yield_paired_budget.R",
                              file.path("Outputs/yield_paired_budget", c("yield_budget.csv", "treatment_labels.csv")),
                              s_B2f, B2_main), format = "file"),
  tar_target(B3_stability, step("B3_SHMI_yield_stability.R",
                                file.path(STAB, c("plotdata_stability.rds", "Fig2c_years.csv", "insurance_by_cause.csv")),
                                s_B3, s_B2, B2_main, B1_spei), format = "file"),
  tar_target(B4_Neq, step("B4_SHMI_N_equivalent.R", file.path(MAIN, c("N_equivalent.csv", "N_equivalent_by_site.csv")),
                          s_B4, B2_main), format = "file"),
  tar_target(B5b_weights, step("B5b_SHMI_yield_weights_sensitivity.R", "Outputs/yield_weights_agcros/weights_slopes.csv",
                               s_B5b, s_B2, B2_main), format = "file"),
  tar_target(B6_Nscen, step("B6_SHMI_N_scenarios.R", file.path(MAIN, c("N_scenarios_all_sites.csv", "N_scenarios_WICST.csv")),
                            s_B6, s_B2, B2_main), format = "file"),
  tar_target(B8_treatments, step_py(c("python/build_study_treatments.py", "--wb-dir", "Inputs/yield_workbooks",
                                      "--units-dir", "Inputs/yield_workbooks", "--old", "data/SHMI_study_treatments_descriptions.xlsx",
                                      "--out", "Outputs/SHMI_study_treatments.xlsx"),
                                    "B8_study_treatments", "Outputs/SHMI_study_treatments.xlsx", s_B8, in_workbooks), format = "file"),
  tar_target(B8b_treatments, step("B8b_study_treatments_with_yield.R", "Outputs/SHMI_study_treatments_with_yield.xlsx",
                                  s_B8b, B8_treatments, B2_main, B2b_soybean, B2f_budget), format = "file"),
  tar_target(B9_coarse, step_slow("B9_SHMI_coarse_records.R",
                                  file.path("Outputs/Coarse", c("coarse_napeshm.csv", "coarse_yield.csv", "coarse_scenarios.csv")),
                                  s_B9, s_B2, B2_main, in_partA, key = "B9",
                                  needs = c("Outputs/Final_plus/final_results.rds", "Outputs/site_climate_monthly.rds"),
                                  vars = list(CFG_COARSE_OVERRIDE = list(python = { py <- find_python(""); if (is.na(py)) "" else as.character(py) },
                                                                         reuse_runs = FALSE))), format = "file"),

  # ---- soil-health analyses (Part C) ---------------------------------------------------------
  tar_target(A4b_parts, step("A4b_Fig2_soil_health.R", file.path(WORK, "Fig2a_panel.rds"), s_A4b, in_partA,
                             needs = c("Outputs/Final_plus/us_oof_predictions.csv", "Outputs/Benchmark/benchmark_summary.csv")),
             format = "file"),
  tar_target(C1d_additivity, step("C1d_practice_additivity.R", file.path(PA, "additivity_effects.csv"), s_C1d, in_napeshm), format = "file"),
  tar_target(C1e_labels, step("C1e_labels_vs_SHMI.R", file.path(PA, "labels_vs_SHMI.csv"), s_C1e, in_napeshm), format = "file"),
  tar_target(C1f_budget, step("C1f_paired_data_budget.R", file.path(PA, c("Fig_data_budget_panels.rds", "data_budget.csv")),
                              s_C1f, in_napeshm), format = "file"),
  tar_target(C1g_label, step("C1g_label_vs_SHMI.R", file.path(PA, c("label_vs_SHMI.csv", "Fig_label_vs_SHMI_panel.rds", "Fig_label_within.png")),
                             s_C1g, in_napeshm), format = "file"),
  tar_target(C1h_SH, step("C1h_SHMI_vs_soil_health.R", file.path(PA, c("SHMI_vs_SH_SHAPE.csv", "Fig_SHMI_vs_SH_panel.rds", "Fig_SHMI_vs_SHAPE.png")),
                          s_C1h, in_napeshm), format = "file"),
  tar_target(C1_fig12, step("C1_Fig1_and_Table1.R",
                            file.path(WORK, c("Fig1.png", "Fig2.png", "ED_Fig_additivity.png", "Table1.csv", "SuppTable6_practice_pairs.csv")),
                            s_C1, in_napeshm, B2_main, C1d_additivity, C1f_budget, C1g_label, C1h_SH, A4b_parts), format = "file"),
  tar_target(C1c_coarse, step("C1c_ED_Fig5.R", file.path(WORK, "ED_Fig5.png"), s_C1c, B9_coarse), format = "file"),

  # ---- figures drawn from plot data (seconds) -----------------------------------------------
  tar_target(F3_yield, figure("F3_yield.R", "Outputs/figures/F3_yield.png", s_fig, B2_main, B2b_soybean, B3_stability), format = "file"),
  tar_target(F4_moderators, figure("F4_yield_moderators.R", "Outputs/figures/F4_yield_moderators.png", s_fig, B2_main), format = "file"),

  # ---- tables and the numbered set ------------------------------------------------------------
  tar_target(C2_tables, step("C2_export_tables.R", file.path(WORK, c("Tables.xlsx", "ED_Table1.csv")),
                             s_C2, B2_main, B2d_parts, B9_coarse, C1_fig12, in_partA, in_napeshm), format = "file"),
  tar_target(C3_numbered, step("C3_collect_figures.R", "Outputs/Manuscript/numbered/manifest_numbered.csv",
                               s_C3, C1_fig12, C1c_coarse, C1g_label, C1h_SH, C2_tables, F3_yield, F4_moderators, B2_main, B2d_parts,
                               B8b_treatments, B4_Neq, B5b_weights, B6_Nscen, C1e_labels), format = "file")
)
