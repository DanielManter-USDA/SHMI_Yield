################################################################################
# _targets.R  (21 Nov 2026, submission)   The SHMI_Yield pipeline
#
#   targets::tar_make()          run every step that is out of date (and only those)
#   targets::tar_visnetwork()    the pipeline as a graph: what depends on what, what is stale
#   targets::tar_outdated()      the steps tar_make() would run now
#   targets::tar_read(B2_main)   the files a step wrote (or any target's value)
#   targets::tar_invalidate(B3)  force a step to rerun next time
#
# Once per project: targets::tar_config_set(reporter_make = "timestamp")
#
# Each step runs one script from R/ (its log goes to Outputs/logs/<run>/) and returns the
# files later steps read. A step reruns when its script, its input files or an upstream
# step's files change. Changing a figure script in R/figures/ redraws only that figure.
#
# Part A (the calibration, run_all.R) is frozen: its outputs are tracked as input files.
# Display items: Figs 1-4 (R/figures/F1-F4), Tables 1-2 (C2), Supplementary Tables 1-6
# (C1m, C2, B8/B8b; combined by C3). Every step here feeds a number, figure or table in the
# paper; scripts for analyses not in the paper are in R/archive and are not run.
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
PANEL_DIR <- "Outputs/panels"
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
  tar_target(in_cache, present("data/cache/nass_site_county_corn_yield.csv", "data/cache/site_soils_ssurgo.csv"), format = "file"),
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
  tar_target(s_B2i, sc("B2i_SHMI_variance_shares.R"), format = "file"),
  tar_target(s_B3,  sc("B3_SHMI_yield_stability.R"), format = "file"),
  tar_target(s_B4,  sc("B4_SHMI_N_equivalent.R"), format = "file"),
  tar_target(s_B6,  sc("B6_SHMI_N_scenarios.R"), format = "file"),
  tar_target(s_B8,  "python/build_study_treatments.py", format = "file"),
  tar_target(s_B8b, sc("B8b_study_treatments_with_yield.R"), format = "file"),
  tar_target(s_C1,  sc("C1_Fig1_path.R"), format = "file"),
  tar_target(s_C1f, sc("C1f_paired_data_budget.R"), format = "file"),
  tar_target(s_C1h, sc("C1h_SHMI_vs_soil_health.R"), format = "file"),
  tar_target(s_C1j, sc("C1j_SHMI_vs_parts.R"), format = "file"),
  tar_target(s_C1k, sc("C1k_ceiling.R"), format = "file"),
  tar_target(s_C1l, sc("C1l_rf_within.R"), format = "file"),
  tar_target(s_C1m, sc("C1m_supp_tables.R"), format = "file"),
  tar_target(s_C2,  sc("C2_export_tables.R"), format = "file"),
  tar_target(s_C3,  sc("C3_collect_figures.R"), format = "file"),
  # figure scripts: shared look and the two yield panel builders (used by F3 and F4)
  tar_target(s_fig, sc(c("figures/fig_common.R", "figures/yield_panels.R", "figures/moderator_panels.R")), format = "file"),
  tar_target(s_F1, sc("figures/F1_soil_health_model.R"), format = "file"),
  tar_target(s_F2, sc("figures/F2_system_score.R"), format = "file"),
  tar_target(s_F3, sc("figures/F3_yield.R"), format = "file"),
  tar_target(s_F4, sc("figures/F4_yield_moderators.R"), format = "file"),

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
  # rolling (current-season) scores used by B2d
  tar_target(B2c_record, step("B2c_SHMI_record_length.R",
                              "Outputs/yield_record_length/rolling_shmi.rds",
                              s_B2c, B2_main), format = "file"),
  # Table 2, corn-yield columns: SHMI against each part
  tar_target(B2d_parts, step("B2d_SHMI_yield_parts.R", file.path("Outputs/yield_parts", c("yield_parts.csv", "yield_parts_loo.csv")),
                             s_B2d, B2_main, B2c_record), format = "file"),
  # share of each experiment's within-site SHMI variance from each part (Fig. 3b, Supplementary Table 2)
  tar_target(B2i_shares, step("B2i_SHMI_variance_shares.R",
                              file.path("Outputs/variance_shares", c("variance_shares.csv", "variance_shares_sites.csv")),
                              s_B2i, B2_main, in_partA, needs = "Outputs/Final_plus/final_results.rds"), format = "file"),
  tar_target(B3_stability, step("B3_SHMI_yield_stability.R",
                                file.path(STAB, c("plotdata_stability.rds", "Fig2c_years.csv", "insurance_by_cause.csv")),
                                s_B3, s_B2, B2_main, B1_spei), format = "file"),
  # fertilizer-N equivalent (Results, Supplementary Methods)
  tar_target(B4_Neq, step("B4_SHMI_N_equivalent.R", file.path(MAIN, c("N_equivalent.csv", "N_equivalent_by_site.csv")),
                          s_B4, B2_main), format = "file"),
  # 20 N-accounting alternatives (Supplementary Table 6)
  tar_target(B6_Nscen, step("B6_SHMI_N_scenarios.R", file.path(MAIN, c("N_scenarios_all_sites.csv", "N_scenarios_WICST.csv")),
                            s_B6, s_B2, B2_main), format = "file"),
  # every treatment (Supplementary Table 3)
  tar_target(B8_treatments, step_py(c("python/build_study_treatments.py", "--wb-dir", "Inputs/yield_workbooks",
                                      "--units-dir", "Inputs/yield_workbooks",
                                      "--old", "data/SHMI_study_treatments_descriptions.xlsx",
                                      "--out", "Outputs/SHMI_study_treatments.xlsx"),
                                    "B8_study_treatments", "Outputs/SHMI_study_treatments.xlsx", s_B8, in_workbooks), format = "file"),
  tar_target(B8b_treatments, step("B8b_study_treatments_with_yield.R", "Outputs/SHMI_study_treatments_with_yield.xlsx",
                                  s_B8b, B8_treatments, B2_main, B2b_soybean), format = "file"),

  # ---- soil-health analyses (Part C) ---------------------------------------------------------
  # Fig. 1b, the path diagram
  tar_target(C1_path, step("C1_Fig1_path.R", c(file.path(WORK, "Fig1_path_coefficients.csv"), file.path(PANEL_DIR, "panel_path.rds")),
                           s_C1, in_napeshm, in_partA), format = "file"),
  # share of systems usable by strict single-practice pairs (Discussion, Supplementary Methods)
  tar_target(C1f_budget, step("C1f_paired_data_budget.R", file.path(PA, c("Fig_data_budget_panels.rds", "data_budget.csv")),
                              s_C1f, in_napeshm), format = "file"),
  # Fig. 2, SHMI against soil health at three scales; SHAPE
  tar_target(C1h_SH, step("C1h_SHMI_vs_soil_health.R", c(file.path(PA, c("SHMI_vs_SH_SHAPE.csv", "Fig_SHMI_vs_SH_panel.rds", "Fig_SHMI_vs_SHAPE.png", "SHMI_parts_by_indicator.csv")),
                                                 file.path(PANEL_DIR, c("panel_SHMI_vs_SH.rds", "panel_SHMI_vs_SHAPE.rds"))),
                          s_C1h, in_napeshm), format = "file"),
  # Table 2, soil-health columns: SHMI against each part alone, ΔAICc and leave-one-site-out R²
  tar_target(C1j_parts, step("C1j_SHMI_vs_parts.R", "Outputs/Parts/SHMI_vs_parts.csv", s_C1j, in_partA, in_napeshm), format = "file"),
  # ceiling on explainable within-site variation, from the field replicates
  tar_target(C1k_ceiling, step("C1k_ceiling.R", file.path("Outputs/Ceiling", c("ceiling.csv", "rep_scores.csv")),
                               s_C1k, in_partA, in_napeshm,
                               needs = c("Inputs/napeshm/soildata.csv", "Outputs/semDF_prepped.rds", "Outputs/sem_noMX_projection.rds")),
             format = "file"),
  # random forests fitted to within-site variation of the same record inputs
  tar_target(C1l_rfwithin, step("C1l_rf_within.R", file.path("Outputs/Benchmark_within", c("rf_within.csv", "rf_within_predictions.csv")),
                                s_C1l, in_partA), format = "file"),
  # Supplementary Tables 1, 4, 5 and 6
  tar_target(C1m_supp, step("C1m_supp_tables.R",
                            file.path(WORK, c("SuppTable_sites.csv", "SuppTable_sem.csv", "SuppTable_commonality.csv", "SuppTable_N_scenarios.csv")),
                            s_C1m, in_partA, in_napeshm, B6_Nscen,
                            needs = c("Outputs/semDF_prepped.rds", "Outputs/sem_noMX_projection.rds")), format = "file"),

  # ---- figures drawn from plot data and saved panels (seconds) ------------------------------
  tar_target(F1_fig, figure("F1_soil_health_model.R", "Outputs/figures/F1_soil_health_model.png", s_fig, s_F1, in_partA, C1_path), format = "file"),
  tar_target(F2_fig, figure("F2_system_score.R", "Outputs/figures/F2_system_score.png", s_fig, s_F2, C1h_SH), format = "file"),
  tar_target(F3_fig, figure("F3_yield.R", "Outputs/figures/F3_yield.png", s_fig, s_F3, B2_main, B2b_soybean), format = "file"),
  tar_target(F4_fig, figure("F4_yield_moderators.R", "Outputs/figures/F4_yield_moderators.png", s_fig, s_F4, B2_main, B3_stability), format = "file"),

  # ---- tables and the numbered set ------------------------------------------------------------
  tar_target(C2_tables, step("C2_export_tables.R", file.path(WORK, c("Tables.xlsx", "ED_Table1.csv", "Table2.csv", "SuppTable1.csv")),
                             s_C2, B2_main, B2d_parts, C1j_parts, C1k_ceiling, C1l_rfwithin, in_partA, in_napeshm), format = "file"),
  tar_target(C3_numbered, step("C3_collect_figures.R", "Outputs/Manuscript/numbered/manifest_numbered.csv",
                               s_C3, C2_tables, C1m_supp, F1_fig, F2_fig, F3_fig, F4_fig, B8b_treatments,
                               B4_Neq, B2i_shares, C1f_budget), format = "file")
)
