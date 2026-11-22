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
# coarse-record test, hours) reruns only when forced: targets::tar_invalidate(B9_coarse);
# Sys.setenv(SHMI_FORCE = "B9"); targets::tar_make(); Sys.unsetenv("SHMI_FORCE").
# B2i (14 Nov 2026) writes the variance shares cited in Supplementary Methods (not the retired B2e).
# Figures (13 Nov 2026, Nature Communications): Figs 1-7 and Supplementary Fig. 1 are drawn by R/figures/F1-F7, S1 from plot
# data and saved panels; the analysis scripts no longer decide the final layout.
# 17 Nov 2026 (streamlined layout): Figs 1-5 = F1 (SEM | path), F2 (SHMI vs lv_SH), F3 (yield; each
# experiment), F4 (N, age, good and bad years), F5 (records); Table 2 = C1j + B2d, Table 3 = C1i + C1f.
# Retired 17 Nov 2026: A4b (parts panel), C1e and C1g (labels vs SHMI), figures F3_practices,
# F4_components, F5_yield, F6_outcome and S1 (scripts moved to R/archive/figures_retired_17Nov).
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
  tar_target(s_B2f, sc("B2f_SHMI_yield_paired_budget.R"), format = "file"),
  tar_target(s_B2i, sc("B2i_SHMI_variance_shares.R"), format = "file"),
  tar_target(s_B3,  sc("B3_SHMI_yield_stability.R"), format = "file"),
  tar_target(s_B4,  sc("B4_SHMI_N_equivalent.R"), format = "file"),
  tar_target(s_B5b, sc("B5b_SHMI_yield_weights_sensitivity.R"), format = "file"),
  tar_target(s_B6,  sc("B6_SHMI_N_scenarios.R"), format = "file"),
  tar_target(s_B8,  "python/build_study_treatments.py", format = "file"),
  tar_target(s_B8b, sc("B8b_study_treatments_with_yield.R"), format = "file"),
  tar_target(s_B9,  c(sc("B9_SHMI_coarse_records.R"), "python/coarse_records/degrade_records.py"), format = "file"),
  tar_target(s_C1d, sc("C1d_practice_additivity.R"), format = "file"),
  tar_target(s_C1f, sc("C1f_paired_data_budget.R"), format = "file"),
  tar_target(s_C1h, sc("C1h_SHMI_vs_soil_health.R"), format = "file"),
  tar_target(s_C1i, sc("C1i_pair_meta.R"), format = "file"),
  tar_target(s_C1k, sc("C1k_ceiling.R"), format = "file"),
  tar_target(s_C1l, sc("C1l_rf_within.R"), format = "file"),
  tar_target(s_C1j, sc("C1j_SHMI_vs_parts.R"), format = "file"),
  tar_target(s_C1,  sc("C1_Fig1_and_Table1.R"), format = "file"),
  tar_target(s_C1c, sc("C1c_ED_Fig5.R"), format = "file"),
  tar_target(s_C2,  sc("C2_export_tables.R"), format = "file"),
  tar_target(s_C3,  sc("C3_collect_figures.R"), format = "file"),
  # figure scripts: shared look and the two yield panel builders (used by F4 and F5)
  tar_target(s_fig, sc(c("figures/fig_common.R", "figures/yield_panels.R", "figures/moderator_panels.R")), format = "file"),
  tar_target(s_F1, sc("figures/F1_soil_health_model.R"), format = "file"),
  tar_target(s_F2, sc("figures/F2_system_score.R"), format = "file"),
  tar_target(s_F3, sc("figures/F3_yield.R"), format = "file"),
  tar_target(s_F4, sc("figures/F4_yield_moderators.R"), format = "file"),
  tar_target(s_F5, sc("figures/F5_records.R"), format = "file"),

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
  # composition of within-site SHMI variation, NAPESHM vs the ten experiments (Supplementary Methods; Fig. 6a)
  tar_target(B2i_shares, step("B2i_SHMI_variance_shares.R",
                              file.path("Outputs/variance_shares", c("variance_shares.csv", "variance_shares_sites.csv")),
                              s_B2i, B2_main, in_partA, needs = "Outputs/Final_plus/final_results.rds"), format = "file"),
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
  tar_target(C1d_additivity, step("C1d_practice_additivity.R", file.path(PA, "additivity_effects.csv"), s_C1d, in_napeshm), format = "file"),
  tar_target(C1f_budget, step("C1f_paired_data_budget.R", file.path(PA, c("Fig_data_budget_panels.rds", "data_budget.csv")),
                              s_C1f, in_napeshm), format = "file"),
  tar_target(C1h_SH, step("C1h_SHMI_vs_soil_health.R", c(file.path(PA, c("SHMI_vs_SH_SHAPE.csv", "Fig_SHMI_vs_SH_panel.rds", "Fig_SHMI_vs_SHAPE.png", "SHMI_parts_by_indicator.csv")),
                                                 file.path(PANEL_DIR, c("panel_SHMI_vs_SH.rds", "panel_SHMI_vs_SHAPE.rds"))),
                          s_C1h, in_napeshm), format = "file"),
  # meta-analysis of the strict matched pairs from field replicates (Table 3; 16 Nov 2026)
  tar_target(C1i_meta, step("C1i_pair_meta.R",
                            c(file.path("Outputs/Pair_meta", c("pooled.csv", "metaregression.csv", "pair_effects.csv", "rep_scores.csv", "plotdata_pair_meta.rds",
                                                              "SuppTable6_practice_pairs.csv", "scores_vs_meta.csv")),
                              file.path(PANEL_DIR, c("panel_meta_forest.rds", "panel_meta_reg.rds"))),
                            s_C1i, in_napeshm, C1_fig12,
                            needs = c("Inputs/napeshm/soildata.csv", "Outputs/semDF_prepped.rds", "Outputs/sem_noMX_projection.rds")),
             format = "file"),
  # Table 2, soil-health columns: SHMI against each part alone, ΔAICc and leave-one-site-out R² (17 Nov 2026)
  tar_target(C1j_parts, step("C1j_SHMI_vs_parts.R", "Outputs/Parts/SHMI_vs_parts.csv", s_C1j, in_partA, in_napeshm), format = "file"),
  # ceiling on explainable within-site variation, from the field replicates (21 Nov 2026)
  tar_target(C1k_ceiling, step("C1k_ceiling.R", file.path("Outputs/Ceiling", c("ceiling.csv", "rep_scores.csv")),
                               s_C1k, in_partA, in_napeshm,
                               needs = c("Inputs/napeshm/soildata.csv", "Outputs/semDF_prepped.rds", "Outputs/sem_noMX_projection.rds")),
             format = "file"),
  # random forests fitted to within-site variation of the same record inputs (21 Nov 2026)
  tar_target(C1l_rfwithin, step("C1l_rf_within.R", file.path("Outputs/Benchmark_within", c("rf_within.csv", "rf_within_predictions.csv")),
                                s_C1l, in_partA), format = "file"),
  tar_target(C1_fig12, step("C1_Fig1_and_Table1.R",
                            c(file.path(WORK, c("Fig1.png", "ED_Fig_additivity.png", "Table1.csv", "SuppTable6_practice_pairs.csv")),
                              file.path(PANEL_DIR, c("panel_path.rds", "panel_additivity.rds"))),
                            s_C1, in_napeshm, B2_main, C1d_additivity), format = "file"),
  tar_target(C1c_coarse, step("C1c_ED_Fig5.R", c(file.path(WORK, "ED_Fig5.png"), file.path(PANEL_DIR, "panel_records.rds")),
                              s_C1c, B9_coarse), format = "file"),

  # ---- figures drawn from plot data and saved panels (seconds; Nature Communications layout) ----
  tar_target(F1_fig, figure("F1_soil_health_model.R", "Outputs/figures/F1_soil_health_model.png", s_fig, s_F1, in_partA, C1_fig12), format = "file"),
  tar_target(F2_fig, figure("F2_system_score.R", "Outputs/figures/F2_system_score.png", s_fig, s_F2, C1h_SH), format = "file"),
  tar_target(F3_fig, figure("F3_yield.R", "Outputs/figures/F3_yield.png", s_fig, s_F3, B2_main, B2b_soybean), format = "file"),
  tar_target(F4_fig, figure("F4_yield_moderators.R", "Outputs/figures/F4_yield_moderators.png", s_fig, s_F4, B2_main, B3_stability), format = "file"),
  tar_target(F5_fig, figure("F5_records.R", "Outputs/figures/F5_records.png", s_fig, s_F5, C1c_coarse), format = "file"),

  # ---- tables and the numbered set ------------------------------------------------------------
  tar_target(C2_tables, step("C2_export_tables.R", file.path(WORK, c("Tables.xlsx", "ED_Table1.csv", "Table2.csv", "Table3.csv")),
                             s_C2, B2_main, B2d_parts, B9_coarse, C1_fig12, C1j_parts, C1k_ceiling, C1l_rfwithin, C1i_meta, C1f_budget, in_partA, in_napeshm), format = "file"),
  tar_target(C3_numbered, step("C3_collect_figures.R", "Outputs/Manuscript/numbered/manifest_numbered.csv",
                               s_C3, C1_fig12, C2_tables, F1_fig, F2_fig, F3_fig, F4_fig, F5_fig, B2_main, B2d_parts,
                               B8b_treatments, B4_Neq, B5b_weights, B6_Nscen, B2i_shares, C1i_meta), format = "file")
)
