#!/usr/bin/env bash
# Run from the SHMI_Yield repository root (Git Bash), BEFORE unzipping the updated files.
# git mv keeps each script's history under its new name.
set -e
git rm -q R/SHMI_benchmark_additions.R
git mv R/Get_site_ssurgo.R R/A1a_Get_site_ssurgo.R
git mv R/Prepare_SEM.R R/A1b_Prepare_SEM.R
git mv R/SH_SEM_refinement_clean.R R/A1c_soil_health_SEM.R
git mv R/Get_site_climate.R R/A1d_Get_site_climate.R
git mv R/SHMI_weight_tests.R R/A2a_SHMI_weight_tests.R
git mv R/Run_calibration_scenarios.R R/A2b_Run_calibration_scenarios.R
git mv R/Project_SH.R R/A2c_Project_SH.R
git mv R/SHMI_calibration_final.R R/A2d_SHMI_calibration_final.R
git mv R/Summarize_results.R R/A2e_Summarize_results.R
git mv R/SHMI_calibration_plus.R R/A3a_SHMI_calibration_plus.R
git mv R/Fig5_Fig6_fit_plus.R R/A3b_calibration_report_figures.R
git mv R/SHMI_benchmark.R R/A4a_SHMI_benchmark.R
git mv R/Fig2_soil_health.R R/A4b_Fig2_soil_health.R
git mv R/Run_within_calibration_scenarios.R R/A5a_Run_within_calibration_scenarios.R
git mv R/ED_Fig7_weights_by_comparison.R R/A5b_SuppFig1_weights_by_comparison.R
git mv R/Prepare_SHAPE.R R/A6a_Prepare_SHAPE.R
git mv R/Merge_SHAPE_SHMI.R R/A6b_Merge_SHAPE_SHMI.R
git mv R/SHMI_SHAPE_age.R R/A6c_SHMI_SHAPE_age.R
git mv R/Check_climate_source.R R/A7_Check_climate_source.R
git mv R/Export_appendix_tables.R R/A8_Export_appendix_tables.R
git mv R/build_napeshm_crosswalk.R R/B0_build_napeshm_crosswalk.R
git mv R/Get_SPEI.R R/B1_Get_SPEI.R
git mv R/SHMI_yield_main.R R/B2_SHMI_yield_main.R
git mv R/SHMI_yield_stability.R R/B3_SHMI_yield_stability.R
git mv R/SHMI_N_equivalent.R R/B4_SHMI_N_equivalent.R
git mv R/SHMI_soil_mediation.R R/B5_SHMI_soil_mediation.R
git mv R/SHMI_yield_weights_sensitivity.R R/B5b_SHMI_yield_weights_sensitivity.R
git mv R/SHMI_N_scenarios.R R/B6_SHMI_N_scenarios.R
git mv R/Assemble_manuscript_figures.R R/B7_Assemble_manuscript_figures.R
git mv R/SHMI_coarse_records.R R/B9_SHMI_coarse_records.R
git mv R/SHMI_calibration_v3.R R/calibration_engine.R
echo "Renamed 31 scripts. Now unzip SHMI_Yield_reorganized.zip over the repository (overwrite), then: git status"
