# SHMI_Yield

Code and derived data for Moore & Manter, *Capturing cropping system complexities at scale*:
- calibration of the Soil Health Management Index (SHMI) against measured soil health;
- its relationship with corn yield in ten long-term experiments.

SHMI scores are computed with the SHMI R package, https://github.com/DanielManter-USDA/SHMI (version 2.0.0 for the paper).

## Running everything

1. Open `SHMI_Yield.Rproj` in RStudio.
2. Run `install_packages.R`, then install the SHMI package. Install python 3.9+ with `pip install pandas openpyxl numpy`.
3. Put the inputs listed in `Inputs/README.md` in `Inputs/`, and the NASS key in `~/.Renviron`.
4. Knit `run_all.Rmd`. To choose parts, use **Knit with Parameters**, or run:

   ```r
   rmarkdown::render("run_all.Rmd", params = list(run_calibration = TRUE, run_site_sensitivity = TRUE))
   ```

`run_all.Rmd` first checks the inputs. It then runs each step in a fresh R session (python for the conversions), logs each step to `Outputs/logs/<time>/`, and ends with the main result tables and figures.

| Parameter | Default | What it runs | Time |
|---|---|---|---|
| `run_napeshm_rebuild` | FALSE | A0: NAPESHM workbook from the raw files | minutes |
| `run_calibration` | FALSE | A1-A8: SEM, SHMI calibration, benchmark, weight scenarios, SHAPE | hours (frozen for the paper) |
| `run_conversion` | FALSE | B0: rebuild the AgCROS workbooks and soil C even when present (they are always built when missing; rebuilding resets their cached site runs) | minutes |
| `run_yield` | TRUE | B1-B5, B7-B8: yield analyses, figures, treatment table | first run about an hour (sites are cached after that) |
| `run_n_scenarios` | TRUE | B5b: SHMI slope under 20 N-accounting assumptions | 30-60 min |
| `run_weights` | TRUE | B6: yield results under three SHMI weight sets | long (rescoring every site twice) |
| `run_site_sensitivity` | FALSE | B9: 7 sites, + INWLTPAC, + MNMOFS, West Lafayette refinements | quick once sites are cached |

## Steps and scripts (`R/`, `python/`)

| Step | Script | Output (in `Outputs/`) |
|---|---|---|
| A0 | `python/napeshm_rebuild/convert.py`, `write_wb.py` | `NAPESHM_rebuilt.xlsx` (before QC) |
| A1 | `Get_site_ssurgo.R`, `Prepare_SEM.R`, `SH_SEM_refinement_clean.R`, `Get_site_climate.R` | SEM data and measured soil health; climate normals |
| A2 | `SHMI_weight_tests.R`, `Run_calibration_scenarios.R` (with `SHMI_calibration_v3.R`), `Project_SH.R`, `SHMI_calibration_final.R`, `Compare_v110_v120.R`, `Summarize_results.R` | calibration scenarios |
| A3 | `SHMI_calibration_plus.R`, `Fig5_Fig6_fit_plus.R` | `Final_plus/`: final weights; Fig. 1-2 panels; Table 10b |
| A4 | `SHMI_benchmark.R` | `Benchmark/` |
| A5 | `Run_within_calibration_scenarios.R`, `ED_Fig7_weights_by_comparison.R` | `Calib_scenarios/`: ED Fig. 7 |
| A6 | `Prepare_SHAPE.R`, then the SHAPE app, `Merge_SHAPE_SHMI.R`, `SHMI_SHAPE_age.R` | `SHAPE/`, `SHAPE_age/` |
| A7-A8 | `Check_climate_source.R`, `Export_appendix_tables.R` | climate source check; appendix tables |
| B0 | `python/agcros2shmi/agcros2shmi.py` (with the supplement and corrections), `python/build_soil_c.py`, `R/build_napeshm_crosswalk.R` | AgCROS workbooks in `Inputs/yield_workbooks/`; `yield/soil_c.csv`; `yield/napeshm_mediators_crosswalk.csv` |
| B1 | `Get_SPEI.R` | `yield/drought_years.csv` |
| B2 | `SHMI_yield_main.R` (with `SHMI_legume_credit.R`, `SHMI_manure_n.R`) | `yield_main_default_total_agcros/`: results_main.csv; Fig. 3-4 panels |
| B3 | `SHMI_yield_stability.R` | `yield_stability_default_total_agcros/`: ED Fig. 8 panels |
| B4 | `SHMI_N_equivalent.R` | N_equivalent.csv |
| B5 | `SHMI_soil_mediation.R` | mediation_results.csv |
| B5b | `SHMI_N_scenarios.R` | N_scenarios_all_sites.csv, N_scenarios_WICST.csv |
| B6 | `SHMI_yield_weights_sensitivity.R` | `yield_weights_agcros/` |
| B7 | `Assemble_manuscript_figures.R` | `manuscript_figures_default_total_agcros/`: Fig3, Fig4, ED_Fig8 |
| B8 | `python/build_study_treatments.py` | `SHMI_study_treatments.xlsx` (Supplementary Table) |
| B9 | `SHMI_yield_main.R` with overrides | `yield_main_<run>/`, summarised in the report |
| (check) | `Validate_converter_SHMI.R`, `python/agcros2shmi/compare_workbooks.py` | converted vs hand-built workbooks |

## Data in this repository (`data/`, `python/agcros2shmi/`)

- **`napeshm_yield_site_units.csv`:** NAPESHM units matched to the yield sites' treatments.
- **`SHAPE_by_unit_with_SHMI.csv`:** SHAPE scores, measured soil health and SHMI by NAPESHM unit (output of A6), used when Part A is not run.
- **`SHMI_study_treatments_descriptions.xlsx`:** treatment descriptions for the hand-built sites.
- **`supplements/NECCIRR_amendments.csv`:** NECCIRR fertilizer and manure, missing from AgCROS (UNL CropWatch 2017).
- **`corrections/NEMLTCRS_corrections.csv`:** two record errors in AgCROS, with their evidence.

## Notes

- Scripts use paths relative to the project root (`here::here()`). Any script can also be run on its own from the project.
- Site runs are cached in `Outputs/yield/site_runs`. The cache key covers each workbook, the period, the settings, the SHMI weights, the SHMI package (version and install time) and the scoring code. Delete the folder to force a full rescore.
