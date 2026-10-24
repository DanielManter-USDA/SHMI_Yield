# SHMI_Yield

Code and derived data for Moore & Manter, *Capturing cropping system complexities at scale*:
- calibration of the Soil Health Management Index (SHMI) against measured soil health;
- its relationship with corn yield in ten long-term experiments.

SHMI scores are computed with the SHMI R package, https://github.com/DanielManter-USDA/SHMI (version 2.0.0 for the paper).

## Running everything

1. Open `SHMI_Yield.Rproj` in RStudio.
2. Run `install_packages.R`, then install the SHMI package. Python 3.9+ (`pip install pandas openpyxl numpy`) is needed only to rebuild workbooks (A0, B0, B8).
3. Download the SPEI grid into `Inputs/climate/` (see `Inputs/README.md`); everything else is in the repository. The raw AgCROS tables are needed only to rebuild the AgCROS workbooks (step B0), which are shipped.
4. Knit `run_all.Rmd`. To choose parts, use **Knit with Parameters**, or run:

   ```r
   rmarkdown::render("run_all.Rmd", params = list(run_calibration = TRUE, run_site_sensitivity = TRUE))
   ```

`run_all.Rmd` first checks the inputs. It then runs each step in a fresh R session (python for the conversions), logs each step to `Outputs/logs/<time>/`, and ends with the main result tables and figures.

| Parameter | Default | What it runs | Time |
|---|---|---|---|
| `run_napeshm_rebuild` | FALSE | A0: NAPESHM workbook from the raw files | minutes |
| `run_calibration` | FALSE | A1-A8: SEM, SHMI calibration, benchmarks, weight scenarios, SHAPE | hours (frozen for the paper) |
| `run_conversion` | FALSE | B0: rebuild the AgCROS workbooks and soil C even when present (they are always built when missing; rebuilding resets their cached site runs) | minutes |
| `run_yield` | TRUE | B1-B5, B7-B8: yield analyses, figures, treatment table | first run about an hour (sites are cached after that) |
| `run_n_scenarios` | TRUE | B6: SHMI slope under 20 N-accounting assumptions | 30-60 min |
| `run_coarse` | FALSE | B9: SHMI from remote-sensing-resolution records (soil-health skill and yield slope) | hours (NAPESHM noise repeats; 18 yield reruns) |

## Steps and scripts (`R/`, `python/`)

Every script below produces something reported in the paper; figure and table numbers are those of the manuscript.

| Step | Script | Output (in `Outputs/`) | In the paper |
|---|---|---|---|
| A0 | `python/napeshm_rebuild/convert.py`, `write_wb.py` | `NAPESHM_rebuilt.xlsx` (before QC) | NAPESHM workbook (Methods) |
| A1 | `Get_site_ssurgo.R`, `Prepare_SEM.R`, `SH_SEM_refinement_clean.R`, `Get_site_climate.R` (with `datasets.R`) | SEM data and measured soil health (lv_SH); climate normals | soil health latent variable; ED Fig. 1 |
| A2 | `SHMI_weight_tests.R`, `Run_calibration_scenarios.R` (with `SHMI_calibration_v3.R`), `Project_SH.R`, `SHMI_calibration_final.R`, `Summarize_results.R` | calibration scenarios; the US model applied to the Mexican sites | Supplementary Data Tables 1-2; Fig. 1b |
| A3 | `SHMI_calibration_plus.R`, `Fig5_Fig6_fit_plus.R` | `Final_plus/`: final weights, out-of-fold predictions | weights; Fig. 2; ED Fig. 2; ED Table 1 |
| A4 | `SHMI_benchmark.R` | `Benchmark/` | ED Fig. 4 |
| A5 | `Run_within_calibration_scenarios.R`, `ED_Fig7_weights_by_comparison.R` | `Calib_scenarios/` | ED Fig. 5 |
| A6 | `Prepare_SHAPE.R`, then the SHAPE app, `Merge_SHAPE_SHMI.R`, `SHMI_SHAPE_age.R` | `SHAPE/`, `SHAPE_age/` | ED Fig. 3 |
| A7-A8 | `Check_climate_source.R`, `Export_appendix_tables.R` | season check; SEM loadings and adjustment slopes | Methods; ED Fig. 1 |
| B0 | `python/agcros2shmi/agcros2shmi.py` (with the supplement and corrections), `python/build_soil_c.py`, `R/build_napeshm_crosswalk.R` | AgCROS workbooks in `Inputs/yield_workbooks/`; `yield/soil_c.csv`; `yield/napeshm_mediators_crosswalk.csv` | yield workbooks; mediation inputs |
| B1 | `Get_SPEI.R` | `yield/drought_years.csv` | drought analysis |
| B2 | `SHMI_yield_main.R` (with `SHMI_legume_credit.R`, `SHMI_manure_n.R`) | `yield_main_default_total_agcros/`: results_main.csv and panels | Fig. 3; Fig. 4a; ED Fig. 6 |
| B3 | `SHMI_yield_stability.R` | `yield_stability_default_total_agcros/` | Fig. 4b; ED Fig. 7 |
| B4 | `SHMI_N_equivalent.R` | N_equivalent.csv | fertilizer-N equivalent (Results) |
| B5 | `SHMI_soil_mediation.R` | mediation_results.csv | soil-health mediation (Results; Supplementary Table) |
| B6 | `SHMI_N_scenarios.R` | N_scenarios_all_sites.csv, N_scenarios_WICST.csv | N-accounting range (Methods) |
| B7 | `Assemble_manuscript_figures.R` | `manuscript_figures_default_total_agcros/`: Fig3, Fig4, ED_Fig8_stability (now ED Fig. 7) | Figs. 3-4; ED Fig. 7 |
| B8 | `python/build_study_treatments.py` | `SHMI_study_treatments.xlsx` | Supplementary Table of sites and treatments |
| B9 | `SHMI_coarse_records.R` with `python/coarse_records/degrade_records.py` | `Coarse/`: coarse_napeshm.csv, coarse_agreement.csv, coarse_yield.csv, coarse_yield_noise.csv, Fig_coarse_records.png | coarse-record test; ED Fig. 8 |

Figures are saved without titles or captions (the legends are in the manuscript); set `options(shmi.ms_figures = FALSE)` before running a script to keep them for working plots.

## Data in this repository (`data/`, `python/agcros2shmi/`)

- **`napeshm_yield_site_units.csv`:** NAPESHM units matched to the yield sites' treatments.
- **`soil_c.csv`:** soil organic carbon by plot, layer and date at the AgCROS experiments (output of `build_soil_c.py`), so the mediation test runs without the AgCROS tables.
- **`cache/nass_site_county_corn_yield.csv`:** USDA-NASS county corn yields for the yield sites, so no NASS key is needed.
- **`cache/site_soils_ssurgo.csv`:** gSSURGO soil groups by NAPESHM site (output of `Get_site_ssurgo.R`), so gSSURGO is not needed.
- **`SHAPE_by_unit_with_SHMI.csv`:** SHAPE scores, measured soil health and SHMI by NAPESHM unit (output of A6), used when Part A is not run.
- **`SHMI_study_treatments_descriptions.xlsx`:** treatment descriptions for the hand-built sites.
- **`supplements/NECCIRR_amendments.csv`:** NECCIRR fertilizer and manure, missing from AgCROS (UNL CropWatch 2017).
- **`corrections/NEMLTCRS_corrections.csv`:** two record errors in AgCROS, with their evidence.

## Notes

- Scripts use paths relative to the project root (`here::here()`). Any script can also be run on its own from the project.
- Site runs are cached in `Outputs/yield/site_runs`. The cache key covers each workbook, the period, the settings, the SHMI weights, the SHMI package (version and install time) and the scoring code. Delete the folder to force a full rescore.
- The converter writes one `CD_term_date` (rule T1) and a zero-mixing placeholder pass at no-till sites (rule D3). Both were workarounds for two issues fixed in the SHMI package (an all-blank required column; an empty sheet). They do not change any score and are kept so the shipped workbooks are reproduced exactly.

