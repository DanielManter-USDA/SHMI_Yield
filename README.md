# SHMI_Yield

Code and derived data for Moore & Manter, *Soil health tracks management systems, not the sum of individual practices* (in review):

- calibration of the Soil Health Management Index (SHMI) against measured soil health on 354 US systems at 74 NAPESHM sites, with the SHAPE benchmark;
- its relationship with corn yield in ten long-term experiments (DRIVES and AgCROS);
- how much management detail a record needs.

SHMI scores are computed with the SHMI R package, version 2.0.0: https://github.com/DanielManter-USDA/SHMI, archived at https://doi.org/10.5281/zenodo.23682170.

## Running everything

1. Open `SHMI_Yield.Rproj` in RStudio.
2. Run `install_packages.R`, then install SHMI 2.0.0 (`remotes::install_github("DanielManter-USDA/SHMI@v2.0.0")`). Python 3.9+ (`pip install pandas openpyxl numpy`) is needed only to rebuild workbooks (A0, B0, B8) and for the coarse-record test (B9).
3. Download the SPEI grid into `Inputs/climate/` (see `Inputs/README.md`); everything else is in the repository. The raw AgCROS tables are needed only to rebuild the AgCROS workbooks (B0), which are shipped.
4. Knit `run_all.Rmd`. To choose parts, use **Knit with Parameters**, or run, for example:

   ```r
   rmarkdown::render("run_all.Rmd", params = list(run_calibration = TRUE, run_weights = TRUE))
   ```

`run_all.Rmd` first checks the inputs. It then runs each step in a fresh R session (python for the conversions), logs each step to `Outputs/logs/<time>/`, and ends with the main result tables and figures.

| Parameter | Default | What it runs | Time |
|---|---|---|---|
| `run_napeshm_rebuild` | FALSE | A0: NAPESHM workbook from the raw files | minutes |
| `run_calibration` | FALSE | A1-A8: soil health model, SHMI calibration, benchmark, weights by comparison, SHAPE | hours (frozen for the paper) |
| `run_conversion` | FALSE | B0: rebuild the AgCROS workbooks and soil C even when present (they are always built when missing) | minutes |
| `run_yield` | TRUE | B1-B5, B7-B8: yield analyses, figures, treatment table | first run about an hour (site runs are cached) |
| `run_weights` | FALSE | B5b: yield slope with weights fitted between or within sites | reruns B2 three times |
| `run_n_scenarios` | TRUE | B6: SHMI slope under 20 N-accounting assumptions | 30-60 min |
| `run_coarse` | FALSE | B9: SHMI from less detailed records | hours |

## Where each result comes from

Figure and table numbers are those of the manuscript. Outputs are in `Outputs/`.

| In the paper | Step and script | Output |
|---|---|---|
| Fig. 1a, path diagram (Diversity → Cover/Roots 0.49; Inverse Disturbance → Cover/Roots 0.08; SHMI → lv_SH 0.41) | coefficients from `data/SHAPE_by_unit_with_SHMI.csv` by least squares; drawing not yet in the repository | |
| Fig. 1b, scores of example systems; Table 1 | SHMI package, example records; not yet in the repository | |
| Fig. 2a-c, out-of-fold fit | A3a `A3a_SHMI_calibration_plus.R` | `Final_plus/us_oof_predictions.csv` |
| Fig. 2d, benchmark; Supplementary Data Table 4 | A4a `A4a_SHMI_benchmark.R`; figure `A4b_Fig2_soil_health.R` | `Benchmark/benchmark_summary.csv`, `benchmark_bootstrap.csv`; `Manuscript/Fig2_soil_health.png` |
| Fig. 3 | B2 `B2_SHMI_yield_main.R`; panels assembled by B7 `B7_Assemble_manuscript_figures.R` | `yield_main_default_total_agcros/results_main.csv`; `manuscript_figures_default_total_agcros/Fig3.png` |
| Fig. 4 | B3 `B3_SHMI_yield_stability.R` (SPEI from B1 `B1_Get_SPEI.R`); assembled by B7 | `yield_stability_default_total_agcros/`; `manuscript_figures_default_total_agcros/Fig4.png` |
| Extended Data Fig. 1, soil health model | A1c `A1c_soil_health_SEM.R`; loadings and adjustment slopes A8 `A8_Export_appendix_tables.R` | SEM outputs |
| Extended Data Fig. 2, weights; Extended Data Table 1 | A3a `A3a_SHMI_calibration_plus.R`; report figures A3b `A3b_calibration_report_figures.R` | `Final_plus/` |
| Extended Data Fig. 3, SHAPE | A6a-A6c (`A6a_Prepare_SHAPE.R`, the SHAPE app, `A6b_Merge_SHAPE_SHMI.R`, `A6c_SHMI_SHAPE_age.R`) | `SHAPE/`, `SHAPE_age/` |
| Extended Data Fig. 4, slope per experiment; Supplementary Data Table 1 | B2 (section 7b: heterogeneity, contrast groups, scale checks); assembled by B7; treatment table B8 | `site_SHMI_ranges.csv`, `SHMI_heterogeneity.csv`, `SHMI_scale_check*.csv`; `ED_Fig4_sites.png`; `SHMI_study_treatments.xlsx` |
| Extended Data Fig. 5, coarse records; Supplementary Data Table 2 | B9 `B9_SHMI_coarse_records.R` with `python/coarse_records/degrade_records.py` | `Coarse/` |
| Supplementary Data Table 3, soil properties as carriers of the slope | B5 `B5_SHMI_soil_mediation.R` (crosswalk from B0 `B0_build_napeshm_crosswalk.R`) | `mediation_results.csv` |
| Supplementary Data Fig. 1, weights by comparison | A5a `A5a_Run_within_calibration_scenarios.R`; figure A5b `A5b_SuppFig1_weights_by_comparison.R` | `Calib_scenarios/` |
| Supplementary Data Fig. 2, effect sizes | B2 | `Fig_effect_sizes.png`, `Fig_effect_sizes.csv` |
| Methods: fertilizer-N equivalent | B4 `B4_SHMI_N_equivalent.R` | `N_equivalent.csv` |
| Methods: yield slope and the choice of weights | B5b `B5b_SHMI_yield_weights_sensitivity.R` | `yield_weights_agcros/weights_slopes.csv` |
| Methods: N-accounting range | B6 `B6_SHMI_N_scenarios.R` | `N_scenarios_all_sites.csv` |
| Methods: climate source check | A7 `A7_Check_climate_source.R` | season check |

## Scripts (`R/`)

Scripts are numbered in the order `run_all.Rmd` runs them: **A** = calibration (Part A), **B** = yield and records (Part B). Every step runs with SHMI 2.0.0. Alternative index structures were compared during development with an earlier build of the package; those scripts are not part of this release (A2a and A2b are therefore absent).

| Step | Script | What it does |
|---|---|---|
| A1a | `A1a_Get_site_ssurgo.R` | gSSURGO soil groups by NAPESHM site (cached in `data/cache/`) |
| A1b | `A1b_Prepare_SEM.R` | SEM data, cropping-system classes (annual, perennial, reference, fallow) |
| A1c | `A1c_soil_health_SEM.R` | climate and texture adjustment; soil health latent variable (lv_SH) |
| A1d | `A1d_Get_site_climate.R` | WorldClim monthly normals by site |
| A2c-A2e | `A2c_Project_SH.R`, `A2d_SHMI_calibration_final.R`, `A2e_Summarize_results.R` | the US model applied to the Mexican sites (a check, not reported); the four-sub-index calibration kept for comparison; summary tables |
| A3a | `A3a_SHMI_calibration_plus.R` | final calibration (three sub-indices plus the amendment bonus); out-of-fold predictions |
| A3b | `A3b_calibration_report_figures.R` | in-sample and out-of-fold report figures and agreement table |
| A4a | `A4a_SHMI_benchmark.R` | SHMI against its parts, reduced versions, fixed published weights and unconstrained models |
| A4b | `A4b_Fig2_soil_health.R` | Fig. 2 |
| A5a-A5b | `A5a_Run_within_calibration_scenarios.R`, `A5b_SuppFig1_weights_by_comparison.R` | weights fitted between sites, to all variation and within sites |
| A6a-A6c | `A6a_Prepare_SHAPE.R`, `A6b_Merge_SHAPE_SHMI.R`, `A6c_SHMI_SHAPE_age.R` | SHAPE benchmark |
| A7 | `A7_Check_climate_source.R` | WorldClim against Daymet growing seasons |
| A8 | `A8_Export_appendix_tables.R` | SEM loadings and adjustment slopes |
| B0 | `B0_build_napeshm_crosswalk.R` | NAPESHM units matched to the yield sites' treatments |
| B1 | `B1_Get_SPEI.R` | growing-season SPEI for every site-year |
| B2 | `B2_SHMI_yield_main.R` | within-site yield model, SHMI x N, SHMI x age, site slopes, heterogeneity, effect sizes |
| B3 | `B3_SHMI_yield_stability.R` | year conditions, drought, worst years |
| B4 | `B4_SHMI_N_equivalent.R` | SHMI in fertilizer-N units |
| B5 | `B5_SHMI_soil_mediation.R` | soil properties as carriers of the slope |
| B5b | `B5b_SHMI_yield_weights_sensitivity.R` | yield slope with between- or within-site weights |
| B6 | `B6_SHMI_N_scenarios.R` | 20 legume-credit x manure-availability assumptions |
| B7 | `B7_Assemble_manuscript_figures.R` | assembles Figs. 3-4 and Extended Data Fig. 4 |
| B9 | `B9_SHMI_coarse_records.R` | coarse-record test |

Helpers sourced by the scripts above: `run_helpers.R` (runs each step and logs it), `datasets.R` (calibration data sets), `SHMI_legume_credit.R` and `SHMI_manure_n.R` (N accounting for B2).

Python (`python/`): `napeshm_rebuild/` (A0), `agcros2shmi/` (B0, with the NECCIRR supplement and NEMLTCRS corrections), `build_soil_c.py` (B0), `build_study_treatments.py` (B8), `coarse_records/` (B9).

## Data in this repository (`data/`, `Inputs/`, `python/agcros2shmi/`)

- **`Inputs/napeshm/`, `Inputs/NAPESHM_rebuilt.xlsx`:** the NAPESHM files (https://doi.org/10.15482/USDA.ADC/25632270) and the quality-checked SHMI workbook built from them.
- **`Inputs/yield_workbooks/`:** SHMI workbooks for the ten experiments (three built by hand from DRIVES, seven converted from AgCROS).
- **`data/napeshm_yield_site_units.csv`:** NAPESHM units matched to the yield sites' treatments.
- **`data/soil_c.csv`:** soil organic carbon by plot, layer and date at the AgCROS experiments (output of `build_soil_c.py`), so the mediation test runs without the AgCROS tables.
- **`data/cache/nass_site_county_corn_yield.csv`:** USDA-NASS county corn yields for the yield sites, so no NASS key is needed.
- **`data/cache/site_soils_ssurgo.csv`:** gSSURGO soil groups by NAPESHM site, so gSSURGO is not needed.
- **`data/SHAPE_by_unit_with_SHMI.csv`:** SHAPE scores, measured soil health and SHMI by NAPESHM unit (output of A6), used when Part A is not run.
- **`data/SHMI_study_treatments_descriptions.xlsx`:** treatment descriptions for the hand-built sites.
- **`python/agcros2shmi/supplements/NECCIRR_amendments.csv`:** NECCIRR fertilizer and manure, missing from AgCROS (UNL CropWatch 2017).
- **`python/agcros2shmi/corrections/NEMLTCRS_corrections.csv`:** two record errors in AgCROS, with their evidence.

The raw AgCROS tables were used with permission and are not redistributed; the full AgCROS database is offline pending relaunch by USDA-ARS.

## Notes

- Scripts use paths relative to the project root (`here::here()`); any script can also be run on its own from the project.
- Site runs are cached in `Outputs/yield/site_runs`. The cache key covers each workbook, the period, the settings, the SHMI weights, the SHMI package (version and install time) and the scoring code. Delete the folder to force a full rescore.
- Figures are saved without titles or captions (the legends are in the manuscript); set `options(shmi.ms_figures = FALSE)` before running a script to keep them for working plots.
- The converter writes one `CD_term_date` (rule T1) and a zero-mixing placeholder pass at no-till sites (rule D3). Both were workarounds for two issues since fixed in the SHMI package. They do not change any score and are kept so the shipped workbooks are reproduced exactly.

## Licence and citation

Public domain (CC0 1.0; see `LICENSE.md`). The derived data keep the terms of their sources (NAPESHM, DRIVES, AgCROS, USDA-NASS, WorldClim, Daymet, SPEIbase).

If you use this code or data, please cite the paper and this archive:

Moore JM, Manter DK. Soil health tracks management systems, not the sum of individual practices. (in review)

Moore JM, Manter DK. (2026). *SHMI_Yield: analysis code and data*, version 1.0-submission. Zenodo. [DOI]
