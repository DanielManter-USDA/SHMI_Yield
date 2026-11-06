# SHMI_Yield

Code and derived data for Moore & Manter, *Soil health tracks management systems, not the sum of individual practices* (in review):

- calibration of the Soil Health Management Index (SHMI) against measured soil health on 354 US systems at 74 NAPESHM sites, with the SHAPE benchmark;
- its relationship with corn yield in ten long-term experiments (DRIVES and AgCROS);
- how much management detail a record needs.

SHMI scores are computed with the SHMI R package, version 2.0.0: https://github.com/DanielManter-USDA/SHMI, archived at https://doi.org/10.5281/zenodo.23682170.

## Running everything

1. Open `SHMI_Yield.Rproj` in RStudio.
2. Run `install_packages.R`, then install SHMI 2.0.0 (`remotes::install_github("DanielManter-USDA/SHMI@v2.0.0")`). Python 3.9+ (`pip install pandas openpyxl numpy`) is needed only to rebuild workbooks (A0, B0, B8) and for the coarse-record test (B9).
3. Everything needed is in the repository. The SPEI grid (`Inputs/climate/spei03.nc`, see `Inputs/README.md`) is optional: without it, the site values in `data/cache/site_spei.csv` are used.
4. Run `source("run_all.R")`. To choose parts, set `RUN_SETTINGS` first (any setting at the top of `run_all.R`), for example:

   ```r
   RUN_SETTINGS <- list(run_calibration = TRUE, run_weights = TRUE)
   source("run_all.R")
   ```

5. To see the results, render the report, which only reads `Outputs/` and takes seconds:

   ```r
   rmarkdown::render("report.Rmd", output_dir = "Outputs")
   ```

`run_all.R` first checks the inputs. It then runs each step in a fresh R session (python for the conversions). The console shows one line per step (start, done or FAILED, minutes); each step's full output goes to its own log in `Outputs/logs/<time>/`, with a table of steps in `steps.csv`. If a step fails, the last lines of its log are printed. On Windows, if a step fails with "callr subprocess failed", rerun with `RUN_SETTINGS <- list(in_session = TRUE)`: each step then runs in the current session, still with its own log. B4 and B5 refit the B2 model in their own sessions, so each takes about as long as B2.

| Setting | Default | What it runs | Time |
|---|---|---|---|
| `run_napeshm_rebuild` | FALSE | A0: NAPESHM workbook from the public NAPESHM tables, to compare with the shipped (QC'd) workbook | minutes |
| `run_calibration` | FALSE | A1-A8: soil health model, SHMI calibration, benchmark, weights by comparison, SHAPE | hours (frozen for the paper) |
| `run_conversion` | FALSE | B0: AgCROS workbooks and soil C from the raw AgCROS tables, which are not redistributed (see below). Shows and checks how the shipped workbooks were made | minutes |
| `run_yield` | TRUE | B1-B5, B7-B8: yield analyses, figures, treatment table | first run about an hour (site runs are cached) |
| `run_weights` | FALSE | B5b: yield slope with weights fitted between or within sites | reruns B2 three times |
| `run_n_scenarios` | TRUE | B6: SHMI slope under 20 N-accounting assumptions | 30-60 min |
| `run_coarse` | FALSE | B9: SHMI from less detailed records | hours |
| `run_record_length` | TRUE | B2c: yield slope by length of record (its current-season score feeds B2d; no longer a figure) | about an hour the first time; cached after |
| `run_yield_parts` | TRUE | B2d: SHMI against its parts for yield (Extended Data Fig. 4b) | 20-40 min |

## Where each result comes from

Every display item ends up in **`Outputs/Manuscript/numbered/`** under the paper's numbering (manuscript v5, 6 Nov 2026, items in order of first citation), with `manifest_numbered.csv` listing each file, its source and whether it is present. Only this folder is shared in git. The scripts that draw items write their own copies, under their earlier names, to `Outputs/Manuscript_work/` (not shared); C3 copies them into `numbered/` and moves anything left directly in `Outputs/Manuscript/` into the working folder.

| In the paper | File in `numbered/` | Computed by | Source file |
|---|---|---|---|
| Fig. 1 (a, path diagram; b, matched practice pairs; c, practices in measured soil health; d, e, how much data matched pairs can use) | `Fig1.png` | C1 `C1_Fig1_and_Table1.R` (c from C1d, d and e from C1f) | `Manuscript_work/Fig1.png` |
| Fig. 2 | `Fig2.png` | A3a, A4a; drawn by A4b `A4b_Fig2_soil_health.R` | `Manuscript_work/Fig2.png` |
| Fig. 3 (a, within sites; b, SHMI × N; c, SHMI × age) | `Fig3.png` | B2 `B2_SHMI_yield_main.R`; assembled by B7 (`FIG3_RECORD_LENGTH <- TRUE` restores the record-length panel) | `manuscript_figures_default_total_agcros/Fig3.png` |
| Fig. 4 | `Fig4.png` | B3 `B3_SHMI_yield_stability.R` (SPEI from B1); assembled by B7 | `manuscript_figures_default_total_agcros/Fig4.png` |
| Table 1 | `Table1.csv` | C1 | `Manuscript_work/Table1.csv` |
| Table 2 | (text) | | |
| Extended Data Fig. 1, soil health model | `ED_Fig1_loadings.csv` (diagram drawn by hand) | A1c, A8 | `sem_noMX_standardized.csv` |
| Extended Data Fig. 2, calibrated weights | `ED_Fig2.png` | A3a; drawn by C1b `C1b_ED_Fig2.R` | `Manuscript_work/ED_Fig2.png` |
| Extended Data Fig. 3, practice labels vs SHMI | `ED_Fig3.png` | C1e `C1e_labels_vs_SHMI.R` | `Practice_additivity/Fig_labels_vs_SHMI.png` |
| Extended Data Fig. 4, SHAPE | `ED_Fig4.png` | A6a-A6c | `SHAPE_age/FigA_SHMI_vs_SHAPE.png` |
| Extended Data Fig. 5, the experiments, their own slopes and SHMI against its parts for yield | `ED_Fig5.png` | B2 (section 7b) and B2d; assembled by B7 | `manuscript_figures_default_total_agcros/ED_Fig4_sites.png` |
| Extended Data Fig. 6, four record types | `ED_Fig6.png` | B9 with `python/coarse_records/degrade_records.py`; figure C1c `C1c_ED_Fig5.R` | `Manuscript_work/ED_Fig5.png` |
| Extended Data Table 1, scores by system type | `ED_Table1.csv` | C2 `C2_export_tables.R` | `Manuscript_work/ED_Table1.csv` |
| Supplementary Data Table 1, SHMI against its variants and parts | `SuppTable1.csv` | A4a; table C2 | `Manuscript_work/SuppTable4.csv` |
| Supplementary Data Table 2, soil properties as carriers of the slope | `SuppTable2.csv` | B5 `B5_SHMI_soil_mediation.R`; table C2 | `Manuscript_work/SuppTable3.csv` |
| Supplementary Data Table 3, the ten experiments | `SuppTable3.csv` | B2 (section 7b), B8; table C2 | `Manuscript_work/SuppTable1.csv` |
| Supplementary Data Table 4, every coarse-record scenario (supplied as a spreadsheet) | `SuppTable4.csv` | B9; table C2 | `Manuscript_work/SuppTable2.csv` |
| Supplementary Data Fig. 1, soybean | `SuppFig1.png` | B2b `B2b_SHMI_soybean.R` | `yield_soybean/Fig_soybean.png` |
| Supplementary Data Fig. 2, matched pairs in the long-term experiments | `SuppFig2.png` | B2f `B2f_SHMI_yield_paired_budget.R` | `yield_paired_budget/Fig_yield_budget.png` |
| Supplementary Data Fig. 3, weights by comparison | `SuppFig3.png` | A5a; figure A5b | `Calib_scenarios/SuppFig2_weights_by_comparison.png` |
| Supplementary Data Fig. 4, effect sizes | `SuppFig4.png` | B2 | `yield_main_default_total_agcros/Fig_effect_sizes.png` |
| Methods: fertilizer-N equivalent | | B4 `B4_SHMI_N_equivalent.R` | `N_equivalent.csv` |
| Methods: yield slope and the choice of weights | | B5b `B5b_SHMI_yield_weights_sensitivity.R` | `yield_weights_agcros/weights_slopes.csv` |
| Methods: N-accounting range | | B6 `B6_SHMI_N_scenarios.R` | `N_scenarios_all_sites.csv` |
| Methods: climate source check | | A7 `A7_Check_climate_source.R` | season check |

Part C runs at the end of every run and takes a minute or two: C1d, C1, C1b, C1c, C2 and C3. When Part A is off, its drawing steps (A4b, A5b, A6c, A8) are rerun from the saved Part A outputs, so every display item is redrawn by the current code. When a part was not run (for example Part A, which is off by default), its items keep the copies already in `Outputs/Manuscript/`, and `manifest.csv` shows when each was written.

## Scripts (`R/`)

Scripts are numbered in the order `run_all.R` runs them: **A** = calibration (Part A), **B** = yield and records (Part B), **C** = display items for the manuscript (Part C). Every step runs with SHMI 2.0.0. Alternative index structures were compared during development with an earlier build of the package; those scripts are not part of this release (A2a and A2b are therefore absent).

| Step | Script | What it does |
|---|---|---|
| A1a | `A1a_Get_site_ssurgo.R` | gSSURGO soil groups by NAPESHM site (cached in `data/cache/`) |
| A1b | `A1b_Prepare_SEM.R` | SEM data, cropping-system classes (annual, perennial, reference, fallow) |
| A1c | `A1c_soil_health_SEM.R` | climate and texture adjustment; soil health latent variable (lv_SH) |
| A1d | `A1d_Get_site_climate.R` | WorldClim monthly normals by site |
| A2c | `A2c_Project_SH.R` | the US soil health model applied to the Mexican sites (the independent test in A3a) |
| A3a | `A3a_SHMI_calibration_plus.R` | final calibration (three sub-indices plus the amendment bonus); out-of-fold predictions |
| A3b | `A3b_calibration_report_figures.R` | in-sample and out-of-fold report figures and agreement table |
| A4a | `A4a_SHMI_benchmark.R` | SHMI against its parts, reduced versions, fixed published weights and unconstrained models |
| A4b | `A4b_Fig2_soil_health.R` | Fig. 2 |
| A5a-A5b | `A5a_Run_within_calibration_scenarios.R`, `A5b_SuppFig2_weights_by_comparison.R` | weights fitted between sites, to all variation and within sites |
| A6a-A6c | `A6a_Prepare_SHAPE.R`, `A6b_Merge_SHAPE_SHMI.R`, `A6c_SHMI_SHAPE_age.R` | SHAPE benchmark |
| A7 | `A7_Check_climate_source.R` | WorldClim against Daymet growing seasons |
| A8 | `A8_Export_appendix_tables.R` | SEM loadings and adjustment slopes |
| B0 | `B0_build_napeshm_crosswalk.R` | NAPESHM units matched to the yield sites' treatments |
| B1 | `B1_Get_SPEI.R` | growing-season SPEI for every site-year |
| B2 | `B2_SHMI_yield_main.R` | within-site yield model, SHMI x N, SHMI x age, site slopes, heterogeneity, effect sizes |
| B2c | `B2c_SHMI_record_length.R` | yield slope with SHMI scored from the last 1-4 crop years before each harvest vs the full record (no longer a figure; its current-season score feeds B2d) |
| B2f | `B2f_SHMI_yield_paired_budget.R` | how much of the yield record matched single-practice pairs can use, against SHMI; paired yield differences against the change in SHMI; the Fig. 3a model on treatments in and outside pairs (`Outputs/yield_paired_budget/`) |
| B2d | `B2d_SHMI_yield_parts.R` | yield: SHMI against each part alone, equal weights and the current season, with each experiment left out in turn; `B2D_ROTATION <- TRUE` adds the preceding legume (Extended Data Fig. 4b) |
| B2e | `B2e_SHMI_yield_adjustment.R` | what the adjustment does: raw yields, site and year removed, then N, other treatments and age (the model of Fig. 3a), on one yield scale |
| B2b | `B2b_SHMI_soybean.R` | the same test for soybean at the experiments that grew it, compared with corn as % of mean yield |
| B3 | `B3_SHMI_yield_stability.R` | year conditions, drought, worst years |
| B4 | `B4_SHMI_N_equivalent.R` | SHMI in fertilizer-N units |
| B5 | `B5_SHMI_soil_mediation.R` | soil properties as carriers of the slope |
| B5b | `B5b_SHMI_yield_weights_sensitivity.R` | yield slope with between- or within-site weights |
| B6 | `B6_SHMI_N_scenarios.R` | 20 legume-credit x manure-availability assumptions |
| B7 | `B7_Assemble_manuscript_figures.R` | assembles Figs. 3-4 and Extended Data Fig. 4 |
| B9 | `B9_SHMI_coarse_records.R` | coarse-record test |
| C1 | `C1_Fig1_and_Table1.R` | Fig. 1 (path coefficients with site-bootstrap intervals; example systems) and Table 1 |
| C1b | `C1b_ED_Fig2.R` | Extended Data Fig. 2 (weights and active-season share), from the saved calibration fits |
| C1c | `C1c_ED_Fig5.R` | Extended Data Fig. 5 (the four record types), from B9's saved results |
| C1d | `C1d_practice_additivity.R` | do practices add up in measured soil health? within-site model of lv_SH on four practice labels and their interactions (Fig. 1c; runs before C1) |
| C1f | `C1f_paired_data_budget.R` | how much of the soil-health data matched practice pairs (Fig. 1b rule) can use, against SHMI; SHMI on the systems in and outside pairs and on identical-checklist pairs (`Outputs/Practice_additivity/data_budget*.csv`, `Fig_data_budget.png`) |
| C1e | `C1e_labels_vs_SHMI.R` | practice labels (with or without interactions) vs SHMI as predictors of soil health at held-out sites, same systems and folds |
| C2 | `C2_export_tables.R` | Extended Data Table 1 and Supplementary Data Tables 1-4; `Tables.xlsx` |
| C3 | `C3_collect_figures.R` | copies every display item into `Outputs/Manuscript/numbered/` under the paper's numbering (`manifest_numbered.csv`); keeps `Outputs/Manuscript/` to that folder |

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
- Random seeds are fixed in every script that resamples (calibration bootstraps and folds, the benchmark, the site bootstraps, random forests), so a rerun reproduces the outputs exactly with the same R and package versions. Each run writes `sessionInfo.txt` to its log folder. When Part A is run, `run_all.R` checks that the recalibrated weights match those built into SHMI 2.0.0 and says so in the console.
- Site runs are cached in `Outputs/yield/site_runs`. The cache key covers each workbook, the period, the settings, the SHMI weights, the SHMI package (version and install time) and the scoring code. Delete the folder to force a full rescore.
- Figures are saved without titles or captions (the legends are in the manuscript); set `options(shmi.ms_figures = FALSE)` before running a script to keep them for working plots.
- The converter writes one `CD_term_date` (rule T1) and a zero-mixing placeholder pass at no-till sites (rule D3). Both were workarounds for two issues since fixed in the SHMI package. They do not change any score and are kept so the shipped workbooks are reproduced exactly.

## Licence and citation

Public domain (CC0 1.0; see `LICENSE.md`). The derived data keep the terms of their sources (NAPESHM, DRIVES, AgCROS, USDA-NASS, WorldClim, Daymet, SPEIbase).

If you use this code or data, please cite the paper and this archive:

Moore JM, Manter DK. Soil health tracks management systems, not the sum of individual practices. (in review)

Moore JM, Manter DK. (2026). *SHMI_Yield: analysis code and data*, version 1.0-submission. Zenodo. [DOI]
