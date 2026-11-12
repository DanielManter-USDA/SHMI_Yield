# SHMI_Yield

Code and derived data for Moore & Manter, *Soil health tracks management systems, not the sum of individual practices* (in review):

- calibration of the Soil Health Management Index (SHMI) against measured soil health on 354 US systems at 74 NAPESHM sites, with the SHAPE benchmark;
- its relationship with corn yield in ten long-term experiments (DRIVES and AgCROS);
- how much management detail a record needs.

SHMI scores are computed with the SHMI R package, version 2.0.0: https://github.com/DanielManter-USDA/SHMI, archived at https://doi.org/10.5281/zenodo.23682170.

## Running everything

1. Open `SHMI_Yield.Rproj` in RStudio.
2. Run `install_packages.R` (includes `targets`), then install SHMI 2.0.0 (`remotes::install_github("DanielManter-USDA/SHMI@v2.0.0")`). Python 3.9+ (`pip install pandas openpyxl numpy`) is needed only to rebuild workbooks (A0, B0, B8) and for the coarse-record test (B9).
3. Everything needed is in the repository. The SPEI grid (`Inputs/climate/spei03.nc`, see `Inputs/README.md`) is optional: without it, the site values in `data/cache/site_spei.csv` are used.
4. **Parts B and C (yield, display items): use the `targets` pipeline.**

   ```r
   targets::tar_make()        # runs only the steps whose code or inputs changed
   targets::tar_outdated()    # what would run next
   targets::tar_visnetwork()  # the dependency graph
   ```

   Analyses save their results and the data behind each figure (`plotdata_*.rds`); the figure scripts in `R/figures/` only read those, so editing a figure reruns that figure and the collection step, nothing else. See `README_pipeline.md` for one-time setup, the slow coarse-record step (B9, `SHMI_FORCE`), accepting existing outputs (`SHMI_KEEP`) and adding steps or figures.
5. **Part A (calibration, frozen for the paper) and the conversions: use `run_all.R`.** The pipeline reads Part A's saved outputs and stops with a message if one is missing. To rerun Part A:

   ```r
   RUN_SETTINGS <- list(run_calibration = TRUE)
   source("run_all.R")
   ```

   `run_all.R` can still run everything end to end (any setting at the top of the file); each step's log goes to `Outputs/logs/<time>/`, with a table of steps in `steps.csv`.
6. To see the results, render the report, which only reads `Outputs/` and takes seconds: `rmarkdown::render("report.Rmd", output_dir = "Outputs")`.

| `run_all.R` setting | Default | What it runs | Time |
|---|---|---|---|
| `run_napeshm_rebuild` | FALSE | A0: NAPESHM workbook from the public NAPESHM tables, to compare with the shipped (QC'd) workbook | minutes |
| `run_calibration` | FALSE | A1-A8: soil health model, SHMI calibration, benchmark, weights by comparison, SHAPE | hours (frozen for the paper) |
| `run_conversion` | FALSE | B0: AgCROS workbooks from the raw AgCROS tables, which are not redistributed (see below). Shows and checks how the shipped workbooks were made | minutes |
| `run_yield` | TRUE | B1-B4, B8: yield analyses, figures, treatment table | first run about an hour (site runs are cached) |
| `run_weights` | TRUE | B5b: yield slope with weights fitted between or within sites | reruns B2 three times |
| `run_n_scenarios` | TRUE | B6: SHMI slope under 20 N-accounting assumptions | 30-60 min |
| `run_coarse` | FALSE | B9: SHMI from less detailed records | hours |
| `run_record_length` | TRUE | B2c: rolling (current-season) scores for B2d; `RL_FULL <- TRUE` adds the slope-by-record-length analysis (not in the paper) | about an hour the first time; cached after |
| `run_yield_parts` | TRUE | B2d: SHMI against its parts for yield (Supplementary Table 1) | 20-40 min (parallel; `B2D_WORKERS`) |

## Where each result comes from

Every display item ends up in **`Outputs/Manuscript/numbered/`** under the paper's numbering (manuscript of 11 Nov 2026), with `manifest_numbered.csv` listing each file, its source and whether it is present. Only this folder is shared in git. C3 copies the items there from where the scripts write them.

| In the paper | File in `numbered/` | Computed by | Drawn by |
|---|---|---|---|
| Fig. 1 (a, SEM paths; b, SHMI against measured soil health) | `Fig1.png` | A1c, A3a; C1h | C1 |
| Fig. 2 (a, practice labels vs SHMI at held-out sites; b, how much data matched pairs can use; c, SHMI against its parts) | `Fig2.png` | C1g (weights refitted in each fold), C1f, A4a/A4b | C1 |
| Fig. 3 (a, corn within experiments; b, soybean; c, poor, dry and good years) | `Fig3.png` | B2, B2b, B3 (SPEI from B1) | `R/figures/F3_yield.R` |
| Fig. 4 (a, SHMI × N; b, SHMI × age; c, each experiment's slope) | `Fig4.png` | B2 | `R/figures/F4_yield_moderators.R` |
| Extended Data Fig. 1, soil health model | `ED_Fig1_loadings.csv` (diagram drawn by hand) | A1c, A8 | |
| Extended Data Fig. 2, practices do not add up | `ED_Fig2.png` | C1d | C1 |
| Extended Data Fig. 3, SHMI within practice-label categories | `ED_Fig3.png` | C1g | C1g |
| Extended Data Fig. 4, SHMI and SHAPE | `ED_Fig4.png` | A6a-A6b, C1h | C1h |
| Extended Data Fig. 5, four record types | `ED_Fig5.png` | B9 (`python/coarse_records/`) | C1c |
| Table 1 | `Table1.csv` | C1 | |
| Extended Data Table 1, scores by system type | `ED_Table1.csv` | C2 | |
| Supplementary Table 1, SHMI against its variants and parts (soil health and corn yield) | `Supplementary_Tables.xlsx` [Table 1] | A4a, B2d; C2, C3 | |
| Supplementary Table 2, the ten experiments | [Table 2] | B2, B8; C2 | |
| Supplementary Table 3, slope by level of the factors SHMI does not score | [Table 3] | B2 | |
| Supplementary Table 4, every treatment | [Table 4] | B8 (`python/build_study_treatments.py`), B8b | |
| Supplementary Table 5, every coarse-record scenario | [Table 5] | B9; C2 | |
| Supplementary Table 6, matched practice pairs | [Table 6] | C1 (pairs from B2f for yield, C1f for soil health) | |
| Methods: fertilizer-N equivalent | | B4 | |
| Methods: N-accounting range | | B6 | |
| Methods: weights by comparison | | A5a, B5b | |
| Methods: current-season score | | B2c → B2d | |
| Methods: four practice labels (r 0.23 vs 0.40) | | C1e | |
| Methods: climate source check | | A7 | |

## Scripts (`R/`)

Scripts are numbered by part: **A** = calibration, **B** = yield and records, **C** = display items and tables; `R/figures/` holds the figure scripts that read saved plot data. Every step runs with SHMI 2.0.0. Analyses tried during the project but not in the paper are not part of this release; alternative index structures compared during development with an earlier build of the package are not part of this release.

| Step | Script | What it does |
|---|---|---|
| A1a | `A1a_Get_site_ssurgo.R` | gSSURGO soil groups by NAPESHM site (cached in `data/cache/`) |
| A1b | `A1b_Prepare_SEM.R` | SEM data, cropping-system classes (annual, perennial, reference, fallow) |
| A1c | `A1c_soil_health_SEM.R` | climate and texture adjustment; soil health latent variable (lv_SH) |
| A1d | `A1d_Get_site_climate.R` | WorldClim monthly normals by site |
| A2c | `A2c_Project_SH.R` | the US soil health model applied to the Mexican sites (the independent test in A3a) |
| A3a | `A3a_SHMI_calibration_plus.R` | final calibration (three sub-indices plus the amendment bonus); out-of-fold predictions |
| A4a | `A4a_SHMI_benchmark.R` | SHMI against its parts, reduced versions, fixed published weights and unconstrained models |
| A4b | `A4b_Fig2_soil_health.R` | SHMI against its parts at held-out sites (Fig. 2c panel) |
| A5a | `A5a_Run_within_calibration_scenarios.R` | weights fitted between sites, to all variation and within sites |
| A6a-A6b | `A6a_Prepare_SHAPE.R`, `A6b_Merge_SHAPE_SHMI.R` | SHAPE benchmark |
| A7 | `A7_Check_climate_source.R` | WorldClim against Daymet growing seasons |
| A8 | `A8_Export_appendix_tables.R` | SEM loadings and adjustment slopes |
| B1 | `B1_Get_SPEI.R` | growing-season SPEI for every site-year |
| B2 | `B2_SHMI_yield_main.R` | within-site yield model, SHMI × N, SHMI × age, site slopes, heterogeneity; saves the model (`B2_model.rds`) and plot data |
| B2b | `B2b_SHMI_soybean.R` | the same test for soybean, compared with corn as % of mean yield |
| B2c | `B2c_SHMI_record_length.R` | rolling (current-season) scores for B2d; `RL_FULL` adds the record-length analysis |
| B2d | `B2d_SHMI_yield_parts.R` | yield: SHMI against each part, equal weights and the current season, each experiment left out in turn |
| B2f | `B2f_SHMI_yield_paired_budget.R` | yield record usable by matched single-practice pairs; paired yield differences against the change in SHMI |
| B3 | `B3_SHMI_yield_stability.R` | year conditions, drought, worst years |
| B4 | `B4_SHMI_N_equivalent.R` | SHMI in fertilizer-N units (reads `B2_model.rds`) |
| B5b | `B5b_SHMI_yield_weights_sensitivity.R` | yield slope with between- or within-site weights |
| B6 | `B6_SHMI_N_scenarios.R` | 20 legume-credit × manure-availability assumptions |
| B8b | `B8b_study_treatments_with_yield.R` | treatment table with SHMI, parts, N and mean yields (after `python/build_study_treatments.py`) |
| B9 | `B9_SHMI_coarse_records.R` | coarse-record test |
| C1d | `C1d_practice_additivity.R` | do practices add up in measured soil health? (runs before C1) |
| C1e | `C1e_labels_vs_SHMI.R` | four practice labels vs SHMI at held-out sites |
| C1f | `C1f_paired_data_budget.R` | how much of the soil-health data matched practice pairs can use |
| C1g | `C1g_label_vs_SHMI.R` | practice labels vs SHMI out of fold, SHMI weights refitted in each fold (`REFIT_IN_FOLD`); SHMI within label categories |
| C1h | `C1h_SHMI_vs_soil_health.R` | SHMI against measured soil health; SHMI and SHAPE |
| C1 | `C1_Fig1_and_Table1.R` | Figs. 1-2, Extended Data Fig. 2, Table 1, Supplementary Table 6 |
| C1c | `C1c_ED_Fig5.R` | Extended Data Fig. 5 (four record types), from B9's saved results |
| C2 | `C2_export_tables.R` | Extended Data Table 1 and the Supplementary Tables |
| C3 | `C3_collect_figures.R` | copies every display item into `Outputs/Manuscript/numbered/` (`manifest_numbered.csv`) |
| C4 | `C4_bundle_for_review.R` | (`run_all.R` only) one zip of the small results files, for checking the text |
| F3, F4 | `figures/F3_yield.R`, `figures/F4_yield_moderators.R` | Figs. 3 and 4 from `plotdata_*.rds`; shared settings in `figures/fig_common.R`; `figures/run_figures.R` redraws both without `targets` |

Pipeline: `_targets.R` (steps and dependencies), `targets_helpers.R` (running and logging each step). Helpers sourced by the scripts: `run_helpers.R`, `datasets.R` (calibration data sets), `SHMI_legume_credit.R` and `SHMI_manure_n.R` (N accounting for B2).

Python (`python/`): `napeshm_rebuild/` (A0), `agcros2shmi/` (B0, with the NECCIRR supplement and NEMLTCRS corrections), `build_study_treatments.py` (B8), `coarse_records/` (B9); `ohio_meadow_species.py` and `wicst_cs5_pea.py` record two one-off edits to the shipped workbooks.

## Data in this repository (`data/`, `Inputs/`, `python/agcros2shmi/`)

- **`Inputs/napeshm/`, `Inputs/NAPESHM_rebuilt.xlsx`:** the NAPESHM files (https://doi.org/10.15482/USDA.ADC/25632270) and the quality-checked SHMI workbook built from them.
- **`Inputs/yield_workbooks/`:** SHMI workbooks for the ten experiments (three built by hand from DRIVES, seven converted from AgCROS).
- **`data/cache/nass_site_county_corn_yield.csv`:** USDA-NASS county corn yields for the yield sites, so no NASS key is needed.
- **`data/cache/site_soils_ssurgo.csv`:** gSSURGO soil groups by NAPESHM site, so gSSURGO is not needed.
- **`data/SHAPE_by_unit_with_SHMI.csv`:** SHAPE scores, measured soil health and SHMI by NAPESHM unit (output of A6a-A6b), used when Part A is not run.
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
