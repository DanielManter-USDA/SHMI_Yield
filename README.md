# SHMI_Yield

[![DOI](https://zenodo.org/badge/1340463740.svg)](https://doi.org/10.5281/zenodo.24432549)

Code and derived data for Moore & Manter, *Scoring management as a system predicts soil health, an emergent property of agricultural soils* (in review):

- a structural equation model of soil health (lv_SH) from eleven indicators measured on 354 US management systems at 74 NAPESHM sites;
- calibration and leave-site-out cross-validation of the Soil Health Management Index (SHMI) against lv_SH, with SHAPE as a second soil-health outcome;
- how the SHMI principles overlap (path model, commonality analysis), SHMI against each single-practice proxy, a replicate-based ceiling on explainable within-site variation, and random forests of the same record inputs as a flexible reference;
- the relationship of SHMI with corn and soybean yield in ten independent long-term experiments (DRIVES and AgCROS), its dependence on nitrogen, time under management and growing conditions, and its fertilizer-N equivalent.

SHMI scores are computed with the SHMI R package, version 2.0.0: https://github.com/DanielManter-USDA/SHMI, archived at https://doi.org/10.5281/zenodo.23682170.

## Running everything

1. Open `SHMI_Yield.Rproj` in RStudio.
2. Run `install_packages.R` (includes `targets`), then install SHMI 2.0.0 (`remotes::install_github("DanielManter-USDA/SHMI@v2.0.0")`). Python 3.9+ (`pip install pandas openpyxl numpy`) is needed only to rebuild workbooks (A0, B0) and the treatment table (B8).
3. Everything needed is in the repository except the raw AgCROS tables (below). The SPEI grid (`Inputs/climate/spei03.nc`, see `Inputs/README.md`) is optional: without it, the site values in `data/cache/site_spei.csv` are used.
4. **Parts B and C (yield, soil-health comparisons, display items): use the `targets` pipeline.**

   ```r
   targets::tar_make()        # runs only the steps whose code or inputs changed
   targets::tar_outdated()    # what would run next
   targets::tar_visnetwork()  # the dependency graph
   ```

   Analyses save their results and the data behind each figure (`plotdata_*.rds`, `Outputs/panels/`); the figure scripts in `R/figures/` only read those, so editing a figure reruns that figure and the collection step, nothing else. See `README_pipeline.md` for one-time setup, accepting existing outputs (`SHMI_KEEP`) and adding steps or figures.
5. **Part A (calibration, frozen for the paper) and the conversions: use `run_all.R`.** The pipeline reads Part A's saved outputs and stops with a message if one is missing. To rerun Part A:

   ```r
   RUN_SETTINGS <- list(run_calibration = TRUE)
   source("run_all.R")
   ```

   `run_all.R` can also run everything end to end; each step's log goes to `Outputs/logs/<time>/`, with a table of steps in `steps.csv`.

| `run_all.R` setting | Default | What it runs | Time |
|---|---|---|---|
| `run_napeshm_rebuild` | FALSE | A0: NAPESHM workbook from the public NAPESHM tables, to compare with the shipped (QC'd) workbook | minutes |
| `run_calibration` | FALSE | A1-A8: soil-health model, SHMI calibration and cross-validation, benchmark, fixed-value grid, SHAPE | hours (frozen for the paper) |
| `run_conversion` | FALSE | B0: AgCROS workbooks from the raw AgCROS tables, which are not redistributed (see below). Shows and checks how the shipped workbooks were made | minutes |
| `run_yield` | TRUE | B1-B4, B8: yield analyses and the treatment table | first run about an hour (site runs are cached) |
| `run_n_scenarios` | TRUE | B6: SHMI slope under 20 N-accounting assumptions | 30-60 min |
| `run_record_length` | TRUE | B2c: rolling (current-season) scores for B2d | about an hour the first time; cached after |
| `run_yield_parts` | TRUE | B2d: SHMI against its parts for corn yield (Table 2) | 20-40 min (parallel; `B2D_WORKERS`) |

## Where each result comes from

Every display item ends up in **`Outputs/Manuscript/numbered/`** under the paper's numbering, with `manifest_numbered.csv` listing each file, its source and whether it is present. C3 copies the items there from where the scripts write them; the Supplementary Tables are combined into one workbook, `Supplementary_Tables.xlsx`, with a Contents sheet giving each title and its notes.

| In the paper | File in `numbered/` | Computed by | Drawn or assembled by |
|---|---|---|---|
| Fig. 1 (a, soil-health model; b, how SHMI scores a system: path diagram) | `Fig1.png` | A1c, A8; A3a, A4a, C1 | `R/figures/F1_soil_health_model.R` |
| Fig. 2 (SHMI against lv_SH: all systems, between sites, within sites) | `Fig2.png` | C1h | `R/figures/F2_system_score.R` |
| Fig. 3 (a, corn and soybean; b, each experiment's slope) | `Fig3.png` | B2, B2b, B2i | `R/figures/F3_yield.R` |
| Fig. 4 (SHMI × N; SHMI × years under management; poor, typical and good years) | `Fig4.png` | B2, B3 (SPEI from B1) | `R/figures/F4_yield_moderators.R` |
| Table 1 (scores by system type) | `Table1.csv` | C2 | |
| Table 2 (SHMI against single-practice proxies: lv_SH, SHAPE, corn yield) | `Table2.csv` | C1j, B2d; C2 | |
| Supplementary Table 1, the 74 NAPESHM sites | `Supplementary_Tables.xlsx` [Table 1] | C1m | C3 |
| Supplementary Table 2, the ten experiments | [Table 2] | B2, B2i; C2 | C3 |
| Supplementary Table 3, every treatment | [Table 3] | B8 (`python/build_study_treatments.py`), B8b | C3 |
| Supplementary Table 4, the structural equation model | [Table 4] | A1c; C1m | C3 |
| Supplementary Table 5, commonality analysis | [Table 5] | C1m | C3 |
| Supplementary Table 6, N-accounting alternatives | [Table 6] | B6; C1m | C3 |
| Results: out-of-fold R² (weights refitted in each fold), sites ranked consistently | | A3a, A4a | |
| Results: replicate ceiling on explainable within-site variation | | C1k | |
| Results and Supplementary Methods: within-site random forests | | C1l (and A4a for the unconstrained forest) | |
| Results: fertilizer-N equivalent | | B4 | |
| Discussion: share of systems usable by strict single-practice pairs | | C1f | |
| Methods: fixed calibration values (Diversity weight, livestock and amendment shares) | | A5a | |
| Supplementary Methods: climate source check | | A7 | |

## Scripts (`R/`)

Scripts are numbered by part: **A** = calibration, **B** = yield, **C** = soil-health comparisons, tables and collection; `R/figures/` holds the figure scripts that read saved plot data and panels. Every step runs with SHMI 2.0.0.

| Step | Script | What it does |
|---|---|---|
| A1a | `A1a_Get_site_ssurgo.R` | gSSURGO soil groups by NAPESHM site (cached in `data/cache/`) |
| A1b | `A1b_Prepare_SEM.R` | SEM data, cropping-system classes (annual, perennial, reference, fallow) |
| A1c | `A1c_soil_health_SEM.R` | outlier screening, climate and texture adjustment, the soil-health model (lv_SH), and the projection used to score field replicates |
| A1d | `A1d_Get_site_climate.R` | WorldClim monthly normals by site |
| A2c | `A2c_Project_SH.R` | the US soil-health model applied to the Mexican sites (read by A3a; not in the paper) |
| A3a | `A3a_SHMI_calibration_plus.R` | final calibration (weights, growing-season share); out-of-fold predictions |
| A4a | `A4a_SHMI_benchmark.R` | out-of-fold performance with the weights refitted in each fold; SHMI against its parts; unconstrained reference models |
| A5a | `A5a_Run_within_calibration_scenarios.R` | the fixed-value grid (Diversity weight, livestock and amendment shares) |
| A6a-A6b | `A6a_Prepare_SHAPE.R`, `A6b_Merge_SHAPE_SHMI.R` | SHAPE scores by system |
| A7 | `A7_Check_climate_source.R` | WorldClim against Daymet growing seasons |
| A8 | `A8_Export_appendix_tables.R` | SEM standardized loadings for Fig. 1a; adjustment slopes |
| B1 | `B1_Get_SPEI.R` | growing-season SPEI for every site-year |
| B2 | `B2_SHMI_yield_main.R` | within-site corn-yield model, SHMI × N, SHMI × age, site slopes, residual screening; saves the model (`B2_model.rds`) and plot data |
| B2b | `B2b_SHMI_soybean.R` | the same test for soybean, compared with corn as % of mean yield |
| B2c | `B2c_SHMI_record_length.R` | rolling (current-season) scores for B2d |
| B2d | `B2d_SHMI_yield_parts.R` | corn yield: SHMI against each part, each experiment left out in turn (Table 2) |
| B2i | `B2i_SHMI_variance_shares.R` | share of each experiment's within-site SHMI variance from each part |
| B3 | `B3_SHMI_yield_stability.R` | county-yield years, SPEI, crop-insurance years, worst years |
| B4 | `B4_SHMI_N_equivalent.R` | SHMI in fertilizer-N units (reads `B2_model.rds`) |
| B6 | `B6_SHMI_N_scenarios.R` | 20 legume-credit × manure-availability assumptions |
| B8b | `B8b_study_treatments_with_yield.R` | treatment table with SHMI, parts, N and mean yields (after `python/build_study_treatments.py`) |
| C1 | `C1_Fig1_path.R` | Fig. 1b: path diagram among the sub-indices, with site-bootstrap intervals |
| C1f | `C1f_paired_data_budget.R` | share of systems usable by strict single-practice pairs |
| C1h | `C1h_SHMI_vs_soil_health.R` | SHMI against lv_SH at three scales; SHMI and SHAPE |
| C1j | `C1j_SHMI_vs_parts.R` | Table 2, soil health: ΔAICc and leave-one-site-out within-site R² for SHMI and each part |
| C1k | `C1k_ceiling.R` | ceiling on explainable within-site variation from the field replicates |
| C1l | `C1l_rf_within.R` | random forests fitted to within-site variation of the same record inputs, each site left out in turn |
| C1m | `C1m_supp_tables.R` | Supplementary Tables 1, 4, 5 and 6 |
| C2 | `C2_export_tables.R` | Tables 1-2 and Supplementary Table 2 |
| C3 | `C3_collect_figures.R` | copies every display item into `Outputs/Manuscript/numbered/` and builds `Supplementary_Tables.xlsx` |
| C4 | `C4_bundle_for_review.R` | (`run_all.R` only) one zip of the small results files, for checking the text |
| F1–F4 | `figures/F1_soil_health_model.R` … `figures/F4_yield_moderators.R` | Figs 1–4 from plot data and panels; `figures/yield_panels.R` and `figures/moderator_panels.R` build the yield panels; shared settings, including text size (`FONT_SCALE`), in `figures/fig_common.R`; `figures/run_figures.R` redraws all |

Pipeline: `_targets.R` (steps and dependencies), `targets_helpers.R` (running and logging each step). Helpers sourced by the scripts: `run_helpers.R`, `datasets.R` (calibration data sets), `SHMI_legume_credit.R` and `SHMI_manure_n.R` (N accounting for B2).

Python (`python/`): `napeshm_rebuild/` (A0), `agcros2shmi/` (B0, with the NECCIRR and IAAM7071 supplements and the NEMLTCRS corrections), `build_study_treatments.py` (B8); `ohio_meadow_species.py` and `wicst_cs5_pea.py` record two one-off edits to the shipped workbooks.

Only the analyses in the paper are included. Scripts for analyses tried during the project but not reported (practice labels against SHMI, practice additivity, a strict-pair meta-analysis, a coarse-record test, weights fitted to yield comparisons, yield pairs) were removed from this version; they remain in the repository's git history.

## Data in this repository (`data/`, `Inputs/`, `python/agcros2shmi/`)

- **`Inputs/napeshm/`, `Inputs/NAPESHM_rebuilt.xlsx`:** the NAPESHM files (https://doi.org/10.15482/USDA.ADC/25632270), including the field-replicate soil measurements, and the quality-checked SHMI workbook built from them.
- **`Inputs/yield_workbooks/`:** SHMI workbooks for the ten experiments (three built by hand from DRIVES, seven converted from AgCROS).
- **`data/cache/nass_site_county_corn_yield.csv`:** USDA-NASS county corn yields for the yield sites, so no NASS key is needed.
- **`data/cache/site_soils_ssurgo.csv`:** gSSURGO soil groups by NAPESHM site, so gSSURGO is not needed.
- **`data/cache/site_spei.csv`:** SPEI-3 at each yield site, so the SPEI grid is not needed.
- **`data/SHAPE_by_unit_with_SHMI.csv`:** SHAPE scores, measured soil health and SHMI by NAPESHM system (output of A6a-A6b), used when Part A is not run.
- **`python/agcros2shmi/supplements/`:** NECCIRR fertilizer and manure (missing from AgCROS; UNL CropWatch 2017) and IAAM7071 cover-crop dates (supplied by the site).
- **`python/agcros2shmi/corrections/NEMLTCRS_corrections.csv`:** two record errors in AgCROS, with their evidence.

The raw AgCROS tables were used with permission and are not redistributed; the full AgCROS database is offline pending relaunch by USDA-ARS.

## Notes

- Scripts use paths relative to the project root (`here::here()`); any script can also be run on its own from the project.
- Random seeds are fixed in every script that resamples (calibration folds and bootstraps, the benchmark, the site bootstraps, random forests), so a rerun reproduces the outputs with the same R and package versions (Supplementary Methods). Each run writes `sessionInfo.txt` to its log folder. When Part A is run, `run_all.R` checks that the recalibrated weights match those built into SHMI 2.0.0.
- Site runs are cached in `Outputs/yield/site_runs`. The cache key covers each workbook, the period, the settings, the SHMI weights, the SHMI package (version and install time) and the scoring code. Delete the folder to force a full rescore.
- Figures are saved without titles or captions (the legends are in the manuscript).
- The converter writes one `CD_term_date` (rule T1) and a zero-mixing placeholder pass at no-till sites (rule D3). Both were workarounds for two issues since fixed in the SHMI package. They do not change any score and are kept so the shipped workbooks are reproduced exactly.

## Licence and citation

Public domain (CC0 1.0; see `LICENSE.md`). The derived data keep the terms of their sources (NAPESHM, DRIVES, AgCROS, USDA-NASS, WorldClim, Daymet, SPEIbase).

If you use this code or data, please cite the paper and this archive (`CITATION.cff`). Archived versions: https://doi.org/10.5281/zenodo.24432549 (all versions).

Moore JM, Manter DK. Scoring management as a system predicts soil health, an emergent property of agricultural soils. (in review)

Moore JM, Manter DK. SHMI_Yield: analysis code and data. Zenodo. https://doi.org/10.5281/zenodo.24432549
