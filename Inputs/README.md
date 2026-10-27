# Inputs

Files the analysis reads but does not create. `run_all.R` checks this list before it starts and reports anything missing for the parts you chose to run.

**What is in the repository.** Everything needed to rerun the paper except the raw AgCROS tables:
- `yield_workbooks/`: the SHMI workbooks of all yield experiments (the AgCROS ones are derived from the AgCROS records);
- `napeshm/`: the NAPESHM tables, copied from the public NAPESHM release on Ag Data Commons (https://doi.org/10.15482/USDA.ADC/25632270);
- `NAPESHM_rebuilt.xlsx`: the quality-checked SHMI workbook built from them;
- `SHAPE/`: the SHAPE app's scores for the NAPESHM units.

**Not in the repository.**
- The raw AgCROS table export (`agcros/`): used with permission and not redistributed; the AgCROS database is offline pending relaunch by USDA-ARS. The analyses use the workbooks and soil carbon table made from it (`yield_workbooks/`, `data/soil_c.csv`). The conversion (`run_conversion = TRUE`, step B0) is kept to show how they were made, and reruns for anyone with the tables.
- Large public files, downloaded if wanted: the SPEI grid (`climate/spei03.nc`; its site values are shipped in `data/cache/site_spei.csv`) and gSSURGO.

Derived files that would otherwise need restricted, very large or keyed sources are shipped in `data/` and copied into `Outputs/` by `run_all.R` when missing: NASS county yields (`data/cache/`, so no NASS key is needed), gSSURGO soil groups (`data/cache/site_soils_ssurgo.csv`), AgCROS soil carbon (`data/soil_c.csv`) and the SHAPE scores by unit.

## Part B: corn yield (the default run)

| Path | What | Source |
|---|---|---|
| `agcros/extracted_tables/*.csv` (only for `run_conversion = TRUE`) | AgCROS (GRACEnet/REAP) table export, one CSV per table. Used: `meassoilchem`, `measyieldnutuptake`, `measharvestfraction`, `measresiduemgnt`, `mgtplanting`, `mgttillage`, `mgtamendments`, `mgtgrazing`, `metatreatments`, `metaunits`, `metasites` | USDA-ARS AgCROS; snapshot exported 2026-04-21 |
| `yield_workbooks/SHMI_OHTVD_Ho.xlsx` | Hoytville, OH: SHMI workbook built by hand from DRIVES, tillage corrected with the Ohio management logs | DRIVES |
| `yield_workbooks/SHMI_OHTVD_Wo.xlsx` | Wooster, OH: SHMI workbook built by hand from DRIVES | DRIVES |
| `yield_workbooks/SHMI_WICST_ARS.xlsx` | Arlington, WI (WICST): SHMI workbook built by hand from DRIVES and the station records | DRIVES; WICST data request (6-5-2026) |
| `climate/spei03.nc` | SPEIbase 3-month SPEI (global netCDF). Optional: without it, B1 uses the site values in `data/cache/site_spei.csv`, extracted from it on an earlier run | https://spei.csic.es/database.html |
| `NASS_API_KEY` in `~/.Renviron` | Key for USDA-NASS Quick Stats (county corn yields). Not a file. Needed only if `data/cache/nass_site_county_corn_yield.csv` is removed or new sites are added. | https://quickstats.nass.usda.gov/api |

**Made from the AgCROS tables.** `agcros2shmi.py` (step B0) wrote the seven AgCROS workbooks in `yield_workbooks/`:
- `SHMI_<site>.xlsx` for NELITCSE, NEMLTCRS, NEMEIRR, NECCIRR, INWLACRE, IAAM7071 and PAUP;
- `SHMI_<site>_units.csv` for each site, which `B2_SHMI_yield_main.R` reads for INWLACRE, IAAM7071 and PAUP.

They are in the repository and used as they are. With `run_conversion = TRUE` and the tables, B0 writes them again over the shipped copies; `git status` then shows whether anything changed.

**Cached** (in `Outputs/yield/`): NASS county yields, copied from `data/cache/` (or downloaded when that copy is missing).

## Part A: calibration (`run_calibration = TRUE`)

| Path | What | Source |
|---|---|---|
| `NAPESHM_rebuilt.xlsx` | NAPESHM management records in the SHMI template, after QC. This is the workbook used in the paper. | built by `python/napeshm_rebuild` from the raw NAPESHM files, then QC |
| `napeshm/sites.csv`, `treatments.csv`, `soildata.csv`, `plfa_ng.csv` | NAPESHM site (climate, descriptors, coordinates), treatment, soil and PLFA tables | NAPESHM, https://doi.org/10.15482/USDA.ADC/25632270 |
| `napeshm/sites_AI_startdate.csv` (optional) | updated treatment start dates, read by `A6c_SHMI_SHAPE_age.R` | NAPESHM site contacts |
| `gSSURGO_CONUS.gpkg` (optional) | gSSURGO CONUS geopackage, for SHAPE soil groups. Large; set `GSSURGO_GPKG` to its path instead of copying it. Without it, put `site_soils_ssurgo.csv` in `Outputs/`. | USDA-NRCS |
| `SHAPE/batch_Combined_scores.csv` | Output of the SHAPE app for `Outputs/SHAPE/SHAPE_multi_indicator_input.csv` | SHAPE app (run outside R) |
| `worldclim/` | WorldClim 2.1 normals, downloaded by `A1d_Get_site_climate.R` on the first run | https://www.worldclim.org |

## Part A0: NAPESHM workbook from the raw files (`run_napeshm_rebuild = TRUE`)

Optional; everything it needs is in the repository. Uses the same `napeshm/` folder: `mgmtdata.csv` (the management event log), `treatments.csv` and `sites.csv`, with `NAPESHM_test.xlsx` as the template layout.

The output, `Outputs/NAPESHM_rebuilt.xlsx`, is before QC. Compare it with `Inputs/NAPESHM_rebuilt.xlsx`.
