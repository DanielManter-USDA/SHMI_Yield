# Inputs

Files the analysis reads but does not create. `run_all.Rmd` checks this list before it starts and reports anything missing for the parts you chose to run.

Everything in this folder except this file is ignored by git. Commit a file only when its terms allow redistribution.

## Part B: corn yield (the default run)

| Path | What | Source |
|---|---|---|
| `agcros/extracted_tables/*.csv` | AgCROS (GRACEnet/REAP) table export, one CSV per table. Used: `meassoilchem`, `measyieldnutuptake`, `measharvestfraction`, `measresiduemgnt`, `mgtplanting`, `mgttillage`, `mgtamendments`, `mgtgrazing`, `metatreatments`, `metaunits`, `metasites` | USDA-ARS AgCROS; snapshot exported 2026-04-21 |
| `yield_workbooks/SHMI_OHTVD_Ho.xlsx` | Hoytville, OH: SHMI workbook built by hand from DRIVES, tillage corrected with the Ohio management logs | DRIVES |
| `yield_workbooks/SHMI_OHTVD_Wo.xlsx` | Wooster, OH: SHMI workbook built by hand from DRIVES | DRIVES |
| `yield_workbooks/SHMI_WICST_ARS.xlsx` | Arlington, WI (WICST): SHMI workbook built by hand from DRIVES and the station records | DRIVES; WICST data request (6-5-2026) |
| `climate/spei03.nc` | SPEIbase 3-month SPEI (global netCDF) | https://spei.csic.es/database.html |
| `NASS_API_KEY` in `~/.Renviron` | Key for USDA-NASS Quick Stats (county corn yields). Not a file. | https://quickstats.nass.usda.gov/api |

**Created by the run (not inputs).** When any of them is missing (or with `run_conversion = TRUE`), `agcros2shmi.py` writes the seven AgCROS workbooks into `yield_workbooks/`:
- `SHMI_<site>.xlsx` for NELITCSE, NEMLTCRS, NEMEIRR, NECCIRR, INWLACRE, IAAM7071 and PAUP;
- `SHMI_<site>_units.csv` for each site, which `SHMI_yield_main.R` reads for INWLACRE, IAAM7071 and PAUP.

If all seven are already here (for example copied from an earlier run), they are used as they are.

**Downloaded or cached on the first run** (in `Outputs/yield/`): NASS county yields.

## Part A: calibration (`run_calibration = TRUE`)

| Path | What | Source |
|---|---|---|
| `NAPESHM_rebuilt.xlsx` | NAPESHM management records in the SHMI template, after QC. This is the workbook used in the paper. | built by `python/napeshm_rebuild` from the raw NAPESHM files, then QC |
| `sites.csv` | NAPESHM site climate and descriptors | NAPESHM, https://doi.org/10.15482/USDA.ADC/25632270 |
| `CSV Files/sites.csv`, `CSV Files/treatments.csv`, `CSV Files/soildata.csv`, `CSV Files/plfa_ng.csv` | NAPESHM site, treatment, soil and PLFA tables | NAPESHM |
| `gSSURGO_CONUS.gpkg` (optional) | gSSURGO CONUS geopackage, for SHAPE soil groups. Large; set `GSSURGO_GPKG` to its path instead of copying it. Without it, put `site_soils_ssurgo.csv` in `Outputs/`. | USDA-NRCS |
| `SHAPE/batch_Combined_scores.csv` | Output of the SHAPE app for `Outputs/SHAPE/SHAPE_multi_indicator_input.csv` | SHAPE app (run outside R) |
| `worldclim/` | WorldClim 2.1 normals, downloaded by `Get_site_climate.R` on the first run | https://www.worldclim.org |

## Part A0: NAPESHM workbook from the raw files (`run_napeshm_rebuild = TRUE`)

`napeshm_raw/` holds:
- `mgmtdata.csv`, `treatments.csv` and `sites.csv` (raw NAPESHM);
- `NAPESHM_test.xlsx`, used as the template layout.

The output, `Outputs/NAPESHM_rebuilt.xlsx`, is before QC. Compare it with `Inputs/NAPESHM_rebuilt.xlsx`.
