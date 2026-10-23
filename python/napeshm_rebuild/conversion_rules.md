# Rebuilding NAPESHM_test.xlsx from the raw NAPESHM files

**Inputs:** `mgmtdata.csv` (event log), `treatments.csv`, `sites.csv`, and the existing `NAPESHM_test.xlsx`. The existing workbook supplied the layout, User_Info and Soil_Test, and the spelling of crop names.

**Outputs:**

- `NAPESHM_rebuilt.xlsx`: same sheets, columns and formats as the existing workbook, plus a `Conversion_Log` sheet with one row per decision.
- `NAPESHM_rebuilt_vs_existing.xlsx`: unit-by-unit comparison.
- `convert.py` / `write_wb.py`: the conversion code.

**General principle:** a date that is not in the raw records is left blank, exactly as the template allows, so the SHMI package's own imputation rules fill it (including the intensive-tillage end rule added to `prepare_shmi_inputs()` in 1.2.0). The rebuilt file therefore contains only what the raw data say, plus the explicit rules below. Every rule application is logged with its rule number.

## Rules

| Rule | What it does |
|---|---|
| **R0** | Records whose date is missing or `None` are skipped. |
| **R1** | **Mgt_Unit:** one row per treatment in `treatments.csv` that has records in `mgmtdata.csv`; the 18 treatments without any records are left out of every sheet, including Soil_Test. Field = site code, lat/lon/state from `sites.csv`, notes = study name, tile drain = "Unsure", study = "NAPESHM". Names are written as values, not formulas, so R reads them directly. |
| **R2** | **Crop plantings:** each `plant cash crop`, `plant cover crop` or `planting` record with a species (not measured in STIR) starts an episode. Species planted on the same date in the same category form one episode (same `CD_seq_num`, `CD_mix = yes`). Planting records that name equipment or fertilizer instead of a species are skipped. `CD_cat`: cash crop → `cash`, cover crop → `cover`, `planting` → `cash`, then R5. |
| **R2b** | A cropped (row-crop) treatment with no crop records at all is a bare fallow: one `Fallow` row from the record's start to its end (end in `CD_harv_date`, following the existing convention). Applies to USMI03_297 and 298. |
| **R2c** | A perennial, non-farmed or woody treatment with no plantings becomes one stand spanning the record, from the "start" row to the "end" row; a missing end date defaults to 31 December of the last event year. **Species:** harvested species if named; otherwise `N-species` with N = `plant_species_ct` from `treatments.csv`; otherwise `Hay` if the harvests say hay; otherwise the perennial species recorded at the same site's other treatments; woody species from the treatment description. Unassigned harvests are attached to the stand. |
| **R2d** | A crop planted while another cash crop is still growing (interseeded, broadcast or frost-seeded into it) is a relay cover crop, whatever its raw category. The companion's harvest doesn't end it; it runs until its own harvest or termination (or the package's imputation). If it is later harvested itself, it is coded `cash`. |
| **R3** | A harvest of a named species with no planting in the record creates a harvest-only episode (plant date blank): the crop was planted before the record starts. |
| **R3c** | A harvest-only crop takes its planting date from the last planting pass (STIR record in a planting category) that has no species and falls after the previous harvest. |
| **R3b** | A `cover crop termination` before any planting in the record creates a pre-record cover episode whose species are taken from the unit's next cover-crop planting. **Inferred; flagged in CD_notes.** |
| **R4** | **Harvests** end the episode of the species named. A generic harvest (equipment name, "hauling") is assigned to an unharvested crop. A generic harvest or residue operation *after* the crop was harvested is ignored, and harvest records describing residue removal (stover, straw, raking, chopping) never end a crop. When the same species has more than one active stand, the harvest goes to the newest. A second harvest of an annual crop within 60 days is ignored (grain then stover). Perennial or forage species get one extra row per harvest (plant date blank), as in the existing workbook. `CD_per_res` comes from the harvest quantity when its unit is % residue removed. |
| **R4** | **Terminations** (`cover crop termination`, `termination`, `pre-plant burndown`, `burndown`, burning) end the active cover or perennial crops. Several products applied on the same day count as one termination. Cash crops are ended only by harvest (or `termination`). |
| **R4b** | A non-selective herbicide (glyphosate, paraquat, glufosinate) recorded as `weed management` ends an active cover crop, provided no unharvested cash crop is in the ground. **Inferred.** |
| **R4c** | `mowing` of a named crop ("hay alfalfa", "sorghum sudangrass") counts as a harvest (a cut); the species is read from the description, or else from the method. |
| **R5** | A perennial species (alfalfa, orchardgrass, bahiagrass, tall fescue, timothy, red clover, gamagrass, prairie mixes) is coded `perennial` when the stand lasts more than 365 days or is harvested in more than one year. Poplar is `woody perennial`. |
| **R5b** | Repeated plantings of the same perennial species, with no termination or other crop in between, are merged into one stand (later plantings are overseeding). |
| **R6** | Raw names are mapped to template names (`corn` → Corn, `cereal rye (secale cereale l)` → Rye, Cereal …), following the existing workbook's choices. **Exceptions:** berseem clover → Clover, Berseem (existing: Clover, Balansa); generic "cover crop" and "native warm season grasses" → 8-species (as existing); "triticale or oats" → Triticale (as existing); "winter legume cover crop" → Legume (as existing). |
| **R7** | **Soil_Disturbance:** every record measured in STIR, whatever its category, including planting passes and harvest traffic. `SD_phase` = raw category; blank and misspelled categories become `disturbance`. |
| **R8** | **Soil_Amendments:** organic amendments only (manure, compost, litter, biosolids, …), matching the existing workbook, which has no inorganic fertilizer. Includes 8 records where the amendment name sits in the category column, and organic materials recorded with a blank category. At USID01_246/247 the record names only the spreader; it is entered as solid manure, as confirmed by the data owner. |
| **R9** | **Animal_Diversity:** grazing periods from "begins"/"ends" pairs. A "stocking rate" record starts a period when there is no "begins". "herd N" plus "acres used for herd N" makes a period with count = head/acres. `AD_intensity` is left blank: its rule could not be inferred, and SHMI uses only presence. |
| **R9b** | Cattle trampling recorded as a disturbance (method "cattle") = a 7-day grazing period, as in the existing workbook, unless the date already falls inside a recorded grazing period. |

**Carried over unchanged:** User_Info and Soil_Test. No soil-test file was supplied, and SHMI doesn't use Soil_Test.

**Not represented, because the template has no field for them:** the raw "start"/"end" record rows (used only in R2b/R2c), soil-sampling dates, inorganic fertilizer, and pesticide applications apart from R4b.

## Comparison with the existing workbook

| Sheet | Existing rows | Rebuilt rows | Units compared | Units identical |
|---|---|---|---|---|
| Mgt_Unit | 490 | 490 | 490 | all |
| Crop_Diversity | 1,694 | 1,808 | 490 | 348 (identical crop records) |
| Soil_Disturbance | 3,293 | 3,324 | 448 | 443 (same count) |
| Soil_Amendments | 108 | 118 | 52 | 48 |
| Animal_Diversity | 124 | 105 | 28 | 21 |

The 18 treatments without records in `mgmtdata.csv` are excluded from both sides of the comparison.

### Findings that need your input

1. **18 treatments have no records in `mgmtdata.csv`** (USCO01_185–187, USCO02_197–199, USCO03_209–211, USCO05_225, USKS02_265, USNE01_387, USNY01_395, USNY02_400, USNY03_403, USNY04_406, USOK01_433, USWY01_508). They are removed from the rebuilt workbook.
2. **The existing workbook contains information that isn't in the raw log.** Examples:
   - alfalfa at USAR01 (the raw says only "hay")
   - some pre-record crops
   - cover-crop termination dates filled in from the next planting
   - winter wheat harvests coded as terminations (USAL01_137–141)
3. **Species counts look swapped between two units.** USWI01_501 is "4-species" in the existing workbook, but `plant_species_ct` = 8. USWV01_505 is "8-species", but its count = 4.
4. **Existing dates that differ from the raw log.** Examples:
   - **USMI01_282–284:** raw rye+vetch planted 14 Sep 2015 and flail-mowed 24 May 2016; the existing workbook has 15 Sep 2015 and 2 May 2016, with similar offsets in later years.
   - **USKY02_273/274:** soybean year 2020 vs 2021.
   - **USNE04_392:** soybean harvest 17 vs 7 Oct.
5. **Records after soil sampling** (e.g. USPA03_458 through 2025) appear in both files. That needs the window-cap decision (`end_date_max`).

### Suggested check in R

Run `prepare_shmi_inputs()` and `build_shmi()` on both workbooks and compare SHMI and each pillar unit by unit. Crop differences affect Cover and Diversity; animal differences affect OrgInput. Units where the rebuilt file leaves an end date blank are scored with the package's imputation. Comparing those with the existing file's hand-filled dates is exactly the "your rules vs mine" test.
