# AgCROS → SHMI template converter

`agcros2shmi.py` fills the SHMI input template from the AgCROS (GRACEnet/REAP) table export (`extracted_tables/`).

```
python3 agcros2shmi.py extracted_tables MNMOFS INWLACRE ... --template SHMI_template.xlsx --out new_sites
```

Each location produces three files:

- `SHMI_<loc>.xlsx`: the filled template.
- `_log.csv`: one row per rule applied.
- `_units.csv`: each unit with its AgCROS treatment descriptors and record period.

Dates that are not recorded are left blank, so the package's own imputation fills them, as in the NAPESHM rebuild. The rules are listed at the top of the script:

- **A1, A1b: units.** A plot that changes treatment is split into one unit per treatment. A few records carrying a different treatment label inside the main treatment's period are relabelled to the main treatment.
- **C1–C4: crops.** (C3c: an annual crop with no harvest record takes the harvest date most often recorded for that crop and year on the site's other plots, the rule used by hand for NECCIRR; it reproduces all 52 hand-filled NECCIRR dates.)
   Species come from the AgCROS crop name, and mixtures are split into species. The category (cash, cover, perennial) is inferred. The harvest date and grain yield come from the yield tables. Stover and biomass records of the same harvest are merged.
- **D1–D2: tillage.** Each AgCROS tillage event is mapped to a T-DISC implement, using the recorded depth. Planter passes are left out by default, matching the hand-built workbooks.
- **M1, G1: amendments and grazing.** Organic amendments and N fertilizer (kg N/ha) go into Soil_Amendments; grazing records go into Animal_Diversity.

## Records missing from AgCROS: supplements

`--supplement FILE` (rule S1) adds records that AgCROS lacks. Each row gives the plot, sheet, date, template values, notes and source.

`supplements/NECCIRR_amendments.csv` (992 rows) holds NECCIRR's fertilizer N and manure, which AgCROS does not record:
- **Fertilizer:** the design rate from the treatment code, at corn planting.
- **Manure:** manure plots in the fall after harvest every two years (UNL CropWatch 2017).

It was written from the hand-built workbook by `make_supplement.py`. **Rebuild it from your current DRIVES/final copy if that one has changed:**

```
python3 make_supplement.py DRIVES/final/SHMI_NECCIRR.xlsx NECCIRR supplements/NECCIRR_amendments.csv --source "..."
```

Step B0 of `run_all.Rmd` rebuilds every workbook used in the paper (the commands are in that chunk).

## Corrections to AgCROS records

`--corrections FILE` (rule K1) fixes errors found in AgCROS. Each row carries its evidence, which is written to the log.

`corrections/NEMLTCRS_corrections.csv` fixes two plot pairs:
- **115S / 115M (continuous corn, 0 vs 90 N):** the N records are swapped between the plots.
- **104N / 104S (rotation 7, low vs 0 N):** the planting, tillage and fertilizer tables label the plots the opposite way from the yield table. The DRIVES codes and the yields support the yield table.

A scan of every site for plots whose N disagrees with their treatment found no other cases. In the rotations, the differences between plots are crop phases, not errors.

## Validation against the hand-built workbooks (`validation/`)

The comparison is field by field over the same years, using `compare_workbooks.py`. Each cell is the % of records that match.

| Site | Crop plantings | Harvests | Tillage passes | N fertilizer |
|---|---|---|---|---|
| NELITCSE | 100 | 100 | 100 | 100 |
| NEMEIRR | 100 | 100 | 100 | 100 |
| NEMLTCRS | 100* | 100 | 100 | 83† |
| NECCIRR | 100 | 100 | 100 | 100 with the supplement‡ |

\* Decision: two species. The oats + clover underseeding is written as two species (sweet clover and red clover, as AgCROS records it). The hand-built workbook has one "clover mix", so the converted rows include 990 extra clover rows.

† After the corrections, the only remaining N differences are AgCROS amounts (67 kg) versus the design rate in the hand-built workbook (69 kg), in 65 plots. Zero-N plots have no AgCROS amendment rows; the yield script already treats missing N as 0 at those sites.

‡ AgCROS has no amendment records for NECCIRR. Its fertilizer and manure come from the supplement, and match the hand-built workbook record for record (864 N and 128 manure rows).

**SHMI scores (computed in R from both versions):** NEMEIRR is identical, NELITCSE within 1.2 points and NEMLTCRS within 2.4 (the clover split). NECCIRR differed by up to 12.5 points for two reasons:
- missing harvest dates, now fixed by rule C3c;
- manure plots with no amendment records in AgCROS, now supplied by `supplements/NECCIRR_amendments.csv`.

**Final check:** SHMI was computed in R from both versions of each workbook and compared field by field (a one-off check; the script is not part of the pipeline).

## New candidate sites (`new_sites/`, `candidate_summary.csv`)

| Site | Corn years | Corn units | Contrasts that move SHMI | Notes |
|---|---|---|---|---|
| MNMOFS (Morris Farming Systems) | 8 (2002–09) | 192 | 4-yr vs 2-yr rotation × strip vs conventional till × organic | Strongest addition. Some plantings missing (harvest-only rows, plant date imputed). |
| IAAM7071 (Ames) | 5 (2008–12) | 88 | rye cover, clover living mulch, conservation vs no-till | Good |
| INWLACRE (West Lafayette) | 9 (2003–11) | 56 | no-till vs conventional; ryegrass cover | Good |
| INWLTPAC (West Lafayette) | 8 (2003–11) | 48 | no-till vs conventional; rye cover | Good |
| PAUP (University Park) | 5 (2008–12) | 60 | rye cover; continuous corn vs corn–soybean | No tillage records (no-till study) |
| MNSPReap (Lamberton) | 4 (2008–11) | 64 | strip vs conventional till | Short |
| COFOARD1 (Fort Collins) | 7 (2002–14) | 12 | no-till vs conventional, irrigated | One plot per treatment; weak |
| SDBRREAP (Brookings) | 6 (2008–13) | 4 | legume cover | Too few plots |

Screened out:

- **Residue-removal-only studies:** NEMERREM, IAAMBRUN, INWLREAP, SCFLSGI. SHMI does not score residue, so they have no SHMI contrast.
- **No corn treatment contrasts:** PAHAW, NDMABCS.

The full screen of all 43 locations is in `agcros_location_screen.csv`.

## To add a site to the yield analysis

1. Check the converted workbook: open the log and look at the harvest-only rows (C4) and the units file.
2. Add the workbook to the `workbooks` list in `B2_SHMI_yield_main.R`, along with the site's:
   - N rule (`zero_if_missing` / `optimum_if_missing`)
   - `trt_regex`
   - county for NASS yields
3. Add the site to `B1_Get_SPEI.R`.
4. Its soil C is already in `meassoilchem` for the mediation test.
