#!/usr/bin/env python3
"""make_supplement.py WORKBOOK.xlsx LOCATION OUT.csv [--sheet Soil_Amendments] [--source TEXT]

Writes the records of one sheet of an SHMI workbook as a supplement file for
agcros2shmi.py (--supplement). Used once to carry over records that AgCROS lacks
(e.g. NECCIRR's fertilizer design rates and manure, entered by hand with their
assumptions in SA_notes). Plots are identified by MGT_field, so the supplement does
not depend on how the workbook named its units."""
import argparse, pandas as pd, openpyxl
ap = argparse.ArgumentParser()
ap.add_argument("workbook"); ap.add_argument("location"); ap.add_argument("out")
ap.add_argument("--sheet", default="Soil_Amendments"); ap.add_argument("--source", default="")
a = ap.parse_args()
wb = openpyxl.load_workbook(a.workbook, read_only=True, data_only=True)
def sheet(n):
    r = list(wb[n].iter_rows(values_only=True)); d = pd.DataFrame(r[4:], columns=r[3])
    return d[[c for c in d.columns if isinstance(c, str) and c]].dropna(how="all")
mu = sheet("Mgt_Unit").dropna(subset=["MGT_field"])
mu["MGT_combo"] = [f"{s}_{f}_{fl}_{t}" for s, f, fl, t in zip(mu.MGT_study, mu.MGT_farm, mu.MGT_field, mu.MGT_trt)]
field = dict(zip(mu.MGT_combo, mu.MGT_field.astype(str)))
d = sheet(a.sheet).dropna(subset=["MGT_combo"])
datecol = {"Soil_Amendments": "SA_date", "Soil_Disturbance": "SD_date", "Crop_Diversity": "CD_plant_date"}[a.sheet]
out = pd.DataFrame({"location_id": a.location, "unit_id": d.MGT_combo.map(field), "sheet": a.sheet,
                    "date": pd.to_datetime(d[datecol]).dt.strftime("%Y-%m-%d")})
for c in d.columns:
    if c not in ("MGT_combo", datecol) and not c.endswith("_notes") and not c.startswith("SD_tdisc"):
        out[c] = d[c].values
notes = [c for c in d.columns if c.endswith("_notes")]
out["notes"] = d[notes[0]].values if notes else ""
out["source"] = a.source or f"carried over from {a.workbook.split('/')[-1]}"
miss = out.unit_id.isna().sum()
out = out.dropna(subset=["unit_id"])
out.to_csv(a.out, index=False)
print(f"{len(out)} rows written to {a.out}" + (f"; {miss} rows without a matching plot dropped" if miss else ""))
