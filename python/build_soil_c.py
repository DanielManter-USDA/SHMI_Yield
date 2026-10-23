"""
build_soil_c.py
Soil carbon by plot, date and layer for the AgCROS yield sites, from the AgCROS
table export (meassoilchem.csv), for SHMI_soil_mediation.R.

    python build_soil_c.py <extracted_tables> <out_dir> [--sites NEMLTCRS NELITCSE ...]

Writes <out_dir>/soil_c.csv (analysis years) and soil_c_all_years.csv (all years).

Rules
  carbon      organic C (organic_c_gc_kg) where recorded, otherwise total C
              (tsc_gc_kg); values <= 0 are treated as missing (AgCROS stores
              unmeasured samples as 0, e.g. IAAM7071 2008-2009, 2012 total C)
  layers      a layer that spans two or more other layers sampled in the same
              site-year is a composite and is dropped (NEMLTCRS 1998 0-30 cm)
  depth       layers ending at or above 30 cm (SHMI_soil_mediation.R uses 0-30 cm)
  years       soil_c.csv drops baseline years (on or before the year the treatments
              began, BASELINE below) and years sampling fewer than 80% of the
              site's units (the most sampled year = 100%)
"""
import argparse, os
import pandas as pd

SITES = ["NEMLTCRS", "NELITCSE", "NEMEIRR", "NECCIRR", "IAAM7071", "INWLACRE", "PAUP"]
COVERAGE = 0.8
MAX_DEPTH = 30
# last pre-treatment (baseline) sampling year at each site
BASELINE = {"NEMLTCRS": 1984, "NELITCSE": 1982, "NECCIRR": 2010, "NEMEIRR": 2001,
            "INWLACRE": 2003, "IAAM7071": 2007, "PAUP": 2008}


def composite(g):
    lay = set(zip(g.top_cm, g.bottom_cm))
    def spans(t, b):
        return sum(1 for (t2, b2) in lay if (t2, b2) != (t, b) and t2 >= t and b2 <= b) >= 2
    return pd.Series([spans(t, b) for t, b in zip(g.top_cm, g.bottom_cm)], index=g.index)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("tables"); ap.add_argument("out_dir")
    ap.add_argument("--sites", nargs="+", default=SITES)
    a = ap.parse_args()

    d = pd.read_csv(os.path.join(a.tables, "meassoilchem.csv"), low_memory=False)
    d = d[d.location_id.isin(a.sites)].copy()
    for c in ["organic_c_gc_kg", "tsc_gc_kg", "upper_cm", "lower_cm"]:
        d[c] = pd.to_numeric(d[c], errors="coerce")
    oc = d.organic_c_gc_kg.where(d.organic_c_gc_kg > 0)
    tc = d.tsc_gc_kg.where(d.tsc_gc_kg > 0)
    d["soc_gkg"] = oc.fillna(tc)
    d["c_source"] = oc.notna().map({True: "organic C", False: "total C"})
    n0 = len(d)
    d = d[d.soc_gkg.notna()].copy()
    d["year"] = pd.to_datetime(d.date, errors="coerce").dt.year.astype(int)
    out = pd.DataFrame({
        "site": d.location_id, "year": d.year,
        "unit": [u[len(s) + 1:] if u.startswith(s + "_") else u for u, s in zip(d.unit_id, d.location_id)],
        "agcros_treatment": d.treatment_id, "top_cm": d.upper_cm, "bottom_cm": d.lower_cm,
        "soc_gkg": d.soc_gkg.round(4), "bd_gcm3": pd.NA, "c_source": d.c_source})
    comp = out.groupby(["site", "year"], group_keys=False).apply(composite)
    print(f"{n0 - len(d)} records without carbon (missing or 0) dropped; {int(comp.sum())} composite-layer records dropped")
    out = out[~comp.reindex(out.index).fillna(False) & (out.bottom_cm <= MAX_DEPTH)]
    out = out.sort_values(["site", "year", "unit", "top_cm"])
    os.makedirs(a.out_dir, exist_ok=True)
    out.to_csv(os.path.join(a.out_dir, "soil_c_all_years.csv"), index=False)

    first = BASELINE
    keep = []
    for (s, y), g in out.groupby(["site", "year"]):
        units = g.unit.nunique()
        top = out[out.site == s].groupby("year").unit.nunique().max()
        ok = y > first.get(s, -1) and units >= COVERAGE * top
        print(f"  {s} {y}: {units} units{'' if ok else '  (dropped: ' + ('baseline' if y <= first.get(s, -1) else 'coverage') + ')'}")
        if ok:
            keep.append((s, y))
    res = out[[(s, y) in keep for s, y in zip(out.site, out.year)]]
    res.to_csv(os.path.join(a.out_dir, "soil_c.csv"), index=False)
    print("Saved soil_c.csv:", len(res), "rows;", res.groupby("site").size().to_dict())


if __name__ == "__main__":
    main()
