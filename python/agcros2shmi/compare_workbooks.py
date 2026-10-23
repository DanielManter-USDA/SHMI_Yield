#!/usr/bin/env python3
"""compare_workbooks.py HANDBUILT.xlsx CONVERTED.xlsx OUT.csv
Compares two SHMI workbooks field by field (MGT_field), over the years both cover:
crop plantings (year, species), harvest dates, tillage passes (date, implement) and
amendments (date, category, N). Prints a summary and writes per-field agreement."""
import sys, re, collections
import pandas as pd, openpyxl

def load(path):
    wb = openpyxl.load_workbook(path, read_only=True, data_only=True)
    def sheet(name):
        rows = list(wb[name].iter_rows(values_only=True))
        hdr = rows[3]; d = pd.DataFrame(rows[4:], columns=hdr)
        return d[[c for c in d.columns if c]].dropna(how="all")
    mu = sheet("Mgt_Unit")
    mu = mu[mu.MGT_field.notna()].copy()
    # MGT_combo may be a formula in hand-built files: rebuild it
    mu["MGT_combo"] = [f"{a}_{b}_{c}_{d}" for a, b, c, d in zip(mu.MGT_study, mu.MGT_farm, mu.MGT_field, mu.MGT_trt)]
    f = dict(zip(mu.MGT_combo, mu.MGT_field.map(lambda x: str(x).lstrip("0") or "0")))
    out = {}
    for nm in ["Crop_Diversity", "Soil_Disturbance", "Soil_Amendments"]:
        d = sheet(nm); d = d[d.MGT_combo.notna()].copy(); d["field"] = d.MGT_combo.map(f); out[nm] = d
    return out

norm = lambda s: re.sub(r"[^a-z]", "", str(s).lower().replace("cereal", "").replace("spring", "").replace("winter", ""))
def species_key(s):
    s = norm(s)
    for k in ["corn", "soybean", "oats", "sorghum", "wheat", "alfalfa", "rye", "vetch", "turnip", "clover", "radish"]:
        if k in s: return k
    return s

A, B = load(sys.argv[1]), load(sys.argv[2])
rows = []
for kind in ["crop", "harvest", "tillage", "amend_org", "amend_N"]:
    def keys(W):
        S = collections.defaultdict(collections.Counter)
        if kind in ("crop", "harvest"):
            d = W["Crop_Diversity"]
            col = "CD_plant_date" if kind == "crop" else None
            for _, r in d.iterrows():
                if kind == "crop":
                    if pd.notna(r.CD_plant_date): S[r.field][(pd.Timestamp(r.CD_plant_date).year, species_key(r.CD_name))] += 1
                else:
                    end = r.CD_harv_date if pd.notna(r.CD_harv_date) else r.CD_term_date
                    if pd.notna(end) and str(r.CD_cat).lower() != "cover": S[r.field][(pd.Timestamp(end).date(), species_key(r.CD_name))] += 1
        elif kind == "tillage":
            for _, r in W["Soil_Disturbance"].iterrows():
                if pd.notna(r.SD_date) and pd.notna(r.SD_equip) and not str(r.SD_equip).startswith("PLANTER"):
                    S[r.field][(pd.Timestamp(r.SD_date).date(), str(r.SD_equip))] += 1
        else:
            for _, r in W["Soil_Amendments"].iterrows():
                if pd.isna(r.SA_date): continue
                if kind == "amend_org" and str(r.SA_cat) == "Organic": S[r.field][(pd.Timestamp(r.SA_date).date(),)] += 1
                if kind == "amend_N" and str(r.SA_cat) == "Fertilizer" and pd.notna(r.SA_N):
                    S[r.field][(pd.Timestamp(r.SA_date).year, round(float(r.SA_N)))] += 1
        return S
    KA, KB = keys(A), keys(B)
    for fld in sorted(set(KA) | set(KB)):
        a, b = KA.get(fld, collections.Counter()), KB.get(fld, collections.Counter())
        ya = {k[0] if isinstance(k[0], int) else k[0].year for k in a}; yb = {k[0] if isinstance(k[0], int) else k[0].year for k in b}
        yrs = ya & yb if ya and yb else ya | yb
        yr = lambda k: k[0] if isinstance(k[0], int) else k[0].year
        a = collections.Counter({k: v for k, v in a.items() if yr(k) in yrs}); b = collections.Counter({k: v for k, v in b.items() if yr(k) in yrs})
        both = sum((a & b).values())
        rows.append(dict(kind=kind, field=fld, handbuilt=sum(a.values()), converted=sum(b.values()), matched=both,
                         only_handbuilt=sum((a - b).values()), only_converted=sum((b - a).values()),
                         ex_only_hand="; ".join(map(str, list((a - b).keys())[:3])), ex_only_conv="; ".join(map(str, list((b - a).keys())[:3]))))
R = pd.DataFrame(rows)
R.to_csv(sys.argv[3], index=False)
S = R.groupby("kind")[["handbuilt", "converted", "matched", "only_handbuilt", "only_converted"]].sum()
S["agreement_%"] = (100 * S.matched / S[["handbuilt", "converted"]].max(axis=1)).round(1)
print(S.to_string())
