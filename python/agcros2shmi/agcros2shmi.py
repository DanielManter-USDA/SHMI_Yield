#!/usr/bin/env python3
"""
agcros2shmi.py - fill the SHMI input template from AgCROS (GRACEnet/REAP) tables.

    python3 agcros2shmi.py TABLES_DIR LOCATION_ID [LOCATION_ID ...] \
        --template SHMI_template.xlsx --out OUTDIR [--study AgCROS] [--planter-passes]

TABLES_DIR holds the AgCROS table exports (extracted_tables/: metaunits.csv,
metatreatments.csv, metasites.csv, mgtplanting.csv, mgttillage.csv, mgtamendments.csv,
mgtgrazing.csv, mgtresidue.csv, measresiduemgnt.csv, measharvestfraction.csv,
measyieldnutuptake.csv). For each location it writes
    OUTDIR/SHMI_<LOCATION>.xlsx          the filled template (values only)
    OUTDIR/SHMI_<LOCATION>_log.csv       one row per rule application
    OUTDIR/SHMI_<LOCATION>_units.csv     unit x treatment x period, with descriptors

Principle (as in the NAPESHM rebuild): write only what the records say plus the explicit,
logged rules below; dates that are not recorded stay blank so the SHMI package's own
imputation fills them. Rules:

  A1  One management unit (MGT_combo) per AgCROS unit x treatment. A plot that changes
      treatment (e.g. a baseline period) becomes two units, split by the treatment_id
      on each record. MGT_study = --study, MGT_farm = location, MGT_field = unit_id
      without the location prefix, MGT_trt = treatment_id without the prefix.
  A2  Coordinates, state and county from metasites (unit coordinates are rarely given).
      MGT_irr_cat = "Sprinkler, Pivot" when the treatment says irrigation = Yes
      (method not recorded in AgCROS; only presence matters to SHMI), else "None".
  C1  Each planting record is one crop row. Species = common name in parentheses
      ("Zea mays (Corn)" -> Corn), mapped to the template's crop list where possible.
      Mixtures ("Oats+[80% Melilotus ... 20% Trifolium pratense(Clover)]") become one
      row per named species planted the same day (same CD_plant_date).
  C1b A crop planted twice before one harvest (a replant) keeps only the later planting.
  C2  CD_cat: Perennial for perennial species (alfalfa, perennial grasses, clovers
      that persist > 365 days before the next planting); Cover for cover-crop species
      (rye, ryegrass, vetch, turnip, radish, clovers, sweet clover) when not harvested
      for grain; otherwise Annual.
  C3  Harvest date and grain yield (dry matter, kg/ha) from measresiduemgnt,
      measharvestfraction or measyieldnutuptake (grain fractions) for the same unit and
      crop: the first record with a grain yield 30-400 days after planting, else the first
      record of any kind (7 Nov 2026: at IAAM7071 a mid-season above-ground sampling,
      about 20 June, had been taken as the harvest and the grain harvest left orphaned). Without one, the date is blank
      (package imputes). Perennials get CD_harv_date only on the first harvest.
  C3c An annual crop without a harvest record takes the harvest date most often recorded
      for the same crop planted that year on the site's other plots (yield stays blank).
  C4  Harvest-only records with a yield (a planting missed in the log) create a crop
      row; its planting date is the date most often recorded for the same crop and year
      on the site's other plots, if any and if the plot has no other annual crop planted
      that year (7 Nov 2026: INWLACRE's modified-planter plots in 2004), else blank. Harvest-only records without a yield (samplings) are ignored.
  D1  Each tillage record is one pass, mapped by tillage_event (and the method text) to
      a T-DISC implement (TILLAGE_MAP). Recorded depth (cm) is written as SD_depth
      (inches); blank or 0 depth leaves the T-DISC default. "(Twice)" = two passes.
      Events that do not move soil (flail mowing, stalk shredding) are not written.
  D2  Planting passes are not written unless --planter-passes is given (the
      hand-built DRIVES/AgCROS workbooks did not include them).
  D3  A site with no tillage record at all (no-till) gets one zero-mixing placeholder pass
      (SPRAYER/AERIAL/BROADCAST APPLICATION, mixing 0) so Soil_Disturbance is not empty.
  D4  Tillage records that the site's own publications say did not happen are not written
      (TILLAGE_NOT_DONE; NELITCSE no-till plots: row cultivation 1981-2000 and 2013, which
      Sindelar et al. 2015 rule out). --keep-documented-errors keeps them (sensitivity).
  T1  If no crop has a termination date, the first annual crop with a harvest date
      gets that date as CD_term_date too (same day, so SHMI is unchanged): an all-blank
      CD_term_date column stops SHMI <= 1.2.x (SHMI_package_todo.md). Remove once fixed.
  M2  Biochar (and charcoal) records are not written: biochar is not scored as an organic
      amendment (stable pyrolysed carbon; the pairs analysis already excludes it).
  M1  Amendments: organic materials (manure, compost, litter, ...) ->
      SA_cat Organic with rate and N; N fertilizers -> SA_cat Fertilizer with SA_N
      (kg N/ha). Lime, P/K-only products, pesticides and irrigation rows are skipped.
  G1  Grazing periods from mgtgrazing (start, end, species).
  K1  --corrections FILE fixes AgCROS records, with the evidence for each fix logged:
      - table = a table name: move that table's records from unit_id to new_unit_id
        (rows of one table are applied together, so a swap is two rows);
      - table = "*" with new_treatment_id: set the plot's treatment in every table.
  S1  --supplement FILE adds records that AgCROS lacks (e.g. NECCIRR's fertilizer and
      manure). One row per record: location_id, unit_id (plot), sheet
      (Soil_Amendments, Soil_Disturbance or Crop_Diversity), date, the template columns
      of that sheet (SA_cat, SA_source, SA_N, ...; SD_equip, SD_depth; CD_cat, CD_name,
      CD_harv_date, CD_term_date, CD_yield), source and notes. The row joins the plot's unit whose record
      period contains the date; source and notes go into the notes column.
"""
import argparse, os, re, sys, copy
from datetime import timedelta
import pandas as pd
import openpyxl

# ------------------------------------------------------------------ dictionaries
CROP_NAME = {  # lower-case common name (from parentheses) -> template name
    "corn": "Corn", "soybean": "Soybean", "oats": "Oats", "wheat": "Wheat", "winter wheat": "Wheat, Winter",
    "spring wheat": "Wheat, Spring", "barley": "Barley", "rye": "Rye, Cereal", "winter rye": "Rye, Cereal",
    "cereal rye": "Rye, Cereal", "winter ryegrass": "Ryegrass, Annual", "annual ryegrass": "Ryegrass, Annual",
    "ryegrass": "Ryegrass", "alfalfa": "Alfalfa", "grain sorghum": "Sorghum, Grain", "sorghum": "Sorghum, Grain",
    "sudangrass": "Sudangrass", "hairy vetch": "Vetch, Hairy", "turnip": "Turnip", "radish": "Radish",
    "red clover": "Clover, Red", "clover": "Clover", "sweet clover": "Clover, Sweet", "white clover": "Clover, White",
    "switchgrass": "Switchgrass", "pea": "Pea", "canola": "Canola", "sunflower": "Sunflower", "triticale": "Triticale",
    "millet": "Millet", "camelina": "Camelina",
}
SCI_NAME = {"melilotus officinalis": "Clover, Sweet", "trifolium pratense": "Clover, Red",
            "trifolium repens": "Clover, White", "secale cereale": "Rye, Cereal", "vicia villosa": "Vetch, Hairy"}
# Tillage records AgCROS holds that the site's own publications say did not happen. Each entry:
# (location, regex on treatment_id, regex on "event / method", evidence). Rule D4; skipped with
# --keep-documented-errors (sensitivity).
TILLAGE_NOT_DONE = [
    ("NELITCSE", r"^NELITCSE_NT", r"cultivator",
     "Sindelar et al. 2015, Agron. J. 107:2241: cultivation 'in all tillage treatments except NT'"),
]
PERENNIAL = {"Alfalfa", "Switchgrass", "Ryegrass, Perennial", "Fescue, Tall", "Reed Canary Grass", "Timothy Grass"}
COVER = {"Rye, Cereal", "Ryegrass, Annual", "Ryegrass", "Vetch, Hairy", "Turnip", "Radish", "Clover, Red",
         "Clover", "Clover, Sweet", "Clover, White", "Clover, Crimson", "Triticale"}
GRAIN = re.compile(r"grain", re.I)

# (pattern on event, pattern on method or None) -> T-DISC implement; first match wins
TILLAGE_MAP = [
    (r"moldboard", None, "PLOW, MOLDBOARD"),
    (r"chisel", r"disk chisel|mulch", "HARROW, DISK, CHISEL (MULCH TILLER)"),
    (r"chisel", None, "PLOW, CHISEL"),
    (r"subsoil|disk rip|ripper", None, "RIPPER/SUBSOILER"),
    (r"strip till", None, "NUTRIENT MGMT: FERT APPL, STRIP TILL"),
    (r"^disk", r"offset", "HARROW, DISK, OFFSET, HEAVYDUTY"),
    (r"^disk", None, "HARROW, DISK, TANDEM, HEAVYDUTY"),
    (r"field cultivator", r"row|ridger|buffalo", "CULTIVATOR, ROW"),
    (r"field cultivator", r"field cultivator", "CULTIVATOR, FIELD"),
    (r"field cultivator", None, "SEASON"),          # NA method: in-season (Jun-Aug) = row, else field
    (r"row.?crop cultivat|cultivat.*row", None, "CULTIVATOR, ROW"),
    (r"cultimulcher", None, "ROLLER, CULTI-MULCHER"),
    (r"field roll|roller|packer", None, "ROLLER, PACKER"),
    (r"rotary hoe", None, "CULTIVATOR, ROTARY HOE"),
    (r"harrow", r"spring|spike", "HARROW, SPIKE TOOTH OR SPRING TOOTH"),
    (r"harrow", None, "HARROW, FINISHING"),
    (r"land level|ridge till", None, "LEVELING, LEVELER"),
    (r"ridge", None, "CULTIVATOR, ROW"),
    (r"rototill|rotary till", None, "ROTO-TILLER, ROTO-TILLER DIGGER, ROTERRA"),
    (r"flail|mow|shred", None, None),               # no soil movement
]
ORGANIC = re.compile(r"manure|compost|litter|biosolid|slurry|sludge|guano|meal", re.I)
# Biochar is not scored as an organic amendment (decided 9 Nov 2026): it is stable, pyrolysed
# carbon rather than an input that feeds soil biology, and the long-term pairs analysis (B2f)
# already excludes it. Its records are not written (rule M2).
NOT_SCORED = re.compile(r"biochar|charcoal", re.I)
NOT_N = re.compile(r"lime|^0-|potash|k2o|super ?phosphate|^K20$|gypsum", re.I)

# ------------------------------------------------------------------ helpers
def strip_prefix(x, loc):
    x = str(x)
    return x[len(loc) + 1:] if x.startswith(loc + "_") else x

def species_from(crop):
    """'Zea mays (Corn)' -> ['Corn']; mixtures -> several names."""
    out = []
    s = str(crop)
    parts = re.split(r"\+|\band\b", s)
    for p in parts:
        p = p.strip(" []")
        if not p: continue
        common = re.findall(r"\(([^()]*)\)", p)
        sci = re.sub(r"\(.*?\)|\d+%|[\[\]]", "", p).strip().lower()
        if sci in SCI_NAME: out.append(SCI_NAME[sci]); continue
        nm = common[-1].strip().lower() if common else sci
        out.append(CROP_NAME.get(nm, nm.title()))
    seen = []
    for o in out:
        if o and o not in seen: seen.append(o)
    return seen

def tdisc_for(event, method, date):
    e = str(event or "").lower(); m = str(method or "").lower()
    if e in ("", "na", "nan") and m in ("", "na", "nan"): return None, "no event"
    for pe, pm, imp in TILLAGE_MAP:
        if re.search(pe, e) or (pe and re.search(pe, m) and e in ("", "na", "nan")):
            if pm is not None and not re.search(pm, m): continue
            if imp == "SEASON":
                return ("CULTIVATOR, ROW" if date.month in (6, 7, 8) else "CULTIVATOR, FIELD"), "field cultivator by season"
            return imp, None
    return "UNMAPPED", "no rule"

def _merge(cl):
    g = [r for r in cl if pd.notna(r.y)]
    r0 = (g or cl)[0].copy()
    return r0

def read(tables, name):
    p = os.path.join(tables, name + ".csv")
    return pd.read_csv(p, low_memory=False, na_values=["NA", ""], keep_default_na=True) if os.path.exists(p) else pd.DataFrame()

def to_date(s):
    return pd.to_datetime(s, errors="coerce", format="mixed")

# ------------------------------------------------------------------ conversion
def convert(tables, loc, template, outdir, study="AgCROS", planter=False, supplements=(), corrections=(),
            keep_documented_errors=False):
    log = []
    L = lambda rule, unit, msg: log.append({"rule": rule, "unit": unit, "message": msg})
    T = {k: read(tables, k) for k in ["metaunits", "metatreatments", "metasites", "mgtplanting", "mgttillage",
                                      "mgtamendments", "mgtgrazing", "measresiduemgnt", "measharvestfraction",
                                      "measyieldnutuptake"]}
    for k, d in T.items():
        if len(d) and "location_id" in d: T[k] = d[d.location_id == loc].copy()
    # K1: corrections to AgCROS records (applied before anything else): move a table's
    # records from one plot to another; listed pairs are applied together, so a swap is
    # two rows. Each row carries its evidence.
    for f in corrections:
        C = pd.read_csv(f, dtype=str, keep_default_na=False)
        C = C[C.location_id == loc]
        for tbl, g in C[C.table != "*"].groupby("table"):
            if tbl not in T or not len(T[tbl]): L("K1", "", f"correction for unknown or empty table {tbl} skipped"); continue
            d = T[tbl]
            dd = to_date(d.start_date if "start_date" in d else d.date)
            new_unit = d.unit_id.copy()
            for _, c in g.iterrows():
                m = d.unit_id == c.unit_id
                if c.get("from_date", ""): m &= dd >= to_date(c.from_date)
                if c.get("to_date", ""): m &= dd <= to_date(c.to_date)
                new_unit[m] = c.new_unit_id
                L("K1", c.unit_id, f"{int(m.sum())} {tbl} records moved to {c.new_unit_id}: {c.evidence}")
            # the treatment follows the receiving plot's own label in the other tables
            t_of = {}
            for k2 in ["mgtplanting", "mgttillage", "measresiduemgnt"]:
                if len(T[k2]): t_of.update(T[k2].groupby("unit_id").treatment_id.agg(lambda s: s.mode().iat[0]).to_dict())
            moved = new_unit != d.unit_id
            d = d.copy(); d.loc[moved, "treatment_id"] = new_unit[moved].map(t_of).fillna(d.loc[moved, "treatment_id"])
            d["unit_id"] = new_unit; T[tbl] = d
        # relabel: table "*" with new_treatment_id sets the plot's treatment in every table
        R_ = C[(C.table == "*") & (C.get("new_treatment_id", "") != "")] if "new_treatment_id" in C else C.iloc[0:0]
        for _, c in R_.iterrows():
            n_ = 0
            for k2 in T:
                if len(T[k2]) and "unit_id" in T[k2] and "treatment_id" in T[k2]:
                    m = T[k2].unit_id == c.unit_id
                    n_ += int((T[k2].loc[m, "treatment_id"] != c.new_treatment_id).sum())
                    T[k2].loc[m, "treatment_id"] = c.new_treatment_id
            L("K1", c.unit_id, f"treatment set to {c.new_treatment_id} in all tables ({n_} records relabelled): {c.evidence}")
    site = T["metasites"].iloc[0] if len(T["metasites"]) else pd.Series(dtype=object)
    trt = T["metatreatments"].drop_duplicates("treatment_id").set_index("treatment_id") if len(T["metatreatments"]) else pd.DataFrame()

    # unit x treatment combinations from every management record
    recs = []
    for k in ["mgtplanting", "mgttillage", "mgtamendments", "mgtgrazing"]:
        d = T[k]
        if len(d): recs.append(d[["unit_id", "treatment_id", "start_date"]])
    for k, dc in [("measresiduemgnt", "date"), ("measharvestfraction", "sampling_date"), ("measyieldnutuptake", "sampling_date")]:
        d = T[k]
        if len(d): recs.append(d[["unit_id", "treatment_id", dc]].rename(columns={dc: "start_date"}))
    R = pd.concat(recs).dropna(subset=["unit_id", "treatment_id"])
    R["d"] = to_date(R.start_date)
    units = R.groupby(["unit_id", "treatment_id"]).d.agg(["min", "max", "size"]).reset_index()
    # A1b: a minor treatment label whose records fall inside the period of the plot's
    # main treatment is a labelling error (e.g. one table carrying the wrong treatment
    # for a few years): its records are reassigned to the main treatment.
    remap = {}
    for uid, g in units.groupby("unit_id"):
        if len(g) < 2: continue
        g = g.sort_values("size", ascending=False)
        main = g.iloc[0]
        for _, m in g.iloc[1:].iterrows():
            inside = m["min"] >= main["min"] - timedelta(days=60) and m["max"] <= main["max"] + timedelta(days=60)
            if inside:
                remap[(uid, m.treatment_id)] = main.treatment_id
                L("A1b", f"{uid}", f"{m['size']} records labelled {m.treatment_id} ({m['min'].date()} to {m['max'].date()}) "
                                    f"fall inside {main.treatment_id}'s period; reassigned to {main.treatment_id}")
    if remap:
        fix = lambda d: d.assign(treatment_id=[remap.get((u, t_), t_) for u, t_ in zip(d.unit_id, d.treatment_id)])
        for k in T:
            if len(T[k]) and "unit_id" in T[k] and "treatment_id" in T[k]: T[k] = fix(T[k])
        R = fix(R)
        units = R.groupby(["unit_id", "treatment_id"]).d.agg(["min", "max", "size"]).reset_index()
    combo = {}
    for _, u in units.iterrows():
        f = strip_prefix(u.unit_id, loc); t = strip_prefix(u.treatment_id, loc)
        combo[(u.unit_id, u.treatment_id)] = f"{study}_{loc}_{f}_{t}"
    if units.groupby("unit_id").size().max() > 1:
        L("A1", "", f"{(units.groupby('unit_id').size() > 1).sum()} plots carry more than one treatment over time; split into separate units")

    # --- Mgt_Unit
    irr = lambda t: ("Sprinkler, Pivot" if str(trt.loc[t, "irrigation"]).strip().lower() == "yes" else "None") \
        if t in trt.index and "irrigation" in trt else "None"
    mgt_rows, unit_rows = [], []
    for _, u in units.iterrows():
        mc = combo[(u.unit_id, u.treatment_id)]
        tinfo = trt.loc[u.treatment_id] if u.treatment_id in trt.index else pd.Series(dtype=object)
        mgt_rows.append([mc, study, loc, strip_prefix(u.unit_id, loc), strip_prefix(u.treatment_id, loc), irr(u.treatment_id),
                         site.get("lat_dec_deg"), site.get("long_dec_deg"), STATE.get(str(site.get("state_province")), site.get("state_province")),
                         site.get("county"), None, None, None, str(tinfo.get("treatment_descriptor", ""))[:250]])
        unit_rows.append({"MGT_combo": mc, "unit_id": u.unit_id, "treatment_id": u.treatment_id,
                          "first_record": u["min"].date() if pd.notna(u["min"]) else None,
                          "last_record": u["max"].date() if pd.notna(u["max"]) else None, "records": u["size"],
                          **{k: tinfo.get(k) for k in ["treatment_descriptor", "rotation_descriptor", "tillage_descriptor",
                                                       "n_treatment_descriptor", "cover_crop", "residue_removal",
                                                       "irrigation", "organic_management"]}})

    # --- harvests (grain yields and dates)
    hv = []
    rm = T["measresiduemgnt"]
    if len(rm):
        g = rm[rm.grain_dry_matt_kg_ha.notna() | rm.harvested_frac.astype(str).str.contains("grain", case=False)]
        hv.append(g.assign(d=to_date(g.date), y=g.grain_dry_matt_kg_ha)[["unit_id", "treatment_id", "crop", "d", "y"]])
    for k, col in [("measharvestfraction", "frac_dry_matt_kg_ha"), ("measyieldnutuptake", "frac_dry_matt_kg_ha_ynu")]:
        d = T[k]
        if len(d):
            g = d[d.plant_fraction.astype(str).str.contains("grain|aboveground", case=False)]
            g = g.assign(d=to_date(g.sampling_date), y=g[col].where(g.plant_fraction.astype(str).str.contains("grain", case=False)))
            hv.append(g[["unit_id", "treatment_id", "crop", "d", "y"]])
    hv = pd.concat(hv) if hv else pd.DataFrame(columns=["unit_id", "treatment_id", "crop", "d", "y"])
    hv = hv.dropna(subset=["d"])
    hv["sp"] = hv.crop.map(lambda c: species_from(c)[0] if species_from(c) else None)
    hv = hv.groupby(["unit_id", "treatment_id", "sp", "d"], as_index=False).y.mean()
    # C3b: records of the same annual crop within 60 days (grain, stover, biomass or the
    # same harvest in two tables) are one harvest: date of the grain-yield record if any
    hv = hv.sort_values(["unit_id", "sp", "d"])
    out_h, n_merged = [], 0
    for (u_, s_), g in hv.groupby(["unit_id", "sp"]):
        if s_ in PERENNIAL:
            out_h.append(g); continue
        cl = []
        for _, r in g.iterrows():
            if cl and (r.d - cl[0].d).days <= 60: cl.append(r); continue
            if cl: out_h.append(pd.DataFrame([_merge(cl)])); n_merged += len(cl) - 1
            cl = [r]
        if cl: out_h.append(pd.DataFrame([_merge(cl)])); n_merged += len(cl) - 1
    hv = pd.concat(out_h, ignore_index=True) if out_h else hv
    # ... and a later record without grain yield up to 150 days after a grain harvest
    # (stover baling, residue sampling) is not a new crop either
    hv = hv.sort_values(["unit_id", "sp", "d"]).reset_index(drop=True)
    drop = []
    for (u_, s_), g in hv.groupby(["unit_id", "sp"]):
        if s_ in PERENNIAL: continue
        lg = None
        for i, r in g.iterrows():
            if pd.notna(r.y): lg = r.d
            elif lg is not None and (r.d - lg).days <= 150: drop.append(i)
    if drop: L("C3b", "", f"{len(drop)} later non-grain harvest records (<=150 days after grain harvest) ignored")
    hv = hv.drop(index=drop)
    if n_merged: L("C3b", "", f"{n_merged} extra harvest records of annual crops within 60 days merged into one harvest")
    used = set()

    # --- Crop_Diversity
    crop_rows = []
    pl = T["mgtplanting"].copy()
    pl["d"] = to_date(pl.start_date)
    pl = pl.dropna(subset=["d"]).sort_values(["unit_id", "d"])
    for (uid, tid), g in pl.groupby(["unit_id", "treatment_id"]):
        mc = combo.get((uid, tid))
        if mc is None: continue
        H = hv[(hv.unit_id == uid)]
        plant_dates = sorted(g.d.unique())
        for _, r in g.iterrows():
            sps = species_from(r.crop)
            if len(sps) > 1: L("C1", mc, f"{r.d.date()} mixture {r.crop} -> {' + '.join(sps)}")
            for sp in sps:
                cand = H[(H.sp == sp) & (H.d >= r.d + timedelta(days=30)) & (H.d <= r.d + timedelta(days=400))]
                harv, y = (None, None)
                if len(cand):
                    cg = cand[cand.y.notna()]                      # C3: a grain harvest first
                    c0 = (cg if len(cg) else cand).sort_values("d").iloc[0]; harv, y = c0.d, c0.y
                    used.add((uid, sp, c0.d))
                nxt = [p for p in plant_dates if p > r.d]
                if sp in PERENNIAL: cat = "Perennial"
                elif sp in COVER and (y is None or pd.isna(y)):
                    cat = "Cover"
                    if sp.startswith("Clover") and nxt and (nxt[0] - r.d).days > 365: cat = "Perennial"
                else: cat = "Annual"
                if cat == "Annual" and harv is None:
                    L("C3", mc, f"{r.d.date()} {sp}: no harvest record; harvest date left blank")
                crop_rows.append([mc, cat, None, sp, r.d.to_pydatetime(), harv.to_pydatetime() if harv is not None else None, None,
                                  None if y is None or pd.isna(y) else round(float(y), 1),
                                  None if y is None or pd.isna(y) else "kgs/hectare", None, f"AgCROS planting: {r.crop}"])
    # C1b: a crop planted twice before one harvest (a replant) keeps the later planting
    seen, keep_rows = {}, []
    for row in sorted(crop_rows, key=lambda x: (x[0], x[3], str(x[5]), str(x[4])), reverse=True):
        k = (row[0], row[3], row[5])
        if row[5] is not None and k in seen:
            L("C1b", row[0], f"{row[4].date() if row[4] else ''} {row[3]}: replanted on {seen[k].date()}; earlier planting dropped")
            continue
        if row[5] is not None: seen[k] = row[4]
        keep_rows.append(row)
    crop_rows = sorted(keep_rows, key=lambda x: (x[0], str(x[4] or x[5])))
    # C3c: an annual crop with no harvest record takes the harvest date most often recorded
    # for the same crop, planted in the same year, on the site's other plots (as done by
    # hand for NECCIRR); the yield stays blank
    from collections import Counter
    hd = Counter()
    for row in crop_rows:
        if row[1] == "Annual" and row[4] is not None and row[5] is not None:
            hd[(row[3], row[4].year, row[5])] += 1
    best = {}
    for (sp, yr, d), n in hd.items():
        if (sp, yr) not in best or n > best[(sp, yr)][1]: best[(sp, yr)] = (d, n)
    for row in crop_rows:
        if row[1] == "Annual" and row[4] is not None and row[5] is None and (row[3], row[4].year) in best:
            d = best[(row[3], row[4].year)][0]
            if d > row[4]:
                row[5] = d
                row[10] = (row[10] or "") + "; harvest date from the site's other plots (no record for this plot)"
                L("C3c", row[0], f"{row[4].date()} {row[3]}: no harvest record; site's most common harvest date {d.date()} used")
    # C4: harvests with no planting. Without a yield (a sampling) it is ignored; with one,
    # the planting date is the site's most common planting date for that crop and year
    pdate = Counter()
    for row in crop_rows:
        if row[1] == "Annual" and row[4] is not None: pdate[(row[3], row[4].year, row[4])] += 1
    best_p = {}
    for (sp, yr, d), n in pdate.items():
        if (sp, yr) not in best_p or n > best_p[(sp, yr)][1]: best_p[(sp, yr)] = (d, n)
    planted = {(row[0], row[4].year) for row in crop_rows if row[1] == "Annual" and row[4] is not None}
    n_samp = 0
    for _, h in hv.iterrows():
        if (h.unit_id, h.sp, h.d) in used or h.sp is None: continue
        mc = combo.get((h.unit_id, h.treatment_id))
        if mc is None: continue
        if pd.isna(h.y):
            n_samp += 1; continue
        cat = "Perennial" if h.sp in PERENNIAL else "Annual"
        # not when the plot already has an annual crop planted that year (e.g. NEMLTCRS 1993:
        # an oat harvest on plots that grew corn), which would add a crop the plot did not grow
        p = best_p.get((h.sp, h.d.year)) if cat == "Annual" and (mc, h.d.year) not in planted else None
        p = p[0] if p is not None and p[0] < h.d.to_pydatetime() else None
        crop_rows.append([mc, cat, None, h.sp, p, h.d.to_pydatetime(), None,
                          round(float(h.y), 1), "kgs/hectare", None,
                          "harvest record without planting record" +
                          (f"; planting date from the site's other plots ({p.date()})" if p is not None else "")])
        L("C4", mc, f"{h.d.date()} {h.sp}: harvest without planting; " +
          (f"site's most common planting date {p.date()} used" if p is not None else "plant date blank"))
    if n_samp: L("C4", "", f"{n_samp} harvest-type records without a yield and not paired with a planting (samplings) ignored")

    # --- Soil_Disturbance
    sd_rows = []
    ti = T["mgttillage"].copy()
    if len(ti):
        ti["d"] = to_date(ti.start_date)
        for _, r in ti.dropna(subset=["d"]).iterrows():
            mc = combo.get((r.unit_id, r.treatment_id))
            if mc is None: continue
            hit = None if keep_documented_errors else next(
                (ev for lc, tre, ere, ev in TILLAGE_NOT_DONE if lc == loc and re.search(tre, str(r.treatment_id))
                 and re.search(ere, f"{r.tillage_event} / {r.tillage_event_method}", re.I)), None)
            if hit:
                L("D4", mc, f"{r.d.date()} {r.tillage_event} / {r.tillage_event_method}: not written ({hit})"); continue
            imp, note = tdisc_for(r.tillage_event, r.tillage_event_method, r.d)
            if imp is None:
                L("D1", mc, f"{r.d.date()} {r.tillage_event} / {r.tillage_event_method}: no soil movement, not written"); continue
            if imp == "UNMAPPED":
                L("D1", mc, f"{r.d.date()} UNMAPPED tillage: {r.tillage_event} / {r.tillage_event_method}")
            dep = pd.to_numeric(r.get("tillage_event_depth_cm"), errors="coerce")
            dep_in = round(dep / 2.54, 1) if pd.notna(dep) and dep > 0 else None
            n = 2 if re.search(r"twice", f"{r.tillage_event} {r.tillage_event_method}", re.I) else 1
            for _ in range(n):
                sd_rows.append([mc, r.d.to_pydatetime(), imp, None, None, None, dep_in, None,
                                f"AgCROS: {r.tillage_event} / {r.tillage_event_method}, {dep if pd.notna(dep) else 'NA'} cm" +
                                (f"; {note}" if note else "")])
    if planter:
        for _, r in pl.iterrows():
            mc = combo.get((r.unit_id, r.treatment_id))
            if mc is None: continue
            no_till = "no till" in str(trt.loc[r.treatment_id, "tillage_descriptor"]).lower() if r.treatment_id in trt.index else False
            imp = "PLANTER, DRILL, NO-TILL, MIN-TILL" if str(r.planting_method).lower() == "drill" and no_till else \
                  "PLANTER, DRILL, DISCGRAIN" if str(r.planting_method).lower() == "drill" else \
                  "PLANTER, NO-TILL" if no_till else "PLANTER, REGULAR"
            if str(r.planting_method).lower() == "broadcast": continue
            sd_rows.append([mc, r.d.to_pydatetime(), imp, None, None, None, None, None, f"planting pass ({r.planting_method})"])

    # D3: if the whole site has no disturbance record (a no-till study), one zero-mixing
    # pass is written for its first unit, so the Soil_Disturbance sheet is not empty (an
    # empty sheet stops SHMI <= 1.2.x; see SHMI_package_todo.md). Mixing efficiency 0:
    # it does not change Inverse disturbance. Remove once the package is fixed.
    if not sd_rows and crop_rows:
        first = min((r for r in crop_rows if r[4] is not None), key=lambda r: r[4], default=None)
        if first is not None:
            sd_rows.append([first[0], first[4], "SPRAYER/AERIAL/BROADCAST APPLICATION", None, None, None, None, None,
                            "no tillage recorded at this site (no-till); zero-disturbance placeholder so the sheet is not empty"])
            L("D3", first[0], "site has no tillage records; one zero-disturbance placeholder pass written")

    # --- Soil_Amendments
    sa_rows = []
    am = T["mgtamendments"].copy()
    if len(am):
        am["d"] = to_date(am.start_date)
        for _, r in am.dropna(subset=["d"]).iterrows():
            mc = combo.get((r.unit_id, r.treatment_id))
            if mc is None: continue
            typ = str(r.amend_type) if pd.notna(r.amend_type) else ""
            N = pd.to_numeric(r.total_n_amount_kgn_ha, errors="coerce")
            rate = pd.to_numeric(r.total_amend_amount_kg_ha, errors="coerce")
            if NOT_SCORED.search(typ):
                L("M2", mc, f"{r.d.date()} {typ} {'' if pd.isna(rate) else rate} kg/ha: not scored as an organic amendment, not written"); continue
            if ORGANIC.search(typ):
                sa_rows.append([mc, r.d.to_pydatetime(), "Organic", typ, None if pd.isna(rate) else float(rate), "kgs/hectare",
                                None if pd.isna(N) else float(N), None, None, None, None, "AgCROS amendment"])
            elif pd.notna(N) and N >= 0 and not NOT_N.search(typ):
                sa_rows.append([mc, r.d.to_pydatetime(), "Fertilizer", typ or "N fertilizer (type not recorded)",
                                None if pd.isna(rate) else float(rate), "kgs/hectare", float(N), None, None, None, None, "AgCROS amendment"])
                if not typ: L("M1", mc, f"{r.d.date()} N {N} kg/ha with no amendment type: kept as fertilizer")

    # --- Animal_Diversity
    ad_rows = []
    gz = T["mgtgrazing"]
    if len(gz):
        for _, r in gz.iterrows():
            mc = combo.get((r.unit_id, r.treatment_id))
            if mc is None: continue
            ad_rows.append([mc, str(r.animal_species), to_date(r.start_date).to_pydatetime(), to_date(r.end_date).to_pydatetime(),
                            None, r.stocking_rate_noanimals_ha if pd.notna(r.stocking_rate_noanimals_ha) else None, None, "AgCROS grazing"])

    # --- S1: supplement records (records missing from AgCROS, each with its source)
    for f in supplements:
        S = pd.read_csv(f, dtype=str, keep_default_na=False)
        S = S[S.location_id == loc]
        n_added = 0
        for _, s in S.iterrows():
            uid = s.unit_id if s.unit_id.startswith(loc + "_") else f"{loc}_{s.unit_id}"
            d = to_date(s.date)
            if pd.isna(d):
                L("S1", uid, f"supplement row without a valid date skipped ({f})"); continue
            cand = units[units.unit_id == uid]
            if not len(cand):
                L("S1", uid, f"supplement row for a unit not in AgCROS skipped ({f})"); continue
            inside = cand[(cand["min"] - timedelta(days=366) <= d) & (d <= cand["max"] + timedelta(days=366))]
            u = (inside if len(inside) else cand).sort_values("size", ascending=False).iloc[0]
            mc = combo[(u.unit_id, u.treatment_id)]
            num = lambda k: (float(s[k]) if k in s and str(s[k]).strip() not in ("", "nan", "None") else None)
            txt = lambda k: (s[k] if k in s and str(s[k]).strip() not in ("", "nan", "None") else None)
            note = "; ".join(x for x in [txt("notes"), f"supplement: {txt('source')}" if txt("source") else "supplement"] if x)
            sheet = s.sheet
            if sheet == "Soil_Amendments":
                sa_rows.append([mc, d.to_pydatetime(), txt("SA_cat"), txt("SA_source"), num("SA_rate"), txt("SA_units"),
                                num("SA_N"), num("SA_P2O5"), num("SA_K2O"), num("SA_H2O_per"), num("SA_CN"), note])
            elif sheet == "Soil_Disturbance":
                sd_rows.append([mc, d.to_pydatetime(), txt("SD_equip"), None, None, None, num("SD_depth"), None, note])
            elif sheet == "Crop_Diversity":
                hd = to_date(txt("CD_harv_date")) if txt("CD_harv_date") else None
                td = to_date(txt("CD_term_date")) if txt("CD_term_date") else None
                crop_rows.append([mc, txt("CD_cat"), None, txt("CD_name"), d.to_pydatetime(),
                                  hd.to_pydatetime() if hd is not None and pd.notna(hd) else None,
                                  td.to_pydatetime() if td is not None and pd.notna(td) else None,
                                  num("CD_yield"), txt("CD_yield_units"), None, note])
            else:
                L("S1", uid, f"unknown sheet '{sheet}' in supplement; row skipped"); continue
            n_added += 1
        if n_added: L("S1", "", f"{n_added} supplement rows added from {os.path.basename(f)}")
        print(f"{loc}: {n_added} supplement rows from {os.path.basename(f)}")

    # --- write workbook
    # T1: one CD_term_date so the column is not entirely blank (see the rule list)
    if crop_rows and all(r[6] is None for r in crop_rows):
        first = next((r for r in crop_rows if r[5] is not None and str(r[1]).lower() == "annual"), None)
        if first is not None:
            first[6] = first[5]
            L("T1", first[0], "no termination dates at this site; harvest date copied to CD_term_date of one crop")

    wb = openpyxl.load_workbook(template)
    def put(sheet, rows, first_col=2):
        ws = wb[sheet]
        for i, row in enumerate(rows):
            for j, v in enumerate(row):
                c = ws.cell(row=5 + i, column=first_col + j)
                c.value = None if (isinstance(v, float) and pd.isna(v)) else v
                if hasattr(v, "year") and not isinstance(v, str): c.number_format = "mm/dd/yyyy"
        # clear template formula rows beyond the data in the key column
        for r in range(5 + len(rows), ws.max_row + 1):
            if isinstance(ws.cell(row=r, column=first_col).value, str) and ws.cell(row=r, column=first_col).value.startswith("="):
                ws.cell(row=r, column=first_col).value = None
    put("Mgt_Unit", mgt_rows)
    put("Crop_Diversity", crop_rows)
    # SD autofill formula columns (E,F) are left to Excel; the package looks up SD_equip itself
    ws = wb["Soil_Disturbance"]
    for i, row in enumerate(sd_rows):
        r = 5 + i
        vals = [row[0], row[1], row[2]]
        for j, v in enumerate(vals): ws.cell(row=r, column=2 + j).value = v
        ws.cell(row=r, column=3).number_format = "mm/dd/yyyy"
        ws.cell(row=r, column=5).value = f'=IF(D{r}="","",IFERROR(VLOOKUP(D{r},TDISC_Lookup,2,FALSE),"operation name"))'
        ws.cell(row=r, column=6).value = f'=IF(D{r}="","",IFERROR(ROUND(VLOOKUP(D{r},TDISC_Lookup,4,FALSE),1),""))'
        ws.cell(row=r, column=8).value = row[6]
        ws.cell(row=r, column=10).value = row[8]
    put("Soil_Amendments", sa_rows)
    put("Animal_Diversity", ad_rows)
    os.makedirs(outdir, exist_ok=True)
    out = os.path.join(outdir, f"SHMI_{loc}.xlsx")
    wb.save(out)
    pd.DataFrame(log).to_csv(os.path.join(outdir, f"SHMI_{loc}_log.csv"), index=False)
    pd.DataFrame(unit_rows).to_csv(os.path.join(outdir, f"SHMI_{loc}_units.csv"), index=False)
    unm = sum(1 for x in log if "UNMAPPED" in x["message"])
    print(f"{loc}: {len(mgt_rows)} units, {len(crop_rows)} crop rows, {len(sd_rows)} passes, {len(sa_rows)} amendments, "
          f"{len(ad_rows)} grazing rows; {len(log)} log entries ({unm} unmapped tillage) -> {out}")
    return out

STATE = {"NE": "Nebraska", "MN": "Minnesota", "IA": "Iowa", "IN": "Indiana", "CO": "Colorado", "PA": "Pennsylvania",
         "SD": "South Dakota", "ND": "North Dakota", "SC": "South Carolina", "KY": "Kentucky", "AL": "Alabama",
         "WI": "Wisconsin", "OH": "Ohio", "MD": "Maryland", "TX": "Texas", "GA": "Georgia", "MT": "Montana",
         "WY": "Wyoming", "ID": "Idaho", "WA": "Washington", "OR": "Oregon", "IL": "Illinois", "MO": "Missouri",
         "KS": "Kansas", "MI": "Michigan", "NY": "New York"}

if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("tables"); ap.add_argument("locations", nargs="+")
    ap.add_argument("--template", required=True); ap.add_argument("--out", default=".")
    ap.add_argument("--study", default="AgCROS"); ap.add_argument("--planter-passes", action="store_true")
    ap.add_argument("--corrections", action="append", default=[],
                    help="CSV of corrections to AgCROS records (repeatable): location_id, table, unit_id, new_unit_id, new_treatment_id, from_date, to_date, evidence")
    ap.add_argument("--supplement", action="append", default=[],
                    help="CSV of records missing from AgCROS (repeatable); see make_supplement.py")
    ap.add_argument("--keep-documented-errors", action="store_true",
                    help="keep AgCROS tillage records that the site's publications say did not happen (rule D4)")
    a = ap.parse_args()
    for loc in a.locations:
        convert(a.tables, loc, a.template, a.out, a.study, a.planter_passes, a.supplement, a.corrections,
                a.keep_documented_errors)
