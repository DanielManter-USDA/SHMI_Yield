"""
build_study_treatments.py
Builds SHMI_study_treatments.xlsx (Supplementary table of sites and treatments) from the
SHMI study workbooks used by SHMI_yield_main.R, over each site's analysis period.

    python3 build_study_treatments.py --wb-dir DRIVES/final --units-dir DRIVES/final \
        --old SHMI_study_treatments_7sites.xlsx --out SHMI_study_treatments.xlsx

--old: the earlier table, whose descriptive text (rotation, tillage, N, cover, organic,
irrigation, other factors, notes) is reused for sites whose codes did not change
(Ohio, WICST, NEMEIRR, NECCIRR). Counts and tillage passes are always recomputed.
"""
import argparse
import sys, os, re, warnings
import pandas as pd

warnings.filterwarnings("ignore")

# ---- sites: period, exclusions and treatment code, as in SHMI_yield_main.R ------------
SITES = {
    "NELITCSE":  dict(start=1986, end=2013, loc="Lincoln, Nebraska"),
    "NEMLTCRS":  dict(start=1983, end=2006, loc="Mead, Nebraska"),
    "OHTVD_Ho":  dict(start=1964, end=2012, loc="Hoytville, Ohio"),
    "OHTVD_Wo":  dict(start=1962, end=2021, loc="Wooster, Ohio"),
    "WICST_ARS": dict(start=1990, end=2020, loc="Arlington, Wisconsin"),
    "NEMEIRR":   dict(start=2001, end=2015, loc="Mead, Nebraska"),
    "NECCIRR":   dict(start=2010, end=2018, loc="Clay Center, Nebraska"),
    "INWLACRE":  dict(start=2004, end=2009, loc="West Lafayette, Indiana (ACRE)", agcros=True),
    "IAAM7071":  dict(start=2008, end=2012, loc="Ames, Iowa", agcros=True),
    "PAUP":      dict(start=2008, end=2012, loc="University Park, Pennsylvania", agcros=True),
}

DESIGN = {
    "NELITCSE": "3 rotations x 6 tillage systems, optimum N; 600 block excluded (no harvests 1986-1995)",
    "NEMLTCRS": "7 rotations x 3 N rates, conventional tillage (first era, 1983-2006; no-till from 2007 not used)",
    "OHTVD_Ho": "3 rotations x 3 tillage systems, optimum N",
    "OHTVD_Wo": "3 rotations x 3 tillage systems, optimum N",
    "WICST_ARS": "5 cropping systems (conventional and organic grain and forage)",
    "NEMEIRR": "Irrigated continuous corn: no-till vs disk x 0, 50, 100% stover removal",
    "NECCIRR": "Irrigated no-till continuous corn: full vs deficit irrigation x amelioration (none, rye cover, manure) x N x stover removal",
    "INWLACRE": "Corn-soybean: conventional till vs no-till, no-till with ryegrass cover, split N, modified planter; perennial plots (no corn)",
    "INWLTPAC": "Corn-soybean: conventional till vs no-till, no-till with rye cover, split N, modified planter; sorghum-sudangrass plots (no corn). Ends before the 2008 redesign",
    "IAAM7071": "Continuous corn: chisel vs no-till x stover removal (0, medium, high) x N management; rye cover; clover living mulch; biochar",
    "PAUP": "No-till continuous corn and corn-soybean x stover removal (0, 50, 100%) x rye cover; twin rows",
}
RECORDS = {
    "NELITCSE": "AgCROS (GRACEnet), converted with agcros2shmi.py",
    "NEMLTCRS": "AgCROS, converted with agcros2shmi.py; 2 corrections (115S/115M N, 104N/104S labels)",
    "NEMEIRR": "AgCROS, converted with agcros2shmi.py",
    "NECCIRR": "AgCROS, converted with agcros2shmi.py; fertilizer and manure from a supplement (UNL CropWatch 2017)",
    "INWLACRE": "AgCROS, converted with agcros2shmi.py",
    "INWLTPAC": "AgCROS, converted with agcros2shmi.py",
    "IAAM7071": "AgCROS, converted with agcros2shmi.py",
    "PAUP": "AgCROS, converted with agcros2shmi.py (no tillage records: no-till)",
}

# ---- T-DISC implement names, shortened ----------------------------------------------
IMPL = {
    "BED SHAPING, HIPPER": "hipper", "CULTIVATOR, DUCKFOOT": "duckfoot cultivator",
    "CULTIVATOR, FIELD": "field cultivator", "CULTIVATOR, ROTARY HOE": "rotary hoe",
    "CULTIVATOR, ROTOVATOR-BEDDER": "rotovator", "CULTIVATOR, ROW": "row cultivator",
    "HARROW, DISK, CHISEL (MULCH TILLER)": "disk chisel", "HARROW, DISK, TANDEM, HEAVYDUTY": "tandem disk",
    "HARROW, FLEX-TINE": "flex-tine harrow", "HARROW, SPIKE TOOTH OR SPRING TOOTH": "spike-tooth harrow",
    "NUTRIENT MGMT: FERT APPL, STRIP TILL": "strip-till applicator", "PLOW, CHISEL": "chisel plow",
    "PLOW, COULTER-CHISEL": "coulter chisel", "PLOW, DISK": "disk plow", "PLOW, MOLDBOARD": "moldboard plow",
    "RESIDUE MGMT: STALK PULLER CHOPPER": "stalk chopper", "RIPPER/SUBSOILER": "subsoiler",
    "ROLLER, CULTI-PACKER/PULVERIZER": "cultipacker", "ROTO-TILLER, ROTO-TILLER DIGGER, ROTERRA": "rototiller",
    "Reel disk vertical tiller": "vertical tiller",
}
PLACEHOLDER = "SPRAYER/AERIAL/BROADCAST APPLICATION"   # agcros2shmi rule D3 (no-till site)

# ---- codes and descriptions for the sites whose codes are AgCROS codes ----------------
NEL_TILL = {"NT": "no-till", "CH": "chisel plow (fall, 20 cm)", "DK": "tandem disk (twice)",
            "RT": "ridge-till", "MP": "moldboard plow (20 cm)", "SS": "subsoil (25 cm)"}
NEL_ROT = {"CC": "continuous corn", "CS": "corn-soybean", "SB": "continuous soybean"}
NEM_ROT = {"1": "continuous corn", "2": "continuous soybean", "3": "continuous sorghum",
           "4": "corn-soybean", "5": "sorghum-soybean",
           "6": "corn - oats/clover - sorghum - soybean", "7": "corn - soybean - sorghum - oats/clover"}
NEM_N = {"0": "0 kg N/ha", "1": "low: 90 kg N/ha corn and sorghum; 34 soybean and oats",
         "2": "high: 180 kg N/ha corn and sorghum; 69 soybean and oats"}
# "other" factors as in the analysis as run (B2, AGCROS_REFINE = FALSE): the modified planter
# at INWLACRE, plus stover removal at all three; biochar and twin rows are not modelled
AGCROS_OTHER = {"INWLACRE": {"modified planter": "no-till planter with row cleaners"}}


def code_for(site, trt):
    trt = str(trt)
    if site == "NELITCSE":
        m = re.match(r"^(CH|DK|MP|NT|RT|SS)(CC|CSB|SBC|SB)$", trt)
        return None if not m else m.group(1) + ("CS" if m.group(2) in ("CSB", "SBC") else m.group(2))
    if site in ("OHTVD_Ho", "OHTVD_Wo", "WICST_ARS"):
        return re.sub(r"P[0-9]+$", "", trt)
    return trt


def read(wb, sheet):
    d = pd.read_excel(wb, sheet_name=sheet, skiprows=3)
    return d.drop(columns=[c for c in d.columns if str(c).startswith("Unnamed")])


def year_of(df, cols):
    y = pd.Series(pd.NA, index=df.index, dtype="Int64")
    for c in cols:
        if c in df:
            v = pd.to_datetime(df[c], errors="coerce").dt.year.astype("Int64")
            y = y.fillna(v)
    return y


def site_rows(site, cfg, wb_dir, units_dir, old):
    wb = os.path.join(wb_dir, f"SHMI_{site}.xlsx")
    s0, s1 = cfg["start"], cfg["end"]
    mgt = read(wb, "Mgt_Unit"); cd = read(wb, "Crop_Diversity")
    sd = read(wb, "Soil_Disturbance"); sa = read(wb, "Soil_Amendments")
    mgt["field"] = mgt.MGT_field.astype(str).str.replace(r"\.0$", "", regex=True)
    if site == "NELITCSE":   # EXCLUDE_NEL_600
        mgt = mgt[~pd.to_numeric(mgt.field, errors="coerce").between(601, 624)]

    if cfg.get("agcros"):
        u = pd.read_csv(os.path.join(units_dir, f"SHMI_{site}_units.csv")).fillna("")
        u["system"] = (u.treatment_descriptor + " | " + u.rotation_descriptor + " | " + u.tillage_descriptor + " | " +
                       u.n_treatment_descriptor + " | " + u.cover_crop.replace("", "no cover") + " | residue removal " +
                       u.residue_removal + " | " + u.organic_management.map({"Yes": "organic"}).fillna("conventional"))
        mgt = mgt.merge(u[["MGT_combo", "treatment_id", "system", "treatment_descriptor", "rotation_descriptor",
                           "tillage_descriptor", "n_treatment_descriptor", "cover_crop", "residue_removal",
                           "organic_management", "irrigation"]], on="MGT_combo", how="left")
        mgt["code"] = mgt.system
    else:
        mgt["code"] = [code_for(site, t) for t in mgt.MGT_trt]
    mgt = mgt[mgt.code.notna()]

    cd["year"] = year_of(cd, ["CD_harv_date", "CD_term_date", "CD_plant_date"])
    cd = cd[cd.MGT_combo.isin(mgt.MGT_combo) & cd.year.between(s0, s1)]
    active = mgt[mgt.MGT_combo.isin(cd.MGT_combo)]
    corn = cd[cd.CD_name.astype(str).str.strip().str.lower().eq("corn") & cd.CD_yield.notna()]
    corn = corn.drop_duplicates(["MGT_combo", "year"]).merge(active[["MGT_combo", "code", "field"]], on="MGT_combo")

    sd["year"] = year_of(sd, ["SD_date"])
    sd = sd[sd.MGT_combo.isin(active.MGT_combo) & sd.year.between(s0, s1) & sd.SD_equip.ne(PLACEHOLDER)]
    sd = sd.merge(active[["MGT_combo", "code"]], on="MGT_combo")
    sa["year"] = year_of(sa, ["SA_date"])
    fert = sa[sa.SA_N.notna() & sa.SA_cat.astype(str).str.lower().eq("fertilizer")].groupby(["MGT_combo", "year"]).SA_N.sum()
    corn["N_fert"] = [fert.get((m, y), float("nan")) for m, y in zip(corn.MGT_combo, corn.year)]

    nyears = s1 - s0 + 1
    rows = []
    for code, g in active.groupby("code", sort=True):
        plots = g.field.nunique()
        c = corn[corn.code == code]
        t = sd[sd.code == code]
        # passes per plot-year: plot-years = each unit's years with crops in the period
        py = cd[cd.MGT_combo.isin(g.MGT_combo)].groupby("MGT_combo").year.agg(lambda y: y.max() - y.min() + 1).sum()
        py = max(py, 1)
        passes = []
        for eq, e in t.groupby("SD_equip"):
            txt = f"{IMPL.get(eq, eq.lower())} {len(e) / py:.2f}"
            if "SD_depth" in e and e.SD_depth.notna().any():
                txt += f" (recorded depth {e.SD_depth.mean():.1f} in)"
            passes.append((len(e), txt))
        passes = "; ".join(x for _, x in sorted(passes, reverse=True)) or "none"
        nrec = c.N_fert.notna().mean() if len(c) else float("nan")
        n_txt = f"{c.N_fert.mean():.0f} kg N/ha (fertilizer on corn, records; {nrec:.0%} of corn plot-years)" \
            if len(c) and c.N_fert.notna().any() else ""
        r = dict(Site=site, **{"Treatment code": code, "Key (trt in the R outputs)": f"{site} {code}"}, Plots=plots,
                 **{"Corn plot-years with yield": len(c), "Years with corn yield": c.year.nunique(),
                    "Tillage passes recorded (implement: passes per plot-year)": passes, "_n": n_txt})
        r.update(describe(site, code, g, old))
        if len(c) == 0:
            r["Notes"] = "; ".join(x for x in [r.get("Notes", ""), "No corn: not in the corn-yield analysis."] if x)
        rows.append(r)
    return rows


def describe(site, code, g, old):
    d = {}
    if site in ("OHTVD_Ho", "OHTVD_Wo", "WICST_ARS", "NEMEIRR", "NECCIRR"):
        # hand-built (DRIVES) and supplemented sites: descriptions come only from --old
        if old is None:
            sys.exit("build_study_treatments.py: --old (data/SHMI_study_treatments_descriptions.xlsx) is needed "
                     f"for the descriptions of {site}")
        o = old[(old.Site == site) & (old["Treatment code"] == code)]
        if len(o):
            o = o.iloc[0]
            for k in ["Rotation", "Tillage treatment", "N fertility", "Cover crop", "Organic inputs",
                      "Irrigation", "Other factors (not scored by SHMI)", "Tillage source", "Notes"]:
                d[k] = "" if pd.isna(o[k]) else o[k]
            return d
    base = {"Cover crop": "none", "Organic inputs": "none", "Irrigation": "none",
            "Other factors (not scored by SHMI)": "none", "Notes": ""}
    if site == "NELITCSE":
        base.update({"Rotation": NEL_ROT[code[2:]], "Tillage treatment": NEL_TILL[code[:2]],
                     "N fertility": "optimum (site practice)", "Tillage source": "records (AgCROS)"})
        return base
    if site == "NEMLTCRS":
        r, n = code[3], code[4]
        base.update({"Rotation": NEM_ROT[r], "Tillage treatment": "conventional (tandem disk, row cultivation)",
                     "N fertility": NEM_N[n], "Tillage source": "records (AgCROS)",
                     "Cover crop": "sweet and red clover underseeded in oats" if r in "67" else "none"})
        return base
    # AgCROS sites: from the AgCROS descriptors
    x = g.iloc[0]
    till = x.tillage_descriptor
    if x.treatment_descriptor and x.treatment_descriptor.lower() not in ("no-till", "no till"):
        till = f"{till} ({x.treatment_descriptor.rstrip('.').replace('chisel.spring', 'chisel, spring')})"
    other = []
    if x.residue_removal in ("Partial", "Yes"):
        other.append("stover removal: " + x.residue_removal.lower())
    for k, lab in AGCROS_OTHER.get(site, {}).items():
        if k in x.system.lower():
            other.append(lab)
    if x.organic_management == "Yes":
        other.append("organic management")
    base.update({"Rotation": x.rotation_descriptor, "Tillage treatment": till,
                 "N fertility": x.n_treatment_descriptor or "not stated",
                 "Cover crop": x.cover_crop or ("perennial cover crop (living mulch), per the treatment description"
                                                if "perennial cover crop" in str(x.treatment_descriptor).lower() else "none"),
                 "Irrigation": "none" if x.irrigation in ("", "No") else x.irrigation,
                 "Other factors (not scored by SHMI)": "; ".join(other) or "none",
                 "Tillage source": "records (AgCROS)" if site != "PAUP" else "no tillage (no-till study)",
                 "Notes": "AgCROS treatment(s): " + ", ".join(sorted(g.treatment_id.unique(), key=lambda s: (len(s), s)))})
    if site == "IAAM7071" and "charcoal" in x.system.lower():
        base["Notes"] += "; biochar, not scored by SHMI (not an organic amendment; no amendment bonus) and not modelled separately"
    if site == "PAUP" and "twin row" in x.system.lower():
        base["Notes"] += "; twin-row planting (not modelled)"
    return base


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--wb-dir", required=True, help="folder with every SHMI_<site>.xlsx")
    ap.add_argument("--units-dir", required=True, help="folder with SHMI_<site>_units.csv for the AgCROS sites")
    ap.add_argument("--old", help="earlier SHMI_study_treatments.xlsx (descriptions reused)")
    ap.add_argument("--out", default="SHMI_study_treatments.xlsx")
    a = ap.parse_args()
    old = pd.read_excel(a.old, sheet_name="Treatments") if a.old else None
    old_sites = pd.read_excel(a.old, sheet_name="Sites") if a.old else None

    rows = []
    for site, cfg in SITES.items():
        rows += site_rows(site, cfg, a.wb_dir, a.units_dir, old)
        print(f"{site}: {sum(r['Site'] == site for r in rows)} treatments")
    t = pd.DataFrame(rows)
    # recorded N appended to the design N text
    t["N fertility"] = [f"{n}; {r}" if r and not n.startswith("0 kg") else n for n, r in zip(t["N fertility"], t["_n"])]
    t = t.drop(columns="_n")
    # AgCROS systems: short codes instead of the long system key
    ag = t.Site.isin([s for s, c in SITES.items() if c.get("agcros")])
    t.loc[ag, "Treatment code"] = t.loc[ag, "Notes"].str.extract(r"AgCROS treatment\(s\): ([^;]+)")[0] \
        .str.replace(r"\b[A-Z0-9]+_", "", regex=True)
    cols = ["Site", "Treatment code", "Key (trt in the R outputs)", "Rotation", "Tillage treatment",
            "Tillage passes recorded (implement: passes per plot-year)", "N fertility", "Cover crop",
            "Organic inputs", "Irrigation", "Other factors (not scored by SHMI)", "Plots",
            "Corn plot-years with yield", "Years with corn yield", "Tillage source", "Notes"]
    t = t[cols].sort_values(["Site", "Treatment code"], key=lambda s: s.map(list(SITES).index) if s.name == "Site" else s)

    srow = []
    for site, cfg in SITES.items():
        g = t[t.Site == site]
        rec = RECORDS.get(site)
        if rec is None and old_sites is not None:
            o = old_sites[old_sites.Site == site]
            rec = o["Records used"].iloc[0] if len(o) else ""
        srow.append({"Site": site, "Location": cfg["loc"], "Analysis period": f"{cfg['start']}-{cfg['end']}",
                     "Design": DESIGN[site], "Records used": rec, "Treatments": len(g), "Plots": int(g.Plots.sum()),
                     "Corn plot-years with yield": int(g["Corn plot-years with yield"].sum()),
                     "Treatments with corn": int((g["Corn plot-years with yield"] > 0).sum())})
    s = pd.DataFrame(srow)
    s = pd.concat([s, pd.DataFrame([{"Site": "Total", "Treatments": s.Treatments.sum(), "Plots": s.Plots.sum(),
                                     "Corn plot-years with yield": s["Corn plot-years with yield"].sum(),
                                     "Treatments with corn": s["Treatments with corn"].sum()}])], ignore_index=True)

    notes = [
        "Treatment summary for the SHMI corn-yield analysis (ten long-term experiments)", "",
        "Counts and tillage passes are computed from the SHMI study workbooks (SHMI_<site>.xlsx) over each site's analysis period (build_study_treatments.py).",
        "Key: the treatment as named in the R outputs (trt = site + MGT_trt_base; e.g. treatment_labels.csv, N_by_treatment.csv), for joining. AgCROS sites: the system string B2 builds from the AgCROS descriptors.",
        "Treatment code: the treatment part of MGT_combo without the rotation phase. AgCROS sites (Indiana, Iowa, Pennsylvania): the AgCROS treatment number(s) of one system; treatments that differ only in rotation phase are listed together.",
        "Plots: plots with crops in the analysis period. Corn plot-years: corn crops with a yield in the period (the yield model may use slightly fewer).",
        "Tillage passes: passes per plot per year in the analysis period, by implement (T-DISC implement names shortened). Recorded depth: mean working depth where the records give one.",
        "N fertility: from the treatment design; the mean recorded fertilizer N on corn is added where records exist (manure N is computed by the SHMI package and not shown).",
        "Other factors: treatment factors SHMI does not score, entered as an 'other' covariate in the yield models.",
        "", "Records and assumptions",
        "Nebraska (NELITCSE, NEMLTCRS, NEMEIRR, NECCIRR): workbooks built from the AgCROS tables by agcros2shmi.py; they match the hand-built workbooks field by field (see agcros2shmi/README.md).",
        "NEMLTCRS: first era only (1983-2006, conventional tillage). Two AgCROS errors corrected (115S/115M N swapped; 104N/104S labels).",
        "NECCIRR: no amendment records in AgCROS; design N at corn planting and manure every other fall after harvest (UNL CropWatch, 2017) from a supplement.",
        "OHTVD_Ho: DRIVES tillage corrected with the Ohio management logs (duplicates and phase misassignment removed).",
        "OHTVD_Wo: DRIVES records only the tillage class; passes assumed from the Hoytville protocol. Coordinates approximate.",
        "WICST_ARS: implements, passes and depths from the station tillage log; all passes match DRIVES by plot and date.",
        "INWLACRE: 2004-2009, including the 2008 redesign (split N became reduced N, 175 vs 220 lbs/A; the modified planter became intermittent tillage). N is the recorded amount; the modified planter enters as an 'other' factor. A sensitivity version (AGCROS_REFINE in B2) ends in 2007 and sets every corn plot to the site-year median N.",
        "IAAM7071: stover removal enters as an 'other' factor. Biochar (Charcoal1, Charcoal2) is not scored by SHMI: it is stable pyrolysed carbon, not an organic amendment, so those plots score as their no-biochar controls. Clover living-mulch plots (21, 22): perennial cover 2008-2010, then corn-rye-soybean (21) and rye-soybean-corn (22) from 2011 (Ag Data Commons 10.15482/USDA.ADC/1528303).",
        "PAUP: no tillage records (no-till study); stover removal enters as an 'other' factor; twin-row planting is not modelled.",
        "", "Screened out",
        "MNMOFS (Morris, MN): fertilizer N missing in 2008-09 and only starter N recorded on many plots; organic plots have no manure N in AgCROS. Run as a sensitivity analysis only.",
        "INWLTPAC (West Lafayette, TPAC): split-N plots became reduced-N plots in 2008 and N records are incomplete (split N records one of two applications; modified planter about half the design N).",
        "Other AgCROS locations were screened out for having no SHMI contrast (residue-removal-only studies), no corn treatment contrasts, too few plots, or locations needing to be split (agcros2shmi/agcros_location_screen.csv).",
        "", "Flags",
        "OHTVD_Ho chisel plots (T2) record moldboard plowing in many years: worth confirming with the site.",
        "IAAM7071: AgCROS logs a mid-season sampling (about 20 June) like a harvest and has no cover-crop plantings. The converter (agcros2shmi.py, 7 Nov 2026) pairs each planting with its grain harvest and ignores samplings; the covers come from supplements/IAAM7071_covers.csv: cereal rye for treatments 19 and 20 with the planting and harvest dates supplied by the site (2008-09 to 2012-13), and a perennial white-clover living mulch assumed for treatments 21 and 22, planted with the first corn in 2008 and kept for the whole experiment (to be confirmed with the site).",
        "INWLACRE: AgCROS has no 2004 corn planting record for the four modified-planter plots; the converter gives them the date of the site's other 2004 corn plantings (19 April 2004).",
        "NELITCSE: AgCROS has no 2014 corn yields (the DRIVES export has them), so the converted workbook's corn yields end in 2013 (60 fewer corn plot-years than the hand-built workbook).",
        "INWLACRE: AgCROS records less N on some plots than the design (modified planter about 106 kg N/ha; split N one of its two applications); the analysis as run uses the recorded amounts.",
        "Sources: AgCROS (GRACEnet/REAP) database backup; DRIVES database export (yield_plant_harv_trt); WICST data request (6-5-2026); Ohio management logs; UNL CropWatch (2017).",
    ]
    with pd.ExcelWriter(a.out, engine="openpyxl") as w:
        s.to_excel(w, sheet_name="Sites", index=False)
        t.to_excel(w, sheet_name="Treatments", index=False)
        pd.DataFrame({notes[0]: notes[1:]}).to_excel(w, sheet_name="Notes", index=False)
        fmt(w.book)
    print("Saved", a.out)


def fmt(book):
    from openpyxl.styles import Alignment, Font, PatternFill
    widths = {"Sites": [11, 26, 10, 60, 50, 11, 8, 12, 12],
              "Treatments": [11, 14, 30, 26, 36, 48, 42, 24, 18, 14, 30, 7, 11, 10, 18, 50],
              "Notes": [160]}
    for ws in book.worksheets:
        for i, wdt in enumerate(widths[ws.title], start=1):
            ws.column_dimensions[ws.cell(1, i).column_letter].width = wdt
        for row in ws.iter_rows():
            for c in row:
                c.alignment = Alignment(wrap_text=True, vertical="top")
        if ws.title != "Notes":
            for c in ws[1]:
                c.font = Font(bold=True); c.fill = PatternFill("solid", fgColor="DDEBF7")
            ws.freeze_panes = "C2" if ws.title == "Treatments" else "B2"
        else:
            ws["A1"].font = Font(bold=True, size=12)
            for row in ws.iter_rows(min_row=2):
                if row[0].value in ("Records and assumptions", "Screened out", "Flags"):
                    row[0].font = Font(bold=True)
    sites = book["Sites"]
    for c in sites[sites.max_row]:
        c.font = Font(bold=True)


if __name__ == "__main__":
    main()
