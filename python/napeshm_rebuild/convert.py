"""
convert.py - rebuild NAPESHM_test.xlsx from the raw NAPESHM files.

Inputs : mgmtdata.csv (event log), treatments.csv, sites.csv,
         the existing NAPESHM_test.xlsx (template layout, User_Info, Soil_Test)
Output : NAPESHM_rebuilt.xlsx + conversion log

Every rule is labelled R1..Rn; the same labels appear in the log and the
rules document. Dates that are not in the raw records are left blank, so
the SHMI package's own imputation decides them (as in the template).
"""
import csv, re, datetime as dt, collections, json, sys, os

# usage: python convert.py <folder with mgmtdata.csv, treatments.csv, sites.csv> [converted.json]
UP  = os.path.join(sys.argv[1], "") if len(sys.argv) > 1 else ""
JSON_OUT = sys.argv[2] if len(sys.argv) > 2 else "converted.json"
RAW = UP + "mgmtdata.csv"
TRT = UP + "treatments.csv"
SIT = UP + "sites.csv"

LOG = []          # (rule, unit, date, detail)
def log(rule, unit, date, detail):
    LOG.append((rule, unit, date, detail))

def D(s):
    try:
        return dt.date.fromisoformat(s.strip())
    except Exception:
        return None

def num(s):
    try:
        return float(str(s).strip())
    except Exception:
        return None

# ---------------------------------------------------------------------------
# R6  Species names: raw planting/harvest descriptions -> template CD_name
#     and CD_group. Names follow the existing workbook's choices where those
#     were consistent; exceptions are marked in SPECIES_NOTES.
# ---------------------------------------------------------------------------
G = {  # template name -> group
 "Corn":"Grass_warm","Sorghum":"Grass_warm","Sorghum, Sudan":"Grass_warm","Bahiagrass":"Grass_warm",
 "Millet, Proso":"Grass_warm","Gamagrass":"Grass_warm",
 "Wheat":"Grass_cool","Wheat, Winter":"Grass_cool","Wheat, Spring":"Grass_cool","Barley":"Grass_cool",
 "Oats":"Grass_cool","Rye":"Grass_cool","Rye, Cereal":"Grass_cool","Rye, Winter":"Grass_cool",
 "Triticale":"Grass_cool","Ryegrass, Annual":"Grass_cool","Orchardgrass":"Grass_cool",
 "Fescue, Tall":"Grass_cool","Timothy Grass":"Grass_cool",
 "Soybean":"Legume_warm","Peanut":"Legume_warm","Cowpea":"Legume_warm","Sunnhemp":"Legume_warm",
 "Alfalfa":"Legume_cool","Clover":"Legume_cool","Clover, Red":"Legume_cool","Clover, Crimson":"Legume_cool",
 "Clover, Berseem":"Legume_cool","Vetch":"Legume_cool","Vetch, Common":"Legume_cool","Vetch, Hairy":"Legume_cool",
 "Vetch, Lana":"Legume_cool","Vetch, Purple":"Legume_cool","Pea":"Legume_cool","Pea, Field":"Legume_cool",
 "Pea, Dry":"Legume_cool","Pea, Winter":"Legume_cool","Pea, Forage":"Legume_cool","Lentils":"Legume_cool",
 "Chickpea":"Legume_cool","Beans, Dry":"Legume_cool","Beans, Bell":"Legume_cool","Beans, Snap":"Legume_cool",
 "Legume":"Legume_cool",
 "Canola":"Brassica","Canola, Winter":"Brassica","Radish":"Brassica","Rapeseed":"Brassica","Broccoli":"Brassica",
 "Turnip, Pasja":"Brassica",
 "Cotton":"Broadleaf_warm","Safflower":"Broadleaf_warm","Sunflower":"Broadleaf_warm","Buckwheat":"Broadleaf_warm",
 "Tomato":"Broadleaf_warm",
 "Potato":"Broadleaf_cool","Sugarbeet":"Broadleaf_cool","Squash, Acorn":"Broadleaf_cool","Squash, Winter":"Broadleaf_cool",
 "Poplar":"Horticulture","8-species":"Mix","Fallow":"Fallow","Hay":"Grass_warm",
}
SPECIES = {  # raw (lower-case) -> list of template names
 "acorn squash":["Squash, Acorn"], "alfalfa":["Alfalfa"], "annual rye grass":["Ryegrass, Annual"],
 "bahiagrass":["Bahiagrass"], "barley":["Barley"], "bean":["Beans, Dry"], "bell bean":["Beans, Bell"],
 "berseem clover":["Clover, Berseem"], "broccoli":["Broccoli"], "buckwheat":["Buckwheat"], "canola":["Canola"],
 "cereal rye":["Rye, Cereal"], "cereal rye (secale cereale l)":["Rye, Cereal"], "chickpea":["Chickpea"],
 "clover":["Clover"], "clover (trifolium incarnatum)":["Clover"], "common vetch":["Vetch, Common"],
 "common vetch (vicia sativa)":["Vetch, Common"], "corn":["Corn"], "cotton":["Cotton"],
 "cover crop":["8-species"], "cow pea":["Cowpea"], "cow pea and buckwheat":["Cowpea","Buckwheat"],
 "crimson clover":["Clover, Crimson"], "crotalaria":["Sunnhemp"], "field pea":["Pea, Field"],
 "forage oat (25%), winter pea (50%), lentils (25%)":["Oats","Pea, Winter","Lentils"],
 "gamagrass":["Gamagrass"], "garbanzo bean":["Chickpea"], "hairy vetch":["Vetch, Hairy"],
 "hard red spring wheat":["Wheat, Spring"], "hard red winter wheat":["Wheat, Winter"], "lana vetch":["Vetch, Lana"],
 "native warm season grasses":["8-species"],
 "nm-6 hybrid poplar tree (populus nigra x populus maximowiczii)":["Poplar"],
 "oats":["Oats"], "orchardgrass":["Orchardgrass"], "pea":["Pea"], "pea, dry edible (pisum sativum)":["Pea, Dry"],
 "peanut":["Peanut"], "potato":["Potato"], "purple vetch":["Vetch, Purple"], "raddish":["Radish"],
 "radish":["Radish"], "radish (raphanus sativus)":["Radish"],
 "radish (20%) buckwheat (40%) forage sorghum (20%)":["Radish","Buckwheat","Sorghum"],
 "rapeseed":["Rapeseed"], "red clover":["Clover, Red"], "rye":["Rye"], "rye + winter wheat":["Rye","Wheat, Winter"],
 "safflower":["Safflower"], "snap beans":["Beans, Snap"], "sorghum":["Sorghum"],
 "sorghum sudangrass":["Sorghum, Sudan"], "sorghum sudangrass (cover crop)":["Sorghum, Sudan"],
 "soybean":["Soybean"], "spring wheat":["Wheat, Spring"], "sugar beet":["Sugarbeet"], "sunflower":["Sunflower"],
 "tall fescue":["Fescue, Tall"], "timothy":["Timothy Grass"], "tomato":["Tomato"], "triticale":["Triticale"],
 "triticale  + common vetch":["Triticale","Vetch, Common"], "triticale + common vetch":["Triticale","Vetch, Common"],
 "triticale (triticosecale wittm.)":["Triticale"], "triticale or oats":["Triticale"], "vetch":["Vetch"],
 "warm season legume (soybean or chickpea)":["Soybean"], "wheat":["Wheat"], "winter legume cover crop":["Legume"],
 "winter pea":["Pea, Winter"], "winter rye":["Rye, Winter"], "winter squash":["Squash, Winter"],
 "winter wheat":["Wheat, Winter"], "winter wheat (cover crop)":["Wheat, Winter"],
 "cover crop mix (tritcale, proso millet, sunflower, winter canola, soybean, forage pea, pasja turnip)":
     ["Triticale","Millet, Proso","Sunflower","Canola, Winter","Soybean","Pea, Forage","Turnip, Pasja"],
 "cover crop mix":["Triticale","Millet, Proso","Sunflower","Canola, Winter","Soybean","Pea, Forage"],
 # harvest-only names
 "sugar beet (roots)":["Sugarbeet"], "native grasses or crabgrass":["8-species"],
 "prarie grass":["8-species"], "prairie grass":["8-species"],
}
SPECIES_NOTES = {
 "cover crop": "generic 'cover crop' coded as '8-species' mix, following the existing workbook",
 "native warm season grasses": "coded as '8-species' mix, following the existing workbook",
 "native grasses or crabgrass": "coded as '8-species' mix (harvest-only record)",
 "cover crop mix": "generic mix; species list taken from the existing workbook",
 "triticale or oats": "'triticale or oats' coded as Triticale, following the existing workbook",
 "warm season legume (soybean or chickpea)": "'soybean or chickpea' coded as Soybean",
 "winter legume cover crop": "coded as generic 'Legume', following the existing workbook",
 "berseem clover": "coded Clover, Berseem (existing workbook used Clover, Balansa)",
 "gamagrass": "not in the template list; coded 'Gamagrass' (Grass_warm)",
}
NOT_SPECIES = re.compile(r"planter|drill|opener|seeder|\bdisc\b|\bdisk\b|fertilizer|npk|phosphate|"
                         r"ammonium|urea|potash|nitrate|insecticide|herbicide|fungicide|row planter", re.I)

# Species kept as perennial stands when they persist beyond one season (R5)
PERENNIAL_SPECIES = {"Alfalfa","Orchardgrass","Bahiagrass","Fescue, Tall","Timothy Grass","Clover, Red",
                     "Gamagrass","8-species"}
WOODY = {"Poplar"}
REGROWTH = {"Sorghum, Sudan"}        # annual forages that are cut repeatedly

def species_of(desc):
    k = re.sub(r"\s+", " ", desc.strip().lower())
    if k in SPECIES: return SPECIES[k], SPECIES_NOTES.get(k)
    k2 = desc.strip().lower()
    if k2 in SPECIES: return SPECIES[k2], SPECIES_NOTES.get(k2)
    return None, None

def harvest_species(desc):
    """Species named in a harvest description, or [] if generic."""
    s = desc.strip().lower()
    s = re.sub(r"\(.*?\)", "", s)                     # drop '(grain)', '(biomass)', ...
    s = re.sub(r"\b(grain|biomass|seed|ears|straw|tubers|roots|cobs|fruit|silage)\b", "", s)
    s2 = re.sub(r"\bhay\b", "", s).strip()
    if s2 and s2 != s.strip():
        s = s2
    parts = re.split(r"\s*\+\s*|\s+or\s+|,\s*", s)
    out = []
    for p in parts:
        p = p.strip()
        if not p: continue
        sp, _ = species_of(p)
        if sp: out += sp
    return out

def base(name):  # 'Wheat, Winter' -> 'wheat'; used to match harvests to plantings
    return name.split(",")[0].strip().lower()

RESIDUE = re.compile(r"stover|straw|residue|\brake\b|chopp|shred", re.I)
BURNDOWN = re.compile(r"glyphosate|paraquat|glufosinate|gramoxone|roundup", re.I)
TERM_CATS = {"cover crop termination","termination","pre-plant burndown","burndown","burning","controlled burn"}
PLANT_CATS = {"plant cash crop","plant cover crop","planting"}
ORGANIC = re.compile(r"manure|compost|litter|biosolid|biochar|digestate|slurry|sludge|bedding|guano|"
                     r"vermicompost|bokashi|blood meal|bone meal|feather meal|alfalfa meal|fish emulsion", re.I)

# ---------------------------------------------------------------------------
# Read raw files
# ---------------------------------------------------------------------------
raw = list(csv.DictReader(open(RAW, encoding="utf-8-sig")))
for r in raw:
    r["key"] = f"NAPESHM_{r['sitecode'].strip()}_{r['treatmentid'].strip()}"
    r["cat"] = (r["category"] or "").strip().lower()
    r["unit"] = (r["unit"] or "").strip()
    r["d"] = D(r["date"] or "")
    r["q"] = num(r["quantity"])
    r["desc"] = (r["eventdesc"] or "").strip()
treatments = list(csv.DictReader(open(TRT, encoding="utf-8-sig")))
sites = {s["sitecode"]: s for s in csv.DictReader(open(SIT, encoding="utf-8-sig"))}

by_unit = collections.defaultdict(list)
for r in raw:
    if r["d"] is None:
        log("R0", r["key"], r["date"], f"record without a valid date skipped: {r['cat']} / {r['desc']}")
        continue
    by_unit[r["key"]].append(r)

# ---------------------------------------------------------------------------
# R1  Mgt_Unit: one row per treatment (treatments.csv + sites.csv)
# ---------------------------------------------------------------------------
mgt_rows = []
trt_info = {}
REMOVED_UNITS = []
REFERENCE_UNITS = []
# Vegetation of reference soils whose single species is not named elsewhere
REFERENCE_SPECIES = {"USNY": ("Alfalfa", "Legume_cool")}   # NY benchmarks: alfalfa (data owner)
for t in treatments:
    sc = t["sitecode"].strip(); key = f"NAPESHM_{sc}_{t['treatmentid'].strip()}"
    s = sites.get(sc, {})
    trt_info[key] = t
    if key not in by_unit:
        # R2e: treatments without management records are reference soils
        # (sod, remnant prairie, exclosures, benchmarks); kept, with a stand
        # row added below
        REFERENCE_UNITS.append(key)
        log("R2e", key, None, "no records in mgmtdata.csv: reference soil, entered by the reference-soil convention")
    mgt_rows.append(dict(MGT_combo=key, MGT_study="NAPESHM", MGT_farm=None, MGT_field=sc,
                         MGT_lat=num(s.get("site_latitude")), MGT_lon=num(s.get("site_longitude")),
                         MGT_state=s.get("state_prov"), MGT_county=None, MGT_polygon=None,
                         MGT_trt=int(t["treatmentid"]), MGT_tile_drain="Unsure", MGT_tile_years="Unsure",
                         MGT_notes=s.get("studyname")))
    if not s: log("R1", key, None, f"site {sc} not in sites.csv; lat/lon/state blank")
for k in by_unit:
    if k not in trt_info: log("R1", k, None, "unit has management records but no treatment row")

# ---------------------------------------------------------------------------
# R2-R5  Crop_Diversity
# ---------------------------------------------------------------------------
def build_crops(key, events):
    evs = sorted(events, key=lambda r: (r["d"], 0 if r["cat"] == "harvest" else
                                         1 if r["cat"] in TERM_CATS else 2))
    start = next((r["d"] for r in evs if r["cat"] == "reserved index row" and r["desc"].lower() == "start"), None)
    end   = next((r["d"] for r in evs if r["cat"] == "reserved index row" and r["desc"].lower() == "end"), None)
    if end is None and evs:
        end = dt.date(max(r["d"] for r in evs).year, 12, 31)
        log("R2c", key, end, "record 'end' row has no date; record end set to 31 Dec of the last event year")
    if start is None and evs:
        start = dt.date(min(r["d"] for r in evs).year, 1, 1)
        log("R2c", key, start, "record 'start' row missing; record start set to 1 Jan of the first event year")
    orphans = []      # harvests that could not be assigned (used by perennial stands, R2c)
    episodes, active = [], []
    first_plant = next((r["d"] for r in evs if r["cat"] in PLANT_CATS and r["unit"] != "STIR"
                        and species_of(r["desc"])[0]), None)

    def new_ep(species, cat, plant, note=None, src=""):
        ep = dict(species=list(dict.fromkeys(species)), cat=cat, plant=plant, harvests=[], term=None,
                  note=note, src=src)
        episodes.append(ep); active.append(ep); return ep

    # planting events grouped by (date, category) -> one episode (R2, R3 mixes)
    plant_groups = collections.OrderedDict()
    for r in evs:
        if r["cat"] in PLANT_CATS and r["unit"] != "STIR":
            sp, note = species_of(r["desc"])
            if sp is None:
                if NOT_SPECIES.search(r["desc"]):
                    log("R2", key, r["d"], f"planting record is equipment/fertilizer, not a species: '{r['desc']}'")
                else:
                    log("R6", key, r["d"], f"unrecognized species '{r['desc']}' - planting skipped")
                continue
            if note: log("R6", key, r["d"], f"'{r['desc']}': {note}")
            plant_groups.setdefault((r["d"], r["cat"]), []).extend(sp)
    plant_done = set()
    species_plant_dates = {pd for (pd, pc) in plant_groups}
    passes = sorted({r["d"] for r in evs if r["cat"] in PLANT_CATS and r["unit"] == "STIR"
                     and r["d"] not in species_plant_dates})
    last_h = [start or dt.date.min]

    for r in evs:
        d = r["d"]
        # ---- harvest (R4)
        cut_sp = (harvest_species(r["desc"]) or harvest_species(r["method"])) if r["cat"] == "mowing" else []
        is_cut = bool(cut_sp)
        if r["cat"] == "harvest" or is_cut:
            sp = cut_sp if is_cut else harvest_species(r["desc"])
            if is_cut:
                log("R4c", key, d, f"mowing '{r['desc']}' of a named crop treated as a harvest (hay cut)")
            pres = r["q"] if ("resid" in r["unit"] or "removed" in r["unit"]) else None
            targets = [e for e in active if sp and any(base(s) == base(x) for s in e["species"] for x in sp)]
            if len(targets) > 1:                      # same species planted more than once: newest stand
                targets = [max(targets, key=lambda e: e["plant"] or dt.date.min)]
            if not sp and RESIDUE.search(r["desc"]):
                log("R4", key, d, f"residue operation '{r['desc']}' recorded as harvest - not a crop harvest, ignored")
                continue
            if not targets and not sp:
                cand = [e for e in active if e["cat"] != "cover" and not e["harvests"]]
                if cand:
                    targets = [cand[-1]]
                    log("R4", key, d, f"generic harvest '{r['desc']}' assigned to unharvested crop {cand[-1]['species']}")
                elif any(e["harvests"] for e in active):
                    log("R4", key, d, f"generic harvest/residue operation '{r['desc']}' after the crop was harvested - ignored")
                    continue
            if not targets:
                if sp:
                    cand = [p for p in passes if last_h[0] < p < d]
                    if cand:
                        ep = new_ep(sp, "cash", cand[-1], "planting date from planting pass without species", "harvest")
                        passes.remove(cand[-1])
                        log("R3c", key, d, f"harvest of {sp} without a species planting: planting date {cand[-1]} "
                                           f"taken from the last unassigned planting pass")
                    else:
                        ep = new_ep(sp, "cash", None, "planted before the record starts", "harvest")
                        log("R3", key, d, f"harvest of {sp} with no planting in the record: harvest-only episode")
                    targets = [ep]
                else:
                    orphans.append((d, pres, r["desc"]))
                    continue
            for e in targets:
                if e["harvests"]:
                    stand = any(x in PERENNIAL_SPECIES or x in REGROWTH for x in e["species"])
                    gap = (d - e["harvests"][-1][0]).days
                    if not stand and gap <= 60:
                        log("R4", key, d, f"second harvest of annual {e['species']} {gap} days after the first - ignored")
                        continue
                e["harvests"].append((d, pres))
            last_h[0] = max(last_h[0], d)
            continue
        # ---- termination (R4)
        if r["cat"] in TERM_CATS:
            if any(e["term"] == d for e in episodes):
                continue                     # several products applied the same day = one termination
            if r["cat"] == "cover crop termination":
                tg = [e for e in active if e["cat"] == "cover" and e["term"] is None]
                if not tg:
                    tg = [e for e in active if not e["harvests"] and e["term"] is None]
            else:
                tg = [e for e in active if (e["cat"] != "cash" or any(x in PERENNIAL_SPECIES for x in e["species"]))
                      and e["term"] is None]
                tg += [e for e in active if e["cat"] == "cash" and not e["harvests"] and e["term"] is None
                       and r["cat"] == "termination"]
            if not tg and r["cat"] == "cover crop termination" and (first_plant is None or d < first_plant):
                # R3b: termination before any planting in the record -> infer the crop from the
                # next cover-crop planting in this unit
                nxt = next((sp for (pd, pc), sp in plant_groups.items() if pc == "plant cover crop" and pd > d), None)
                if nxt:
                    ep = new_ep(nxt, "cover", None, "species inferred from next cover planting", "inferred")
                    tg = [ep]
                    log("R3b", key, d, f"termination before any planting; crop inferred as {nxt} (next cover planting)")
            if not tg:
                log("R4", key, d, f"termination '{r['cat']}' with no active crop to end - ignored")
                continue
            for e in tg:
                e["term"] = d
            active[:] = [e for e in active if e["term"] is None]
            continue
        # ---- non-selective herbicide on an active cover crop (R4b)
        if r["cat"] == "weed management" and BURNDOWN.search(r["desc"]):
            live_cash = [e for e in active if e["cat"] == "cash" and not e["harvests"] and e["plant"]]
            tg = [e for e in active if e["cat"] == "cover" and e["term"] is None]
            if tg and not live_cash and not any(e["term"] == d for e in episodes):
                for e in tg:
                    e["term"] = d
                active[:] = [e for e in active if e["term"] is None]
                log("R4b", key, d, f"cover crop {tg[0]['species']} ended by non-selective herbicide '{r['desc']}' (weed management)")
            continue
        # ---- planting (R2)
        if r["cat"] in PLANT_CATS and (d, r["cat"]) in plant_groups and (d, r["cat"]) not in plant_done:
            plant_done.add((d, r["cat"]))
            sp = plant_groups[(d, r["cat"])]
            cat = "cover" if r["cat"] == "plant cover crop" else "cash"
            growing = [e for e in active if e["cat"] == "cash" and e["plant"] and e["plant"] < d
                       and not e["harvests"] and e["term"] is None
                       and not set(e["species"]) & set(sp)]
            relay = bool(growing)
            if relay and cat == "cash":
                cat = "cover"
                log("R2d", key, d, f"{sp} planted into growing {growing[-1]['species']}: coded as a relay/interseeded "
                                   f"cover crop (raw category '{r['cat']}'); it is not ended by the companion's harvest")
            # a new cash planting ends harvested crops that are still listed as active
            if cat == "cash":
                active[:] = [e for e in active if not e["harvests"]]
            ep = new_ep(sp, cat, d, src="planting")
            ep["relay"] = relay

    # R2c: perennial / non-farmed / woody treatments without plantings: one stand spanning the record
    t = trt_info.get(key, {})
    ccat = (t.get("eu_crop_cat") or "").strip().lower()
    if ccat in ("perennial", "non_farmed", "woody_perennial") and not any(e["src"] == "planting" for e in episodes):
        hsp = [e for e in episodes if e["src"] in ("harvest", "inferred")]
        if hsp:
            for e in hsp:
                if e["plant"] is None and start:
                    e["plant"] = start
                if e["term"] is None and end:
                    e["term"] = end
                e["cat"] = "perennial"
            log("R2c", key, None, f"perennial treatment without plantings: harvested stand {hsp[0]['species']} set to span the record {start}..{end}")
        elif start and end:
            n = int(t.get("plant_species_ct") or 0)
            woody = re.search(r"\(([^)]+)\)", t.get("eu_treat_desc", ""))
            if ccat == "woody_perennial":
                name = woody.group(1).strip().capitalize() if woody else "Forest"
                sp, cat = [name], "woody perennial"
            elif n >= 2:
                sp, cat = [f"{n}-species"], "perennial"
            elif any(o[2].strip().lower() == "hay" for o in orphans):
                sp, cat = ["Hay"], "perennial"
            else:
                sp, cat = ["BORROW"], "perennial"
            episodes.append(dict(species=sp, cat=cat, plant=start, harvests=sorted((o[0], o[1]) for o in orphans),
                                 term=end, note="stand spanning the record (no crop events)", src="stand"))
            if orphans:
                log("R2c", key, None, f"{len(orphans)} unassigned harvest(s) attached to the stand")
            orphans = []
            log("R2c", key, None, f"perennial treatment without crop events: {sp[0]} stand {start}..{end} "
                                  f"(species count from treatments.csv plant_species_ct = {n})")

    for o in orphans:
        log("R4", key, o[0], f"harvest '{o[2]}' with no crop to assign - ignored")

    # R2b: units with no crop at all in a cropped (non-perennial) treatment = bare fallow
    if not episodes:
        t = trt_info.get(key, {})
        desc = " ".join([t.get("eu_treat_desc",""), t.get("eu_crop_cat",""), t.get("eu_other_trt","")]).lower()
        if ("fallow" in desc or ccat == "row_crop") and start and end:
            episodes.append(dict(species=["Fallow"], cat="fallow", plant=start, harvests=[(end, 0.0)],
                                 term=None, note="fallow treatment", src="fallow"))
            log("R2b", key, None, f"no crops recorded; fallow treatment -> Fallow row {start}..{end}")
        else:
            log("R2b", key, None, "no crop records in the event log")

    # R5b: repeated plantings of the same perennial species with nothing else in between = one stand
    episodes.sort(key=lambda e: e["plant"] or (e["harvests"][0][0] if e["harvests"] else e["term"]) or dt.date.max)
    merged = []
    for e in episodes:
        prev = merged[-1] if merged else None
        if (prev and set(e["species"]) == set(prev["species"]) and prev["term"] is None
                and all(x in PERENNIAL_SPECIES for x in e["species"]) and e["src"] == "planting"):
            prev["harvests"] += e["harvests"]
            prev["term"] = e["term"]
            log("R5b", key, e["plant"], f"repeated planting of {e['species']} merged into the existing stand (overseeding)")
            continue
        merged.append(e)
    episodes[:] = merged

    # R5 categories: perennial if a perennial species persists beyond one season
    rows = []
    episodes.sort(key=lambda e: e["plant"] or (e["harvests"][0][0] if e["harvests"] else e["term"]) or dt.date.max)
    for seq, e in enumerate(episodes, 1):
        dates = [x for x in [e["plant"], e["term"]] + [h[0] for h in e["harvests"]] if x]
        span = (max(dates) - min(dates)).days if dates else 0
        years = {h[0].year for h in e["harvests"]}
        cat = e["cat"]
        if e.get("relay") and e["harvests"] and cat == "cover":
            cat = "cash"
            log("R2d", key, e["plant"], f"relay crop {e['species']} was harvested: coded cash")
        if any(s in WOODY for s in e["species"]):
            cat = "woody perennial"
        elif cat in ("cash", "cover") and any(s in PERENNIAL_SPECIES for s in e["species"]) \
                and (span > 365 or len(years) > 1 or (e["plant"] is None and len(e["harvests"]) > 1)):
            cat = "perennial"
            log("R5", key, e["plant"], f"{e['species']} coded perennial (stand spans {span} days, harvests in {len(years)} year(s))")
        mix = "yes" if len(e["species"]) > 1 else "no"
        hs = e["harvests"]
        for s in e["species"]:
            first_h = hs[0] if hs else (None, None)
            term_main = e["term"] if len(hs) <= 1 else None
            grp = "Mix" if re.match(r"\d+-species$", s) else ("Horticulture" if cat == "woody perennial" else G.get(s, "Unknown"))
            rows.append(dict(MGT_combo=key, CD_seq_num=seq, CD_mix=mix, CD_cat=cat, CD_group=grp,
                             CD_name=s, CD_plant_date=e["plant"], CD_harv_date=first_h[0], CD_term_date=term_main,
                             CD_per_res=first_h[1], CD_notes=e["note"]))
            for i, (hd, pr) in enumerate(hs[1:], 1):
                last = i == len(hs) - 1
                rows.append(dict(MGT_combo=key, CD_seq_num=seq, CD_mix=mix, CD_cat=cat, CD_group=grp,
                                 CD_name=s, CD_plant_date=None, CD_harv_date=hd,
                                 CD_term_date=e["term"] if last else None, CD_per_res=pr, CD_notes=None))
    return rows

crop_rows = []
for key in sorted(by_unit):
    crop_rows += build_crops(key, by_unit[key])

# R2c (cont.): stands with unknown species borrow the perennial species recorded
# at other treatments of the same site
for row in crop_rows:
    if row["CD_name"] == "BORROW":
        site = row["MGT_combo"].rsplit("_", 1)[0]
        donors = collections.Counter(r["CD_name"] for r in crop_rows
                                     if r["MGT_combo"].startswith(site + "_") and r["CD_name"] != "BORROW"
                                     and r["CD_cat"] == "perennial")
        name = donors.most_common(1)[0][0] if donors else "Unknown"
        row["CD_name"], row["CD_group"] = name, G.get(name, "Unknown")
        log("R2c", row["MGT_combo"], None, f"stand species unknown; borrowed '{name}' from other treatments at the site")

# R2e (cont.): reference soils without records get one perennial stand row
# spanning the calendar years recorded at the same site; no disturbance,
# amendments or animals. Species: "N-species" from plant_species_ct, or the
# named species for single-species references (REFERENCE_SPECIES).
site_years = collections.defaultdict(list)
for key, evs in by_unit.items():
    site_years[key.rsplit("_", 1)[0]] += [r["d"].year for r in evs]
for key in REFERENCE_UNITS:
    site = key.rsplit("_", 1)[0]
    yrs = site_years.get(site)
    if not yrs:
        log("R2e", key, None, "no records at the site either: no evaluation window, reference soil left without crops")
        continue
    t = trt_info[key]
    n = int(t.get("plant_species_ct") or 0)
    sc = site.replace("NAPESHM_", "")
    named = next((v for k, v in REFERENCE_SPECIES.items() if sc.startswith(k)), None)
    if n <= 1 and named:
        name, grp = named
    elif n >= 2:
        name, grp = f"{n}-species", "Mix"
    else:
        name, grp = "Unknown", "Unknown"
        log("R2e", key, None, "single-species reference soil with no named species: 'Unknown'")
    start = dt.date(min(yrs), 1, 1); end = dt.date(max(yrs), 12, 31)
    crop_rows.append(dict(MGT_combo=key, CD_seq_num=1, CD_mix="no", CD_cat="perennial", CD_group=grp,
                          CD_name=name, CD_plant_date=start, CD_harv_date=None, CD_term_date=end,
                          CD_per_res=None, CD_notes="reference soil (no management records): stand spanning the site's record years"))
    log("R2e", key, None, f"reference stand '{name}' {start}..{end} (site record years)")

# ---------------------------------------------------------------------------
# R7  Soil_Disturbance: every record measured in STIR, whatever its category
# ---------------------------------------------------------------------------
dist_rows = []
for key in sorted(by_unit):
    for r in sorted(by_unit[key], key=lambda r: r["d"]):
        if r["unit"] == "STIR":
            if r["q"] is None:
                log("R7", key, r["d"], f"STIR record without a value skipped: '{r['desc']}'"); continue
            phase = r["cat"] if r["cat"] and r["cat"] != "disturbnace" else "disturbance"
            if r["cat"] in ("", "disturbnace"):
                log("R7", key, r["d"], f"category '{r['cat']}' recorded as 'disturbance'")
            dist_rows.append(dict(MGT_combo=key, SD_phase=phase, SD_date=r["d"], SD_equip=r["desc"], SD_mixeff=r["q"]))

# ---------------------------------------------------------------------------
# R8  Soil_Amendments: organic amendments only
# ---------------------------------------------------------------------------
amend_rows = []
for key in sorted(by_unit):
    for r in sorted(by_unit[key], key=lambda r: r["d"]):
        is_amend = (r["cat"] == "amendment" or r["cat"].startswith("composted")
                    or (r["cat"] == "" and r["unit"] != "STIR"))
        if is_amend and ORGANIC.search(r["desc"] + " " + r["cat"]):
            src = r["desc"] if r["cat"] in ("amendment", "") else r["category"].strip()
            if src.lower().startswith("solid manure spreader"):
                src = "solid manure (material not named; raw: '" + r["desc"] + "')"
                log("R8", key, r["d"], "amendment named only by the spreader; recorded as solid manure (confirmed by data owner)")
            amend_rows.append(dict(MGT_combo=key, SA_cat="Organic", SA_source=src, SA_date=r["d"],
                                   SA_rate=r["q"], SA_units=r["unit"] or None))
            if r["cat"] == "":
                log("R8", key, r["d"], f"organic amendment with a blank category: '{r['desc']}'")
            elif r["cat"] != "amendment":
                log("R8", key, r["d"], f"amendment recorded in the category field: '{r['category'].strip()}'")

# ---------------------------------------------------------------------------
# R9  Animal_Diversity: grazing 'begins'/'ends' pairs
# ---------------------------------------------------------------------------
def animal_type(r):
    t = (r["desc"] + " " + r["method"] + " " + r["unit"]).lower()
    for k in ("sheep", "steers", "cow/calf", "cattle", "cow"):
        if k in t: return k
    return "animal"

animal_rows = []
for key in sorted(by_unit):
    evs = sorted(by_unit[key], key=lambda r: r["d"])
    gz = [r for r in evs if r["cat"] == "grazing"]
    begins_dates = {r["d"] for r in gz if "begin" in (r["desc"] + " " + r["method"]).lower()}
    open_ = None
    used = set()
    for r in gz:
        txt = (r["desc"] + " " + r["method"]).lower()
        herd = re.search(r"^herd (\d+)", r["desc"].lower())
        acres = re.search(r"acres used for herd (\d+)", r["desc"].lower())
        is_begin = "begin" in txt or herd or ("stocking" in txt and r["d"] not in begins_dates)
        is_end = re.search(r"\bends?\b|end grazing", txt) or acres
        if is_begin:
            if open_:
                animal_rows.append(open_); log("R9", key, open_["AD_start_date"], "grazing start without an end")
            open_ = dict(MGT_combo=key, AD_type=animal_type(r), AD_intensity=None, AD_count=r["q"],
                         AD_start_date=r["d"], AD_end_date=None, AD_notes=r["unit"] or None)
            if "stocking" in txt and "begin" not in txt:
                log("R9", key, r["d"], "stocking-rate record used as the start of a grazing period")
            if herd:
                open_["_herd_count"] = r["q"]
            used.add(id(r))
        elif is_end:
            used.add(id(r))
            if open_:
                open_["AD_end_date"] = r["d"]
                if acres and open_.get("_herd_count") and r["q"]:
                    open_["AD_count"] = open_["_herd_count"] / r["q"]; open_["AD_notes"] = "cattle ac-1"
                    log("R9", key, open_["AD_start_date"], "herd record = start, 'acres used for herd' = end; count = head / acres")
                open_.pop("_herd_count", None)
                animal_rows.append(open_); open_ = None
            else:
                animal_rows.append(dict(MGT_combo=key, AD_type=animal_type(r), AD_intensity=None, AD_count=None,
                                        AD_start_date=None, AD_end_date=r["d"], AD_notes=None))
                log("R9", key, r["d"], "grazing end without a start")
    if open_:
        open_.pop("_herd_count", None)
        animal_rows.append(open_); log("R9", key, open_["AD_start_date"], "grazing start without an end")
    others = [r for r in gz if id(r) not in used]
    if others:
        log("R9", key, None, f"{len(others)} grazing records not used as period boundaries (e.g. '{others[0]['desc']}')")
    # R9b: grazing recorded only as animal trampling in Soil_Disturbance (method 'cattle')
    for r in evs:
        periods = [(a["AD_start_date"], a["AD_end_date"]) for a in animal_rows if a["MGT_combo"] == key]
        inside = any(s0 and e0 and s0 <= r["d"] <= e0 for s0, e0 in periods)
        if r["cat"] == "disturbance" and "cattle" in r["method"].lower() and not inside:
            animal_rows.append(dict(MGT_combo=key, AD_type="cattle", AD_intensity=None, AD_count=None,
                                    AD_start_date=r["d"], AD_end_date=r["d"] + dt.timedelta(days=7), AD_notes=None))
            log("R9b", key, r["d"], "cattle trampling (disturbance record) = grazing period of 7 days")

json.dump(dict(removed=REMOVED_UNITS, reference=REFERENCE_UNITS, mgt=mgt_rows, crop=crop_rows, dist=dist_rows, amend=amend_rows, animal=animal_rows,
               log=LOG), open(JSON_OUT, "w"), default=str)
print({k: len(v) for k, v in dict(mgt=mgt_rows, crop=crop_rows, dist=dist_rows, amend=amend_rows,
                                   animal=animal_rows, log=LOG).items()})
print(collections.Counter(l[0] for l in LOG))
