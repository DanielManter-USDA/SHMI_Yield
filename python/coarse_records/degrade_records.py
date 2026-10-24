#!/usr/bin/env python3
"""
degrade_records.py: reduce an SHMI workbook to coarser, remote-sensing-like records.

Writes a copy of an SHMI workbook in which the management records keep only the
detail that public or remote-sensing products provide, so SHMI can be rescored
with the published weights and compared with the full-record score
(R/SHMI_coarse_records.R runs the scenarios).

Operations (--ops, comma-separated; always applied in this order):
  dates      every crop's planting, harvest and termination dates set to the
             site's typical day of year for that crop and date (median over the
             site's records; the year is kept): programs know what was grown,
             not when.
  cover      cover crops reduced to present/absent: all cover rows between two
             cash crops become one generic cover crop (the workbook's most common
             cover species) spanning the group's dates (crops without an end date
             end at the next planting), within the gap between the cash crops.
  maincrop   one crop at a time: a crop that overlaps an earlier-listed crop of
             the same unit for more than half its duration is dropped
             (mixtures, companion crops, intercrops); covers are kept.
  tillage    each crop interval keeps only its tillage class (no tillage,
             no-till, reduced, conventional: T-DISC designations): its passes are
             replaced by one pass at the typical intensity of that class in
             this workbook (median of the intervals in the class).
  noamend    organic amendments (SA_cat "Organic") and grazing removed (amendments
             become zero fertilizer rows); fertilizer records are kept.
  last:N     only the last N calendar years of each unit's record are kept.
  jitter:N   every crop and tillage date shifted by a random 0-N days, either
             way, within its calendar year (a row keeps its original dates if the
             shift breaks its order).
  sjitter:N  as jitter, but the shift is shared: every crop and tillage date that
             falls on the same day at a site moves by the same amount (replicate
             plots, treatments planted the same day, and a pass and the planting it
             goes with stay together), as when a field's dates are known only to
             within N days.
  misclass:P with cover and/or tillage: each cover-crop gap is flipped
             (present <-> absent) and each tilled interval's class moved one
             step (no-till <-> reduced <-> conventional) with probability P.

Tillage intensity follows SHMI's compute_disturbance(): crop intervals from one
cash harvest to the next, five windows around planting, implements mixed within
a window from shallow to deep (EPA, depth capped at 30 cm, TI = mixed depth / 30)
or, with --tillage stir, window STIR / 135; the interval takes its most intense
window. NAPESHM records carry STIR (the calibration scores them with STIR);
the long-term experiments carry implements (EPA).

Usage:
  python degrade_records.py --in SHMI_X.xlsx --out SHMI_X_D3.xlsx --ops dates,cover,tillage \
      [--tillage epa|stir] [--seed 1] [--site-regex '^(NAPESHM_[^_]+)_'] [--log log.csv]
"""
import argparse, csv, datetime as dt, os, random, re, statistics, sys
from collections import Counter, defaultdict

import warnings
import openpyxl
warnings.filterwarnings("ignore", message="Data Validation extension")

HERE = os.path.dirname(os.path.abspath(__file__))
STIR_MAX = 135.0
HDR_ROW = 4          # row 4 holds the column names SHMI reads; data start in row 5
COVER_GAP_DAYS = 14

# ----------------------------------------------------------------------------- workbook I/O
class Sheet:
    def __init__(self, ws):
        self.ws = ws
        self.cols = {}
        for c in range(1, ws.max_column + 1):
            v = ws.cell(HDR_ROW, c).value
            if isinstance(v, str) and v.strip():
                self.cols[v.strip()] = c
        self.ncol = ws.max_column
        key = self.cols.get("MGT_combo")
        self.rows = []
        last = HDR_ROW
        for r in range(HDR_ROW + 1, ws.max_row + 1):
            vals = [ws.cell(r, c).value for c in range(1, self.ncol + 1)]
            if key and vals[key - 1] not in (None, ""):
                self.rows.append({"_row": r, "_vals": vals})
                last = r
        self.last = last

    def get(self, row, name):
        c = self.cols.get(name)
        return row["_vals"][c - 1] if c else None

    def set(self, row, name, value):
        c = self.cols.get(name)
        if c: row["_vals"][c - 1] = value

    def new_row(self, template):
        return {"_row": template["_row"], "_vals": list(template["_vals"])}

    def write(self):
        ws = self.ws
        end = max(self.last, HDR_ROW + len(self.rows))
        for r in range(HDR_ROW + 1, end + 1):
            for c in range(1, self.ncol + 1):
                ws.cell(r, c).value = None
        for i, row in enumerate(self.rows):
            r = HDR_ROW + 1 + i
            for c, v in enumerate(row["_vals"], start=1):
                if isinstance(v, str) and v.startswith("="):      # formulas: move row references
                    v = re.sub(r"(\$?[A-Z]{1,3}\$?)%d\b" % row["_row"], lambda m: m.group(1) + str(r), v)
                ws.cell(r, c).value = v


def as_date(v):
    if isinstance(v, dt.datetime): return v.date()
    if isinstance(v, dt.date): return v
    if isinstance(v, str):
        for fmt in ("%m/%d/%Y", "%Y-%m-%d"):
            try: return dt.datetime.strptime(v.strip(), fmt).date()
            except ValueError: pass
    return None


def to_cell(d):
    return dt.datetime(d.year, d.month, d.day) if d else None


# ----------------------------------------------------------------------------- crop helpers
def cat(sh, row):
    return str(sh.get(row, "CD_cat") or "").strip().lower()

def is_cash(sh, row): return cat(sh, row) in ("cash", "annual")
def is_cover(sh, row): return cat(sh, row) == "cover"

def crop_dates(sh, row):
    return (as_date(sh.get(row, "CD_plant_date")), as_date(sh.get(row, "CD_harv_date")),
            as_date(sh.get(row, "CD_term_date")))

def crop_span(sh, row):
    p, h, t = crop_dates(sh, row)
    end = max([d for d in (h, t) if d], default=None)
    return p, end

def spans(sh, unit_rows):
    """(start, end) per row; a crop without an end date ends the day before the unit's
    next planting, or on 31 December (as prepare_shmi_inputs() closes open crops)."""
    plants = sorted({crop_dates(sh, r)[0] for r in unit_rows if crop_dates(sh, r)[0]})
    out = {}
    for r in unit_rows:
        s, e = crop_span(sh, r)
        if s and not e:
            nxt = [p for p in plants if p > s]
            e = nxt[0] - dt.timedelta(days=1) if nxt else dt.date(s.year, 12, 31)
        out[id(r)] = (s, e)
    return out

def cash_list(sh, unit_rows):
    """Cash crops of a unit as (plant, harvest) sorted by harvest."""
    out = []
    for r in unit_rows:
        if is_cash(sh, r):
            p, h, t = crop_dates(sh, r)
            if h: out.append((p, h))
    out.sort(key=lambda x: x[1])
    return out

def by_unit(sh):
    g = defaultdict(list)
    for r in sh.rows: g[sh.get(r, "MGT_combo")].append(r)
    return g


# ----------------------------------------------------------------------------- operations
def op_dates(sh, site_of, log):
    doy = defaultdict(list)
    for r in sh.rows:
        site = site_of(sh.get(r, "MGT_combo")); name = sh.get(r, "CD_name")
        for f in ("CD_plant_date", "CD_harv_date", "CD_term_date"):
            d = as_date(sh.get(r, f))
            if d: doy[(site, name, f)].append(d.timetuple().tm_yday)
    med = {k: int(round(statistics.median(v))) for k, v in doy.items()}
    changed = reverted = 0
    for r in sh.rows:
        site = site_of(sh.get(r, "MGT_combo")); name = sh.get(r, "CD_name")
        new = {}
        for f in ("CD_plant_date", "CD_harv_date", "CD_term_date"):
            d = as_date(sh.get(r, f))
            if d:
                y0 = dt.date(d.year, 1, 1)
                new[f] = y0 + dt.timedelta(days=min(med[(site, name, f)], 366 if d.year % 4 == 0 else 365) - 1)
        h0, t0 = as_date(sh.get(r, "CD_harv_date")), as_date(sh.get(r, "CD_term_date"))
        if h0 and t0 and h0 == t0: new["CD_term_date"] = new["CD_harv_date"]   # term copied from harvest (rule T1)
        p = new.get("CD_plant_date"); e = [new.get("CD_harv_date"), new.get("CD_term_date")]
        if p and any(x and x <= p for x in e):
            reverted += 1; continue
        for f, d in new.items(): sh.set(r, f, to_cell(d))
        changed += 1
    log.append(("dates", f"{changed} crop rows set to site-typical dates; {reverted} kept (order would break)"))


def gaps_of(cash):
    """Gaps between consecutive cash crops: list of (prev_harvest or None, next_plant or None)."""
    out = [(None, cash[0][0] if cash else None)]
    for k in range(len(cash)):
        nxt = cash[k + 1][0] if k + 1 < len(cash) else None
        out.append((cash[k][1], nxt))
    return out

def gap_index(cash, d):
    return sum(1 for (_, h) in cash if h < d)

def op_cover(sh, rng, p_flip, log):
    names = Counter(sh.get(r, "CD_name") for r in sh.rows if is_cover(sh, r))
    if not names:
        log.append(("cover", "no cover crops in this workbook")); return
    generic = names.most_common(1)[0][0]
    template = next(r for r in sh.rows if is_cover(sh, r))
    units = by_unit(sh)
    keep, made, flipped = [], 0, 0
    covered_gaps = {}
    for u, rows in units.items():
        cash = cash_list(sh, rows)
        sp = spans(sh, rows)
        groups = defaultdict(list)
        for r in rows:
            if not is_cover(sh, r): keep.append(r); continue
            s, e = sp[id(r)]
            mid = s + (e - s) / 2 if s and e else (s or e)
            groups[gap_index(cash, mid) if mid else -1].append(r)
        gaps = gaps_of(cash)
        present = {g for g in groups if g >= 0}
        if p_flip > 0:                                   # misclassification of cover presence
            for g in range(len(gaps)):
                if rng.random() < p_flip:
                    prev_h, next_p = gaps[g]
                    if g in present:
                        present.discard(g); groups.pop(g, None); flipped += 1
                    elif prev_h and next_p and (next_p - prev_h).days >= 60:
                        present.add(g); groups[g] = []; flipped += 1
        for g, grp in groups.items():
            if g < 0:
                keep.extend(grp); continue
            if g not in present: continue
            prev_h, next_p = gaps[g] if g < len(gaps) else (None, None)
            gs = [sp[id(r)] for r in grp]
            start = min([s for s, _ in gs if s], default=None)
            end = max([e for _, e in gs if e], default=None)
            if not grp:                                    # misclassified as present: standard window
                start = prev_h + dt.timedelta(days=COVER_GAP_DAYS) if prev_h else None
                end = next_p - dt.timedelta(days=COVER_GAP_DAYS) if next_p else None
            if start and prev_h and start <= prev_h: start = prev_h + dt.timedelta(days=1)
            if end and next_p and end >= next_p: end = next_p - dt.timedelta(days=1)
            if not start or not end or end <= start: continue
            nr = sh.new_row(grp[0] if grp else template)
            sh.set(nr, "MGT_combo", u)
            sh.set(nr, "CD_name", generic); sh.set(nr, "CD_cat", sh.get(template, "CD_cat"))
            sh.set(nr, "CD_plant_date", to_cell(start)); sh.set(nr, "CD_term_date", to_cell(end))
            sh.set(nr, "CD_harv_date", None)
            for f in ("CD_mix", "CD_group", "CD_seq_num", "CD_yield", "CD_yield_units"):
                if f in sh.cols and f != "CD_seq_num": sh.set(nr, f, None if f != "CD_mix" else "no")
            sh.set(nr, "CD_notes", "coarse record: cover crop present")
            keep.append(nr); made += 1
    order = {id(r): i for i, r in enumerate(sh.rows)}
    keep.sort(key=lambda r: (str(sh.get(r, "MGT_combo")), order.get(id(r), 10**9)))
    n0 = sum(1 for r in sh.rows if is_cover(sh, r))
    sh.rows = keep
    log.append(("cover", f"{n0} cover rows -> {made} generic '{generic}' covers; {flipped} gaps flipped"))


def op_maincrop(sh, log):
    keep, dropped = [], 0
    for u, rows in by_unit(sh).items():
        kept_spans, kept_open = [], set()
        sp = spans(sh, rows)
        for r in rows:
            if is_cover(sh, r): keep.append(r); continue
            s, e = sp[id(r)]
            if s and e and e > s:
                dur = (e - s).days
                ov = max([(min(e, e2) - max(s, s2)).days for s2, e2 in kept_spans] + [0])
                if ov > 0.5 * dur:
                    dropped += 1; continue
                kept_spans.append((s, e))
            elif s or e:                 # one date only (e.g. planted before the record): same date = same stand
                k = (s, e)
                if k in kept_open:
                    dropped += 1; continue
                kept_open.add(k)
            keep.append(r)
    sh.rows = keep
    log.append(("maincrop", f"{dropped} overlapping crop rows dropped"))


# tillage ----------------------------------------------------------------------
def load_tdisc():
    imp = {}
    with open(os.path.join(HERE, "tdisc_implements.csv"), newline="", encoding="utf-8") as f:
        for r in csv.DictReader(f):
            imp[key(r["implement"])] = (float(r["mixing_efficiency"]), float(r["depth_cm"]))
    mp = defaultdict(list)
    with open(os.path.join(HERE, "tdisc_mapping.csv"), newline="", encoding="utf-8") as f:
        for r in csv.DictReader(f):
            mp[key(r["operation"])].append(key(r["implement"]))
    return imp, mp

def key(s): return re.sub(r"\s+", " ", str(s or "").strip()).upper()

def num(v):
    try: return float(v)
    except (TypeError, ValueError): return None

def pass_values(sd, row, method, imp, mp):
    """List of (me, depth_cm) or [stir] for one pass; [] if it cannot be resolved."""
    if method == "stir":
        st = num(sd.get(row, "SD_stir"))
        return [st] if st is not None else []
    me_o = num(sd.get(row, "SD_mixeff")); dep_o = num(sd.get(row, "SD_depth"))
    dep_o = dep_o * 2.54 if dep_o is not None else None
    if me_o is not None and dep_o is not None: return [(me_o, dep_o)]
    k = key(sd.get(row, "SD_equip"))
    if k in imp:
        me, dep = imp[k]
        return [(me_o if me_o is not None else me, dep_o if dep_o is not None else dep)]
    if k in mp and all(i in imp for i in mp[k]):
        return [imp[i] for i in mp[k]]
    return []

def window_ti(vals, method):
    if method == "stir": return min(sum(vals) / STIR_MAX, 1.0)
    S = 0.0
    for me, d in sorted(vals, key=lambda x: (min(x[1], 30), x[0])):
        S += me * max(min(d, 30) - S, 0)
    return min(S / 30, 1.0)

def tdisc_window(d, plant, harvests):
    if plant is None: return "single"
    if d in harvests: return "harvest"
    if d > plant: return "after_planting"
    k = (plant - d).days
    return "planting" if k <= 6 else ("before_planting" if k <= 55 else "field_preparation")

CLASSES = ["Z", "NT", "RT", "CT"]
REP_IMPLEMENT = {"NT": "PLANTER, NO-TILL", "RT": "CULTIVATOR, FIELD", "CT": "PLOW, MOLDBOARD"}  # names only; intensity set below
def designation(ti):
    if ti < 0.001: return "Z"
    return "NT" if ti <= 0.075 else ("RT" if ti <= 0.252 else "CT")

def intervals(cd, unit_rows, first, last):
    """(start, end, plant, harvest days) per cash crop, then calendar years for the rest."""
    cash = defaultdict(list)
    plants = defaultdict(list)
    for r in unit_rows:
        p, h, t = crop_dates(cd, r)
        if is_cash(cd, r) and h:
            cash[(key(cd.get(r, "CD_name")), p)].append(h)
        if p: plants[key(cd.get(r, "CD_name"))].append(p)
    groups = sorted(((max(hs), name, p, sorted(set(hs))) for (name, p), hs in cash.items()), key=lambda x: x[0])
    out, start = [], first
    for e, name, p, hs in groups:
        if e < start or e > last: continue
        if p is None:
            ps = [x for x in plants[name] if x <= hs[0]]
            p = max(ps) if ps else None
        out.append((start, e, p, set(hs))); start = e + dt.timedelta(days=1)
    while start <= last:
        e = min(dt.date(start.year, 12, 31), last)
        out.append((start, e, None, set())); start = e + dt.timedelta(days=1)
    return out

def op_tillage(sd, cd, method, rng, p_flip, log):
    imp, mp = load_tdisc()
    crop_units = by_unit(cd); sd_units = by_unit(sd)
    unresolved = 0
    plan = {}                       # unit -> list of (class, date, ti)
    tis = defaultdict(list)
    for u in set(sd_units) | set(crop_units):
        passes = []
        for r in sd_units.get(u, []):
            d = as_date(sd.get(r, "SD_date"))
            if not d: continue
            v = pass_values(sd, r, method, imp, mp)
            if not v: unresolved += 1; continue
            passes.append((d, v))
        if not passes: continue
        dates = [d for d, _ in passes] + [x for r in crop_units.get(u, []) for x in crop_dates(cd, r) if x]
        first = dt.date(min(dates).year, 1, 1); last = dt.date(max(dates).year, 12, 31)
        rows = []
        for s, e, p, hs in intervals(cd, crop_units.get(u, []), first, last):
            pk = [(d, v) for d, v in passes if s <= d <= e]
            if not pk: continue
            win = defaultdict(list)
            for d, v in pk: win[tdisc_window(d, p, hs)].extend(v)
            ti = max(window_ti(v, method) for v in win.values())
            cls = designation(ti)
            if cls != "Z": tis[cls].append(ti)
            top = max(pk, key=lambda x: window_ti(x[1], method))[0]
            when = p if p and s <= p <= e else top
            rows.append([cls, when, ti])
        plan[u] = rows
    typical = {c: statistics.median(v) for c, v in tis.items() if v}
    flipped = 0
    if p_flip > 0:
        for rows in plan.values():
            for x in rows:
                if x[0] != "Z" and rng.random() < p_flip:
                    i = CLASSES.index(x[0]); j = i + rng.choice([-1, 1])
                    j = min(max(j, 1), 3)
                    if j != i: x[0] = CLASSES[j]; flipped += 1
    # replace every pass by one pass per tilled interval at the class's typical intensity
    template = sd.rows[0] if sd.rows else None
    keep = [r for r in sd.rows if sd.get(r, "MGT_combo") not in plan]
    made = 0
    for u, rows in plan.items():
        for cls, when, _ in rows:
            if cls == "Z" or cls not in typical: continue
            t = typical[cls]
            nr = sd.new_row(template)
            sd.set(nr, "MGT_combo", u); sd.set(nr, "SD_date", to_cell(when))
            for f in ("SD_id", "SD_phase", "SD_goal", "SD_name", "SD_equip_cat"):
                if f in sd.cols: sd.set(nr, f, None)
            sd.set(nr, "SD_equip", REP_IMPLEMENT[cls])
            for f in ("SD_tdisc_mixeff", "SD_tdisc_depth"):          # display-only lookups
                if f in sd.cols: sd.set(nr, f, None)
            if method == "stir":
                sd.set(nr, "SD_stir", round(t * STIR_MAX, 3))
                if "SD_depth" in sd.cols: sd.set(nr, "SD_depth", None)
            else:
                sd.set(nr, "SD_mixeff", 1.0)                      # with SD_depth, overrides the implement; sd.set(nr, "SD_depth", round(t * 30 / 2.54, 4))
                if "SD_stir" in sd.cols: sd.set(nr, "SD_stir", None)
            sd.set(nr, "SD_notes", f"coarse record: tillage class {cls}")
            keep.append(nr); made += 1
    if not keep and template is not None:
        # every interval was zero-disturbance (e.g. a no-till study): keep one zero-mixing
        # placeholder pass, as agcros2shmi.py does, since an empty sheet stops
        # prepare_shmi_inputs() (package to-do item 2)
        nr = sd.new_row(template)
        for f in ("SD_id", "SD_phase", "SD_goal", "SD_name", "SD_equip_cat", "SD_tdisc_mixeff", "SD_tdisc_depth", "SD_depth"):
            if f in sd.cols: sd.set(nr, f, None)
        sd.set(nr, "SD_equip", "SPRAYER/AERIAL/BROADCAST APPLICATION")
        if method == "stir": sd.set(nr, "SD_stir", 0)
        elif "SD_mixeff" in sd.cols: sd.set(nr, "SD_mixeff", 0)
        if "SD_notes" in sd.cols: sd.set(nr, "SD_notes", "coarse record: no tillage; zero-disturbance placeholder so the sheet is not empty")
        keep.append(nr)
    sd.rows = keep
    n_iv = Counter(x[0] for rows in plan.values() for x in rows)
    log.append(("tillage", f"intervals by class {dict(n_iv)}; typical TI {({k: round(v, 3) for k, v in typical.items()})}; "
                           f"{made} passes written; {flipped} classes flipped; {unresolved} passes unresolved (ignored)"))


def op_noamend(wb, log):
    msg = []
    if "Soil_Amendments" in wb.sheetnames:
        # organic amendments are re-entered as zero fertilizer rows rather than deleted: SHMI
        # credits only "Organic", and an emptied sheet would trip the empty-sheet problem in
        # prepare_shmi_inputs() (package to-do item 2)
        sa = Sheet(wb["Soil_Amendments"]); n = 0
        for r in sa.rows:
            if str(sa.get(r, "SA_cat") or "").strip().lower() == "organic":
                sa.set(r, "SA_cat", "Fertilizer"); sa.set(r, "SA_source", "coarse record: organic amendment not observed")
                for f in ("SA_rate", "SA_N", "SA_P2O5", "SA_K2O"): sa.set(r, f, 0)
                n += 1
        sa.write(); msg.append(f"{n} organic amendments removed (kept as zero fertilizer rows)")
    if "Animal_Diversity" in wb.sheetnames:
        ad = Sheet(wb["Animal_Diversity"]); n0 = len(ad.rows); ad.rows = []
        ad.write(); msg.append(f"{n0} grazing records removed")
    log.append(("noamend", "; ".join(msg)))


DATE_FIELDS = {"Crop_Diversity": ("CD_plant_date", "CD_harv_date", "CD_term_date"),
               "Soil_Disturbance": ("SD_date",), "Soil_Amendments": ("SA_date",),
               "Animal_Diversity": ("AD_start_date", "AD_end_date")}

def op_last(wb, n, log):
    last = defaultdict(int)
    for name, fields in DATE_FIELDS.items():
        if name not in wb.sheetnames: continue
        sh = Sheet(wb[name])
        for r in sh.rows:
            for f in fields:
                d = as_date(sh.get(r, f))
                if d: last[sh.get(r, "MGT_combo")] = max(last[sh.get(r, "MGT_combo")], d.year)
    dropped = 0
    for name, fields in DATE_FIELDS.items():
        if name not in wb.sheetnames: continue
        sh = Sheet(wb[name]); keep = []
        for r in sh.rows:
            ds = [as_date(sh.get(r, f)) for f in fields]
            ds = [d for d in ds if d]
            cut = last[sh.get(r, "MGT_combo")] - n + 1
            if ds and max(ds).year < cut: dropped += 1; continue
            keep.append(r)
        sh.rows = keep; sh.write()
    log.append((f"last:{n}", f"{dropped} records before each unit's last {n} years removed"))


def shift_in_year(d, n, rng):
    """d moved by a random 0-n days either way, reflected back if it would leave its
    calendar year: the scoring window runs from 1 January of the first year with a record
    to 31 December of the last, so a date pushed across New Year would add an empty year."""
    k = rng.randint(-n, n)
    x = d + dt.timedelta(days=k)
    if x.year != d.year: x = d - dt.timedelta(days=k)
    return x if x.year == d.year else d

def op_jitter(cd, sd, n, rng, log):
    moved = kept = 0
    for r in cd.rows:
        p, h, t = crop_dates(cd, r)
        new = [shift_in_year(d, n, rng) if d else None for d in (p, h, t)]
        if new[0] and any(x and x <= new[0] for x in new[1:]): kept += 1; continue
        for f, d in zip(("CD_plant_date", "CD_harv_date", "CD_term_date"), new):
            if d: cd.set(r, f, to_cell(d))
        moved += 1
    for r in sd.rows:
        d = as_date(sd.get(r, "SD_date"))
        if d: sd.set(r, "SD_date", to_cell(shift_in_year(d, n, rng)))
    log.append((f"jitter:{n}", f"{moved} crop rows and all tillage dates shifted; {kept} crop rows kept"))


def op_sjitter(cd, sd, n, rng, site_of, log):
    shift = {}
    def move(u, d):
        key = (site_of(u), d)
        if key not in shift: shift[key] = rng.randint(-n, n)
        x = d + dt.timedelta(days=shift[key])
        if x.year != d.year: x = d - dt.timedelta(days=shift[key])
        return x if x.year == d.year else d
    moved = kept = 0
    for r in cd.rows:
        u = cd.get(r, "MGT_combo"); p, h, t = crop_dates(cd, r)
        new = [move(u, d) if d else None for d in (p, h, t)]
        if new[0] and any(x and x <= new[0] for x in new[1:]): kept += 1; continue
        for f, d in zip(("CD_plant_date", "CD_harv_date", "CD_term_date"), new):
            if d: cd.set(r, f, to_cell(d))
        moved += 1
    for r in sd.rows:
        d = as_date(sd.get(r, "SD_date"))
        if d: sd.set(r, "SD_date", to_cell(move(sd.get(r, "MGT_combo"), d)))
    log.append((f"sjitter:{n}", f"{moved} crop rows and all tillage dates shifted by date (shared at each site, "
                                f"{len(shift)} distinct dates); {kept} crop rows kept"))


def freeze_formulas(wb, path):
    """Replace every formula by the value Excel last calculated for it. openpyxl writes
    formulas without their cached results, and readxl (prepare_shmi_inputs) reads only
    cached results, so a workbook whose MGT_combo is a formula (the hand-built DRIVES
    workbooks) would otherwise come back with blank unit IDs."""
    vals = openpyxl.load_workbook(path, data_only=True)
    n = 0
    for ws in wb.worksheets:
        wv = vals[ws.title]
        for row in ws.iter_rows():
            for c in row:
                if isinstance(c.value, str) and c.value.startswith("="):
                    c.value = wv[c.coordinate].value; n += 1
    return n


# ----------------------------------------------------------------------------- main
def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--in", dest="inp", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--ops", required=True, help="e.g. dates,cover,maincrop,tillage,noamend,last:2,jitter:14,misclass:0.1")
    ap.add_argument("--tillage", choices=["epa", "stir"], default="epa")
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--site-regex", default=r"^(NAPESHM_[^_]+)_", help="site = first group of MGT_combo; NONE or no match = the whole workbook is one site")
    ap.add_argument("--log", default=None)
    a = ap.parse_args(argv)

    ops = [o.strip() for o in a.ops.split(",") if o.strip()]
    known = {"dates", "cover", "maincrop", "tillage", "noamend"}
    par = {}
    for o in ops:
        k, _, v = o.partition(":")
        if k in ("last", "jitter", "sjitter"): par[k] = int(v)
        elif k == "misclass": par[k] = float(v)
        elif k not in known: sys.exit(f"unknown operation: {o}")
    flags = {o.partition(":")[0] for o in ops}
    rng = random.Random(a.seed)
    rx = None if a.site_regex in ("", "NONE") else re.compile(a.site_regex)
    site_of = lambda u: (rx.match(str(u)).group(1) if rx and rx.match(str(u)) else "ALL")

    wb = openpyxl.load_workbook(a.inp)
    freeze_formulas(wb, a.inp)
    log = []
    cd = Sheet(wb["Crop_Diversity"]); sd = Sheet(wb["Soil_Disturbance"])
    p_flip = par.get("misclass", 0.0)
    if "dates" in flags: op_dates(cd, site_of, log)
    if "cover" in flags: op_cover(cd, rng, p_flip, log)
    if "maincrop" in flags: op_maincrop(cd, log)
    if "jitter" in flags: op_jitter(cd, sd, par["jitter"], rng, log)
    if "sjitter" in flags: op_sjitter(cd, sd, par["sjitter"], rng, site_of, log)
    if "tillage" in flags: op_tillage(sd, cd, a.tillage, rng, p_flip, log)
    cd.write(); sd.write()
    if "noamend" in flags: op_noamend(wb, log)
    if "last" in flags: op_last(wb, par["last"], log)
    os.makedirs(os.path.dirname(os.path.abspath(a.out)), exist_ok=True)
    wb.save(a.out)
    for k, m in log: print(f"  {k}: {m}")
    if a.log:
        new = not os.path.exists(a.log)
        with open(a.log, "a", newline="", encoding="utf-8") as f:
            w = csv.writer(f)
            if new: w.writerow(["workbook", "ops", "step", "detail"])
            for k, m in log: w.writerow([os.path.basename(a.inp), a.ops, k, m])


if __name__ == "__main__":
    main()
