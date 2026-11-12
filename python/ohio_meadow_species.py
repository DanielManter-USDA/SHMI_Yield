"""Ohio (Triplett-Van Doren) meadow phase as a grass-legume mixture (documented assumption).

The corn-oats-meadow rotation at Wooster and Hoytville seeded a forage stand with oats as a
nurse crop. de Camargo Santos et al. (2025, Field Crops Res. 331, 109993) describe it as "a
grass-legume mix including species such as orchard grass, tall fescue, alfalfa, red clover
and white clover"; Dick & Van Doren (1985, Agron. J. 77, 459-465) as "alfalfa or mixed grass
meadow". Which species were sown in which years is not reported. The DRIVES workbooks record
one perennial per stand ("perennial mix" at Wooster, "alfalfa" at Hoytville).

Rule (option B, agreed 9 Nov 2026): each stand is one legume plus one grass, with the same
dates as the recorded stand. The pair changes every BLOCK years of seeding date (default 9,
i.e. three 3-year cycles), counted from the workbook's first stand:
    legume  Alfalfa -> Clover, Red -> Clover, White -> Alfalfa ...
    grass   Orchardgrass -> Fescue, Tall -> Orchardgrass ...
so every plot seeded in the same calendar year gets the same pair. Sensitivity (option A:
Alfalfa + Orchardgrass in every cycle): --block 0.

The workbook is edited at the XML level, and only its Crop_Diversity sheet: every other part
of the file (Mgt_Unit and Soil_Disturbance formulas with their calculated values, dropdowns,
lookups) is copied byte for byte. (Saving with openpyxl drops the calculated values of the
formulas, so R then reads MGT_combo in Mgt_Unit as empty.)

Hoytville has field notes (Dick, Triplett & Van Doren 2013, OSU Knowledge Bank, "Hoytville
Field Notes"), so its stands follow the record instead of the rule: alfalfa with grass seeded
after the oats through 1982 ("No grass seeded after oats this year", 1979; NT hay "mostly
grass"), and from 1983 mammoth red clover seeded with the oats "instead of seeding of grass"
(note of 29 Jul 1982; "oat-clover plots" from 1985). --documented ho applies it (default for
workbooks whose name contains OHTVD_Ho).
Wooster follows Dick, Van Doren, Triplett & Henry (1986, OARDC Research Bulletin 1181): alfalfa
1962-1967, alfalfa-orchardgrass 1968-1975, Kentucky 31 tall fescue 1976, perennial ryegrass
1977-1984; after 1984 not documented, red clover assumed as at Hoytville (--documented wo,
default for OHTVD_Wo). --documented none with --block 9 gives option B as a sensitivity.

Usage:  python python/ohio_meadow_species.py Inputs/yield_workbooks/SHMI_OHTVD_Wo.xlsx [out.xlsx] [--block 9]
Idempotent: a converted workbook (its notes carry the tag below) is left unchanged.
"""
import argparse, copy, datetime as dt, re, shutil, sys, tempfile, zipfile, os
import openpyxl
from lxml import etree

LEGUMES = ["Alfalfa", "Clover, Red", "Clover, White"]
GRASSES = ["Orchardgrass", "Fescue, Tall"]
STAND_NAMES = {"perennial mix", "alfalfa"}
TAG = "meadow as grass-legume mix (assumed; de Camargo Santos et al. 2025)"
NS = "http://schemas.openxmlformats.org/spreadsheetml/2006/main"
RNS = "http://schemas.openxmlformats.org/officeDocument/2006/relationships"
q = lambda t: "{%s}%s" % (NS, t)


def sheet_part(z, name):
    wbx = etree.fromstring(z.read("xl/workbook.xml"))
    rid = next(s.get("{%s}id" % RNS) for s in wbx.iter(q("sheet")) if s.get("name") == name)
    rels = etree.fromstring(z.read("xl/_rels/workbook.xml.rels"))
    tgt = next(r.get("Target") for r in rels if r.get("Id") == rid)
    return tgt.lstrip("/") if tgt.startswith("/") else "xl/" + tgt


def set_str(row, ref, text):
    c = next((c for c in row.findall(q("c")) if c.get("r") == ref), None)
    if c is None:                                              # cell absent: insert in column order
        c = etree.Element(q("c")); c.set("r", ref)
        col = re.match(r"[A-Z]+", ref).group(0)
        key = lambda r_: (len(r_), r_)
        after = [x for x in row.findall(q("c")) if key(re.match(r"[A-Z]+", x.get("r")).group(0)) < key(col)]
        (after[-1].addnext(c) if after else row.insert(0, c))
    for ch in list(c): c.remove(ch)
    c.set("t", "inlineStr")
    is_ = etree.SubElement(c, q("is")); t = etree.SubElement(is_, q("t")); t.text = text
    t.set("{http://www.w3.org/XML/1998/namespace}space", "preserve")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("workbook"); ap.add_argument("out", nargs="?")
    ap.add_argument("--block", type=int, default=9, help="years per species pair; 0 = one pair throughout")
    ap.add_argument("--documented", choices=["auto", "ho", "wo", "none"], default="auto",
                    help="documented species: ho = Hoytville field notes, wo = Wooster bulletin "
                         "(auto: by the file name, OHTVD_Ho or OHTVD_Wo)")
    a = ap.parse_args()
    base = os.path.basename(a.workbook)
    ho = a.documented == "ho" or (a.documented == "auto" and "OHTVD_Ho" in base)
    wo = a.documented == "wo" or (a.documented == "auto" and "OHTVD_Wo" in base)
    # read: header, columns and the stands (values only)
    wb = openpyxl.load_workbook(a.workbook, read_only=True, data_only=True)
    ws = wb["Crop_Diversity"]
    rows = list(ws.iter_rows(values_only=True))
    hdr = next(i for i, r in enumerate(rows[:10]) if len(r) > 1 and r[1] == "MGT_combo")
    col = {v: i for i, v in enumerate(rows[hdr]) if v}
    letter = {k: openpyxl.utils.get_column_letter(i + 1) for k, i in col.items()}
    data = [(i + 1, r) for i, r in enumerate(rows) if i > hdr and r[col["MGT_combo"]]]
    if any(TAG in str(r[col["CD_notes"]] or "") for _, r in data):
        print("already converted:", a.workbook); return
    stands = [(rn, r) for rn, r in data if str(r[col["CD_cat"]]).strip() == "Perennial"
              and str(r[col["CD_name"]]).strip().lower() in STAND_NAMES]
    if not stands: sys.exit("no meadow stands found in " + a.workbook)
    yr = lambda v: v.year if isinstance(v, (dt.date, dt.datetime)) else int(str(v)[:4])
    y0 = min(yr(r[col["CD_plant_date"]]) for _, r in stands)
    last = max(rn for rn, _ in data)
    wb.close()

    # edit Crop_Diversity's XML only
    with zipfile.ZipFile(a.workbook) as z:
        part = sheet_part(z, "Crop_Diversity")
        sx = etree.fromstring(z.read(part))
        sd = sx.find(q("sheetData"))
        byrow = {int(r.get("r")): r for r in sd.findall(q("row"))}
        nxt = last + 1; counts = {}
        for rn, r in stands:
            old = r[col["CD_name"]]; prev = r[col["CD_notes"]]
            if wo:
                y = yr(r[col["CD_plant_date"]])
                sp = (["Alfalfa"] if y <= 1967 else ["Alfalfa", "Orchardgrass"] if y <= 1975 else
                      ["Fescue, Tall"] if y == 1976 else ["Ryegrass, Perennial"] if y <= 1984 else ["Clover, Red"])
                leg, gr = sp[0], (sp[1] if len(sp) > 1 else None)
                note = (f"{TAG}: recorded as '{old}'; Wooster: Dick et al. 1986, OARDC Res. Bull. 1181 (alfalfa 1962-67, "
                        f"alfalfa-orchardgrass 1968-75, Kentucky 31 tall fescue 1976, perennial ryegrass 1977-84)"
                        + ("" if y <= 1984 else "; after 1984 not documented: red clover ASSUMED, as at Hoytville from 1983"))
            elif ho:
                y = yr(r[col["CD_plant_date"]])
                leg, gr = ("Alfalfa", "Orchardgrass") if y <= 1982 else ("Clover, Red", None)
                src = ("Hoytville field notes" if "OHTVD_Ho" in os.path.basename(a.workbook)
                       else "assumed from the Hoytville field notes (same experiment and protocol; no Wooster notes found)")
                note = (f"{TAG}: recorded as '{old}'; {src}: alfalfa with grass through 1982 "
                        f"(grass species not named: orchardgrass assumed), mammoth red clover alone from 1983")
            else:
                k = 0 if a.block <= 0 else (yr(r[col["CD_plant_date"]]) - y0) // a.block
                leg, gr = LEGUMES[k % len(LEGUMES)], GRASSES[k % len(GRASSES)]
                note = (f"{TAG}: recorded as '{old}'; pair {k + 1} ({leg} + {gr}), species switched every {a.block} years"
                        if a.block > 0 else f"{TAG}: recorded as '{old}'; {leg} + {gr} throughout")
            row = byrow[rn]
            new = copy.deepcopy(row)
            set_str(row, f"{letter['CD_name']}{rn}", leg)
            if "CD_group" in letter:
                set_str(row, f"{letter['CD_group']}{rn}", "Legume_cool" if leg in ("Alfalfa", "Clover, Red", "Clover, White") else "Grass_cool")
            set_str(row, f"{letter['CD_notes']}{rn}", f"{prev}; {note}" if prev else note)
            if gr is None:                                        # legume alone: no companion row
                counts[(leg, "-")] = counts.get((leg, "-"), 0) + 1; continue
            new.set("r", str(nxt))
            for c in new.findall(q("c")):
                c.set("r", re.sub(r"\d+$", str(nxt), c.get("r")))
                if c.find(q("f")) is not None: c.remove(c.find(q("f")))
            set_str(new, f"{letter['CD_name']}{nxt}", gr)
            if "CD_group" in letter: set_str(new, f"{letter['CD_group']}{nxt}", "Grass_cool")
            set_str(new, f"{letter['CD_notes']}{nxt}", f"{prev}; {note}" if prev else note)
            if nxt in byrow: sd.replace(byrow[nxt], new)
            else:
                later = [byrow[k_] for k_ in sorted(byrow) if k_ > nxt]
                (later[0].addprevious(new) if later else sd.append(new))
            byrow[nxt] = new
            nxt += 1; counts[(leg, gr)] = counts.get((leg, gr), 0) + 1
        dim = sx.find(q("dimension"))
        if dim is not None:
            ref = dim.get("ref"); m = re.match(r"([A-Z]+\d+):([A-Z]+)(\d+)", ref or "")
            if m and int(m.group(3)) < nxt - 1: dim.set("ref", f"{m.group(1)}:{m.group(2)}{nxt - 1}")
        new_xml = etree.tostring(sx, xml_declaration=True, encoding="UTF-8", standalone=True)
        out = a.out or a.workbook
        fd, tmp = tempfile.mkstemp(suffix=".xlsx"); os.close(fd)
        with zipfile.ZipFile(tmp, "w", zipfile.ZIP_DEFLATED) as zo:
            for it in z.infolist():
                zo.writestr(it, new_xml if it.filename == part else z.read(it.filename))
    shutil.move(tmp, out)
    print(f"{a.workbook}: {len(stands)} stands from {y0}; pairs: " +
          "; ".join(f"{l} + {g}: {n}" for (l, g), n in counts.items()) + f" -> {out}")

if __name__ == "__main__":
    main()
