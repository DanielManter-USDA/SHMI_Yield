"""WICST CS5: field pea in the oat companion crop (documented correction).

Posner, Baldock & Hedtcke (2008, Agron. J. 100:253-260, Table 1): in the organic forage
system CS5 (corn - oat/alfalfa - alfalfa), "the alfalfa in Phase 1 of CS5 was companion
seeded with the oat-pea mix". The DRIVES workbook records the oats and the alfalfa but not
the pea. This adds a 'Pea, Field' row with the oat row's dates to every CS5 oat seeding
(--to YEAR limits it; Posner et al. document 1990-2002, later years assume the design was
unchanged). Only the Crop_Diversity sheet's XML is edited (see ohio_meadow_species.py).

Usage: python python/wicst_cs5_pea.py Inputs/yield_workbooks/SHMI_WICST_ARS.xlsx [out.xlsx] [--to 2020]
Idempotent (rows carry the tag below).
"""
import argparse, copy, datetime as dt, os, re, shutil, sys, tempfile, zipfile
import openpyxl
from lxml import etree
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ohio_meadow_species import sheet_part, set_str, q   # noqa: E402

TAG = "field pea in the CS5 oat companion (Posner et al. 2008, Table 1)"

def main():
    ap = argparse.ArgumentParser(); ap.add_argument("workbook"); ap.add_argument("out", nargs="?")
    ap.add_argument("--to", type=int, default=9999); a = ap.parse_args()
    wb = openpyxl.load_workbook(a.workbook, read_only=True, data_only=True)
    rows = list(wb["Crop_Diversity"].iter_rows(values_only=True))
    hdr = next(i for i, r in enumerate(rows[:10]) if len(r) > 1 and "MGT_combo" in r)
    col = {v: i for i, v in enumerate(rows[hdr]) if v}
    letter = {k: openpyxl.utils.get_column_letter(i + 1) for k, i in col.items()}
    data = [(i + 1, r) for i, r in enumerate(rows) if i > hdr and r[col["MGT_combo"]]]
    if any(TAG in str(r[col["CD_notes"]] or "") for _, r in data): print("already converted:", a.workbook); return
    yr = lambda v: v.year if isinstance(v, (dt.date, dt.datetime)) else int(str(v)[:4])
    oats = [(rn, r) for rn, r in data if "CS5" in str(r[col["MGT_combo"]]) and str(r[col["CD_name"]]).strip().lower() == "spring oats"
            and r[col["CD_plant_date"]] and yr(r[col["CD_plant_date"]]) <= a.to]
    last = max(rn for rn, _ in data); wb.close()
    with zipfile.ZipFile(a.workbook) as z:
        part = sheet_part(z, "Crop_Diversity"); sx = etree.fromstring(z.read(part)); sd = sx.find(q("sheetData"))
        byrow = {int(r.get("r")): r for r in sd.findall(q("row"))}; nxt = last + 1
        for rn, r in oats:
            new = copy.deepcopy(byrow[rn]); new.set("r", str(nxt))
            for c in new.findall(q("c")):
                c.set("r", re.sub(r"\d+$", str(nxt), c.get("r")))
                if c.find(q("f")) is not None: c.remove(c.find(q("f")))
            set_str(new, f"{letter['CD_name']}{nxt}", "Pea, Field")
            if "CD_group" in letter: set_str(new, f"{letter['CD_group']}{nxt}", "Legume_cool")
            y = yr(r[col["CD_plant_date"]])
            set_str(new, f"{letter['CD_notes']}{nxt}", TAG + ("" if y <= 2002 else "; after 2002 ASSUMED (design unchanged)"))
            for k in ("CD_yield", "CD_yield_units"):             # the pea has no yield of its own: drop those cells
                if k in letter:
                    c = next((c for c in new.findall(q("c")) if c.get("r") == f"{letter[k]}{nxt}"), None)
                    if c is not None: new.remove(c)
            if nxt in byrow: sd.replace(byrow[nxt], new)
            else:
                later = [byrow[k] for k in sorted(byrow) if k > nxt]
                (later[0].addprevious(new) if later else sd.append(new))
            byrow[nxt] = new; nxt += 1
        dim = sx.find(q("dimension"))
        if dim is not None:
            m = re.match(r"([A-Z]+\d+):([A-Z]+)(\d+)", dim.get("ref") or "")
            if m and int(m.group(3)) < nxt - 1: dim.set("ref", f"{m.group(1)}:{m.group(2)}{nxt - 1}")
        xml = etree.tostring(sx, xml_declaration=True, encoding="UTF-8", standalone=True)
        fd, tmp = tempfile.mkstemp(suffix=".xlsx"); os.close(fd)
        with zipfile.ZipFile(tmp, "w", zipfile.ZIP_DEFLATED) as zo:
            for it in z.infolist(): zo.writestr(it, xml if it.filename == part else z.read(it.filename))
    shutil.move(tmp, a.out or a.workbook)
    print(f"{a.workbook}: field pea added to {len(oats)} CS5 oat seedings -> {a.out or a.workbook}")

if __name__ == "__main__":
    main()
