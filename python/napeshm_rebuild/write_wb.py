import json, openpyxl, datetime as dt, sys
from openpyxl.styles import Font, Alignment
# usage: python write_wb.py <NAPESHM_test.xlsx (layout)> <converted.json> <NAPESHM_rebuilt.xlsx>
SRC = sys.argv[1] if len(sys.argv) > 1 else "NAPESHM_test.xlsx"
JSN = sys.argv[2] if len(sys.argv) > 2 else "converted.json"
OUT = sys.argv[3] if len(sys.argv) > 3 else "NAPESHM_rebuilt.xlsx"
j = json.load(open(JSN))
USER = "jennifer.moore"

wb  = openpyxl.load_workbook(SRC)
wbv = openpyxl.load_workbook(SRC, data_only=True)

def to_date(v):
    return dt.date.fromisoformat(v) if isinstance(v, str) and len(v) == 10 and v[4] == "-" else v

def rewrite(sheet, rows, keycol_name="MGT_combo"):
    ws = ws_ = wb[sheet]
    hdr = {c.value: c.column for c in ws[4] if c.value}
    date_fmt = "mm-dd-yy"
    old_max = ws.max_row
    for r in ws.iter_rows(min_row=5, max_row=old_max):
        for c in r:
            c.value = None
    for i, row in enumerate(rows, start=5):
        ws.cell(i, hdr["user_name"]).value = USER
        for k, v in row.items():
            if k not in hdr:
                continue
            v = to_date(v) if ("date" in k.lower()) else v
            cell = ws.cell(i, hdr[k]); cell.value = v
            if isinstance(v, dt.date):
                cell.number_format = date_fmt
    last = 4 + len(rows)
    if old_max > last:
        ws.delete_rows(last + 1, old_max - last)
    return last

# Replace formulas by their cached values in sheets we keep unchanged, so R reads them
for sheet in ["User_Info", "Soil_Test"]:
    for r_f, r_v in zip(wb[sheet].iter_rows(), wbv[sheet].iter_rows()):
        for cf, cv in zip(r_f, r_v):
            if isinstance(cf.value, str) and cf.value.startswith("="):
                cf.value = cv.value

# Remove treatments without management records from Soil_Test as well
removed = set(j.get("removed", []))
ws = wb["Soil_Test"]
hdr = {c.value: c.column for c in ws[4] if c.value}
for r in range(ws.max_row, 4, -1):
    if ws.cell(r, hdr["MGT_combo"]).value in removed:
        ws.delete_rows(r, 1)

n = {}
n["Mgt_Unit"]         = rewrite("Mgt_Unit", j["mgt"])
n["Crop_Diversity"]   = rewrite("Crop_Diversity", j["crop"])
n["Soil_Disturbance"] = rewrite("Soil_Disturbance", j["dist"])
n["Soil_Amendments"]  = rewrite("Soil_Amendments", j["amend"])
n["Animal_Diversity"] = rewrite("Animal_Diversity", j["animal"])

# Conversion log sheet
ws = wb.create_sheet("Conversion_Log")
ws.append(["rule", "MGT_combo", "date", "detail"])
for c in ws[1]: c.font = Font(bold=True)
for rule, unit, date, detail in j["log"]:
    ws.append([rule, unit, to_date(date) if date not in (None, "None") else None, detail])
for col, w in zip("ABCD", [7, 22, 11, 110]): ws.column_dimensions[col].width = w

wb.calculation.fullCalcOnLoad = True
wb.save(OUT)
print(n)
