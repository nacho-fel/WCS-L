"""Extract WCS specs from the GeoVision workbooks (Sandia, Lowry et al. 2017).

Reads either a folder of .xlsx files or a JSON dump of their cells, and writes
``data/wcs/geovision_specs.json.gz``: for every workbook, the inputs and formula
variants as a WCS spec, plus the total the spreadsheet itself computes (cell
'Well Cost Summary'!C29). The tests use this file to check that ``wcsl.wcs``
reproduces every workbook.

Usage:
    python tools/extract_geovision_specs.py --xlsx-dir PATH/TO/GeoVision_files
    python tools/extract_geovision_specs.py --json-dump all.json
"""
from __future__ import annotations

import argparse
import gzip
import json
import re
from pathlib import Path

COLS = ["C", "F", "I", "L", "O", "R"]    # cost columns of the six intervals
TCOLS = ["B", "E", "H", "K", "N", "Q"]   # time columns
DV_KEYS = {7: "bha", 8: "csg_price", 9: "cem_price", 10: "run", 11: "cem_h", 14: "trip_speed", 15: "bha_h",
           16: "log_speed", 17: "log_svc", 18: "wh_h", 19: "mud0", 20: "mud_day", 21: "dir_eng",
           22: "dir_tool", 23: "extra"}


def load_xlsx(path: Path) -> dict:
    import openpyxl
    wf = openpyxl.load_workbook(path, data_only=False)
    wv = openpyxl.load_workbook(path, data_only=True)
    d = {"sheets": {}}
    for ws in wf.worksheets:
        vs = wv[ws.title]
        cells = {}
        for row in ws.iter_rows():
            for c in row:
                if c.value is None:
                    continue
                f = c.value if isinstance(c.value, (int, float, str)) else str(c.value)
                cells[c.coordinate] = [f, vs[c.coordinate].value]
        d["sheets"][ws.title] = {"cells": cells}
    return d


def spec_from_workbook(d: dict) -> tuple[dict, float]:
    W = d["sheets"]["Well Description"]["cells"]
    IC = d["sheets"]["Interval Costs"]["cells"]
    DV = d["sheets"]["Default Values"]["cells"]
    WS = d["sheets"]["Well Cost Summary"]["cells"]

    def v(S, a, default=0):
        x = S.get(a, [None, None])[1]
        return default if x is None else x

    def f(a):
        x = IC.get(a, [""])[0]
        return x if isinstance(x, str) else ""

    hor = W.get("M2", [None])[0] == "Measured Depth"
    params = {name: v(DV, f"I{k}") for k, name in DV_KEYS.items()}
    params.update(rig=v(IC, "E5"), cont=v(W, ("R" if hor else "Q") + "19"),
                  mob=v(WS, "C24"), site=v(WS, "C25"), eng=v(WS, "C26"), wellhead=v(WS, "C27"))
    intervals = []
    for i, c in enumerate(COLS):
        t = TCOLS[i]
        g = dict(end=v(IC, c + "9"), D=v(IC, c + "10"), rop=v(IC, c + "11"), bit=v(IC, c + "12"),
                 Dc=v(IC, c + "13"), w=v(IC, c + "14"), CL=v(IC, c + "15"))
        g["mult"] = {k: v(IC, f"{c}{r}", 1) for k, r in zip(
            ["bit", "bha", "csg", "cem", "mud", "dir", "logt"], range(35, 42))}
        if i > 0:
            g["n"] = v(IC, c + "44", 1)
        g["trip_mult"] = v(IC, c + "45", 1)
        wd_row = 3 + i
        g["log"] = int(v(W, ("W" if hor else "V") + str(wd_row)) == 1)
        g["dir"] = int(v(W, ("X" if hor else "W") + str(wd_row)) == 1)
        # drilling-time formula variant
        if i > 0:
            fe = f(t + "17")
            if fe.endswith("/" + c + "11") or re.search("/" + c + r"11\)?$", fe):
                g["drill_mode"] = "plain"
            elif "IF(" in fe:
                ref = re.search(r"'Well Description'!(\w+)", fe).group(1)
                g["drill_mode"] = "switch"
                g["drill_switch"] = 0 if v(W, ref) == 0 else 1
            else:
                g["drill_mode"] = "eps"
                g["drill_eps"] = 0.001 if "+0.001)" in fe else 0.0001
        # bit price formula variant
        fb = f(c + "19")
        if "$I$3" in fb:
            g["bit_a"], g["bit_b"] = v(IC, "I3"), v(IC, "I4")
        g["bit_eps"] = 0.001 if "+0.001)" in fb else (0.0001 if "+0.0001)" in fb else 0.0)
        g["bit_n"] = (c + "44") in fb
        g["csg_n"] = (c + "44") in f(c + "22")
        g["cem_n"] = (c + "44") in f(c + "24")
        g["cem_pi"] = "PI()" in f(c + "24")
        ft = f(t + "23")
        g["pump"] = (t + "51") in ft
        g["pump_rate"] = v(IC, t + "49", 5) or 5
        g["cemt_oh"] = (c + "43") in ft
        g["logt_n"] = (c + "44") in f(t + "28")
        g["logc_n"] = (c + "44") in f(c + "28")
        fw = f(t + "30")
        g["wh_oh"] = (c + "43") in fw
        g["wh_n"] = (c + "44") in fw
        mm = re.search(r"IF\('Well Description'!(\w+)=0", fw)
        g["wh_switch"] = (0 if v(W, mm.group(1)) == 0 else 1) if mm else 1
        g["trouble_h"], g["trouble_c"] = v(IC, t + "25"), v(IC, c + "26")
        hh = re.search(r"\$I22\*(\d+)\)", f(c + "27"))
        g["dir_hours"] = float(hh.group(1)) if hh else None
        intervals.append(g)
    expected = v(WS, "C29")
    return dict(params=params, intervals=intervals), expected


def describe(name: str) -> dict:
    parts = name.split("/")
    m = re.search(r"(\d+)\.xlsx$", name)
    return dict(file=name, layout=parts[0] if len(parts) > 1 else "", scenario=parts[1] if len(parts) > 2 else "",
                diameter=parts[2] if len(parts) > 3 else "", depth_m=int(m.group(1)) if m else None)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--xlsx-dir")
    ap.add_argument("--json-dump")
    ap.add_argument("--out", default=str(Path(__file__).resolve().parents[1] / "data/wcs/geovision_specs.json.gz"))
    a = ap.parse_args()
    books: dict[str, dict] = {}
    if a.json_dump:
        books = json.load(open(a.json_dump))
    elif a.xlsx_dir:
        root = Path(a.xlsx_dir)
        for p in sorted(root.rglob("*.xlsx")):
            books[str(p.relative_to(root))] = load_xlsx(p)
    out = []
    for name, d in books.items():
        try:
            spec, expected = spec_from_workbook(d)
        except (KeyError, AttributeError):
            continue  # not a GeoVision cost workbook (e.g. the broken example file)
        if not expected:
            continue
        out.append(dict(**describe(name), expected_total=expected, spec=spec))
    with gzip.open(a.out, "wt") as fh:
        json.dump(out, fh)
    print(f"wrote {len(out)} specs to {a.out}")


if __name__ == "__main__":
    main()
