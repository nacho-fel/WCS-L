"""Assemble the table of points for the MILP and the data for the figures.

Reads physics.json, drilling_table.json and grid.json (reports/depth_table) and writes:
  table_points.csv   one row per (depth, producers): temperature, plant, MW, wells, drilling P10/P50/P90,
                     plant, stimulation, total capex, LCOE (GEOPHIRES with WCS-L drilling cost)
  point100.csv       the 100 MW point of every depth, rule plant and the other plant, ROP sensitivities
  figdata.json       everything the page draws
"""
from __future__ import annotations
import csv, json, math, sys
from pathlib import Path
import numpy as np
ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))
from wcsl import rop  # noqa: E402
OUT = ROOT / "reports" / "depth_table"


def n_wells(p):
    return p + math.ceil(0.666 * p)


def main():
    phys = json.loads((OUT / "physics.json").read_text())
    drill = json.loads((OUT / "drilling_table.json").read_text())
    g = json.loads((OUT / "grid.json").read_text())
    meta, runs, tsw = g["meta"], g["runs"], g["t_switch"]
    dc = {(r["profile"], r["depth_km"], r["n"]): r for r in drill}
    depths = sorted(float(d) for d in meta)
    rows = []
    for k, r in runs.items():
        _, d, pl, p, prof = k.split("_")
        d, p = float(d[1:]), int(p[1:])
        N = n_wells(p)
        D = dc[(prof, d, N)]
        rows.append(dict(depth_km=d, t_prod_C=meta[str(d)]["t_prod"], plant=pl, rule_plant=meta[str(d)]["rule"],
                         rop_profile=prof, producers=p, wells=N, net_mw=r["net_mw"], gross_mw=r["gross_mw"],
                         drill_per_well_p10=D["p10"] / 1e6, drill_per_well_p50=D["p50"] / 1e6, drill_per_well_p90=D["p90"] / 1e6,
                         first_well=D["first"] / 1e6, last_well=D["last"] / 1e6, drilling_total=r["drill"], stimulation=r["stim"],
                         plant_cost=r["plant"], plant_usd_per_kw_net=r["plant"] / r["net_mw"] * 1000 if r["net_mw"] else None,
                         gathering=r["gath"], exploration=r["expl"], overnight_capex=r["overnight"], lcoe=r["lcoe"], npv=r["npv"],
                         capex_per_kw_net=r["overnight"] / r["net_mw"] * 1000 if r["net_mw"] else None))
    rows.sort(key=lambda x: (x["rop_profile"] != "base", x["depth_km"], x["plant"], x["producers"]))
    keys = list(rows[0])
    with open(OUT / "table_points.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=keys); w.writeheader(); w.writerows(rows)
    # 100 MW points
    pts = []
    for d in depths:
        m = meta[str(d)]
        for pl in ("ORC", "DF"):
            p = m["p100"][pl]
            for prof in ("base", "k0", "k030", "k0656"):
                hit = [r for r in rows if r["depth_km"] == d and r["plant"] == pl and r["producers"] == p and r["rop_profile"] == prof]
                if hit:
                    pts.append(dict(hit[0], is_rule=(pl == m["rule"])))
    with open(OUT / "point100.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(pts[0])); w.writeheader(); w.writerows(pts)
    base = [r for r in drill if r["profile"] == "base"]
    Ns = sorted({r["n"] for r in base})
    fig = dict(
        t_switch=tsw, depths=depths,
        temp=[meta[str(d)]["t_prod"] for d in depths],
        rule=[meta[str(d)]["rule"] for d in depths],
        mwpp={pl: [meta[str(d)]["mw_per_producer"][pl] for d in depths] for pl in ("ORC", "DF")},
        plant_kw_gross={pl: [phys[f"phys_d{d:.2f}_{pl}"]["plant"] / phys[f"phys_d{d:.2f}_{pl}"]["gross_mw"] * 1000 for d in depths] for pl in ("ORC", "DF")},
        Ns=Ns,
        drill_surface=[[dc[("base", d, n)]["p50"] / 1e6 for n in Ns] for d in depths],
        drill_first={prof: [dc[(prof, d, 1)]["first"] / 1e6 for d in depths] for prof in ("k0", "k030", "base", "k0656")},
        drill_mean20={prof: [dc[(prof, d, 20)]["mean"] / 1e6 for d in depths] for prof in ("k0", "k030", "base", "k0656")},
        ts_first=[dc[("base", d, 20)]["ts_first"] for d in depths], ts_last=[dc[("base", d, 20)]["ts_last"] for d in depths],
        ts_first_k0=[dc[("k0", d, 20)]["ts_first"] for d in depths], ts_last_k0=[dc[("k0", d, 20)]["ts_last"] for d in depths],
        saving20=[1 - dc[("base", d, 20)]["total"] / dc[("base", d, 20)]["no_learning_total"] for d in depths],
        saving20_k0=[1 - dc[("k0", d, 20)]["total"] / dc[("k0", d, 20)]["no_learning_total"] for d in depths],
        saving20_usd=[(dc[("base", d, 20)]["no_learning_total"] - dc[("base", d, 20)]["total"]) / 1e6 for d in depths],
        first_days=[dc[("base", d, 1)]["first_days"] for d in depths],
        rop={name: [float(pr.rop(d * rop.FT_PER_KM)) for d in depths] for name, pr in
             dict(k0=rop.CONSTANT, k030=rop.HALF, base=rop.BASE, k0656=rop.PABLO_ALLIN).items()},
        grid=[r for r in rows if r["rop_profile"] == "base"],
        p100=pts,
        fervo_depth_days=[list(map(float, rop.FERVO_DEPTH_DAYS[0])), list(map(float, rop.FERVO_DEPTH_DAYS[1]))],
        fit_allin=rop.fit_allin(), fit_onbottom=rop.fit_onbottom(),
        band=json.loads((ROOT / "reports" / "band" / "band.json").read_text()),
    )
    (OUT / "figdata.json").write_text(json.dumps(fig))
    print("rows", len(rows), "points100", len(pts))
    for p in pts:
        if p["rop_profile"] == "base" and p["is_rule"]:
            print(p["depth_km"], round(p["t_prod_C"]), p["plant"], p["producers"], p["wells"], round(p["net_mw"], 1),
                  round(p["drill_per_well_p50"], 2), round(p["drilling_total"], 1), round(p["plant_cost"], 1), round(p["overnight_capex"], 1), p["lcoe"])


if __name__ == "__main__":
    main()
