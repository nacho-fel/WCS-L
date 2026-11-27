"""Pass 2: GEOPHIRES with WCS-L drilling costs over the (depth, producers) grid, plant chosen by temperature.

Plant rule: supercritical ORC below T_SWITCH (produced temperature), double flash at or above it, both
priced with GEOPHIRES' built-in correlations (Cape's $1,900/kW override removed). Water is assumed available.
Per-well drilling and completion cost = WCS-L P50 campaign average for N = p + ceil(0.666 p) wells
(base ROP profile). Also: the 100 MW point of every depth with the other plant, and with the other ROP profiles.
Usage: python geo_grid.py <GEOPHIRES-X>/tests/examples
"""
from __future__ import annotations
import json, math, sys
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "geophires"))
import sweep, depth  # noqa: E402

sweep.OUT = ROOT / "reports" / "depth_table" / "geophires_runs"
sweep.FIELDS.update(depth.sweep.FIELDS)
OUT = ROOT / "reports" / "depth_table"
DEPTHS = [round(2.0 + 0.25 * i, 2) for i in range(19)]
PRODUCERS = [2, 3, 4, 5, 6, 8, 10, 12, 14, 16]
T_SWITCH = 210.0
TARGET = 100.0
PT = {"ORC": 2, "DF": 4}


def n_wells(p):
    return p + math.ceil(0.666 * p)


def main():
    base = (Path(sys.argv[1]) / "Fervo_Project_Cape-5.txt").read_text()
    phys = json.loads((OUT / "physics.json").read_text())
    missing = [(f"phys_d{d:.2f}_{n}", depth.case(base, d, t, False)) for d in DEPTHS for n, t in PT.items()
               if f"phys_d{d:.2f}_{n}" not in phys]
    with ProcessPoolExecutor(2) as ex:
        phys.update(dict(ex.map(sweep.run_one, missing)))
    (OUT / "physics.json").write_text(json.dumps(phys, indent=1))
    drill = json.loads((OUT / "drilling_table.json").read_text())
    dc = {(r["profile"], r["depth_km"], r["n"]): r for r in drill}
    jobs, meta = [], {}
    for d in DEPTHS:
        T = phys[f"phys_d{d:.2f}_DF"]["t_prod"]
        rule = "ORC" if T < T_SWITCH else "DF"
        plants = ["ORC", "DF"] if d <= 3.5 else [rule]
        mw = {pl: phys[f"phys_d{d:.2f}_{pl}"]["net_mw"] / 12 for pl in PT}
        p100 = {pl: math.ceil(TARGET / mw[pl] - 1e-9) for pl in PT}
        meta[str(d)] = dict(depth_km=d, t_prod=T, rule=rule, mw_per_producer=mw, p100=p100)
        todo = {(pl, p, "base") for pl in plants for p in PRODUCERS}
        todo |= {(pl, p100[pl], "base") for pl in PT}
        todo |= {(rule, p100[rule], prof) for prof in ("k0", "k030", "k0656")}
        for pl, p, prof in sorted(todo):
            c = dc[(prof, d, n_wells(p))]["p50"]
            t = depth.case(base, d, PT[pl], False)
            t = sweep.edit(t, **{"Number of Production Wells": p, "Well Drilling and Completion Capital Cost": round(c / 1e6, 4)})
            jobs.append((f"grid_d{d:.2f}_{pl}_p{p}_{prof}", t))
    print(len(jobs), "runs", flush=True)
    with ProcessPoolExecutor(2) as ex:
        res = dict(ex.map(sweep.run_one, jobs))
    (OUT / "grid.json").write_text(json.dumps(dict(meta=meta, runs=res, t_switch=T_SWITCH), indent=1))
    print("done")


if __name__ == "__main__":
    main()
