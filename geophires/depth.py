"""Deeper wells: GEOPHIRES' Cape case at greater reservoir depths, with each surface plant type.

Fixed: 12 producers (+ 8 injectors), the data-centre size. Varied: reservoir depth (km), Power Plant Type,
and how the plant is priced: the case's fixed $1,900/kW, or GEOPHIRES' built-in correlation for that type.
Drilling here is still GEOPHIRES' flat NREL 2025 curve x 0.9 (it rises with depth); WCS-L replaces it later.

Usage: python geophires/depth.py --examples <GEOPHIRES-X>/tests/examples
"""
from __future__ import annotations
import argparse, json, re, sys
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import sweep  # noqa: E402

sweep.OUT = sweep.ROOT / "reports" / "geophires" / "depth"
sweep.FIELDS.update({
    "t_prod": r"Average Production Temperature:\s+([\d.]+) degC",
    "t_bh": r"Bottom-hole temperature:\s+([\d.]+) degC",
    "gross_mw": r"Average Total Electricity Generation:\s+([\d.]+) MW",
    "pump_mw": r"Average Pumping Power:\s+([\d.]+) MW",
    "eff": r"First law efficiency \(average\):\s+([\d.]+)",
})
PLANT_KEY = "Capital Cost for Power Plant for Electricity Generation"
DEPTHS = [2.68, 3.0, 3.5, 4.0, 4.5]
CONFIGS = {  # name: (plant type, keep $1,900/kW override?)
    "sORC_1900": (2, True), "sORC_corr": (2, False), "sflash_corr": (3, False), "dflash_corr": (4, False),
}


def case(base: str, depth: float, ptype: int, fixed: bool) -> str:
    t = sweep.edit(base, **{"Reservoir Depth": depth, "Power Plant Type": ptype, "Number of Production Wells": 12})
    if not fixed:
        t = "\n".join(l for l in t.splitlines() if not l.startswith(PLANT_KEY + ",")) + "\n"
    return t


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--examples", required=True); ap.add_argument("--workers", type=int, default=4)
    a = ap.parse_args()
    base = (Path(a.examples) / "Fervo_Project_Cape-5.txt").read_text()
    jobs = [(f"d{d:.2f}_{c}", case(base, d, p, f)) for d in DEPTHS for c, (p, f) in CONFIGS.items()]
    with ProcessPoolExecutor(a.workers) as ex:
        res = dict(ex.map(sweep.run_one, jobs))
    sweep.SUMMARY.mkdir(parents=True, exist_ok=True)
    (sweep.SUMMARY / "depth_results.json").write_text(json.dumps(res, indent=1))
    for k, r in res.items():
        print(k, {x: r.get(x) for x in ["t_prod", "gross_mw", "net_mw", "pump_mw", "drill_per_well", "plant", "lcoe", "npv"]})


if __name__ == "__main__":
    main()
