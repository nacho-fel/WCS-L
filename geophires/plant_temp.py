"""Plant cost and output against produced temperature, for each Power Plant Type (GEOPHIRES' own cost correlations).

Cape case, 12 producers, reservoir depth 2.0-5.0 km (the temperature comes from the case's gradients). The plant
cost override ($1,900/kW) is removed, so GEOPHIRES prices each type with its built-in correlation.
"""
from __future__ import annotations
import json, sys
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import sweep, depth  # noqa: E402

sweep.OUT = sweep.ROOT / "reports" / "geophires" / "depth"
sweep.FIELDS.update({"conv_eff": r"Heat to Power Conversion Efficiency:\s+([\d.]+) %"})
DEPTHS = [2.0, 2.3, 2.68, 3.0, 3.5, 4.0, 4.5, 5.0]
TYPES = {1: "subORC", 2: "supORC", 3: "sflash", 4: "dflash"}

if __name__ == "__main__":
    base = (Path(sys.argv[1]) / "Fervo_Project_Cape-5.txt").read_text()
    jobs = [(f"pt_d{d:.2f}_{n}", depth.case(base, d, t, False)) for d in DEPTHS for t, n in TYPES.items()]
    with ProcessPoolExecutor(2) as ex:
        res = dict(ex.map(sweep.run_one, jobs))
    (sweep.SUMMARY / "plant_temp_results.json").write_text(json.dumps(res, indent=1))
    for k, r in res.items():
        g = r["gross_mw"] or 0
        print(f"{k:22s} T={r['t_prod']} gross={g} net={r['net_mw']} plant={r['plant']} $/kWgross={r['plant']/g*1000 if g else None:.0f} "
              f"$/kWnet={r['plant']/r['net_mw']*1000 if r['net_mw'] else 0:.0f} eff={r['conv_eff']} lcoe={r['lcoe']}")
