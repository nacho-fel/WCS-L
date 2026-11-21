"""Deeper wells for a 100 MW data centre: wells needed, WCS-L drilling cost and LCOE by depth.

For each reservoir depth: the fewest producers giving >= 100 MW net (from depth.py's 12-producer runs),
injectors = ceil(0.666 x producers); drilling at the flat NREL curve (GEOPHIRES) or the WCS-L campaign
average (P50) for that many wells, priced for the deeper geometry (vertical depth + lateral; low, central and
high ROP); plant = supercritical ORC at $1,900/kW or at GEOPHIRES' correlation, or double flash at GEOPHIRES' correlation.
"""
from __future__ import annotations
import json, math, sys, warnings
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parent)); sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
warnings.simplefilter("ignore")
import sweep, depth  # noqa: E402

sweep.OUT = sweep.ROOT / "reports" / "geophires" / "depth"
FT = 1 / 0.3048
GEOM = {2.68: (8700, None), 3.0: (3000 * FT, 5000), 3.5: (3500 * FT, 5000), 4.0: (4000 * FT, 2000 * FT), 4.5: (4500 * FT, 2000 * FT)}
ROP = {"low": (50, 35), "central": (75, 50), "high": (100, 60)}
PLANTS = {"sORC_1900": (2, True), "sORC_corr": (2, False), "dflash_corr": (4, False)}
TARGET = 100.0


def wcsl_cost(depth_km, n_wells, rop_case="central", n_draws=3000):
    import wcsl
    from wcsl import wcs
    cal = wcsl.calibrate(theta=0.75); m = wcsl.WCSL(cal.params)
    tvd, lat = GEOM[depth_km]; rop, bit = ROP[rop_case]
    if lat is None:
        lat = wcs.lateral_for_md(14100)
    rigs = sweep.rigs_for(n_wells)
    plans = [wcsl.WellPlan(spec=wcs.horizontal_well(tvd, lat, 10, rop, bit, **wcs.CAPE_LIKE_PARAMS),
                           t_years=(i // rigs) / 1.5 / 12, rig=f"r{i % rigs}", project="p", pad=i // 6, name=f"w{i+1}")
             for i in range(n_wells)]
    plans = [plans[j] for j in np.argsort([p.t_years for p in plans], kind="stable")]
    ex = m.expected(plans); sim = m.simulate(plans, n_draws=n_draws, uncertainty=cal.uncertainty())
    pw = sim["total"] / n_wells
    return dict(rigs=rigs, p10=float(np.percentile(pw, 10)), p50=float(np.percentile(pw, 50)), p90=float(np.percentile(pw, 90)),
                first=ex["wells"][0]["cost"], first_days=ex["wells"][0]["days"], last=ex["wells"][-1]["cost"],
                last_days=ex["wells"][-1]["days"], md=ex["wells"][0]["md_ft"], path=[w["cost"] for w in ex["wells"]])


def main():
    ex_dir = Path(sys.argv[1])
    base = (ex_dir / "Fervo_Project_Cape-5.txt").read_text()
    d12 = json.loads((sweep.SUMMARY / "depth_results.json").read_text())
    jobs, meta = [], {}
    for d in GEOM:
        for pl, (ptype, fixed) in PLANTS.items():
            per = d12[f"d{d:.2f}_{pl}"]["net_mw"] / 12
            p = math.ceil(TARGET / per - 1e-9); n = p + math.ceil(0.666 * p)
            costs = {rc: wcsl_cost(d, n, rc) for rc in ROP}
            meta[f"{d}|{pl}"] = dict(depth=d, plant=pl, producers=p, wells=n, mw_per_producer=per, wcsl=costs)
            t = depth.case(base, d, ptype, fixed)
            t = sweep.edit(t, **{"Number of Production Wells": p})
            jobs.append((f"dc_d{d:.2f}_{pl}_flat", t))
            for rc in ROP:
                jobs.append((f"dc_d{d:.2f}_{pl}_wcsl_{rc}", sweep.edit(t, **{"Well Drilling and Completion Capital Cost": round(costs[rc]["p50"] / 1e6, 3)})))
    with ProcessPoolExecutor(2) as ex:
        res = dict(ex.map(sweep.run_one, jobs))
    (sweep.SUMMARY / "depth_dc_results.json").write_text(json.dumps(dict(meta=meta, runs=res), indent=1))
    for k, v in meta.items():
        d, pl = v["depth"], v["plant"]
        f = res[f"dc_d{d:.2f}_{pl}_flat"]
        w = {rc: res[f"dc_d{d:.2f}_{pl}_wcsl_{rc}"] for rc in ROP}
        print(f"{d:4.2f} {pl:12s} p={v['producers']:2d} wells={v['wells']:2d} net={f['net_mw']:6.1f}  flat/well {f['drill_per_well']:5.2f} "
              f"wcsl/well {v['wcsl']['central']['p50']/1e6:5.2f}  LCOE flat {f['lcoe']:5.2f} wcsl L/C/H {w['low']['lcoe']:5.2f}/{w['central']['lcoe']:5.2f}/{w['high']['lcoe']:5.2f}"
              f"  NPV flat {f['npv']:6.1f} wcsl-c {w['central']['npv']:6.1f}")


if __name__ == "__main__":
    main()
