"""Extend the GEOPHIRES grid to 18-26 producers (rule plant, base ROP), so the 100 MW path at shallow depth
(26 producers at 2 km, 17 at 2.25 km) lies on the capex/LCOE surfaces. Merges into grid.json."""
from __future__ import annotations
import json, sys
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "geophires")); sys.path.insert(0, str(Path(__file__).parent))
import sweep, depth  # noqa: E402
import geo_grid as gg  # noqa: E402

EXTRA = [18, 20, 22, 24, 26]

def main():
    base = (Path(sys.argv[1]) / "Fervo_Project_Cape-5.txt").read_text()
    g = json.loads((gg.OUT / "grid.json").read_text())
    dc = {(r["profile"], r["depth_km"], r["n"]): r for r in json.loads((gg.OUT / "drilling_table.json").read_text())}
    jobs = []
    for d in gg.DEPTHS:
        rule = g["meta"][str(d)]["rule"]
        for p in EXTRA:
            key = f"grid_d{d:.2f}_{rule}_p{p}_base"
            if key in g["runs"]:
                continue
            c = dc[("base", d, gg.n_wells(p))]["p50"]
            t = depth.case(base, d, gg.PT[rule], False)
            t = sweep.edit(t, **{"Number of Production Wells": p, "Well Drilling and Completion Capital Cost": round(c / 1e6, 4)})
            jobs.append((key, t))
    print(len(jobs), "runs", flush=True)
    with ProcessPoolExecutor(2) as ex:
        g["runs"].update(dict(ex.map(sweep.run_one, jobs)))
    (gg.OUT / "grid.json").write_text(json.dumps(g, indent=1))
    print("done")

if __name__ == "__main__":
    main()
