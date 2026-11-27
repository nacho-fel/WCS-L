"""Pass 1: GEOPHIRES physics by depth and plant type (Cape-5 case, 12 producers, plant priced by GEOPHIRES' correlation).

Gives, for each reservoir depth, the produced temperature, gross and net MW per producer and the plant $/kW.
Usage: python geo_physics.py <GEOPHIRES-X>/tests/examples
"""
import json, sys
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "geophires"))
import sweep, depth  # noqa: E402

sweep.OUT = ROOT / "reports" / "depth_table" / "geophires_runs"
sweep.FIELDS.update(depth.sweep.FIELDS)
DEPTHS = [round(2.5 + 0.25 * i, 2) for i in range(17)]          # 2.5 ... 6.5 km
TYPES = {2: "ORC", 4: "DF"}                                       # supercritical ORC, double flash

if __name__ == "__main__":
    base = (Path(sys.argv[1]) / "Fervo_Project_Cape-5.txt").read_text()
    jobs = [(f"phys_d{d:.2f}_{n}", depth.case(base, d, t, False)) for d in DEPTHS for t, n in TYPES.items()]
    with ProcessPoolExecutor(2) as ex:
        res = dict(ex.map(sweep.run_one, jobs))
    out = ROOT / "reports" / "depth_table" / "physics.json"
    out.write_text(json.dumps(res, indent=1))
    for k, r in res.items():
        print(k, r["t_prod"], r["gross_mw"], r["net_mw"], r["plant"], r["lcoe"])
