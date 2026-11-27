"""WCS-L drilling cost by reservoir depth and number of wells, with ROP falling with depth.

For each ROP profile (k = 0, 0.30, 0.609 base, 0.656), WCS-L is recalibrated to Fervo's Cape wells
(calibrate_cape17, so Cape is reproduced whatever k is), then a campaign of N identical horizontal wells
(TVD = reservoir depth, 5,000 ft lateral, Cape casing design) is priced:
rigs: 1 up to 10 wells, 2 up to 30, 3 above; 1.5 wells per rig per month; 6 wells per pad.
Writes reports/depth_table/drilling_table.json (+ .csv).
"""
from __future__ import annotations
import csv, json, sys
from pathlib import Path
import numpy as np
ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))
import wcsl  # noqa: E402
from wcsl import wcs, rop  # noqa: E402

OUT = ROOT / "reports" / "depth_table"
DEPTHS = [round(2.0 + 0.25 * i, 2) for i in range(19)]          # 2.0 ... 6.5 km
NS = list(range(1, 61))
PROFILES = {"k0": rop.CONSTANT, "k030": rop.HALF, "base": rop.BASE, "k0656": rop.PABLO_ALLIN}
LATERAL = 5000.0


def rigs_for(n):
    return 1 if n <= 10 else (2 if n <= 30 else 3)


def campaign(model, cal, tvd_ft, n, profile, n_draws=2000):
    rigs = rigs_for(n)
    plans = [wcsl.WellPlan(spec=wcs.egs_well(tvd_ft, LATERAL, profile), t_years=(i // rigs) / 1.5 / 12,
                           rig=f"r{i % rigs}", project="p", pad=i // 6, name=f"w{i + 1}") for i in range(n)]
    ex = model.expected(plans)
    sim = model.simulate(plans, n_draws=n_draws, uncertainty=cal.uncertainty(), seed=11)
    pw = sim["total"] / n
    w = ex["wells"]
    return dict(n=n, rigs=rigs, total=ex["total"], mean=ex["mean_cost_per_well"],
                p10=float(np.percentile(pw, 10)), p50=float(np.percentile(pw, 50)), p90=float(np.percentile(pw, 90)),
                first=w[0]["cost"], last=w[-1]["cost"], first_days=w[0]["days"], last_days=w[-1]["days"],
                ts_first=w[0]["time_share"], ts_last=w[-1]["time_share"], md=w[0]["md_ft"],
                no_learning_total=n * w[0]["cost"] + sum(ex["campaign_fixed"].values()))


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    rows = []
    for name, prof in PROFILES.items():
        cal = wcsl.calibrate_cape17(spec_builder=lambda md, p=prof: wcs.cape_like_well(md_ft=md, rop_profile=p))
        model = wcsl.WCSL(cal.params)
        for d in DEPTHS:
            tvd = d * rop.FT_PER_KM
            for n in NS:
                r = campaign(model, cal, tvd, n, prof, n_draws=2000 if name == "base" else 500)
                r.update(profile=name, k_per_km=prof.k_per_km, depth_km=d, rop_at_depth=float(prof.rop(tvd)))
                rows.append(r)
            print(name, d, round(rows[-40]["first"] / 1e6, 2), round(rows[-1]["mean"] / 1e6, 2), flush=True)
    (OUT / "drilling_table.json").write_text(json.dumps(rows))
    keys = ["profile", "k_per_km", "depth_km", "n", "rigs", "rop_at_depth", "md", "first", "last", "mean", "p10", "p50",
            "p90", "total", "no_learning_total", "first_days", "last_days", "ts_first", "ts_last"]
    with open(OUT / "drilling_table.csv", "w", newline="") as fh:
        wr = csv.DictWriter(fh, fieldnames=keys, extrasaction="ignore"); wr.writeheader(); wr.writerows(rows)


if __name__ == "__main__":
    main()
