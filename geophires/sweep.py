"""Where does WCS-L matter? GEOPHIRES runs across project sizes, plant costs and stimulation assumptions.

Template: GEOPHIRES' Fervo Cape case (Fervo_Project_Cape-5). For each project size (number of producers;
injectors = ceil(0.666 x producers), as GEOPHIRES does), the per-well drilling and completion cost is set to:
  flat   : the case's own NREL 2025 curve x 0.9 ($4.65M for every well);
  wcsl   : the WCS-L campaign average for that number of wells (P10, P50, P90 and expected);
and optionally the stimulation cost per well is scaled by a stimulation learning factor.
Everything else (reservoir physics, plant at $/kW, exploration $30M, finance) is the GEOPHIRES case.

Usage: python geophires/sweep.py --examples <GEOPHIRES-X>/tests/examples --what size|plant|stim|dc
"""
from __future__ import annotations

import argparse
import json
import math
import re
import subprocess
import sys
import warnings
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
warnings.simplefilter("ignore")
OUT = ROOT / "reports" / "geophires" / "sweep"
SUMMARY = ROOT / "reports" / "geophires" / "sweep_summary"
FIELDS = {
    "net_mw": r"Average Net Electricity Production:\s+([\d.]+) MW",
    "lcoe": r"Electricity breakeven price:\s+([\d.]+) cents/kWh",
    "npv": r"Project NPV:\s+(-?[\d.]+) MUSD",
    "irr": r"After-tax IRR:\s+(-?[\d.]+) %",
    "nprod": r"Number of production wells:\s+(\d+)",
    "ninj": r"Number of injection wells:\s+(\d+)",
    "expl": r"Exploration costs:\s+([\d.]+) MUSD",
    "drill": r"Drilling and completion costs:\s+([\d.]+) MUSD",
    "drill_per_well": r"Drilling and completion costs per well:\s+([\d.]+) MUSD",
    "stim": r"Stimulation costs:\s+([\d.]+) MUSD",
    "plant": r"Surface power plant costs:\s+([\d.]+) MUSD",
    "gath": r"Field gathering system costs:\s+([\d.]+) MUSD",
    "overnight": r"Overnight Capital Cost:\s+([\d.]+) MUSD",
    "capex": r"Total CAPEX:\s+([\d.]+) MUSD",
    "redrill_cost_yr": r"Redrilling costs:\s+([\d.]+) MUSD/yr",
    "om_total_yr": r"Total operating and maintenance costs:\s+([\d.]+) MUSD/yr",
    "redrills": r"Number of times redrilling:\s+(\d+)",
}


def edit(text: str, **kw) -> str:
    """Replace or add GEOPHIRES parameters (name -> value)."""
    lines = text.splitlines()
    for name, val in kw.items():
        lines = [l for l in lines if not l.startswith(name + ",")]
        lines.append(f"{name}, {val}")
    return "\n".join(lines) + "\n"


def run_one(args):
    name, text = args
    OUT.mkdir(parents=True, exist_ok=True)
    inp, out = OUT / f"{name}.txt", OUT / f"{name}.out"
    inp.write_text(text)
    subprocess.run([sys.executable, "-m", "geophires_x", str(inp), str(out)], check=True,
                   stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    s = out.read_text()
    r = {}
    for k, pat in FIELDS.items():
        m = re.search(pat, s)
        r[k] = float(m.group(1)) if m else None
    return name, r


def rigs_for(n_wells: int) -> int:
    return 1 if n_wells <= 10 else (2 if n_wells <= 30 else 3)


def wcsl_campaign(n_wells: int, months_gap_after: int | None = None, gap_months: float = 0.0,
                  n_draws: int = 3000) -> dict:
    """WCS-L campaign for n wells: rigs by size (<=10 wells: 1, <=30: 2, else 3), 1.5 wells per rig per
    month, 6 wells per pad, 14,100 ft measured depth (Cape-5 geometry). Optional drilling pause."""
    import wcsl
    from wcsl import wcs
    cal = wcsl.calibrate(theta=0.75)
    model = wcsl.WCSL(cal.params)
    rigs = rigs_for(n_wells)
    plans = []
    for i in range(n_wells):
        t = (i // rigs) / 1.5 / 12
        if months_gap_after is not None and i >= months_gap_after:
            t += gap_months / 12
        plans.append(wcsl.WellPlan(spec=wcs.cape_like_well(md_ft=14100), t_years=t, rig=f"rig-{i % rigs}",
                                   project="cape", pad=i // 6, name=f"w{i + 1}"))
    order = np.argsort([p.t_years for p in plans], kind="stable")
    plans = [plans[i] for i in order]
    ex = model.expected(plans)
    sim = model.simulate(plans, n_draws=n_draws, uncertainty=cal.uncertainty())
    pw = sim["total"] / n_wells
    return dict(n=n_wells, rigs=rigs, expected=ex["mean_cost_per_well"], p10=float(np.percentile(pw, 10)),
                p50=float(np.percentile(pw, 50)), p90=float(np.percentile(pw, 90)),
                costs=[w["cost"] for w in ex["wells"]], E_last=float(ex["wells"][-1]["E"]))


def stim_factor(n_wells: int, floor: float, lr_k: float) -> float:
    """Average stimulation cost over n wells relative to the first, with a Brett-Millheim decay in the
    well count: floor + (1 - floor) exp(-k (n - 1)); the case's $4M is treated as a mature-well value, so
    the factor is normalised to the floor (=1 for a fully learned well)."""
    n = np.arange(1, n_wells + 1)
    f = floor + (1 - floor) * np.exp(-lr_k * (n - 1))
    return float(f.mean() / floor)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--examples", required=True)
    ap.add_argument("--what", default="size")
    a = ap.parse_args()
    base = (Path(a.examples) / "Fervo_Project_Cape-5.txt").read_text()
    producers = [1, 2, 3, 5, 8, 12, 20, 35, 56]
    jobs, meta = [], {}
    if a.what == "size":
        for p in producers:
            n = p + math.ceil(0.666 * p)
            c = wcsl_campaign(n)
            meta[str(p)] = dict(n_wells=n, campaign=c)
            t = edit(base, **{"Number of Production Wells": p})
            jobs.append((f"size_p{p}_flat", t))
            for q in ("p10", "p50", "p90"):
                jobs.append((f"size_p{p}_wcsl_{q}", edit(t, **{"Well Drilling and Completion Capital Cost": f"{c[q] / 1e6:.4f}"})))
    elif a.what == "stim":
        # stimulation learning scenarios on top of WCS-L P50 drilling
        scen = json.loads((ROOT / "geophires" / "stim_scenarios.json").read_text())
        size = json.loads((SUMMARY / "size_meta.json").read_text())
        for p in producers:
            n = size[str(p)]["n_wells"]
            c = size[str(p)]["campaign"]
            t = edit(base, **{"Number of Production Wells": p, "Well Drilling and Completion Capital Cost": f"{c['p50'] / 1e6:.4f}"})
            for s in scen:
                f = stim_factor(n, s["floor"], s["k"])
                meta[f"{p}_{s['id']}"] = dict(factor=f)
                jobs.append((f"stim_p{p}_{s['id']}", edit(t, **{"Reservoir Stimulation Capital Cost per Injection Well": f"{4 * f:.4f}",
                                                                 "Reservoir Stimulation Capital Cost per Production Well": f"{4 * f:.4f}"})))
    elif a.what == "plant":
        size = json.loads((SUMMARY / "size_meta.json").read_text())
        for p in (3, 12, 56):
            c = size[str(p)]["campaign"]
            for kw in (1400, 1900, 2500, 3200):
                t = edit(base, **{"Number of Production Wells": p, "Capital Cost for Power Plant for Electricity Generation": kw})
                jobs.append((f"plant_p{p}_{kw}_flat", t))
                jobs.append((f"plant_p{p}_{kw}_wcsl", edit(t, **{"Well Drilling and Completion Capital Cost": f"{c['p50'] / 1e6:.4f}"})))
    elif a.what == "dc":
        # data-centre preview: 12 producers (~100 MW net); one continuous campaign vs a 12-month pause after 10 wells
        p = 12
        n = p + math.ceil(0.666 * p)
        for label, gap in (("continuous", None), ("pause12", 10)):
            c = wcsl_campaign(n, months_gap_after=gap, gap_months=12 if gap else 0)
            meta[label] = c
            jobs.append((f"dc_{label}", edit(base, **{"Number of Production Wells": p,
                                                       "Well Drilling and Completion Capital Cost": f"{c['p50'] / 1e6:.4f}"})))
    with ProcessPoolExecutor(max_workers=2) as ex:
        res = dict(ex.map(run_one, jobs))
    SUMMARY.mkdir(parents=True, exist_ok=True)
    (SUMMARY / f"{a.what}_meta.json").write_text(json.dumps(meta, indent=1))
    (SUMMARY / f"{a.what}_results.json").write_text(json.dumps(res, indent=1))
    print(f"{a.what}: {len(res)} runs")


if __name__ == "__main__":
    main()
