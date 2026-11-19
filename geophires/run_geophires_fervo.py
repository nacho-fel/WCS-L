"""GEOPHIRES cross-check: Fervo's Cape and Project Red cases, with WCS-L well costs plugged in.

1. Reproduce GEOPHIRES' own Fervo Cape Station case study (Fervo_Project_Cape-5: 500 MW, 56 producers,
   38 injectors, flat $4.65M per well from the NREL 2025 curves).
2. Run GEOPHIRES' Project Red case (Fervo_Norbeck_Latimer_2023) and compare with Fervo's reported
   3.5 MW at 63 L/s and 191 C.
3. Replace the flat well cost with (a) the observed average of Fervo's six Cape wells, (b) the WCS with
   NREL 2025 parameters, (c) WCS-L calibrated to the six wells, P10 / P50 / P90 of the campaign mean.

Needs GEOPHIRES-X installed from https://github.com/NREL/GEOPHIRES-X (pip install .) and the path to its
tests/examples folder: python geophires/run_geophires_fervo.py --examples PATH
"""
from __future__ import annotations

import argparse
import json
import re
import subprocess
import sys
import tempfile
import warnings
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
import wcsl  # noqa: E402
from wcsl import wcs  # noqa: E402

warnings.simplefilter("ignore")
OUT = ROOT / "reports" / "geophires"

FIELDS = {
    "net_mw": r"Average Net Electricity Production:\s+([\d.]+) MW",
    "lcoe_cents": r"Electricity breakeven price:\s+([\d.]+) cents/kWh",
    "npv_musd": r"Project NPV:\s+(-?[\d.]+) MUSD",
    "irr_pct": r"After-tax IRR:\s+(-?[\d.]+) %",
    "producers": r"Number of production wells:\s+(\d+)",
    "injectors": r"Number of injection wells:\s+(\d+)",
    "drilling_total_musd": r"Drilling and completion costs:\s+([\d.]+) MUSD",
    "drilling_per_well_musd": r"Drilling and completion costs per well:\s+([\d.]+) MUSD",
    "capex_musd": r"Total CAPEX:\s+([\d.]+) MUSD",
    "redrills": r"Number of times redrilling:\s+(\d+)",
    "gross_mw": r"Average Total Electricity Generation:\s+([\d.]+) MW",
    "max_net_mw": r"Maximum Net Electricity Generation:\s+([\d.]+) MW",
    "prod_temp_c": r"Average Production Temperature:\s+([\d.]+) degC",
}


def run(text: str, name: str) -> dict:
    OUT.mkdir(parents=True, exist_ok=True)
    inp = OUT / f"{name}.txt"
    out = OUT / f"{name}.out"
    inp.write_text(text)
    subprocess.run([sys.executable, "-m", "geophires_x", str(inp), str(out)], check=True,
                   stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    s = out.read_text()
    res = {}
    for k, pat in FIELDS.items():
        m = re.search(pat, s)
        res[k] = float(m.group(1)) if m else None
    return res


def with_well_cost(text: str, musd: float) -> str:
    """Override the correlation with a fixed per-well drilling and completion cost."""
    lines = [l for l in text.splitlines() if not l.startswith("Well Drilling and Completion Capital Cost")]
    lines.append(f"Well Drilling and Completion Capital Cost, {musd:.4f}, -- WCS-L cross-check (per well, MUSD)")
    return "\n".join(lines) + "\n"


def wcsl_campaign(n_wells: int, rigs: int = 3, wells_per_rig_month: float = 1.5, wells_per_pad: int = 6,
                  md_ft: float = 14100, n_draws: int = 4000) -> dict:
    """WCS-L cost of an n-well campaign drilled by several rigs in parallel (Fervo: 3 rigs, about 1.5 wells
    per rig per month). Wells are assigned to rigs in turn; experience from another rig counts omega = 0.4."""
    cal = wcsl.calibrate(theta=0.75)
    model = wcsl.WCSL(cal.params)
    plans = []
    for i in range(n_wells):
        rig = i % rigs
        k = i // rigs
        plans.append(wcsl.WellPlan(spec=wcs.cape_like_well(md_ft=md_ft), t_years=k / wells_per_rig_month / 12,
                                   rig=f"rig-{rig}", project="cape", pad=i // wells_per_pad, name=f"w{i + 1}"))
    order = np.argsort([p.t_years for p in plans], kind="stable")
    plans = [plans[i] for i in order]
    ex = model.expected(plans)
    sim = model.simulate(plans, n_draws=n_draws, uncertainty=cal.uncertainty())
    per_well_mean = sim["total"] / n_wells
    return dict(expected_mean=ex["mean_cost_per_well"], p10=float(np.percentile(per_well_mean, 10)),
                p50=float(np.percentile(per_well_mean, 50)), p90=float(np.percentile(per_well_mean, 90)),
                first=ex["wells"][0]["cost"], last=ex["wells"][-1]["cost"], costs=[w["cost"] for w in ex["wells"]])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--examples", required=True, help="GEOPHIRES-X tests/examples folder")
    a = ap.parse_args()
    ex = Path(a.examples)
    cape = (ex / "Fervo_Project_Cape-5.txt").read_text()
    red = (ex / "Fervo_Norbeck_Latimer_2023.txt").read_text()

    results = {"cape_base": run(cape, "cape5_base")}
    reported = {line.split(",")[0]: line for line in cape.splitlines()}
    # Project Red: GEOPHIRES case (average flow 41 kg/s, 650 gpm) and the 63 L/s peak (about 60 kg/s)
    results["red_41kgs"] = run(red, "red_41kgs")
    results["red_60kgs"] = run(red.replace("Production Flow Rate per Well, 41", "Production Flow Rate per Well, 60"), "red_60kgs")

    # Well costs to plug in
    n = int(results["cape_base"]["producers"] + results["cape_base"]["injectors"])
    fervo = [w["cost"] for w in wcsl.load_fervo()]
    wcs_2025 = wcs.well_cost(wcs.cape_like_well(md_ft=14100))["total"]
    camp = wcsl_campaign(n)
    costs = {
        "NREL 2025 curve x 0.9 (GEOPHIRES case)": results["cape_base"]["drilling_per_well_musd"] * 1e6,
        "WCS, NREL 2025 parameters (no learning)": wcs_2025,
        "Fervo's six Cape wells, observed average": float(np.mean(fervo)),
        "Fervo's first Cape well": fervo[0],
        f"WCS-L, {n}-well campaign, P10": camp["p10"],
        f"WCS-L, {n}-well campaign, P50": camp["p50"],
        f"WCS-L, {n}-well campaign, P90": camp["p90"],
    }
    for name, c in costs.items():
        if name.startswith("NREL"):
            results[name] = results["cape_base"]
        else:
            results[name] = run(with_well_cost(cape, c / 1e6), re.sub(r"[^a-z0-9]+", "_", name.lower()).strip("_"))
        results[name]["well_cost_input"] = c
    results["wcsl_campaign"] = camp
    (OUT / "results.json").write_text(json.dumps(results, indent=2))

    b = results["cape_base"]
    L = ["# GEOPHIRES cross-check (generated by geophires/run_geophires_fervo.py)", "",
         "## 1. GEOPHIRES' Fervo Cape case reproduced", "",
         f"- Net generation {b['net_mw']:.1f} MW from {b['producers']:.0f} producers ({b['net_mw'] / b['producers']:.1f} MW per producer) and {b['injectors']:.0f} injectors.",
         f"- LCOE {b['lcoe_cents']:.2f} cents/kWh, NPV {b['npv_musd']:.0f} MUSD, drilling {b['drilling_per_well_musd']:.2f} MUSD per well ({b['drilling_total_musd']:.0f} MUSD), {b['redrills']:.0f} redrilling events.", "",
         "## 2. Project Red (Fervo: 3.5 MW at 63 L/s, 191 C, 30-day test)", "",
         "| GEOPHIRES run | Flow | Gross MW (average) | Net MW (average) | Average production temperature |", "|---|---|---|---|---|"]
    for k, lab in (("red_41kgs", "41 kg/s (650 gpm, the case file)"), ("red_60kgs", "60 kg/s (about 63 L/s)")):
        r = results[k]
        L.append(f"| {lab} | | {r['gross_mw']:.2f} | {r['net_mw']:.2f} | {r['prod_temp_c']:.0f} C |")
    L += ["", "## 3. Cape case with different well costs", "",
          "| Well cost used | $M per well | Drilling total, MUSD | LCOE, cents/kWh | NPV, MUSD | IRR |", "|---|---|---|---|---|---|"]
    for name in costs:
        r = results[name]
        L.append(f"| {name} | {r['well_cost_input'] / 1e6:.2f} | {r['drilling_total_musd']:.0f} | {r['lcoe_cents']:.2f} | {r['npv_musd']:.0f} | {r['irr_pct']:.1f} % |")
    L += ["", f"WCS-L campaign: {n} wells, 3 rigs at 1.5 wells per rig per month, 6 wells per pad, 14,100 ft measured depth. "
          f"First well ${camp['first'] / 1e6:.2f}M, last ${camp['last'] / 1e6:.2f}M, expected mean ${camp['expected_mean'] / 1e6:.2f}M."]
    (OUT / "results.md").write_text("\n".join(L) + "\n")
    print("\n".join(L))


if __name__ == "__main__":
    main()
