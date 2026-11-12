"""Calibrate WCS-L to Fervo's six Cape Station wells and draw the calibration figures.

Writes reports/calibration/{calibration.json, learning_curves.png, spread_rate.png,
cost_path.png, cost_structure.png}.
"""
from __future__ import annotations

import warnings

import numpy as np

from _common import AMBER, GREY, INK, MUTED, ROCK, ROOT, STEEL, WARN, dump, plt, save
import wcsl
from wcsl import validation as V
from wcsl import wcs

warnings.simplefilter("ignore")
OUT = ROOT / "reports" / "calibration"


def main(theta: float = 0.75):
    W = wcsl.load_fervo()
    cal = wcsl.calibrate(W, theta=theta)
    model = wcsl.WCSL(cal.params)
    plans = cal.plans()
    ex = model.expected(plans)
    sim = model.simulate(plans, n_draws=10000, uncertainty=cal.uncertainty())
    Pc, Pd = wcsl.percentiles(sim["cost"]), wcsl.percentiles(sim["days"])
    n = np.array([w["n"] for w in W])
    obs_c = np.array([w["cost"] for w in W])
    obs_d = np.array([w["days"] for w in W])
    td = np.array([w["td_ft"] for w in W])

    # ---------- learning curves on days per 1,000 ft
    y = obs_d / (td / 1000)
    fits = dict(bm=wcsl.fit_brett_millheim(n, y), wright=wcsl.fit_wright(n, y), quad=wcsl.fit_quadratic(n, y))
    x = np.linspace(1, 20, 200)
    fig, ax = plt.subplots(figsize=(7.5, 4.2))
    ax.plot(x, fits["quad"].predict(x), color=WARN, lw=2, label=f"Quadratic (R² {fits['quad'].r2:.2f})")
    ax.plot(x, fits["wright"].predict(x), color=STEEL, lw=2,
            label=f"Wright, {wcsl.learning_rate(fits['wright'].params['b']):.0%} per doubling (R² {fits['wright'].r2:.2f})")
    ax.plot(x, fits["bm"].predict(x), color=AMBER, lw=2.8, label=f"Brett–Millheim, floor {fits['bm'].params['tau_inf']:.2f} (R² {fits['bm'].r2:.2f})")
    ax.scatter(n, y, color=INK, zorder=5, label="Cape wells (observed)")
    ax.axvline(6.5, color=MUTED, ls=":", lw=1)
    ax.text(6.7, 3.7, "extrapolation", color=MUTED, fontsize=8)
    ax.set_ylim(0, 4)
    ax.set_xlabel("Well number at Cape Station")
    ax.set_ylabel("Drilling days per 1,000 ft")
    ax.set_title("Three curves fit six wells; only one stays physical")
    ax.legend(fontsize=8, loc="upper right")
    save(fig, OUT / "learning_curves.png")

    # ---------- spread rate
    sf = cal.spread_fit
    fig, ax = plt.subplots(figsize=(6.5, 3.8))
    ax.scatter(n, np.array(sf["s_obs"]) / 1e3, color=INK, zorder=5, label="(cost − non-time cost) / days")
    xx = np.linspace(1, 12, 100)
    ax.plot(xx, cal.params.spread.rate(xx) / 1e3, color=AMBER, lw=2.5,
            label=f"Fit: ${sf['s_inf'] / 1e3:.0f}k × (1 + {sf['eps']:.2f} e^(−{sf['k_s']:.2f}(E−1)))")
    ax.axhline(33.5, color=WARN, ls="--", lw=1.5)
    ax.text(11.9, 37, "WCS rig rate, $33.5k/day", color=WARN, ha="right", fontsize=8)
    ax.set_ylim(0, 260)
    ax.set_xlabel("Well number")
    ax.set_ylabel("Spread rate, $k per day")
    ax.set_title(f"Implied daily spread at Cape (time share θ = {theta:.0%})")
    ax.legend(fontsize=8, loc="center right")
    save(fig, OUT / "spread_rate.png")

    # ---------- cost path: observed vs WCS vs WCS + learning vs WCS-L
    wb = V.wcs_baseline(W)
    anchored = V.anchored_path(W, V.WCS_THETA)
    fig, ax = plt.subplots(figsize=(7.5, 4.3))
    ax.fill_between(n, Pc[0] / 1e6, Pc[2] / 1e6, color=AMBER, alpha=0.18, lw=0, label="WCS-L P10–P90")
    ax.plot(n, [w["cost"] / 1e6 for w in ex["wells"]], color=AMBER, lw=2.8, marker="o", ms=4, label="WCS-L expected")
    ax.plot(n, anchored / 1e6, color=STEEL, lw=2, ls="--", marker="s", ms=4, label="WCS structure + same learning (anchored on well 1)")
    ax.plot(n, wb["cost"] / 1e6, color=GREY, lw=2, marker="D", ms=4, label="WCS, NREL 2025 parameters")
    ax.scatter(n, obs_c / 1e6, color=INK, s=40, zorder=6, label="Observed (Akindipe & Witter 2025)")
    ax.set_ylim(0, 10)
    ax.set_xlabel("Cape well")
    ax.set_ylabel("Cost per well, $M")
    ax.set_title("Cape Station cost path: data vs models")
    ax.legend(fontsize=8, loc="lower left")
    save(fig, OUT / "cost_path.png")

    # ---------- cost structure per well: WCS vs WCS-L
    fig, ax = plt.subplots(figsize=(7.5, 3.8))
    r = wcs.well_cost(wcs.cape_like_well())
    sp = wcs.cost_split(r)
    labels = ["WCS\n(any well)"] + [f"WCS-L\nwell {i}" for i in n]
    time = [sp["time"] * r["total"]] + [w["time_cost"] + w["skid"] for w in ex["wells"]]
    mat = [sp["materials"] * r["total"]] + [w["materials"] for w in ex["wells"]]
    fix = [sp["fixed"] * r["total"]] + [w["wellhead"] for w in ex["wells"]]
    xs = np.arange(len(labels))
    ax.bar(xs, np.array(time) / 1e6, color=STEEL, label="Time-driven (spread × days)")
    ax.bar(xs, np.array(mat) / 1e6, bottom=np.array(time) / 1e6, color=ROCK, label="Materials")
    ax.bar(xs, np.array(fix) / 1e6, bottom=(np.array(time) + np.array(mat)) / 1e6, color=GREY, label="Fixed / wellhead")
    for xi, t, m_, f in zip(xs, time, mat, fix):
        ax.text(xi, (t + m_ + f) / 1e6 + 0.12, f"{t / (t + m_ + f):.0%}", ha="center", fontsize=8, color=STEEL)
    ax.set_xticks(xs)
    ax.set_xticklabels(labels, fontsize=8)
    ax.set_ylabel("$M per well")
    ax.set_title("Where the money goes (label = time-driven share)")
    ax.legend(fontsize=8, loc="upper right")
    save(fig, OUT / "cost_structure.png")

    out = dict(theta=theta,
               learning=dict(tau_1=cal.bm_fit.params["tau_inf"] + cal.bm_fit.params["excess"], **cal.bm_fit.params,
                             floor=cal.params.floor, r2=cal.bm_fit.r2, note="fitted on planned days per 1,000 ft (expected trouble removed)"),
               curve_comparison={k: dict(params=v.params, r2=v.r2, rmse=v.rmse, aic=v.aic,
                                         at_well_20=float(v.predict(np.array([20.0]))[0])) for k, v in fits.items()},
               time_scale=cal.params.time_scale, materials_mult=cal.materials_mult, nontime_cost_per_well=cal.nontime_cost,
               spread={k: v for k, v in sf.items() if k != "s_obs"}, spread_observed=sf["s_obs"],
               uncertainty=cal.uncertainty().__dict__,
               wells=[dict(well=w["well"], observed_cost=w["cost"], observed_days=w["days"], **{k: e[k] for k in
                           ("E", "days", "planned_days", "npt_days", "spread_rate", "time_cost", "materials", "cost", "time_share")},
                           cost_p10=float(Pc[0, i]), cost_p50=float(Pc[1, i]), cost_p90=float(Pc[2, i]),
                           days_p10=float(Pd[0, i]), days_p90=float(Pd[2, i]))
                      for i, (w, e) in enumerate(zip(W, ex["wells"]))],
               in_sample=dict(mape_cost=V.mape([e["cost"] for e in ex["wells"]], obs_c),
                              mape_days=V.mape([e["days"] for e in ex["wells"]], obs_d),
                              coverage_cost=V.coverage(Pc[0], Pc[2], obs_c), coverage_days=V.coverage(Pd[0], Pd[2], obs_d),
                              wcs_mape_cost=V.mape(wb["cost"], obs_c), anchored_wcs_structure_mape=V.mape(anchored, obs_c)))
    dump(out, OUT / "calibration.json")
    print({k: out[k] for k in ("learning", "time_scale", "materials_mult", "spread", "in_sample")})


if __name__ == "__main__":
    main()
