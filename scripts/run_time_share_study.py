"""Study 1: what share of a well's cost is driven by time?

Reads data/time_share/sources.csv and components.csv, computes each source's share under two
definitions, adds the Cape Station derivation and the WCS model, and writes
reports/time_share/{results.json, results.md, time_share_by_source.png}.

Definitions (shares of total well cost, completion and stimulation excluded where reported):
  spread  : everything billed per unit of time (rig day rate, crew, fuel, water, rentals,
            day-rate services such as directional and mud engineering, supervision) and the
            costs of time-based activities in activity-based models.
  broad   : spread + consumables whose quantity follows drilling time (bits, mud materials).
WCS-L prices time with the spread definition (bits are materials that fall with drilling hours),
so the spread share is the one used to calibrate the model.
"""
from __future__ import annotations

import csv
from collections import defaultdict

import numpy as np

from _common import AMBER, GOOD, GREY, INK, MUTED, ROOT, STEEL, WARN, dump, pct, plt, save
from wcsl import wcs

OUT = ROOT / "reports" / "time_share"
CATEGORY_LABEL = {
    "modern_egs_horizontal": "Modern horizontal EGS (Fervo, 2022-24)",
    "egs_study": "EGS engineering studies (2005-23)",
    "conventional_geothermal": "Conventional geothermal wells (1979-2018)",
    "bound": "Bounds implied by other studies",
    "model_wcs": "WCS model (Sandia 2017, used by NREL 2025)",
}
CATEGORY_COLOR = dict(modern_egs_horizontal=AMBER, egs_study=STEEL, conventional_geothermal=GOOD, bound=GREY,
                      model_wcs=WARN)


def read(name):
    with open(ROOT / "data" / "time_share" / name) as fh:
        return list(csv.DictReader(fh))


def itemised_shares(components):
    by = defaultdict(list)
    for c in components:
        by[c["source"]].append(c)
    out = {}
    for src, rows in by.items():
        tot = sum(float(r["value"]) for r in rows)
        cls = defaultdict(float)
        for r in rows:
            cls[r["class"]] += float(r["value"])
        spread = cls["time"] / tot
        broad = (cls["time"] + cls["rop_consumable"]) / tot
        lab = cls["labour_capital"] / tot
        out[src] = dict(spread=spread, broad=broad, labour_capital=lab,
                        classes={k: v / tot for k, v in cls.items()}, total=tot)
    return out


def cape_derived():
    rows = [r for r in csv.DictReader(open(ROOT / "data" / "fervo" / "fervo_wells.csv")) if r["cost_usd"]]
    cost = np.array([float(r["cost_usd"]) for r in rows])
    res = {}
    for M in (1.2e6, 1.6e6, 2.0e6):
        share = 1 - M / cost
        res[M] = dict(per_well=share.tolist(), campaign=float(1 - M * len(cost) / cost.sum()))
    return res


def wcs_shares():
    cfg = {
        "Vertical 5,000 m, GeoVision base (2016 prices)": wcs.vertical_well(5000),
        "Vertical 5,000 m, NREL 2025 parameters": wcs.with_params(wcs.vertical_well(5000), **wcs.NREL_2025),
        "Cape-like horizontal, NREL 2025 parameters": wcs.cape_like_well(),
    }
    out = {}
    for name, spec in cfg.items():
        r = wcs.well_cost(spec)
        sp = wcs.cost_split(r)
        L = r["lines"]
        bits = (L["bits"] + L["bha"]) * (1 + r["params"]["cont"]) / r["total"]
        out[name] = dict(spread=sp["time"], broad=sp["time"] + bits, materials=sp["materials"], fixed=sp["fixed"],
                         total=r["total"])
    return out


def main():
    sources = {s["id"]: s for s in read("sources.csv")}
    comp = itemised_shares(read("components.csv"))
    cape = cape_derived()
    wm = wcs_shares()

    rows = []   # one row per plotted estimate
    for sid in ("carson_lin_1981", "thorhallsson_2012", "gul_2018"):
        s, c = sources[sid], comp[sid]
        label = dict(carson_lin_1981="Carson & Lin (c. 1981), Baca NM: 15 wells",
                     thorhallsson_2012="Thorhallsson & Sveinbjornsson (2012), Iceland",
                     gul_2018="Gul & Aslanoglu (2018), Turkey")[sid]
        rows.append(dict(id=sid, category=s["category"], label=label,
                         low=c["spread"], central=c["spread"], high=c["broad"], basis=s["basis"],
                         note="point = spread definition; bar extends to the broad definition (+ bits and mud)"))
    c = comp["polsky_2008"]
    rows.append(dict(id="polsky_2008", category="egs_study", label="Polsky et al. (2008), Sandia EGS well model",
                     low=c["spread"], central=c["spread"], high=c["broad"], basis="itemised_model",
                     note="activity costs at a universal $58,130/day spread"))
    c = comp["price_ornl_2023"]
    rows.append(dict(id="price_ornl_2023", category="egs_study", label="Price et al. (2023), ORNL road map",
                     low=c["spread"], central=c["spread"] + c["labour_capital"] / 2, high=c["spread"] + c["labour_capital"],
                     basis="itemised_model", note="labour and cost of capital (25 %) counted as non-time (low) to time (high)"))
    s = sources["bloomfield_laney_2005"]
    rows.append(dict(id=s["id"], category=s["category"], label="Bloomfield & Laney (2005), INL EGS cost study",
                     low=float(s["share_low"]), central=float(s["share_central"]), high=float(s["share_high"]),
                     basis="stated", note="'the remainder of the cost is proportional to the time for drilling'"))
    s = sources["elsadi_2024"]
    rows.append(dict(id=s["id"], category=s["category"], label="El-Sadi et al. (2024), Fervo: 'over 75 %'",
                     low=0.75, central=0.75, high=None, basis="stated", note="lower bound stated by the operator"))
    rows.append(dict(id="cape_derived", category="modern_egs_horizontal",
                     label="Cape Station wells, derived from Akindipe & Witter (2025)",
                     low=cape[2.0e6]["campaign"], central=cape[1.6e6]["campaign"], high=cape[1.2e6]["campaign"],
                     basis="derived", note="1 - non-time cost / cost; non-time cost $1.2M-$2.0M per well (central $1.6M)"))
    rows.append(dict(id="tester_mit_2006", category="bound", label="Tester et al. (2006), MIT: ROP sensitivity",
                     low=0.44, central=None, high=None, basis="lower_bound", note="ROP-linked share >= 44 % (4,000 m well)"))
    rows.append(dict(id="finger_2010", category="bound", label="Finger & Blankenship (2010), Sandia handbook",
                     low=None, central=None, high=0.60, basis="upper_bound", note="casing and cementing > 40 % (conventional)"))
    for name, v in wm.items():
        rows.append(dict(id="wcs_2017", category="model_wcs", label=f"WCS: {name}", low=v["spread"], central=v["spread"],
                         high=v["broad"], basis="model", note="rig time + mud per day + directional; + bits and BHA (broad)"))

    # synthesis by category (central values; spread definition)
    syn = {}
    for cat in ("conventional_geothermal", "egs_study", "modern_egs_horizontal", "model_wcs"):
        vals = [r["central"] for r in rows if r["category"] == cat and r["central"] is not None]
        lows = [r["low"] for r in rows if r["category"] == cat and r["low"] is not None]
        highs = [r["high"] for r in rows if r["category"] == cat and r["high"] is not None]
        syn[cat] = dict(n=len(vals), mean_central=float(np.mean(vals)), min_low=float(min(lows)), max_high=float(max(highs)))
    priors = {
        "modern_egs_horizontal": dict(central=0.75, p05=0.60, p95=0.86,
                                      rationale="Operator statement (>= 0.75) and the Cape derivation (0.66-0.80 campaign, "
                                                "0.76-0.86 on the first wells); lower end allows for a broader definition of materials"),
        "vertical_or_deep_egs": dict(central=0.52, p05=0.45, p95=0.62,
                                     rationale="Itemised conventional wells (0.47-0.49 spread, 0.61-0.64 broad) and EGS studies (0.46-0.60)"),
    }
    # Beta prior matching the modern-EGS central value and 90 % interval
    from scipy import optimize, stats
    def beta_fit(m, lo, hi):
        def f(c):
            a, b = m * c, (1 - m) * c
            return [stats.beta.ppf(0.05, a, b) - lo]
        c = optimize.brentq(lambda c: f(c)[0], 2, 2000)
        return dict(a=m * c, b=(1 - m) * c)
    for k, v in priors.items():
        v["beta"] = beta_fit(v["central"], v["p05"], v["p95"])
        a, b = v["beta"]["a"], v["beta"]["b"]
        v["beta"]["p05"], v["beta"]["p95"] = float(stats.beta.ppf(0.05, a, b)), float(stats.beta.ppf(0.95, a, b))

    # implied Cape spread rate for each theta (average well: $5.85M, 24.7 days)
    fw = [r for r in csv.DictReader(open(ROOT / "data" / "fervo" / "fervo_wells.csv")) if r["cost_usd"]]
    mc, md = np.mean([float(r["cost_usd"]) for r in fw]), np.mean([float(r["drilling_days"]) for r in fw])
    implied = {f"{t:.2f}": dict(spread_per_day=t * mc / md, nontime_per_well=(1 - t) * mc) for t in (0.5, 0.6, 0.75, 0.86)}

    dump(dict(rows=rows, itemised=comp, cape_derived={str(int(k)): v for k, v in cape.items()}, wcs=wm,
              synthesis=syn, priors=priors, implied_cape=implied), OUT / "results.json")

    # ---------------- figure: share of well cost driven by time, by source
    order = ["modern_egs_horizontal", "egs_study", "conventional_geothermal", "bound", "model_wcs"]
    rows_sorted = [r for cat in order for r in rows if r["category"] == cat]
    fig, ax = plt.subplots(figsize=(9.2, 6.4))
    y = np.arange(len(rows_sorted))[::-1]
    ax.axvspan(priors["modern_egs_horizontal"]["p05"], priors["modern_egs_horizontal"]["p95"], color=AMBER, alpha=0.10, lw=0)
    ax.axvspan(priors["vertical_or_deep_egs"]["p05"], priors["vertical_or_deep_egs"]["p95"], color=STEEL, alpha=0.08, lw=0)
    for yi, r in zip(y, rows_sorted):
        col = CATEGORY_COLOR[r["category"]]
        lo, ce, hi = r["low"], r["central"], r["high"]
        if r["basis"] == "lower_bound":
            ax.annotate("", xy=(0.98, yi), xytext=(lo, yi), arrowprops=dict(arrowstyle="->", color=col, lw=1.6))
            ax.plot(lo, yi, "|", color=col, ms=12, mew=2)
        elif r["basis"] == "upper_bound":
            ax.annotate("", xy=(0.02, yi), xytext=(hi, yi), arrowprops=dict(arrowstyle="->", color=col, lw=1.6))
            ax.plot(hi, yi, "|", color=col, ms=12, mew=2)
        elif hi is None:
            ax.annotate("", xy=(0.98, yi), xytext=(lo, yi), arrowprops=dict(arrowstyle="->", color=col, lw=1.6))
            ax.plot(ce, yi, "o", color=col, ms=7)
        else:
            ax.plot([lo, hi], [yi, yi], color=col, lw=2.2, solid_capstyle="round")
            ax.plot(ce, yi, "o", color=col, ms=7, mec="white", mew=1)
    ax.set_yticks(y)
    ax.set_yticklabels([r["label"] for r in rows_sorted], fontsize=8.5)
    ax.set_xlim(0, 1)
    ax.xaxis.set_major_formatter(plt.FuncFormatter(lambda v, _: f"{v:.0%}"))
    ax.set_xlabel("Share of total well cost driven by time")
    ax.grid(axis="y", visible=False)
    ax.set_title("How much of a well's cost depends on time?")
    handles = [plt.Line2D([], [], color=CATEGORY_COLOR[c], marker="o", lw=2, label=CATEGORY_LABEL[c]) for c in order]
    handles.append(plt.Rectangle((0, 0), 1, 1, color=AMBER, alpha=0.2, label="Recommended prior, horizontal EGS (90 %)"))
    handles.append(plt.Rectangle((0, 0), 1, 1, color=STEEL, alpha=0.15, label="Recommended prior, vertical / deep EGS (90 %)"))
    ax.legend(handles=handles, loc="upper center", bbox_to_anchor=(0.42, -0.09), ncol=2, fontsize=8)
    save(fig, OUT / "time_share_by_source.png")

    # ---------------- markdown table
    lines = ["| Source | Category | Basis | Time share (point) | Range | Note |", "|---|---|---|---|---|---|"]
    for r in rows_sorted:
        f = lambda v: "" if v is None else pct(v)
        rng = f"{f(r['low'])} - {f(r['high'])}" if r["low"] is not None and r["high"] is not None else \
            (f"at least {f(r['low'])}" if r["low"] is not None else f"at most {f(r['high'])}")
        lines.append(f"| {r['label']} | {CATEGORY_LABEL[r['category']]} | {r['basis']} | {f(r['central'])} | {rng} | {r['note']} |")
    lines += ["", "Synthesis (spread definition, central values):", ""]
    for cat, v in syn.items():
        lines.append(f"- {CATEGORY_LABEL[cat]}: mean {pct(v['mean_central'])} over {v['n']} estimates "
                     f"(range {pct(v['min_low'])} - {pct(v['max_high'])})")
    (OUT / "results.md").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))
    print(priors)
    print(implied)


if __name__ == "__main__":
    main()
