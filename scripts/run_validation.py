"""Validation suite: is WCS-L working?

V1  Engine fidelity: the WCS engine reproduces all 249 Sandia GeoVision workbooks, and the NREL 2025 update.
V2  Structure test: anchored on the first Cape well, which cost structure reproduces the observed path?
V3  In-sample fit of the full Cape path (WCS-L vs WCS).
V4  Out-of-sample: rolling-origin one-step-ahead forecasts and a 4 -> 2 hold-out, against benchmarks.
V5  Probabilistic calibration: coverage of the P10-P90 band, CRPS and PIT of the out-of-sample forecasts.
V6  Cross-project transfer: can experience from Project Red explain the first Cape well?
V7  Activity check: WCS hours by activity vs Utah FORGE 16A(78)-32 daily reports.
V8  Sensitivity: time share theta and spread-rate learning, on a 10-well campaign.

Writes reports/validation/{results.json, results.md, *.png}.
"""
from __future__ import annotations

import csv
import gzip
import json
import warnings
from collections import defaultdict

import numpy as np

from _common import AMBER, GOOD, GREY, INK, MUTED, ROCK, ROOT, STEEL, WARN, dump, pct, plt, save
import wcsl
from wcsl import validation as V
from wcsl import wcs

warnings.simplefilter("ignore")
OUT = ROOT / "reports" / "validation"
MODELS = [("wcsl", "WCS-L"), ("fervo_structure", "Fervo split + learning (1-line model)"),
          ("persistence", "Persistence (last well, per ft)"), ("wcs", "WCS, NREL 2025 parameters"),
          ("wcs_learning", "WCS split + learning (1-line model)")]


def v1_engine():
    specs = json.load(gzip.open(ROOT / "data" / "wcs" / "geovision_specs.json.gz", "rt"))
    err = np.array([wcs.well_cost(s["spec"])["total"] - s["expected_total"] for s in specs])
    base = wcs.well_cost(wcs.vertical_well(5000))["total"]
    upd = wcs.well_cost(wcs.with_params(wcs.vertical_well(5000), **wcs.NREL_2025))["total"]
    return dict(n_workbooks=len(specs), n_exact=int(np.sum(np.abs(err) < 0.5)), max_abs_error_usd=float(np.abs(err).max()),
                nrel_2025_change_vertical_5000m=upd / base - 1, base_total=base, nrel_total=upd)


def v2_structure(W):
    obs = np.array([w["cost"] for w in W])
    rows = {}
    for th in (0.21, 0.40, 0.60, 0.75, 0.86, 1.00):
        a = V.anchored_path(W, th)
        rows[f"{th:.2f}"] = dict(mape=V.mape(a[1:], obs[1:]), well6_change=float(a[-1] / a[0] - 1), path=a.tolist())
    n = np.arange(1, 7)
    fig, ax = plt.subplots(figsize=(7, 4))
    cols = {"0.21": WARN, "0.60": GREY, "0.75": STEEL, "1.00": GOOD}
    for th, c in cols.items():
        lab = "WCS split (21 %)" if th == "0.21" else f"time share {float(th):.0%}"
        ax.plot(n, np.array(rows[th]["path"]) / 1e6, color=c, lw=2, marker="o", ms=4,
                label=f"{lab}: MAPE {rows[th]['mape']:.0%}")
    ax.scatter(n, obs / 1e6, color=INK, s=45, zorder=6, label="Observed")
    ax.set_ylim(0, 9.5)
    ax.set_xlabel("Cape well")
    ax.set_ylabel("Cost per well, $M")
    ax.set_title("Same learning curve, different cost structures (anchored on well 1)")
    ax.legend(fontsize=8, loc="lower left")
    save(fig, OUT / "v2_structure_test.png")
    return rows


def v4_out_of_sample(W):
    ro = V.rolling_origin(W, first_train=3)
    ho = V.holdout(W, n_train=4)
    res = {}
    for key, _ in MODELS:
        pc = np.array([r[key]["cost"][0] for r in ro])
        pdays = np.array([r[key]["days"][0] for r in ro])
        oc = np.array([r["target"]["cost"] for r in ro])
        od = np.array([r["target"]["days"] for r in ro])
        hc, hoc = ho[key]["cost"], np.array([w["cost"] for w in ho["test"]])
        res[key] = dict(rolling_mape_cost=V.mape(pc, oc), rolling_bias_cost=V.bias(pc, oc),
                        rolling_mape_days=V.mape(pdays, od), holdout_mape_cost=V.mape(hc, hoc),
                        holdout_bias_cost=V.bias(hc, hoc), rolling_pred_cost=pc.tolist(), holdout_pred_cost=list(map(float, hc)),
                        rolling_crps_cost=float(np.mean(np.abs(pc - oc))))
    # probabilistic scores for WCS-L (pooled: rolling forecasts + hold-out wells)
    samples_c = [r["wcsl"]["cost_samples"][:, 0] for r in ro] + [ho["wcsl"]["cost_samples"][:, j] for j in range(len(ho["test"]))]
    samples_d = [r["wcsl"]["days_samples"][:, 0] for r in ro] + [ho["wcsl"]["days_samples"][:, j] for j in range(len(ho["test"]))]
    obs_c = [r["target"]["cost"] for r in ro] + [w["cost"] for w in ho["test"]]
    obs_d = [r["target"]["days"] for r in ro] + [w["days"] for w in ho["test"]]
    Sc, Sd = np.array(samples_c).T, np.array(samples_d).T
    Pc, Pd = np.percentile(Sc, [10, 50, 90], axis=0), np.percentile(Sd, [10, 50, 90], axis=0)
    prob = dict(n_forecasts=len(obs_c), coverage_cost_p10_p90=V.coverage(Pc[0], Pc[2], obs_c),
                coverage_days_p10_p90=V.coverage(Pd[0], Pd[2], obs_d),
                crps_cost_wcsl=float(np.mean(V.crps_samples(Sc, obs_c))), pit_cost=V.pit(Sc, obs_c).tolist(),
                pit_days=V.pit(Sd, obs_d).tolist(), rolling_crps_cost_wcsl=float(np.mean(V.crps_samples(Sc[:, :len(ro)], obs_c[:len(ro)]))))
    # figure: one-step-ahead forecasts
    fig, ax = plt.subplots(figsize=(7.5, 4.2))
    x = np.array([r["target"]["n"] for r in ro], float)
    allW = np.array([w["n"] for w in W])
    ax.scatter(allW, [w["cost"] / 1e6 for w in W], color=INK, s=45, zorder=6, label="Observed")
    q = np.array([r["wcsl"]["cost_pct"][:, 0] for r in ro]) / 1e6
    ax.errorbar(x - 0.12, q[:, 1], yerr=[q[:, 1] - q[:, 0], q[:, 2] - q[:, 1]], fmt="o", color=AMBER, capsize=4, lw=2,
                label="WCS-L forecast, P50 with P10–P90")
    off = dict(fervo_structure=0.0, persistence=0.12, wcs=0.24, wcs_learning=-0.24)
    cols = dict(fervo_structure=STEEL, persistence=GREY, wcs=ROCK, wcs_learning=WARN)
    mk = dict(fervo_structure="s", persistence="^", wcs="D", wcs_learning="v")
    for key, lab in MODELS[1:]:
        ax.scatter(x + off[key], np.array(res[key]["rolling_pred_cost"]) / 1e6, color=cols[key], marker=mk[key], s=30, label=lab)
    ax.set_xticks(allW)
    ax.set_ylim(0, 9.5)
    ax.set_xlabel("Cape well (forecasts for wells 4-6 use only the wells before)")
    ax.set_ylabel("Cost per well, $M")
    ax.set_title("One-step-ahead forecasts")
    ax.legend(fontsize=7.5, loc="upper right")
    save(fig, OUT / "v4_out_of_sample.png")
    return res, prob


def v6_transfer(cal):
    red = [w for w in wcsl.load_fervo(None) if w["project"] == "Project Red"]
    red_kft = [w["days"] / (w["td_ft"] / 1000) for w in red]
    cape1 = cal.wells[0]["days"] / (cal.wells[0]["td_ft"] / 1000)
    out = {}
    for gap in (0.75, 1.0, 1.5):
        E = wcsl.experience([0.0, 0.1, 0.1 + gap], project=["red", "red", "cape"], rig=["r", "r", "r"])
        f = float(wcsl.bm_factor(E[2], cal.params.floor, cal.params.k))
        f_red2 = float(wcsl.bm_factor(E[1], cal.params.floor, cal.params.k))
        out[f"gap_{gap}y"] = dict(E_cape1=float(E[2]), predicted_ratio_to_red1=f, observed_ratio=cape1 / red_kft[0],
                                  implied_technology_step=cape1 / red_kft[0] / f)
    out["red_days_per_kft"] = red_kft
    out["cape1_days_per_kft"] = cape1
    out["red2_vs_red1_observed"] = red_kft[1] / red_kft[0]
    return out


def v7_forge():
    hours = defaultdict(float)
    with open(ROOT / "data" / "forge" / "forge_16A_activity_hours.csv") as fh:
        for r in csv.DictReader(fh):
            hours[r["group"]] += float(r["hours"])
    tot = sum(hours.values())
    forge = {k: v / tot for k, v in hours.items()}
    r = wcs.well_cost(wcs.vertical_well(3350))           # a WCS well of similar depth (about 11,000 ft)
    H = wcsl.activity_hours(r)
    th = sum(H.values())
    wcs_split = {k: v / th for k, v in H.items()}
    fig, ax = plt.subplots(figsize=(7, 3.2))
    groups = ["drilling", "tripping", "casing_cementing", "logging", "flat"]
    labels = ["Drilling on bottom", "Tripping", "Casing & cementing", "Logging", "Flat time"]
    y = np.arange(2)
    left = np.zeros(2)
    cols = [AMBER, STEEL, ROCK, GOOD, GREY]
    for g, lab, c in zip(groups, labels, cols):
        vals = np.array([forge.get(g, 0), wcs_split.get(g, 0)])
        ax.barh(y, vals, left=left, color=c, label=lab)
        for yi, v, l in zip(y, vals, left):
            if v > 0.06:
                ax.text(l + v / 2, yi, f"{v:.0%}", ha="center", va="center", color="white", fontsize=8)
        left += vals
    ax.set_yticks(y)
    ax.set_yticklabels(["Utah FORGE 16A\n(daily reports, 75 days)", "WCS hours\n(vertical, ~11,000 ft)"])
    ax.invert_yaxis()
    ax.set_xlim(0, 1)
    ax.xaxis.set_major_formatter(plt.FuncFormatter(lambda v, _: f"{v:.0%}"))
    ax.set_title("Where rig time goes: real well vs WCS")
    ax.legend(fontsize=7.5, ncol=5, loc="upper center", bbox_to_anchor=(0.5, -0.18))
    ax.grid(axis="y", visible=False)
    save(fig, OUT / "v7_forge_activities.png")
    return dict(forge=forge, forge_total_hours=tot, wcs=wcs_split, wcs_flat_missing=forge.get("flat", 0) - wcs_split.get("flat", 0),
                implied_time_scale=1 / (1 - (forge.get("flat", 0) - wcs_split.get("flat", 0))))


def v8_sensitivity(W):
    out = {}
    for th in (0.60, 0.75, 0.86):
        for eps0 in (False, True):
            cal = wcsl.calibrate(W, theta=th)
            p = cal.params
            if eps0:   # conservative: no spread-rate learning; refit s_inf as the mean implied spread
                s_obs = np.array(cal.spread_fit["s_obs"])
                p.spread = wcsl.SpreadRate(s_inf=float(s_obs.mean()), eps=0.0, k_s=0.5)
            m = wcsl.WCSL(p)
            plans = wcsl.uniform_campaign(10, lambda: wcs.cape_like_well(md_ft=13500), months_between=1, wells_per_pad=4)
            ex = m.expected(plans)
            single = wcsl.WCSL(p).expected(plans[:1])
            obs = np.array([w["cost"] for w in W])
            fit = m.expected(cal.plans())
            # response of a mature well (well 10) to drilling 20 % faster: scale planned time by 0.8
            p_fast = wcsl.WCSLParams(**{**p.__dict__, "time_scale": p.time_scale * 0.8})
            fast = wcsl.WCSL(p_fast).expected(plans)
            out[f"theta={th:.2f}, eps={'0' if eps0 else 'fitted'}"] = dict(
                cost_change_if_20pct_less_time=fast["wells"][-1]["cost"] / ex["wells"][-1]["cost"] - 1,
                mean_cost_10_wells=ex["mean_cost_per_well"], first_well=single["total"],
                saving_vs_single_well=1 - ex["mean_cost_per_well"] / single["total"],
                well10_vs_well1=ex["wells"][-1]["cost"] / ex["wells"][0]["cost"] - 1,
                in_sample_mape=V.mape([w["cost"] for w in fit["wells"]], obs))
    return out


def main():
    W = wcsl.load_fervo()
    cal = wcsl.calibrate(W, theta=0.75)
    calj = json.loads((ROOT / "reports" / "calibration" / "calibration.json").read_text())
    res = dict(v1_engine=v1_engine(), v2_structure=v2_structure(W), v3_in_sample=calj["in_sample"])
    res["v4_out_of_sample"], res["v5_probabilistic"] = v4_out_of_sample(W)
    res["v6_transfer"] = v6_transfer(cal)
    res["v7_forge"] = v7_forge()
    res["v8_sensitivity"] = v8_sensitivity(W)
    dump(res, OUT / "results.json")

    r = res
    L = ["# Validation results (generated by scripts/run_validation.py)", ""]
    v1 = r["v1_engine"]
    L += ["## V1 Engine fidelity", "",
          f"- WCS engine reproduces **{v1['n_exact']} of {v1['n_workbooks']}** GeoVision workbooks to the dollar (max error ${v1['max_abs_error_usd']:.2e}).",
          f"- NREL 2025 parameters on the 5,000 m vertical well: {pct(v1['nrel_2025_change_vertical_5000m'], 1)} (NREL reports about -15 %).", ""]
    L += ["## V2 Structure test (anchored on Cape well 1, same learning curve)", "",
          "| Time-driven share | Cost change, well 6 vs 1 | MAPE wells 2-6 |", "|---|---|---|"]
    for th, v in r["v2_structure"].items():
        L.append(f"| {float(th):.0%} | {pct(v['well6_change'])} | {pct(v['mape'])} |")
    L += ["", "Observed change, well 6 vs 1: -42 %.", ""]
    v3 = r["v3_in_sample"]
    L += ["## V3 In-sample fit (6 Cape wells)", "",
          f"- WCS-L: MAPE cost {pct(v3['mape_cost'], 1)}, MAPE days {pct(v3['mape_days'], 1)}; P10-P90 band contains {pct(v3['coverage_cost'])} of costs and {pct(v3['coverage_days'])} of days.",
          f"- WCS (NREL 2025): MAPE cost {pct(v3['wcs_mape_cost'], 1)}. WCS split + learning anchored on well 1: {pct(v3['anchored_wcs_structure_mape'], 1)}.", ""]
    L += ["## V4 Out-of-sample forecasts", "",
          "| Model | One-step-ahead MAPE, cost (wells 4-6) | Bias | MAPE, days | Hold-out 1-4 -> 5-6 MAPE, cost | Bias |",
          "|---|---|---|---|---|---|"]
    for key, lab in MODELS:
        v = r["v4_out_of_sample"][key]
        L.append(f"| {lab} | {pct(v['rolling_mape_cost'], 1)} | {pct(v['rolling_bias_cost'], 1)} | {pct(v['rolling_mape_days'], 1)} | "
                 f"{pct(v['holdout_mape_cost'], 1)} | {pct(v['holdout_bias_cost'], 1)} |")
    p5 = r["v5_probabilistic"]
    L += ["", "## V5 Probabilistic calibration (WCS-L, pooled out-of-sample forecasts)", "",
          f"- {p5['n_forecasts']} forecasts; P10-P90 band contains {pct(p5['coverage_cost_p10_p90'])} of observed costs and {pct(p5['coverage_days_p10_p90'])} of observed days (nominal 80 %).",
          f"- CRPS (cost, one-step-ahead): ${p5['rolling_crps_cost_wcsl'] / 1e6:.2f}M for WCS-L vs mean absolute error "
          + ", ".join(f"${r['v4_out_of_sample'][k]['rolling_crps_cost'] / 1e6:.2f}M ({lab})" for k, lab in MODELS[1:]) + ".",
          f"- PIT values (cost): {', '.join(f'{x:.2f}' for x in p5['pit_cost'])}.", ""]
    v6 = r["v6_transfer"]
    g = v6["gap_1.0y"]
    L += ["## V6 Cross-project transfer (Project Red -> Cape)", "",
          f"- Red wells: {v6['red_days_per_kft'][0]:.2f} and {v6['red_days_per_kft'][1]:.2f} days per 1,000 ft; first Cape well {v6['cape1_days_per_kft']:.2f}.",
          f"- With a one-year gap and experience from another project at half weight, E = {g['E_cape1']:.2f}: learning alone predicts the first Cape well at {pct(g['predicted_ratio_to_red1'])} of Red-1; observed {pct(g['observed_ratio'])}.",
          f"- The remaining factor ({g['implied_technology_step']:.2f}) is a technology, design or geology step: a lower floor (generation_floor), not learning.", ""]
    v7 = r["v7_forge"]
    L += ["## V7 Activity check (Utah FORGE 16A)", "",
          "| Activity | FORGE 16A share of rig hours | WCS share of hours |", "|---|---|---|"]
    for k in ("drilling", "tripping", "casing_cementing", "logging", "flat"):
        L.append(f"| {k} | {pct(v7['forge'].get(k, 0))} | {pct(v7['wcs'].get(k, 0))} |")
    L += ["", f"- Flat time missing from the WCS: {pct(v7['wcs_flat_missing'])} of rig hours, which implies a time scale of about {v7['implied_time_scale']:.2f}; "
          f"the calibrated time scale at Cape is {calj['time_scale']:.2f}.", ""]
    L += ["## V8 Sensitivity: 10-well campaign (4 per pad, monthly, 13,500 ft)", "",
          "| Case | Single-well campaign | Mean cost per well, 10 wells | Saving vs a single well | Well 10 vs well 1 | Well 10 if 20 % less rig time | In-sample MAPE |", "|---|---|---|---|---|---|---|"]
    for k, v in r["v8_sensitivity"].items():
        L.append(f"| {k} | ${v['first_well'] / 1e6:.2f}M | ${v['mean_cost_10_wells'] / 1e6:.2f}M | {pct(v['saving_vs_single_well'])} | {pct(v['well10_vs_well1'])} | {pct(v['cost_change_if_20pct_less_time'])} | {pct(v['in_sample_mape'], 1)} |")
    (OUT / "results.md").write_text("\n".join(L) + "\n")
    print("\n".join(L))


if __name__ == "__main__":
    main()
