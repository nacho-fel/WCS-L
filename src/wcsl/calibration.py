"""Calibrate WCS-L to an operator's well series (Fervo's Cape Station wells by default).

Steps (each one is reported, so the calibration can be audited):
1. Experience E_n of each well (monthly drilling, one rig, one project).
2. Remove the expected trouble days, then fit the Brett-Millheim curve to planned days per 1,000 ft:
   floor = tau_inf / tau_1 and learning speed k.
3. Time scale: ratio of the observed first-well planned hours to the WCS first-well hours.
4. Cost structure: the non-time cost per well is (1 - theta) x mean observed cost, where theta is
   the time-driven share from the time-share study. Materials are the WCS quantities scaled by one
   multiplier so that they match it.
5. Spread rate: fit s_inf, eps and k_s to the observed (cost - non-time cost) / days.
"""
from __future__ import annotations

import csv
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Callable, Sequence

import numpy as np
from scipy.optimize import least_squares

from . import rop as _rop
from . import wcs
from .experience import experience
from .learning import Fit, bm_factor, fit_brett_millheim
from .model import ACTIVITIES, Uncertainty, WCSLParams, WellPlan, _well_static
from .npt import NPTModel
from .spread import SpreadRate

DATA = Path(__file__).resolve().parents[2] / "data"


def load_fervo(project: str | None = "Cape Station") -> list[dict[str, Any]]:
    """Fervo well data (Akindipe & Witter 2025, Table 1; El-Sadi et al. 2024)."""
    rows = []
    with open(DATA / "fervo" / "fervo_wells.csv") as fh:
        for r in csv.DictReader(fh):
            if project and r["project"] != project:
                continue
            rows.append(dict(well=r["well"], project=r["project"], n=int(r["order_in_project"]),
                             order=int(r["order_overall"]), days=float(r["drilling_days"]), td_ft=float(r["td_ft"]),
                             cost=float(r["cost_usd"]) if r["cost_usd"] else None))
    return rows


def default_spec_builder(md_ft: float) -> dict[str, Any]:
    """First-well WCS geometry for a Cape-like producer of measured depth ``md_ft`` (2025 parameters),
    with the on-bottom ROP falling with TVD (``rop.BASE``; 75 ft/h at the Cape reservoir depth)."""
    return wcs.cape_like_well(md_ft=md_ft, rop_profile=_rop.BASE)


@dataclass
class Calibration:
    params: WCSLParams
    theta: float
    bm_fit: Fit
    spread_fit: dict[str, Any]
    nontime_cost: float
    materials_mult: float
    E: np.ndarray
    wells: list[dict[str, Any]]
    spec_builder: Callable = field(repr=False, default=default_spec_builder)
    months_between: float = 1.0

    def plans(self, wells: Sequence[dict[str, Any]] | None = None) -> list[WellPlan]:
        """Well plans for a list of wells (``n`` = position in the campaign, ``td_ft``)."""
        wells = self.wells if wells is None else wells
        return [WellPlan(spec=self.spec_builder(w["td_ft"]), t_years=(w["n"] - 1) * self.months_between / 12,
                         rig="rig-1", project="cape", pad=1, name=w.get("well", "")) for w in wells]

    time_well_logsd: float = 0.05

    def uncertainty(self) -> Uncertainty:
        """Parameter uncertainty from the fit covariances (with sensible minimums)."""
        p = self.bm_fit.params
        cov = self.bm_fit.cov
        t1 = p["tau_inf"] + p["excess"]
        # floor = tau_inf / (tau_inf + excess): gradient for the delta method
        g = np.array([p["excess"] / t1 ** 2, -p["tau_inf"] / t1 ** 2, 0.0])
        floor_sd = float(np.sqrt(max(g @ cov @ g, 0.0))) if np.all(np.isfinite(cov)) else 0.04
        k_logsd = float(np.sqrt(max(cov[2, 2], 0.0)) / p["k"]) if np.isfinite(cov[2, 2]) else 0.25
        sf = self.spread_fit
        return Uncertainty(floor_sd=float(np.clip(floor_sd, 0.02, 0.10)), k_logsd=float(np.clip(k_logsd, 0.10, 0.50)),
                           s_inf_logsd=float(np.clip(sf.get("s_inf_logsd", 0.10), 0.05, 0.20)),
                           eps_sd=float(np.clip(sf.get("eps_sd", 0.08), 0.04, 0.15)),
                           materials_logsd=0.15, materials_well_logsd=0.05, time_well_logsd=self.time_well_logsd)


def calibrate(wells: Sequence[dict[str, Any]] | None = None, theta: float = 0.75,
              npt: NPTModel | None = None, spec_builder: Callable = default_spec_builder,
              months_between: float = 1.0, k_s: float | None = None,
              spread_wells: Sequence[dict[str, Any]] | None = None) -> Calibration:
    """Calibrate WCS-L to a series of wells drilled by one rig in one project.

    ``wells``: dicts with ``n`` (position in the campaign), ``days``, ``td_ft`` and ``cost``.
    ``theta``: time-driven share of well cost (see studies/time_share). ``k_s`` fixes the speed at
    which the spread-rate premium fades (fitted when None). ``spread_wells``: extra wells with observed
    days and cost (e.g. later wells known only as cost per foot) used in step 5 only, to pin the long-run
    spread rate s_inf.
    """
    wells = list(load_fervo() if wells is None else wells)
    npt = npt or NPTModel()
    n = np.array([w["n"] for w in wells], float)
    days = np.array([w["days"] for w in wells], float)
    td = np.array([w["td_ft"] for w in wells], float)
    cost = np.array([w["cost"] for w in wells], float)
    E = experience((n - 1) * months_between / 12)
    # E for the campaign position (the training wells may not start at n = 1)
    full = experience(np.arange(int(n.max())) * months_between / 12)
    E = full[(n - 1).astype(int)]

    # 2. learning of planned time
    planned = days - npt.expected_days(E)
    bm = fit_brett_millheim(E, planned / (td / 1000))
    tau_inf, excess, k = bm.params["tau_inf"], bm.params["excess"], bm.params["k"]
    tau_1 = tau_inf + excess
    floor = tau_inf / tau_1

    # 3. time scale from the WCS first-well hours per 1,000 ft
    statics = [_well_static(WellPlan(spec=spec_builder(t)), WCSLParams()) for t in td]
    wcs_h_per_kft = np.mean([sum(s["H"].values()) / (s["md"] / 1000) for s in statics])
    time_scale = tau_1 * 24 / wcs_h_per_kft

    # 4. cost structure from theta
    wellhead = WCSLParams().wellhead
    nontime = (1 - theta) * cost.mean()
    fac = bm_factor(E, floor, k)
    base_mat = np.mean([s["mat_fixed"] + s["mat_bits"] * f for s, f in zip(statics, fac)])
    materials_mult = max(nontime - wellhead, 0.0) / base_mat

    # 5. spread rate
    mat = np.array([(s["mat_fixed"] + s["mat_bits"] * f) * materials_mult for s, f in zip(statics, fac)])
    s_obs = (cost - mat - wellhead) / days
    if spread_wells:
        nx = np.array([w["n"] for w in spread_wells], float)
        fullx = experience(np.arange(int(nx.max())) * months_between / 12)
        Ex = fullx[(nx - 1).astype(int)]
        stx = [_well_static(WellPlan(spec=spec_builder(w["td_ft"])), WCSLParams()) for w in spread_wells]
        facx = bm_factor(Ex, floor, k)
        matx = np.array([(q["mat_fixed"] + q["mat_bits"] * f) * materials_mult for q, f in zip(stx, facx)])
        sx = (np.array([w["cost"] for w in spread_wells]) - matx - wellhead) / np.array([w["days"] for w in spread_wells])
        E = np.concatenate([E, Ex])
        s_obs = np.concatenate([s_obs, sx])

    def resid(x):
        s_inf, eps, ks = x[0], x[1], (k_s if k_s is not None else x[2])
        return s_inf * (1 + eps * np.exp(-ks * (E - 1))) - s_obs

    x0 = [s_obs.min(), max(s_obs[0] / s_obs.min() - 1, 0.05), 0.5]
    lo, hi = [1e3, 0.0, 0.05], [1e7, 5.0, 5.0]
    if k_s is not None:
        x0, lo, hi = x0[:2], lo[:2], hi[:2]
    sol = least_squares(resid, x0, bounds=(lo, hi))
    s_inf, eps = sol.x[0], sol.x[1]
    ks = k_s if k_s is not None else sol.x[2]
    J = sol.jac
    dof = max(len(s_obs) - len(sol.x), 1)
    s2 = float(sol.fun @ sol.fun) / dof
    try:
        cov = np.linalg.inv(J.T @ J) * s2
    except np.linalg.LinAlgError:
        cov = np.full((len(sol.x), len(sol.x)), np.nan)
    sf = dict(s_inf=float(s_inf), eps=float(eps), k_s=float(ks), s_obs=s_obs.tolist(),
              s_inf_logsd=float(np.sqrt(abs(cov[0, 0])) / s_inf) if np.isfinite(cov[0, 0]) else 0.1,
              eps_sd=float(np.sqrt(abs(cov[1, 1]))) if np.isfinite(cov[1, 1]) else 0.08,
              r2=float(1 - np.sum(sol.fun ** 2) / np.sum((s_obs - s_obs.mean()) ** 2)) if len(s_obs) > 1 else float("nan"))

    E = E[:len(wells)]
    params = WCSLParams(floor=float(floor), k=float(k), time_scale=float(time_scale), npt=npt,
                        spread=SpreadRate(s_inf=float(s_inf), eps=float(eps), k_s=float(ks)),
                        materials_mult=float(materials_mult))
    return Calibration(params=params, theta=theta, bm_fit=bm, spread_fit=sf, nontime_cost=float(nontime),
                       materials_mult=float(materials_mult), E=E, wells=wells, spec_builder=spec_builder,
                       months_between=months_between)


#: Days of Cape wells 7-17 (Fervo deck, Apr 2025, slide 14, read from the chart) and the late cost per foot.
CAPE_LATE_DAYS = [23.4, 23.7, 21.4, 19.5, 19.8, 17.3, 25.3, 20.6, 22.9, 17.9, 22.1]
CAPE_LATE_COST_PER_FT = 329.0
CAPE_LATE_TD_FT = 13527.0
LAMBDA_INF = 0.35       # long-run trouble events per well (studies/band)
SIGMA_TIME = 0.11       # well-to-well planned-time residual: P10-P90 leaves out only the fastest and slowest Cape wells


def calibrate_cape17(theta: float = 0.75, spec_builder: Callable = default_spec_builder) -> Calibration:
    """Production calibration (24 Nov 2026): learning and materials from Cape 1-6, long-run spread rate
    pinned with the cost per foot of Cape 7-17, band widened so P10-P90 contains every Cape well but the fastest
    (Cape 12) and the slowest (Cape 13), each left out by ~0.3 d (lambda_inf = 0.35, sigma_t = 0.11;
    see studies/band/recalibrate_band.py)."""
    late = [dict(n=7 + i, td_ft=CAPE_LATE_TD_FT, days=d, cost=CAPE_LATE_COST_PER_FT * CAPE_LATE_TD_FT)
            for i, d in enumerate(CAPE_LATE_DAYS)]
    cal = calibrate(theta=theta, npt=NPTModel(lambda_inf=LAMBDA_INF), spread_wells=late, spec_builder=spec_builder)
    cal.time_well_logsd = SIGMA_TIME
    return cal
