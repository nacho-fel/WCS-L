"""WCS-L: the WCS with learning, a spread rate and random trouble time.

For a campaign of wells (in drilling order) the model computes, for each well n:

    E_n  experience                      (experience.py)
    T_n  rig days = planned time of each activity learning toward its floor + trouble days
    s_n  spread rate, $/day              (spread.py)
    M_n  materials: casing, cement, bits, BHA, logging service (bits fall with drilling hours)
    cost_n = s_n * T_n + M_n + wellhead (+ the rig move on a pad)

and adds campaign-level costs once: mobilisation and engineering, site preparation per pad and
the rig move between pads. Geometry, first-well hours and material quantities come from the
validated WCS engine (wcs.py); WCS-L changes how time is priced and how it falls well to well.
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Any, Sequence

import numpy as np

from . import wcs
from .experience import DELTA, GAMMA, OMEGA, experience
from .learning import bm_factor
from .npt import NPTModel
from .spread import SpreadRate

ACTIVITIES = ("drilling", "tripping", "casing_cementing", "logging", "flat")
_WCS_TO_ACTIVITY = dict(drill="drilling", trip="tripping", cemt="casing_cementing", logt="logging",
                        wht="flat", trouble="flat")


def activity_hours(result: dict[str, Any]) -> dict[str, float]:
    """WCS hours grouped into the five WCS-L activities."""
    out = {a: 0.0 for a in ACTIVITIES}
    for k, v in result["hours"].items():
        out[_WCS_TO_ACTIVITY[k]] += v
    return out


@dataclass
class WellPlan:
    """One well of a campaign. ``spec`` is a WCS spec at first-well technology."""
    spec: dict[str, Any]
    t_years: float = 0.0
    rig: str = "rig-1"
    project: str = "project-1"
    pad: Any = 1
    generation: int = 0
    name: str = ""


@dataclass
class WCSLParams:
    # learning of planned time (Brett-Millheim), common to all activities unless overridden
    floor: float = 0.633
    k: float = 0.994
    activity_floor: dict[str, float] = field(default_factory=dict)
    activity_k: dict[str, float] = field(default_factory=dict)
    generation_floor: dict[int, float] = field(default_factory=dict)  # multiplies the floor (technology step)
    time_scale: float = 1.0      # multiplies the WCS first-well hours (calibration)
    flat_hours: float = 0.0      # first-well flat time not in the WCS (rig up, circulating, BOP...), h
    # experience
    delta: float = DELTA
    omega: float = OMEGA
    gamma: float = GAMMA
    # trouble time and daily cost
    npt: NPTModel = field(default_factory=NPTModel)
    spread: SpreadRate = field(default_factory=SpreadRate)
    # materials
    materials_mult: float = 1.0
    price_index: dict[str, float] = field(default_factory=lambda: dict(steel=1.0, cement=1.0, bits=1.0))
    wellhead: float = 80000.0
    skid_days: float = 1.0       # rig move between wells on one pad, days at the spread rate
    # campaign
    mob: float = 260000.0        # rig mobilisation + engineering, once per campaign
    site: float = 250000.0       # site preparation per pad
    pad_move_frac: float = 0.4   # rig move between pads, as a fraction of a full mobilisation
    mob_rig: float = 250000.0

    def floor_for(self, activity: str, generation: int = 0) -> float:
        f = self.activity_floor.get(activity, self.floor)
        return min(1.0, f * self.generation_floor.get(generation, 1.0))

    def k_for(self, activity: str) -> float:
        return self.activity_k.get(activity, self.k)


@dataclass
class Uncertainty:
    """Parameter uncertainty for the Monte Carlo (one draw per simulated campaign)."""
    floor_sd: float = 0.04
    k_logsd: float = 0.25
    s_inf_logsd: float = 0.10
    eps_sd: float = 0.08
    materials_logsd: float = 0.15       # common to the campaign
    materials_well_logsd: float = 0.05  # well to well
    time_well_logsd: float = 0.0        # extra well-to-well noise on planned time


def _well_static(plan: WellPlan, p: WCSLParams) -> dict[str, Any]:
    r = wcs.well_cost(plan.spec)
    H = activity_hours(r)
    H["flat"] += p.flat_hours / max(p.time_scale, 1e-12)
    L = r["lines"]
    idx = p.price_index
    mat_fixed = (L["casing"] * idx.get("steel", 1.0) + L["cement"] * idx.get("cement", 1.0) + L["logs"]) * p.materials_mult
    mat_bits = (L["bits"] + L["bha"]) * idx.get("bits", 1.0) * p.materials_mult
    return dict(H=H, mat_fixed=mat_fixed, mat_bits=mat_bits, md=r["md"], wcs=r)


class WCSL:
    """The WCS-L model for one campaign."""

    def __init__(self, params: WCSLParams | None = None):
        self.p = params or WCSLParams()

    # ---------------------------------------------------------------- experience and pads
    def experience(self, plans: Sequence[WellPlan]) -> np.ndarray:
        return experience([w.t_years for w in plans], [w.rig for w in plans], [w.project for w in plans],
                          self.p.delta, self.p.omega, self.p.gamma)

    @staticmethod
    def _pad_flags(plans: Sequence[WellPlan]) -> tuple[np.ndarray, int]:
        seen, first = set(), []
        for w in plans:
            first.append(w.pad not in seen)
            seen.add(w.pad)
        return np.array(first), len(seen)

    def campaign_fixed(self, plans: Sequence[WellPlan]) -> dict[str, float]:
        _, n_pads = self._pad_flags(plans)
        p = self.p
        return dict(mobilisation=p.mob if plans else 0.0, sites=n_pads * p.site,
                    pad_moves=max(n_pads - 1, 0) * p.pad_move_frac * p.mob_rig)

    # ---------------------------------------------------------------- expected (deterministic) path
    def expected(self, plans: Sequence[WellPlan]) -> dict[str, Any]:
        """Expected cost of each well (trouble days at their mean) and of the campaign."""
        p = self.p
        E = self.experience(plans)
        first_on_pad, _ = self._pad_flags(plans)
        rows = []
        for i, w in enumerate(plans):
            st = _well_static(w, p)
            fac = {a: float(bm_factor(E[i], p.floor_for(a, w.generation), p.k_for(a))) for a in ACTIVITIES}
            planned_h = p.time_scale * sum(st["H"][a] * fac[a] for a in ACTIVITIES)
            npt_d = float(p.npt.expected_days(E[i]))
            days = planned_h / 24 + npt_d
            s = float(p.spread.rate(E[i]))
            time_cost = s * days
            materials = st["mat_fixed"] + st["mat_bits"] * fac["drilling"]
            skid = 0.0 if first_on_pad[i] else p.skid_days * s
            cost = time_cost + materials + p.wellhead + skid
            rows.append(dict(well=w.name or f"well-{i + 1}", E=float(E[i]), md_ft=st["md"], days=days,
                             planned_days=planned_h / 24, npt_days=npt_d, spread_rate=s, time_cost=time_cost,
                             materials=materials, wellhead=p.wellhead, skid=skid, cost=cost,
                             time_share=(time_cost + skid) / cost, activity_factor=fac))
        fixed = self.campaign_fixed(plans)
        total = sum(r["cost"] for r in rows) + sum(fixed.values())
        return dict(wells=rows, campaign_fixed=fixed, total=total, mean_cost_per_well=total / max(len(rows), 1))

    # ---------------------------------------------------------------- Monte Carlo
    def simulate(self, plans: Sequence[WellPlan], n_draws: int = 5000, seed: int = 42,
                 uncertainty: Uncertainty | None = None) -> dict[str, Any]:
        """Monte Carlo of the campaign: parameter uncertainty plus random trouble events.

        Returns arrays of shape (n_draws, n_wells) for days and cost, and campaign totals.
        """
        p, u = self.p, uncertainty or Uncertainty()
        rng = np.random.default_rng(seed)
        E = self.experience(plans)
        first_on_pad, _ = self._pad_flags(plans)
        n = len(plans)
        st = [_well_static(w, p) for w in plans]
        # one parameter draw per simulated campaign
        floor_shift = rng.normal(0.0, u.floor_sd, size=(n_draws, 1))
        k_mult = np.exp(rng.normal(0.0, u.k_logsd, size=(n_draws, 1)))
        s_inf = p.spread.s_inf * np.exp(rng.normal(0.0, u.s_inf_logsd, size=(n_draws, 1)))
        eps = np.clip(p.spread.eps + rng.normal(0.0, u.eps_sd, size=(n_draws, 1)), 0.0, None)
        mat_mult = np.exp(rng.normal(0.0, u.materials_logsd, size=(n_draws, 1))
                          + rng.normal(0.0, u.materials_well_logsd, size=(n_draws, n)))
        planned_h = np.zeros((n_draws, n))
        drill_fac = np.zeros((n_draws, n))
        for i, w in enumerate(plans):
            for a in ACTIVITIES:
                f = np.clip(p.floor_for(a, w.generation) + floor_shift, 0.05, 1.0)
                fac = f + (1 - f) * np.exp(-p.k_for(a) * k_mult * (E[i] - 1.0))
                planned_h[:, i] += p.time_scale * st[i]["H"][a] * fac[:, 0]
                if a == "drilling":
                    drill_fac[:, i] = fac[:, 0]
        if u.time_well_logsd > 0:
            planned_h *= np.exp(rng.normal(0.0, u.time_well_logsd, size=planned_h.shape))
        npt_days = p.npt.sample(E, rng, n_draws)
        days = planned_h / 24 + npt_days
        s = s_inf * (1.0 + eps * np.exp(-p.spread.k_s * (E[None, :] - 1.0)))
        materials = np.array([[x["mat_fixed"] for x in st]]) + np.array([[x["mat_bits"] for x in st]]) * drill_fac
        materials = materials * mat_mult
        skid = np.where(first_on_pad[None, :], 0.0, p.skid_days * s)
        cost = s * days + materials + p.wellhead + skid
        total = cost.sum(axis=1) + sum(self.campaign_fixed(plans).values())
        return dict(E=E, days=days, cost=cost, npt_days=npt_days, total=total,
                    time_share=(s * days + skid) / cost)


def percentiles(a: np.ndarray, q=(10, 50, 90)) -> np.ndarray:
    """Percentiles along the draws axis: shape (len(q), n_wells) or (len(q),)."""
    return np.percentile(a, q, axis=0)


def uniform_campaign(n_wells: int, spec_builder, months_between: float = 1.0, wells_per_pad: int = 4,
                     rig: str = "rig-1", project: str = "project-1") -> list[WellPlan]:
    """n identical wells, one every ``months_between`` months, ``wells_per_pad`` per pad."""
    return [WellPlan(spec=spec_builder(), t_years=i * months_between / 12, rig=rig, project=project,
                     pad=i // wells_per_pad + 1, name=f"well-{i + 1}") for i in range(n_wells)]
