"""Method K: choose the well design by cost per MW, not cost per well.

Each design option has a WCS geometry (for cost) and an output per well in MW. The MW must come
from GEOPHIRES (or another reservoir model) for that geometry; the values in the examples are
placeholders. For a demand in MW, each option needs ceil(demand / MW per well) wells; WCS-L prices
that campaign with learning, and the option with the lowest $ per MW wins.
"""
from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Any, Callable, Sequence

import numpy as np

from .model import WCSL, Uncertainty, WCSLParams, uniform_campaign


@dataclass
class DesignOption:
    name: str
    spec_builder: Callable[[], dict[str, Any]]
    mw_per_well: float
    wells_per_pad: int = 4


def evaluate(options: Sequence[DesignOption], demand_mw: float, params: WCSLParams,
             months_between: float = 1.0, n_draws: int = 3000, seed: int = 11,
             uncertainty: Uncertainty | None = None) -> list[dict[str, Any]]:
    """Campaign cost and $ per MW (P10/P50/P90) for each design option, sorted by expected $ per MW."""
    model = WCSL(params)
    out = []
    for o in options:
        n = math.ceil(demand_mw / o.mw_per_well)
        plans = uniform_campaign(n, o.spec_builder, months_between=months_between, wells_per_pad=o.wells_per_pad)
        ex = model.expected(plans)
        sim = model.simulate(plans, n_draws=n_draws, seed=seed, uncertainty=uncertainty)
        mw = n * o.mw_per_well
        p10, p50, p90 = np.percentile(sim["total"], [10, 50, 90])
        out.append(dict(design=o.name, wells=n, mw=mw, expected_total=ex["total"], usd_per_mw=ex["total"] / mw,
                        usd_per_mw_p10=p10 / mw, usd_per_mw_p50=p50 / mw, usd_per_mw_p90=p90 / mw,
                        mean_cost_per_well=ex["mean_cost_per_well"]))
    return sorted(out, key=lambda r: r["usd_per_mw"])
