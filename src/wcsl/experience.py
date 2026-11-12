"""Experience E_n: how much relevant experience the team has when it starts well n.

    E_n = 1 + sum_{m < n} delta^(t_n - t_m) * omega^1[different rig] * gamma^1[other project]

Every earlier well adds a weight between 0 and 1. The weight is 1 for a well drilled just
before, by the same rig, in the same project; forgetting (delta), a different rig (omega)
and another project (gamma) each shrink it. Values from Kellogg (2011, QJE): a well drilled
a year ago counts about 44 %; experience from another rig is worth about 40 % (5.2 %/yr vs
13 %/yr learning); experience from elsewhere about half.
"""
from __future__ import annotations

from typing import Sequence

import numpy as np

DELTA = 0.44   # forgetting per year
OMEGA = 0.40   # weight of a well drilled by a different rig
GAMMA = 0.50   # weight of a well drilled in another project


def experience(t_years: Sequence[float], rig: Sequence | None = None, project: Sequence | None = None,
               delta: float = DELTA, omega: float = OMEGA, gamma: float = GAMMA) -> np.ndarray:
    """Experience at the start of each well, in the order given.

    ``t_years`` are the start dates in years (any origin). Wells must be listed in drilling
    order; a well's experience counts only the wells listed before it.
    """
    t = np.asarray(t_years, dtype=float)
    n = len(t)
    rig = np.asarray(rig if rig is not None else [0] * n, dtype=object)
    project = np.asarray(project if project is not None else [0] * n, dtype=object)
    if np.any(np.diff(t) < 0):
        raise ValueError("wells must be listed in drilling order (non-decreasing dates)")
    E = np.ones(n)
    for i in range(n):
        if i == 0:
            continue
        gap = t[i] - t[:i]
        w = delta ** gap
        w = w * np.where(rig[:i] != rig[i], omega, 1.0)
        w = w * np.where(project[:i] != project[i], gamma, 1.0)
        E[i] = 1.0 + w.sum()
    return E


def steady_state_ceiling(interval_years: float, delta: float = DELTA) -> float:
    """Largest experience reachable when wells are drilled every ``interval_years`` forever.

    With r = delta^interval the sum is geometric, so E never exceeds 1 + r / (1 - r).
    """
    r = delta ** interval_years
    return float("inf") if r >= 1 else 1.0 + r / (1.0 - r)
