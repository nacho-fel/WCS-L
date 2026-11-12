"""Spread rate: the all-in cost of one day on the well (rig, services, rentals, fuel, supervision).

    s_n = s_inf * (1 + eps * exp(-k_s (E_n - 1)))

The first wells pay a premium per day (extra services, standby tools, supervision) that fades
as routines settle. ``components`` documents what s_inf contains; it does not change the maths.
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np


@dataclass
class SpreadRate:
    s_inf: float = 170000.0
    eps: float = 0.25
    k_s: float = 0.5
    components: dict[str, float] = field(default_factory=dict)

    def rate(self, E):
        return self.s_inf * (1.0 + self.eps * np.exp(-self.k_s * (np.asarray(E, float) - 1.0)))

    def premium_left(self, E):
        """Fraction of the first-well premium still paid at experience E."""
        return np.exp(-self.k_s * (np.asarray(E, float) - 1.0))
