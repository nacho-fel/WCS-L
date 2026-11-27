"""Rate of penetration that falls with depth.

The WCS (and every GeoVision workbook) uses one ROP for the whole well, so a deeper well only costs more
because it is longer. Field data say ROP also falls with depth (harder, more compacted rock, higher
pressure, more bit wear). Following the oil-and-gas drilling models (Bourgoyne & Young's compaction term
reduces to an exponential in depth), the on-bottom ROP is

    ROP(z) = ROP_ref * exp(-k (z - z_ref))          z = true vertical depth (TVD), ft

anchored at the Cape-like well (ROP_ref = 75 ft/h at z_ref = 8,700 ft TVD, the NREL 2025 value), so the
calibration to Fervo's Cape wells is unchanged in level; k only decides how the cost grows below Cape.

k is fitted to Fervo's days-versus-depth curve (Pablo's workbook, Stanford/FORGE data, Samuel et al. 2022):
- all-in progress rate (every hour counted): ROP = 60 exp(-0.709 z_km) ft/h, close to the workbook's
  50 exp(-0.0002 z_ft) = 50 exp(-0.656 z_km);
- on-bottom ROP with the WCS's own tripping (bit life 50 h, 1,000 ft/h trips, 6 h per bit change):
  k = 0.609 per km. This is the value used here, because the WCS engine adds trips on top of the
  on-bottom ROP; using the all-in k would count the slowdown from tripping twice.

The lateral is drilled at constant TVD, so it is drilled at the ROP of the reservoir depth.
"""
from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np

FT_PER_KM = 1000 / 0.3048

#: Fervo days-versus-depth points from Pablo's workbook (Sheet2, columns L-M): depth ft, cumulative days.
FERVO_DEPTH_DAYS = (np.arange(0, 10000, 1000.0),
                    np.array([0, 0.25, 0.75, 1.5, 4.5, 6.5, 9, 11.5, 15.5, 18.5]))


@dataclass(frozen=True)
class RopProfile:
    rop_ref: float = 75.0          # ft/h at tvd_ref
    tvd_ref_ft: float = 8700.0
    k_per_km: float = 0.609        # decay per km of TVD (0 = the WCS's constant ROP)
    rop_min: float = 3.0           # floor, ft/h, to keep the integral finite in extreme cases

    @property
    def k_per_ft(self) -> float:
        return self.k_per_km / FT_PER_KM

    def rop(self, tvd_ft):
        r = self.rop_ref * np.exp(-self.k_per_ft * (np.asarray(tvd_ft, float) - self.tvd_ref_ft))
        return np.maximum(r, self.rop_min)

    def hours(self, md_top: float, md_end: float, tvd_of_md, n: int = 400) -> float:
        """On-bottom hours to drill from md_top to md_end (trapezoid on 1 / ROP)."""
        if md_end <= md_top:
            return 0.0
        md = np.linspace(md_top, md_end, n)
        inv = 1.0 / self.rop(tvd_of_md(md))
        return float(np.sum((inv[1:] + inv[:-1]) / 2 * np.diff(md)))

    def effective_rop(self, md_top: float, md_end: float, tvd_of_md) -> float:
        """The constant ROP that drills the interval in the same hours (harmonic mean)."""
        h = self.hours(md_top, md_end, tvd_of_md)
        return (md_end - md_top) / h if h > 0 else self.rop_ref


CONSTANT = RopProfile(k_per_km=0.0)        # the WCS assumption
BASE = RopProfile()                         # fitted to Fervo's depth curve (on-bottom, with WCS trips)
PABLO_ALLIN = RopProfile(k_per_km=0.656)    # the workbook's k applied to on-bottom ROP (upper case)
HALF = RopProfile(k_per_km=0.30)            # half the fitted decay (granite flatter than the fit)


def horizontal_tvd_of_md(tvd_ft: float, bur: float):
    """TVD as a function of MD for a vertical-build-lateral well ending horizontal at ``tvd_ft``."""
    R = 100 * 180 / (math.pi * bur)
    kop = tvd_ft - R
    end_build = kop + math.pi / 2 * R

    def f(md):
        md = np.asarray(md, float)
        build = kop + R * np.sin(np.clip((md - kop) / R, 0, math.pi / 2))
        return np.where(md <= kop, md, np.where(md <= end_build, build, tvd_ft))
    return f


def fit_allin(depth_ft=None, days=None) -> dict:
    """Fit days(z) = (exp(k z) - 1) / (k ROP0) / 24 to a cumulative days-versus-depth curve."""
    from scipy.optimize import least_squares
    z, d = FERVO_DEPTH_DAYS if depth_ft is None else (np.asarray(depth_ft, float), np.asarray(days, float))

    def f(p, zz):
        r0, k = p
        return (np.exp(k * zz) - 1) / (k * r0) / 24
    s = least_squares(lambda p: f(p, z) - d, [50, 2e-4], x_scale=[10, 1e-4])
    return dict(rop0=float(s.x[0]), k_per_km=float(s.x[1] * FT_PER_KM), rmse_days=float(np.sqrt(np.mean(s.fun ** 2))))


def fit_onbottom(depth_ft=None, days=None, bit_life_h: float = 50, trip_speed: float = 1000, bha_h: float = 6,
                 profile: RopProfile = BASE) -> dict:
    """Fit k of the on-bottom ROP so that drilling + WCS tripping time reproduces the shape of the curve.

    A free level factor A absorbs flat time (casing, BOP, circulating) that the shape test cannot see.
    """
    from scipy.optimize import least_squares
    z, d = FERVO_DEPTH_DAYS if depth_ft is None else (np.asarray(depth_ft, float), np.asarray(days, float))
    g = np.linspace(0, z.max(), 4001)
    dz = g[1] - g[0]

    def hours(k):
        r = np.maximum(profile.rop_ref * np.exp(-k / FT_PER_KM * (g - profile.tvd_ref_ft)), profile.rop_min)
        rate = 1 / r + (2 * g / trip_speed + bha_h) / (bit_life_h * r)
        return np.interp(z, g, np.concatenate([[0], np.cumsum(rate[:-1] * dz)]))

    s = least_squares(lambda p: p[0] * hours(p[1]) / 24 - d, [1.5, 0.5], bounds=([0.1, 0], [50, 3]))
    return dict(level=float(s.x[0]), k_per_km=float(s.x[1]), rmse_days=float(np.sqrt(np.mean(s.fun ** 2))))
