"""Learning curves and their fitting.

Brett & Millheim (1986): time per unit of work decays exponentially toward a technical limit,

    tau(E) = tau_inf + (tau_1 - tau_inf) * exp(-k (E - 1))

Wright's law (power law, no floor) and a quadratic are provided as comparison models.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Callable

import numpy as np
from scipy.optimize import curve_fit


def bm_factor(E, floor: float, k: float):
    """Time relative to the first well: floor + (1 - floor) exp(-k (E - 1)); 1 at E = 1."""
    E = np.asarray(E, dtype=float)
    return floor + (1.0 - floor) * np.exp(-k * (E - 1.0))


def brett_millheim(E, tau_inf: float, excess: float, k: float):
    return tau_inf + excess * np.exp(-k * (np.asarray(E, dtype=float) - 1.0))


def wright(n, t1: float, b: float):
    return t1 * np.asarray(n, dtype=float) ** (-b)


def quadratic(n, a: float, b: float, c: float):
    n = np.asarray(n, dtype=float)
    return a + b * n + c * n * n


def learning_rate(b: float) -> float:
    """Wright learning rate: fractional fall per doubling, LR = 1 - 2^-b."""
    return 1.0 - 2.0 ** (-b)


@dataclass
class Fit:
    name: str
    params: dict[str, float]
    cov: np.ndarray
    r2: float
    rmse: float
    aic: float
    n_obs: int
    func: Callable = field(repr=False, default=None)

    def predict(self, x):
        return self.func(x, *self.params.values())

    def sample_params(self, rng: np.random.Generator, size: int) -> np.ndarray:
        """Parameter draws from the asymptotic normal approximation of the fit."""
        return rng.multivariate_normal(list(self.params.values()), self.cov, size=size, method="eigh")


def _stats(y, yhat, k):
    y = np.asarray(y, float)
    res = y - yhat
    ss_res = float(np.sum(res ** 2))
    ss_tot = float(np.sum((y - y.mean()) ** 2))
    n = len(y)
    r2 = 1 - ss_res / ss_tot if ss_tot > 0 else float("nan")
    rmse = float(np.sqrt(ss_res / n))
    aic = n * np.log(ss_res / n + 1e-300) + 2 * k
    return r2, rmse, float(aic)


def fit_brett_millheim(x, y, floor_bounds=(0.0, np.inf)) -> Fit:
    """Least-squares fit of tau_inf, excess (= tau_1 - tau_inf) and k."""
    x, y = np.asarray(x, float), np.asarray(y, float)
    p0 = [y.min() * 0.95, max(y[0] - y.min(), 1e-3), 0.7]
    lo = [floor_bounds[0], 0.0, 1e-4]
    hi = [floor_bounds[1], np.inf, 20.0]
    p, cov = curve_fit(brett_millheim, x, y, p0=p0, bounds=(lo, hi), maxfev=20000)
    r2, rmse, aic = _stats(y, brett_millheim(x, *p), 3)
    return Fit("brett_millheim", dict(tau_inf=p[0], excess=p[1], k=p[2]), cov, r2, rmse, aic, len(y), brett_millheim)


def fit_wright(x, y) -> Fit:
    """Wright's law fitted in log space (constant percentage errors)."""
    x, y = np.asarray(x, float), np.asarray(y, float)
    A = np.vstack([np.ones_like(x), -np.log(x)]).T
    coef, *_ = np.linalg.lstsq(A, np.log(y), rcond=None)
    resid = np.log(y) - A @ coef
    s2 = float(resid @ resid / max(len(y) - 2, 1))
    cov_log = s2 * np.linalg.inv(A.T @ A)
    t1, b = float(np.exp(coef[0])), float(coef[1])
    J = np.diag([t1, 1.0])                                   # delta method from (log t1, b) to (t1, b)
    r2, rmse, aic = _stats(y, wright(x, t1, b), 2)
    return Fit("wright", dict(t1=t1, b=b), J @ cov_log @ J.T, r2, rmse, aic, len(y), wright)


def fit_quadratic(x, y) -> Fit:
    x, y = np.asarray(x, float), np.asarray(y, float)
    p, cov = np.polyfit(x, y, 2, cov=True) if len(x) > 3 else (np.polyfit(x, y, 2), np.zeros((3, 3)))
    c, b, a = p
    r2, rmse, aic = _stats(y, quadratic(x, a, b, c), 3)
    return Fit("quadratic", dict(a=a, b=b, c=c), cov[::-1, ::-1], r2, rmse, aic, len(y), quadratic)


def technical_limit_hours_per_kft(best_rop_ft_h: float) -> float:
    """Floor of on-bottom drilling time from the best sustained ROP: 1,000 / ROP hours per 1,000 ft."""
    return 1000.0 / best_rop_ft_h
