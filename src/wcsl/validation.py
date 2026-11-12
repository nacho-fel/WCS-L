"""Validation: metrics, out-of-sample procedures and benchmark models.

Point metrics compare expected values with observations; probabilistic metrics check whether the
Monte Carlo bands are honest (coverage of the P10-P90 band, continuous ranked probability score).
Out-of-sample procedures always calibrate on earlier wells and predict later ones, the way the
model would be used to plan a campaign.
"""
from __future__ import annotations

from typing import Any, Sequence

import numpy as np

from . import wcs
from .calibration import Calibration, calibrate
from .learning import bm_factor, fit_brett_millheim
from .model import WCSL, percentiles

NREL_CAPE = dict(wcs.CAPE_LIKE_PARAMS)


# ------------------------------------------------------------------------------------- metrics
def ape(pred, obs) -> np.ndarray:
    pred, obs = np.asarray(pred, float), np.asarray(obs, float)
    return np.abs(pred - obs) / np.abs(obs)


def mape(pred, obs) -> float:
    return float(np.mean(ape(pred, obs)))


def bias(pred, obs) -> float:
    """Mean signed error as a fraction of the observation (+ = over-prediction)."""
    pred, obs = np.asarray(pred, float), np.asarray(obs, float)
    return float(np.mean((pred - obs) / obs))


def coverage(lo, hi, obs) -> float:
    obs = np.asarray(obs, float)
    return float(np.mean((obs >= np.asarray(lo)) & (obs <= np.asarray(hi))))


def crps_samples(samples: np.ndarray, obs) -> np.ndarray:
    """CRPS of each column of ``samples`` (draws x wells) against ``obs`` (lower is better).

    For a single number (no uncertainty) the CRPS equals the absolute error, so deterministic
    and probabilistic forecasts can be compared on the same scale.
    """
    X = np.asarray(samples, float)
    y = np.asarray(obs, float)
    t1 = np.mean(np.abs(X - y[None, :]), axis=0)
    Xs = np.sort(X, axis=0)
    n = X.shape[0]
    i = np.arange(1, n + 1)[:, None]
    t2 = np.sum((2 * i - n - 1) * Xs, axis=0) / (n * n)      # = 0.5 E|X - X'|
    return t1 - t2


def pit(samples: np.ndarray, obs) -> np.ndarray:
    """Probability integral transform: share of draws below each observation (uniform if calibrated)."""
    return np.mean(np.asarray(samples) <= np.asarray(obs)[None, :], axis=0)


# ------------------------------------------------------------------------------------- predictions
def predict(cal: Calibration, wells: Sequence[dict[str, Any]], n_draws: int = 4000, seed: int = 7) -> dict[str, Any]:
    """Expected values and Monte Carlo bands for ``wells`` with a calibrated model.

    The campaign is simulated from well 1 up to the last well asked for, so experience is counted
    correctly; only the requested wells are returned.
    """
    n_max = max(w["n"] for w in wells)
    td_default = float(np.mean([w["td_ft"] for w in cal.wells]))
    by_n = {w["n"]: w for w in list(cal.wells) + list(wells)}
    full = [dict(n=i, td_ft=by_n.get(i, {}).get("td_ft", td_default), well=by_n.get(i, {}).get("well", f"n{i}"))
            for i in range(1, n_max + 1)]
    plans = cal.plans(full)
    model = WCSL(cal.params)
    ex = model.expected(plans)
    sim = model.simulate(plans, n_draws=n_draws, seed=seed, uncertainty=cal.uncertainty())
    idx = [w["n"] - 1 for w in wells]
    return dict(days=np.array([ex["wells"][i]["days"] for i in idx]),
                cost=np.array([ex["wells"][i]["cost"] for i in idx]),
                time_share=np.array([ex["wells"][i]["time_share"] for i in idx]),
                days_samples=sim["days"][:, idx], cost_samples=sim["cost"][:, idx],
                days_pct=percentiles(sim["days"][:, idx]), cost_pct=percentiles(sim["cost"][:, idx]))


def wcs_baseline(wells: Sequence[dict[str, Any]]) -> dict[str, np.ndarray]:
    """The WCS with NREL 2025 parameters: the same cost and days for every well of a given depth."""
    res = [wcs.well_cost(wcs.cape_like_well(md_ft=w["td_ft"])) for w in wells]
    return dict(cost=np.array([r["total"] for r in res]),
                days=np.array([sum(r["hours"].values()) / 24 for r in res]))


def learning_on_structure(train: Sequence[dict[str, Any]], wells: Sequence[dict[str, Any]],
                          theta: float) -> dict[str, np.ndarray]:
    """A one-line cost model with a given cost structure: cost per 1,000 ft = A [(1 - theta) + theta f(n)].

    f(n) is the Brett-Millheim factor fitted to the training wells' days per 1,000 ft, and A is fitted by
    least squares to the training costs. With theta = 0.21 (the WCS split) this is "the WCS plus a decay
    rate" recalibrated to the operator's own level; with theta = 0.75 it has Fervo's split. Comparing the
    two isolates the effect of the cost structure, everything else being equal.
    """
    n_tr = np.array([w["n"] for w in train], float)
    y = np.array([w["days"] / (w["td_ft"] / 1000) for w in train])
    bm = fit_brett_millheim(n_tr, y)
    t1 = bm.params["tau_inf"] + bm.params["excess"]
    floor, k = bm.params["tau_inf"] / t1, bm.params["k"]
    g_tr = np.array([(1 - theta) + theta * float(bm_factor(w["n"], floor, k)) for w in train])
    c_tr = np.array([w["cost"] / (w["td_ft"] / 1000) for w in train])
    A = float(c_tr @ g_tr / (g_tr @ g_tr))
    g = np.array([(1 - theta) + theta * float(bm_factor(w["n"], floor, k)) for w in wells])
    td = np.array([w["td_ft"] / 1000 for w in wells])
    days = np.array([t1 * float(bm_factor(w["n"], floor, k)) for w in wells]) * td
    return dict(cost=A * g * td, days=days)


WCS_THETA = 0.21   # time-driven share of the WCS Cape-like well (studies/time_share)


def persistence(train: Sequence[dict[str, Any]], wells: Sequence[dict[str, Any]]) -> dict[str, np.ndarray]:
    """Naive forecast: the last training well's cost and days per 1,000 ft, scaled to the target depth."""
    last = max(train, key=lambda w: w["n"])
    c_kft, d_kft = last["cost"] / last["td_ft"], last["days"] / last["td_ft"]
    return dict(cost=np.array([c_kft * w["td_ft"] for w in wells]), days=np.array([d_kft * w["td_ft"] for w in wells]))


# ------------------------------------------------------------------------------------- procedures
def rolling_origin(wells: Sequence[dict[str, Any]], first_train: int = 3, theta: float = 0.75,
                   n_draws: int = 4000) -> list[dict[str, Any]]:
    """Calibrate on wells 1..k, predict well k+1, for k = first_train .. N-1 (one-step-ahead)."""
    wells = sorted(wells, key=lambda w: w["n"])
    out = []
    for k in range(first_train, len(wells)):
        train, test = wells[:k], [wells[k]]
        cal = calibrate(train, theta=theta, k_s=0.5 if k < 5 else None)
        pr = predict(cal, test, n_draws=n_draws)
        out.append(dict(n_train=k, target=test[0], wcsl=pr, wcs=wcs_baseline(test),
                        wcs_learning=learning_on_structure(train, test, WCS_THETA),
                        fervo_structure=learning_on_structure(train, test, theta),
                        persistence=persistence(train, test)))
    return out


def holdout(wells: Sequence[dict[str, Any]], n_train: int = 4, theta: float = 0.75, n_draws: int = 4000) -> dict[str, Any]:
    wells = sorted(wells, key=lambda w: w["n"])
    train, test = wells[:n_train], wells[n_train:]
    cal = calibrate(train, theta=theta, k_s=0.5 if n_train < 5 else None)
    return dict(train=train, test=test, cal=cal, wcsl=predict(cal, test, n_draws=n_draws), wcs=wcs_baseline(test),
                wcs_learning=learning_on_structure(train, test, WCS_THETA),
                fervo_structure=learning_on_structure(train, test, theta), persistence=persistence(train, test))


def anchored_path(wells: Sequence[dict[str, Any]], theta: float) -> np.ndarray:
    """Structure test: anchor a model on the FIRST well's cost only and let learning (fitted on days,
    which every model shares) carry it forward. Returns predicted costs for all wells."""
    wells = sorted(wells, key=lambda w: w["n"])
    n = np.array([w["n"] for w in wells], float)
    y = np.array([w["days"] / (w["td_ft"] / 1000) for w in wells])
    bm = fit_brett_millheim(n, y)
    t1 = bm.params["tau_inf"] + bm.params["excess"]
    floor, k = bm.params["tau_inf"] / t1, bm.params["k"]
    c1 = wells[0]["cost"] / (wells[0]["td_ft"] / 1000)
    return np.array([c1 * ((1 - theta) + theta * float(bm_factor(w["n"], floor, k))) * w["td_ft"] / 1000 for w in wells])
