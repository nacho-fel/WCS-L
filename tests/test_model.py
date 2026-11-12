"""WCS-L campaign model, calibration and validation helpers."""
import numpy as np
import pytest

import wcsl
from wcsl import validation as V
from wcsl import wcs


@pytest.fixture(scope="module")
def cal():
    return wcsl.calibrate(theta=0.75)


def test_calibration_reproduces_the_cape_path(cal):
    ex = wcsl.WCSL(cal.params).expected(cal.plans())
    obs = [w["cost"] for w in cal.wells]
    assert V.mape([w["cost"] for w in ex["wells"]], obs) < 0.08
    assert V.mape([w["days"] for w in ex["wells"]], [w["days"] for w in cal.wells]) < 0.06


def test_calibrated_time_share_matches_theta(cal):
    ex = wcsl.WCSL(cal.params).expected(cal.plans())
    time = sum(w["time_cost"] + w["skid"] for w in ex["wells"])
    assert time / sum(w["cost"] for w in ex["wells"]) == pytest.approx(0.75, abs=0.02)


def test_calibrated_parameters_are_physical(cal):
    p = cal.params
    assert 0.4 < p.floor < 0.9 and p.k > 0
    assert 1.0 < p.time_scale < 2.5          # the WCS misses flat time
    assert 0.2 < cal.materials_mult < 1.0     # Fervo materials below the GeoVision design
    assert 100e3 < p.spread.s_inf < 250e3


def test_costs_fall_then_flatten_and_respect_floor(cal):
    plans = wcsl.uniform_campaign(30, lambda: wcs.cape_like_well(md_ft=13500), wells_per_pad=30)
    days = np.array([w["days"] for w in wcsl.WCSL(cal.params).expected(plans)["wells"]])
    assert days[0] > days[5] > 0
    assert abs(days[-1] - days[-2]) < 0.05
    planned_floor = days[0] * cal.params.floor * 0.9
    assert days.min() > planned_floor


def test_pause_in_drilling_costs_money(cal):
    m = wcsl.WCSL(cal.params)
    w = lambda: wcs.cape_like_well(md_ft=13500)  # noqa: E731
    fast = m.expected(wcsl.uniform_campaign(8, w, months_between=1))["total"]
    slow = m.expected(wcsl.uniform_campaign(8, w, months_between=12))["total"]
    assert slow > fast


def test_campaign_total_is_sum_of_wells_plus_fixed(cal):
    plans = wcsl.uniform_campaign(9, lambda: wcs.cape_like_well(), wells_per_pad=4)
    ex = wcsl.WCSL(cal.params).expected(plans)
    assert ex["total"] == pytest.approx(sum(w["cost"] for w in ex["wells"]) + sum(ex["campaign_fixed"].values()))
    assert ex["campaign_fixed"]["sites"] == 3 * cal.params.site


def test_monte_carlo_is_reproducible_and_ordered(cal):
    plans = cal.plans()
    m = wcsl.WCSL(cal.params)
    a = m.simulate(plans, 2000, seed=3, uncertainty=cal.uncertainty())
    b = m.simulate(plans, 2000, seed=3, uncertainty=cal.uncertainty())
    assert np.array_equal(a["cost"], b["cost"])
    P = wcsl.percentiles(a["cost"])
    assert np.all(P[0] < P[1]) and np.all(P[1] < P[2])
    ex = m.expected(plans)
    assert np.mean(a["cost"][:, 0]) == pytest.approx(ex["wells"][0]["cost"], rel=0.05)


def test_technology_step_lowers_cost(cal):
    p = cal.params
    p2 = wcsl.WCSLParams(**{**p.__dict__, "generation_floor": {1: 0.8}})
    plans = wcsl.uniform_campaign(6, lambda: wcs.cape_like_well())
    plans_new = [wcsl.WellPlan(**{**w.__dict__, "generation": 1}) for w in plans]
    assert wcsl.WCSL(p2).expected(plans_new)["total"] < wcsl.WCSL(p).expected(plans)["total"]


def test_metrics():
    assert V.mape([110, 90], [100, 100]) == pytest.approx(0.10)
    assert V.coverage([0, 0], [1, 1], [0.5, 2]) == 0.5
    pt = np.full((10, 1), 5.0)
    assert V.crps_samples(pt, [3.0])[0] == pytest.approx(2.0)     # point forecast: CRPS = absolute error


def test_structure_test_wcs_split_cannot_reproduce_cape():
    W = wcsl.load_fervo()
    obs = np.array([w["cost"] for w in W])
    wcs_path = V.anchored_path(W, V.WCS_THETA)
    assert wcs_path[-1] / wcs_path[0] - 1 > -0.10       # the WCS split predicts less than 10 % fall
    assert obs[-1] / obs[0] - 1 < -0.40                  # the wells fell 42 %


def test_design_evaluation_returns_sorted_usd_per_mw(cal):
    opts = [wcsl.DesignOption("a", lambda: wcs.cape_like_well(), 5.0),
            wcsl.DesignOption("b", lambda: wcs.horizontal_well(11000, 5000, 10, 75, 50, **wcs.CAPE_LIKE_PARAMS), 7.0)]
    res = wcsl.evaluate_designs(opts, 50, cal.params, n_draws=500)
    assert res[0]["usd_per_mw"] <= res[1]["usd_per_mw"]
    assert all(r["usd_per_mw_p10"] <= r["usd_per_mw_p50"] <= r["usd_per_mw_p90"] for r in res)
