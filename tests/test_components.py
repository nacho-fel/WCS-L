"""Experience, learning curves, trouble time and spread rate."""
import numpy as np
import pytest

import wcsl
from wcsl.npt import NPTModel


def test_experience_continuous_campaign_is_close_to_well_number():
    E = wcsl.experience(np.arange(5) * 0.0)
    assert np.allclose(E, [1, 2, 3, 4, 5])


def test_experience_worked_example_from_formula_guide():
    E = wcsl.experience([0, 1 / 12, 2 / 12, 3 / 12, 4 / 12])
    assert E[3] == pytest.approx(3.62, abs=0.01)
    assert E[4] == pytest.approx(4.38, abs=0.01)
    E_pause = wcsl.experience([0, 1 / 12, 2 / 12, 14 / 12])
    assert E_pause[3] == pytest.approx(2.24, abs=0.01)


def test_experience_discounts_multiply():
    E = wcsl.experience([0, 1], rig=["a", "b"], project=["x", "y"])
    assert E[1] == pytest.approx(1 + 0.44 * 0.4 * 0.5)


def test_experience_requires_order():
    with pytest.raises(ValueError):
        wcsl.experience([1, 0])


def test_steady_state_ceiling():
    assert wcsl.steady_state_ceiling(1 / 12) == pytest.approx(15.2, abs=0.1)
    assert wcsl.steady_state_ceiling(1.0) == pytest.approx(1.79, abs=0.01)


def test_bm_factor_starts_at_one_and_never_goes_below_floor():
    f = wcsl.bm_factor(np.linspace(1, 200, 400), 0.63, 1.0)
    assert f[0] == pytest.approx(1.0)
    assert f.min() >= 0.63 - 1e-12
    assert np.all(np.diff(f) <= 0)


def test_brett_millheim_fit_recovers_parameters():
    n = np.arange(1, 9, dtype=float)
    y = wcsl.brett_millheim(n, 1.6, 0.9, 0.8)
    fit = wcsl.fit_brett_millheim(n, y)
    assert fit.params["tau_inf"] == pytest.approx(1.6, rel=1e-3)
    assert fit.params["k"] == pytest.approx(0.8, rel=1e-3)
    assert fit.r2 == pytest.approx(1.0)


def test_wright_fit_and_learning_rate():
    n = np.arange(1, 9, dtype=float)
    fit = wcsl.fit_wright(n, wcsl.wright(n, 2.5, 0.3))
    assert fit.params["b"] == pytest.approx(0.3, rel=1e-6)
    assert wcsl.learning_rate(0.273) == pytest.approx(0.173, abs=0.001)


def test_npt_mean_and_decay():
    m = NPTModel(lambda0=1.5, k_lambda=0.3, median_days=1.5, sigma=0.8)
    rng = np.random.default_rng(0)
    s = m.sample(np.array([1.0, 5.0]), rng, 200000)
    assert s.mean(axis=0) == pytest.approx(m.expected_days(np.array([1.0, 5.0])), rel=0.03)
    assert m.rate(5.0) < m.rate(1.0)
    assert (s >= 0).all()


def test_spread_rate_tends_to_long_run_value():
    s = wcsl.SpreadRate(170000, 0.25, 0.5)
    assert s.rate(1.0) == pytest.approx(212500)
    assert s.rate(100.0) == pytest.approx(170000, rel=1e-6)
