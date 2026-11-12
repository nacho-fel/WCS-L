"""The baseline WCS engine reproduces the Sandia spreadsheets."""
import gzip
import json
from pathlib import Path

import pytest

from wcsl import wcs

SPECS = json.load(gzip.open(Path(__file__).resolve().parents[1] / "data/wcs/geovision_specs.json.gz", "rt"))


def test_all_geovision_workbooks_reproduced_to_the_dollar():
    assert len(SPECS) == 249
    for s in SPECS:
        assert wcs.well_cost(s["spec"])["total"] == pytest.approx(s["expected_total"], abs=0.5), s["file"]


def test_known_values():
    assert wcs.well_cost(wcs.vertical_well(5000))["total"] == pytest.approx(16.40e6, rel=1e-3)
    r = wcs.well_cost(wcs.cape_like_well())
    assert r["total"] == pytest.approx(4.571e6, rel=1e-3)
    assert r["days"] == pytest.approx(24.8, abs=0.1)


def test_nrel_2025_update_is_about_minus_15_percent():
    base = wcs.well_cost(wcs.vertical_well(5000))["total"]
    new = wcs.well_cost(wcs.with_params(wcs.vertical_well(5000), **wcs.NREL_2025))["total"]
    assert new / base - 1 == pytest.approx(-0.154, abs=0.002)


def test_cost_split_sums_to_one_and_wcs_time_share_is_low():
    s = wcs.cost_split(wcs.well_cost(wcs.cape_like_well()))
    assert sum(s.values()) == pytest.approx(1.0)
    assert 0.18 < s["time"] < 0.24


def test_cost_falls_with_rop_and_rises_with_depth():
    c = [wcs.well_cost(wcs.with_params(wcs.cape_like_well(), rop=r))["total"] for r in (25, 75, 150)]
    assert c[0] > c[1] > c[2]
    assert wcs.well_cost(wcs.horizontal_well(10000))["total"] > wcs.well_cost(wcs.horizontal_well(8000))["total"]


def test_lateral_for_md_round_trip():
    spec = wcs.cape_like_well(md_ft=13500)
    assert spec["intervals"][-1]["end"] == pytest.approx(13500)
