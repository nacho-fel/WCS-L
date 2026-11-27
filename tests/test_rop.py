import numpy as np
from wcsl import wcs, rop, calibrate_cape17


def test_constant_profile_reproduces_wcs():
    a = wcs.well_cost(wcs.cape_like_well(md_ft=13500))
    b = wcs.well_cost(wcs.cape_like_well(md_ft=13500, rop_profile=rop.CONSTANT))
    assert abs(a["total"] - b["total"]) < 1.0


def test_lateral_at_reference_depth_drills_at_reference_rop():
    s = wcs.egs_well(8700, 5000, rop.BASE)
    assert abs(s["intervals"][-1]["rop"] - 75.0) < 0.1


def test_deeper_is_slower_and_costlier():
    c = [wcs.well_cost(wcs.egs_well(km * rop.FT_PER_KM, 5000, rop.BASE))["total"] for km in (3, 4, 5)]
    assert c[0] < c[1] < c[2] and (c[2] - c[1]) > (c[1] - c[0])      # increasing and convex in depth


def test_fit_matches_workbook():
    assert abs(rop.fit_allin()["k_per_km"] - 0.709) < 0.01
    assert abs(rop.fit_onbottom()["k_per_km"] - 0.609) < 0.01


def test_band_calibration_parameters():
    cal = calibrate_cape17()
    assert cal.params.npt.lambda_inf == 0.35 and cal.time_well_logsd == 0.11
    assert 1.0e5 < cal.params.spread.s_inf < 1.6e5
