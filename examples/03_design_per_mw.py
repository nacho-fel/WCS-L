"""Method K: pick the well design by $ per MW for a 100 MW demand.

The MW per well below are PLACEHOLDERS. Replace them with GEOPHIRES results for each geometry.
"""
import sys
import warnings
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
warnings.simplefilter("ignore")
import wcsl  # noqa: E402
from wcsl import wcs  # noqa: E402

cal = wcsl.calibrate(theta=0.75)
options = [
    wcsl.DesignOption("Cape-like: 8,700 ft TVD, 5,000 ft lateral", lambda: wcs.cape_like_well(), mw_per_well=5.0),
    wcsl.DesignOption("Longer lateral: 8,700 ft TVD, 7,500 ft lateral",
                      lambda: wcs.horizontal_well(8700, 7500, 10, 75, 50, **wcs.CAPE_LIKE_PARAMS), mw_per_well=6.5),
    wcsl.DesignOption("Deeper: 11,000 ft TVD, 5,000 ft lateral",
                      lambda: wcs.horizontal_well(11000, 5000, 10, 75, 50, **wcs.CAPE_LIKE_PARAMS), mw_per_well=7.0),
]
for r in wcsl.evaluate_designs(options, demand_mw=100, params=cal.params, uncertainty=cal.uncertainty()):
    print(f"{r['design']:48s} {r['wells']:3d} wells  ${r['usd_per_mw'] / 1e6:5.2f}M/MW "
          f"(P10 {r['usd_per_mw_p10'] / 1e6:.2f} - P90 {r['usd_per_mw_p90'] / 1e6:.2f}), ${r['mean_cost_per_well'] / 1e6:.2f}M per well")
