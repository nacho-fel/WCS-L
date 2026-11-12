"""A 10-well campaign with the Cape-calibrated WCS-L: continuity, pads and rig changes."""
import sys
import warnings
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
warnings.simplefilter("ignore")
import numpy as np  # noqa: E402

import wcsl  # noqa: E402
from wcsl import wcs  # noqa: E402

cal = wcsl.calibrate(theta=0.75)                     # Fervo Cape Station, time share 75 %
model = wcsl.WCSL(cal.params)
well = lambda: wcs.cape_like_well(md_ft=13500)        # noqa: E731

cases = {
    "monthly, one rig, 4 wells per pad": wcsl.uniform_campaign(10, well, months_between=1, wells_per_pad=4),
    "one well every 6 months": wcsl.uniform_campaign(10, well, months_between=6, wells_per_pad=4),
    "monthly, new rig at every pad": [wcsl.WellPlan(spec=well(), t_years=i / 12, rig=f"rig-{i // 4}", pad=i // 4)
                                      for i in range(10)],
}
for name, plans in cases.items():
    ex = model.expected(plans)
    sim = model.simulate(plans, n_draws=4000, uncertainty=cal.uncertainty())
    p10, p50, p90 = np.percentile(sim["total"], [10, 50, 90]) / 1e6
    costs = " ".join(f"{w['cost'] / 1e6:4.1f}" for w in ex["wells"])
    print(f"\n{name}\n  cost per well ($M): {costs}\n  campaign ${ex['total'] / 1e6:.1f}M "
          f"(P10 {p10:.1f}, P50 {p50:.1f}, P90 {p90:.1f}); mean per well ${ex['mean_cost_per_well'] / 1e6:.2f}M")
