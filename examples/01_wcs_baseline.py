"""The baseline WCS engine: one well, its cost lines and its cost structure."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from wcsl import wcs  # noqa: E402

for name, spec in {
    "Vertical 5,000 m, GeoVision base (2016)": wcs.vertical_well(5000),
    "Vertical 5,000 m, NREL 2025 parameters": wcs.with_params(wcs.vertical_well(5000), **wcs.NREL_2025),
    "Cape-like horizontal (8,700 ft TVD, 5,000 ft lateral), NREL 2025": wcs.cape_like_well(),
}.items():
    r = wcs.well_cost(spec)
    s = wcs.cost_split(r)
    print(f"\n{name}\n  total ${r['total'] / 1e6:.2f}M, {r['days']:.1f} days, measured depth {r['md']:,.0f} ft")
    print(f"  time-driven {s['time']:.0%}, materials {s['materials']:.0%}, fixed {s['fixed']:.0%}")
    for k, v in sorted(r["lines"].items(), key=lambda kv: -kv[1]):
        if v > 0:
            print(f"    {k:10s} ${v / 1e3:>9,.0f}k")
