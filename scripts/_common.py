"""Shared plotting style and helpers for the report scripts."""
from __future__ import annotations

import json
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

INK, MUTED, RULE = "#1a1f22", "#5b6468", "#d6dad8"
STEEL, AMBER, ROCK, GOOD, WARN = "#3f5f73", "#b5651d", "#b8a888", "#2f7d55", "#b0413e"
GREY = "#9aa3a6"

plt.rcParams.update({
    "figure.dpi": 150, "savefig.dpi": 200, "savefig.bbox": "tight", "font.size": 10,
    "axes.spines.top": False, "axes.spines.right": False, "axes.edgecolor": MUTED, "axes.labelcolor": INK,
    "xtick.color": MUTED, "ytick.color": MUTED, "axes.grid": True, "grid.color": RULE, "grid.linewidth": 0.6,
    "legend.frameon": False, "axes.axisbelow": True, "axes.titleweight": "bold", "axes.titlesize": 11,
})


def save(fig, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path)
    plt.close(fig)


def dump(obj, path: Path) -> None:
    def conv(o):
        if isinstance(o, (np.floating, np.integer)):
            return o.item()
        if isinstance(o, np.ndarray):
            return o.tolist()
        raise TypeError(type(o))
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(obj, indent=2, default=conv))


def pct(x: float, d: int = 0) -> str:
    return f"{100 * x:.{d}f} %"
