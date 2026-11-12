"""Baseline engine: the Sandia Well Cost Simplified (WCS) model.

A line-by-line re-implementation of the GeoVision WCS spreadsheet
(Lowry et al. 2017, SAND2017-9977, https://doi.org/10.2172/1394062).
It reproduces all 249 GeoVision workbooks to the dollar
(see ``tests/test_wcs_reproduction.py``).

A well is a *spec*: a dict with ``params`` (prices, rates, fixed costs) and a
list of ``intervals`` (one per casing section). Every interval key has a
default, so a new well only needs geometry, ROP and bit life. The optional
keys reproduce the formula variants found in the original workbooks.

Units: feet, inches, hours, US dollars.
"""
from __future__ import annotations

import copy
import math
from typing import Any

#: Default Values sheet of the GeoVision workbooks (2016 prices).
DEFAULT_PARAMS: dict[str, float] = dict(
    bha=0.5,            # BHA cost as a fraction of bit cost
    csg_price=1.5,      # $ per lb of casing steel
    cem_price=125.0,    # cement price factor ($ per ft3 equivalent in the WCS formula)
    run=300.0,          # casing running speed, ft/h
    cem_h=5.0,          # fixed cementing time per string, h
    trip_speed=1000.0,  # tripping speed, ft/h
    bha_h=6.0,          # BHA handling time per trip, h
    log_speed=300.0,    # logging speed, ft/h
    log_svc=50000.0,    # logging service charge per run, $
    wh_h=40.0,          # wellhead installation time, h
    mud0=10000.0,       # initial mud cost, $
    mud_day=4000.0,     # mud cost per day, $
    dir_eng=75.0,       # directional engineering, $/h (all hours of a directional interval)
    dir_tool=500.0,     # directional tools, $/h (drilling hours)
    extra=0.0,          # extra cost per interval, $
    bit_a=3775.0,       # bit price = bit_a * diameter + bit_b
    bit_b=-21750.0,
    bit_mult=1.0,       # bit price multiplier (2.0 in NREL 2025)
    rig=38400.0,        # rig rate, $/day
    cont=0.15,          # contingency, fraction of interval cost
    mob=250000.0, site=250000.0, eng=10000.0, wellhead=80000.0,
    mob_days=3.0, site_days=2.0,
)

#: GeoVision technology scenarios (ROP ft/h, bit life h, build rate deg/100 ft, contingency, rig $/day).
SCENARIOS: dict[str, dict[str, float]] = {
    "base":  dict(rop=25, bit=50, bur=3, cont=0.15, rig=38400),
    "int1":  dict(rop=50, bit=100, bur=10, cont=0.10, rig=33500),
    "int2":  dict(rop=75, bit=150, bur=20, cont=0.05, rig=28500),
    "ideal": dict(rop=100, bit=200, bur=30, cont=0.00, rig=28500),
}

#: NREL 2025 update (Akindipe & Witter 2025): parameter changes on top of the GeoVision base.
NREL_2025: dict[str, float] = dict(rop=75, bit_mult=2.0, run=800.0, cont=0.10)

TIME_LINES = ("drill", "trip", "cemt", "logt", "wht", "trouble")   # charged at the rig rate
TIME_LINKED = ("mud", "dir")                                       # also scale with hours
MATERIAL_LINES = ("bits", "bha", "casing", "cement", "logs", "trouble_c", "extra")


def _iv(g: dict[str, Any], key: str, default: Any) -> Any:
    v = g.get(key, default)
    return default if v is None else v


def well_cost(spec: dict[str, Any]) -> dict[str, Any]:
    """Cost of one well with the WCS formulas.

    Returns a dict with ``total``, ``interval_cost``, ``cont``, ``fixed``, ``lines``
    (cost by line item, $), ``hours`` (by activity), ``days`` and per-interval detail.
    """
    p = dict(DEFAULT_PARAMS)
    p.update(spec.get("params", {}))
    h = p["rig"] / 24.0
    lines = {k: 0.0 for k in TIME_LINES + TIME_LINKED + MATERIAL_LINES}
    hours = {k: 0.0 for k in TIME_LINES}
    detail = []
    top = 0.0
    for i, g in enumerate(spec["intervals"]):
        end = float(g["end"])
        D, Dc, w = float(g["D"]), float(_iv(g, "Dc", 0.0)), float(_iv(g, "w", 0.0))
        CL, rop, bl = int(_iv(g, "CL", 1)), float(g["rop"]), float(g["bit"])
        m = {k: 1.0 for k in ("bit", "bha", "csg", "cem", "mud", "dir", "logt")}
        m.update(_iv(g, "mult", {}))
        n = 1.0 if i == 0 else float(_iv(g, "n", 1.0))
        oh = 1.0 if i == 0 else float(_iv(g, "oh", 0.0 if w == 0 else 1.0))
        # drilling time
        if i == 0:
            drill = end / rop if rop else 0.0
        else:
            mode = _iv(g, "drill_mode", "plain")
            if mode == "plain":
                drill = (end - top) / rop if rop else 0.0
            elif mode == "switch":
                drill = (end - top) / rop * float(_iv(g, "drill_switch", 1))
            else:  # "eps"
                drill = (end - top) / (rop + float(_iv(g, "drill_eps", 0.0))) * n
        # tripping time
        if i == 0:
            trip = 2 * end / p["trip_speed"] + p["bha_h"]
        else:
            nb = drill / (bl + float(_iv(g, "trip_eps", 0.01)))
            trip = ((nb * (end + top) if nb > 1 else end) / p["trip_speed"]
                    + (nb if nb > 1 else 1) * p["bha_h"]
                    + float(_iv(g, "trip_mult", 1.0)) * (top + end) / p["trip_speed"]) * n
        # bits and BHA
        a = float(_iv(g, "bit_a", p["bit_a"]))
        b = float(_iv(g, "bit_b", p["bit_b"]))
        be = float(_iv(g, "bit_eps", 0.0))
        bits = (drill / (bl + be) if (bl + be) else 0.0) * (a * D + b) * p["bit_mult"] * m["bit"] \
            * (n if _iv(g, "bit_n", False) else 1.0)
        bha = bits * p["bha"] * m["bha"]
        # casing and cement
        L = end if (i == 0 or CL == 1) else end - top
        casing = w * p["csg_price"] * m["csg"] * oh * L * (n if _iv(g, "csg_n", False) else 1.0)
        pi_ = math.pi if _iv(g, "cem_pi", False) else 3.14
        cement = (pi_ / 576) * (D ** 2 - Dc ** 2) * L * p["cem_price"] * m["cem"] * oh \
            * (n if _iv(g, "cem_n", False) else 1.0)
        pump = 0.0
        if _iv(g, "pump", False) and cement > 0:
            vol = 0.25 * math.pi / 144 / 5.615 * (D ** 2 - Dc ** 2) * end
            pump = vol / float(_iv(g, "pump_rate", 5.0)) / 60
        cemt = (end / p["run"] + p["cem_h"] + pump) * (oh if _iv(g, "cemt_oh", True) else 1.0)
        # logging
        logt = (2 + end / p["log_speed"] * m["logt"]) if _iv(g, "log", 0) else 0.0
        logt_time = logt * (n if _iv(g, "logt_n", False) else 1.0)
        logt_cost = logt_time * h * (n if _iv(g, "logc_n", False) else 1.0)
        logs = 0.0 if logt_time == 0 else p["log_svc"]
        # wellhead time
        wht = (p["wh_h"] if CL == 1 else 0.0) * (oh if _iv(g, "wh_oh", True) else 1.0) \
            * (n if _iv(g, "wh_n", False) else 1.0) * float(_iv(g, "wh_switch", 1))
        trouble = float(_iv(g, "trouble_h", 0.0))
        trouble_c = float(_iv(g, "trouble_c", 0.0))
        T = drill + trip + cemt + trouble + logt_time + wht
        # directional and mud
        dirc = 0.0
        if _iv(g, "dir", 0):
            tool_h = _iv(g, "dir_hours", None)
            dirc = (p["dir_eng"] * T + p["dir_tool"] * (drill if tool_h is None else float(tool_h))) * m["dir"]
        mud = (p["mud0"] * m["mud"] if i == 0 else 0.0) + p["mud_day"] / 24 * T * m["mud"]
        row = dict(drill=drill * h, trip=trip * h, cemt=cemt * h, logt=logt_cost, wht=wht * h,
                   trouble=trouble * h, mud=mud, dir=dirc, bits=bits, bha=bha, casing=casing,
                   cement=cement, logs=logs, trouble_c=trouble_c, extra=p["extra"])
        for k, v in row.items():
            lines[k] += v
        for k, v in dict(drill=drill, trip=trip, cemt=cemt, logt=logt_time, wht=wht, trouble=trouble).items():
            hours[k] += v
        detail.append(dict(cost=sum(row.values()), hours=T, **row))
        top = end
    interval_cost = sum(lines.values())
    cont = interval_cost * p["cont"]
    fixed = {k: p[k] for k in ("mob", "site", "eng", "wellhead")}
    total = interval_cost + cont + sum(fixed.values())
    T = sum(hours.values())
    return dict(total=total, interval_cost=interval_cost, cont=cont, fixed=fixed, lines=lines,
                hours=hours, days=T / 24 + p["mob_days"] + p["site_days"], md=top,
                intervals=detail, params=p)


def cost_split(result: dict[str, Any], contingency_to: str = "pro_rata") -> dict[str, float]:
    """Split a WCS result into time-driven, materials and fixed shares of the total.

    Time-driven = everything charged per hour (rig time, mud per day, directional).
    The contingency is allocated pro rata to the interval lines by default.
    """
    L = result["lines"]
    time = sum(L[k] for k in TIME_LINES + TIME_LINKED)
    mat = sum(L[k] for k in MATERIAL_LINES)
    fixed = sum(result["fixed"].values())
    ic = result["interval_cost"]
    c = result["cont"]
    if contingency_to == "pro_rata" and ic > 0:
        time += c * time / ic
        mat += c * mat / ic
        cont = 0.0
    else:
        cont = c
    tot = result["total"]
    return dict(time=time / tot, materials=mat / tot, fixed=fixed / tot, contingency=cont / tot)


# ----------------------------------------------------------------------------------------------
# Geometry builders (the GeoVision layouts, plus a horizontal well with a free lateral length)
# ----------------------------------------------------------------------------------------------

def vertical_well(target_m: float = 5000, rop: float = 25, bit: float = 50, **params: Any) -> dict[str, Any]:
    """GeoVision vertical well, base scenario, large diameter (5 intervals, GETEM casing rules).

    The production interval is 1,000 ft of open hole below 500 ft past the target depth.
    """
    tgt = target_m / 0.3048
    bot = round(tgt) + 500
    prod_top = bot - 1000
    s = (prod_top - 1000) / 3
    iv = [dict(end=1000, D=36, Dc=30, w=215.4, CL=1, log=0, dir=0),
          dict(end=1000 + s, D=28, Dc=24, w=216.3, CL=1, log=1, dir=0),
          dict(end=1000 + 2 * s, D=22, Dc=18.625, w=131.8, CL=1, log=1, dir=0),
          dict(end=prod_top, D=17, Dc=13.625, w=78.5, CL=1, log=1, dir=0),
          dict(end=bot, D=12.25, Dc=0, w=0, CL=1, log=1, dir=0)]
    for g in iv:
        g.update(rop=rop, bit=bit)
    return dict(params=dict(params), intervals=iv)


def horizontal_well(tvd_ft: float = 8700, lateral_ft: float = 5000, bur: float = 10, rop: float = 25,
                    bit: float = 50, **params: Any) -> dict[str, Any]:
    """Horizontal EGS producer: surface casing, production casing through the build, liner in the lateral.

    The build ends horizontal at ``tvd_ft``; its measured length is ``tvd + (pi/2 - 1) R`` with
    ``R = 5730 / bur``. The lateral liner (9-5/8 in.) is not cemented, as in the GeoVision files.
    """
    R = 100 * 180 / (math.pi * bur)
    end_build = tvd_ft + (math.pi / 2 - 1) * R
    iv = [dict(end=1000, D=22, Dc=18.625, w=131.8, CL=1, log=0, dir=0),
          dict(end=end_build, D=17, Dc=13.625, w=78.5, CL=1, log=0, dir=1),
          dict(end=end_build + lateral_ft, D=12.25, Dc=9.625, w=39.4, CL=0, log=1, dir=1, mult=dict(cem=0))]
    for g in iv:
        g.update(rop=rop, bit=bit)
    return dict(params=dict(params), intervals=iv)


def lateral_for_md(md_ft: float, tvd_ft: float = 8700, bur: float = 10) -> float:
    """Lateral length that gives a horizontal well of measured depth ``md_ft``."""
    R = 100 * 180 / (math.pi * bur)
    return md_ft - (tvd_ft + (math.pi / 2 - 1) * R)


def with_params(spec: dict[str, Any], **params: Any) -> dict[str, Any]:
    """Copy of a spec with some parameters changed (``rop``/``bit`` apply to every interval)."""
    s = copy.deepcopy(spec)
    for k in ("rop", "bit"):
        if k in params:
            v = params.pop(k)
            for g in s["intervals"]:
                g[k] = v
    s.setdefault("params", {}).update(params)
    return s


#: The Cape-like well used throughout Beyond WCS: 2025 NREL parameters, Fervo Phase I geometry.
CAPE_LIKE_PARAMS = dict(bit_mult=2.0, run=800.0, cont=0.10, rig=33500.0, mud0=5000.0)


def cape_like_well(md_ft: float | None = None, rop: float = 75, bit: float = 50) -> dict[str, Any]:
    lateral = 5000 if md_ft is None else lateral_for_md(md_ft)
    return horizontal_well(8700, lateral, 10, rop, bit, **CAPE_LIKE_PARAMS)
