# GEOPHIRES cross-check: production and cost

**Result.**
- GEOPHIRES' own Fervo cases reproduce, and GEOPHIRES lands within 4 % of Project Red's reported 3.5 MW, read as gross output at peak flow.
- Replacing the flat $4.65M per well with WCS-L's learning campaign (94 wells) raises the average well cost to $5.16M (P10–P90: $4.72–5.65M).
- That moves the LCOE from 8.50 to 8.67 cents/kWh (8.52–8.83) and lowers NPV by 13 % (P50).
- The flat NREL number matches Fervo's mature wells but misses the expensive first wells of a campaign.

Run: `python geophires/run_geophires_fervo.py --examples <GEOPHIRES-X>/tests/examples`. It needs GEOPHIRES-X installed from https://github.com/NREL/GEOPHIRES-X. Results are in [`reports/geophires/results.md`](../reports/geophires/results.md); the input and output files of every run are in `reports/geophires/`.

## 1. Fervo Cape case (GEOPHIRES `Fervo_Project_Cape-5`) reproduced

| | Value |
|---|---|
| Net generation | 510 MW from 56 producers and 38 injectors |
| Per producer | 9.1 MW |
| Flow per well | 107 kg/s |
| LCOE | 8.50 cents/kWh |
| NPV | $199M |
| After-tax IRR | 22.4 % |
| Drilling cost | $4.65M per well (NREL 2025 curve × 0.9), $437M in total |
| Redrilling events | 3 |

## 2. Project Red (GEOPHIRES `Fervo_Norbeck_Latimer_2023`)

Fervo reported 3.5 MW from one well pair at 63 L/s and 191 °C over a 30-day test.

| GEOPHIRES run | Gross MW | Net MW | Average production temperature (10-year life) |
|---|---|---|---|
| 41 kg/s (650 gpm, as in the case file) | 2.90 | 2.26 | 174 °C |
| 60 kg/s (≈ 63 L/s) | 3.37 | 1.87 | 157 °C |

At the peak flow, GEOPHIRES' gross output is within 4 % of Fervo's 3.5 MW. Net output is lower because GEOPHIRES subtracts pumping power, and at higher flow the reservoir cools faster over its 10-year life. Fervo's release does not say whether 3.5 MW is gross or net. It was a 30-day test, so the averages over the plant's life are not directly comparable.

## 3. Cape case with different well costs

| Well cost used | $M per well | LCOE, ¢/kWh | NPV, $M | IRR |
|---|---|---|---|---|
| NREL 2025 curve × 0.9 (GEOPHIRES case) | 4.65 | 8.50 | 199 | 22.4 % |
| WCS, NREL 2025 parameters (no learning) | 4.58 | 8.48 | 202 | 22.5 % |
| **WCS-L, 94-well campaign, P50** | **5.16** | **8.67** | **173** | **21.4 %** |
| WCS-L, P10 – P90 | 4.72 – 5.65 | 8.52 – 8.83 | 196 – 148 | 22.3 – 20.6 % |
| Fervo's six Cape wells, observed average | 5.85 | 8.90 | 138 | 20.2 % |
| Fervo's first Cape well, for every well | 8.44 | 9.76 | 7 | 15.3 % |

WCS-L campaign assumptions:
- 94 wells, drilled by 3 rigs at about 1.5 wells per rig per month (Fervo's stated steady-state pace);
- 6 wells per pad;
- 14,100 ft measured depth;
- experience from another rig counts at 0.4.

The first well costs $8.8M and the last $5.0M.

## How to read it

- **The flat NREL value is a mature-well cost.** $4.65M sits close to WCS-L's floor ($5.0M) and Fervo's best wells ($4.6–4.9M). A real campaign starts expensive, so its average is about 11 % higher.
- **The effect on LCOE is small but not negligible: +2 %, P10–P90 +0.2 % to +4 %.** Drilling is about one sixth of this project's capex; the surface plant dominates. The effect is larger on NPV (−13 %) and much larger for smaller projects with fewer wells, where the expensive first wells weigh more. That case matters for the data-centre MILP.
- **Independent support for the floor.** GEOPHIRES' newest Cape case (`Fervo_Project_Cape-7`, 2026) has added a redrilling learning rate. It states that geothermal well costs "flatten at 70 to 75% of the first well's", citing El-Sadi et al. (2024), Akindipe & Witter (2025), Lowry et al. (2017) and Kellogg. WCS-L's fitted floor for planned time is 66 %, and its last-to-first cost ratio in the campaign above is 56 %.

## Caveats

- GEOPHIRES applies one per-well cost to every well, redrills included. In WCS-L, redrills come late in the experience curve, so the campaign average slightly overstates their cost.
- The exploration cost ($30M) was derived in the case file from the $4.65M well and is not rescaled here.
- The Cape-5 geometry (2.68 km TVD, 5,000 ft lateral) is Phase I. Cape-7 moves to the Phase II design (3.06 km, 7,500 ft lateral, $8.5M per well all-in). Repeating the cross-check on Cape-7 is the next step.
