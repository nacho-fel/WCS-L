# Data

| File | What it is | Source |
|---|---|---|
| `data/wcs/geovision_specs.json.gz` | Inputs and formula variants of the 249 GeoVision WCS workbooks, plus each workbook's own computed total, as test fixtures | Extracted with `tools/extract_geovision_specs.py` from the Sandia WCS files (Lowry et al. 2017) shared by T. Lowry. The original spreadsheets are **not** in this repository. |
| `data/fervo/fervo_wells.csv` | Fervo's eight horizontal EGS wells: drilling days, total depth, cost per ft, cost | Cape wells: Akindipe & Witter (2025), Table 1. Project Red days and notes: El-Sadi et al. (2024). |
| `data/forge/forge_16A_activity_hours.csv` | 609 time-log rows from the Utah FORGE 16A(78)-32 daily drilling reports (Oct 2020 to Jan 2021), each with an activity code mapped to the five WCS-L activities | Utah FORGE, [GDR submission 1283](https://gdr.openei.org/submissions/1283) (public). 13 rows had no parsable code and are labelled `unparsed`. |
| `data/time_share/sources.csv`, `components.csv` | The time-share evidence base: one row per source, plus itemised cost components with their class | See `docs/01_time_share_study.md` |

**Before making the repository public:** the GeoVision fixture is derived from files Sandia shared on request. The numbers are inputs and totals of the published GeoVision curves, but check with T. Lowry before publishing it. The engine tests need it, and every other part of the repository works without it.

To regenerate the fixture from the original workbooks:

```bash
pip install openpyxl
python tools/extract_geovision_specs.py --xlsx-dir "path/to/GeoVision WCS files"
```
