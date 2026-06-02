# MCM\_Hyetographs v1.0.0

**Synthetic Design Hyetographs Under Non-Stationary Climate Conditions**  
*Microcanonical Multiplicative Cascade · GEV-GAMLSS · ENSO/ONI Covariate · HEC-HMS Ready*

[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)
[![R ≥ 4.0.0](https://img.shields.io/badge/R-%E2%89%A54.0.0-276DC3?logo=r)](https://www.r-project.org/)
[![DOI preprint](https://img.shields.io/badge/DOI-10.31224%2F7062-orange)](https://doi.org/10.31224/7062)
[![ORCID](https://img.shields.io/badge/ORCID-0009--0003--4328--5691-A6CE39?logo=orcid)](https://orcid.org/0009-0003-4328-5691)

---

## Overview

`MCM_Hyetographs_v1.0.0.R` is a self-contained R script that generates
**synthetic sub-hourly design hyetographs** from daily maximum precipitation
records. It integrates a full non-stationarity diagnostic pipeline with
stochastic disaggregation via a **Microcanonical Multiplicative Cascade (MCM)**,
and exports results directly in **HEC-HMS Specified Hyetograph format**.

The script is documented in the companion preprint:

> Victoria Niño, M.J. (2026). *Synthetic Design Hyetographs Under
> Non-Stationary Climate Conditions: A Microcanonical Multiplicative Cascade
> Model with GEV-GAMLSS Frequency Analysis and ENSO Covariate for the Colombian
> Andes.* EngrXiv (preprint). DOI: [10.31224/7062](https://doi.org/10.31224/7062)

---

## Key features

- **Six-test stationarity ensemble** — majority vote (≥ 3/6) selects the appropriate frequency model automatically:
  Mann-Kendall, Pettitt change-point, Sequential Sneyers, Moving-Window MK, White heteroscedasticity, Augmented Dickey-Fuller.
- **Dual frequency path** — stationary GEV (MLE via `evd`) or non-stationary GEV-GAMLSS with the Oceanic Niño Index (ONI) as a linear covariate in the location parameter.
- **Bell (1969) IDF curves** anchored on station-specific GEV quantiles; the only regional parameter is `ratio_60_1440` = P(1 h)/P(24 h).
- **MCM disaggregation** — binary multiplicative cascade with Beta(α, β) multipliers; temporal structure imposed via Huff (1967) mass curves (Q1–Q4); quartile selectable automatically from the coefficient of advance *r* or set manually.
- **Monte Carlo ensemble** — 500 realisations per return period; Q10/Q50/Q90 percentile hyetographs with proportional mass correction.
- **Nine return periods** — T = 2, 2.33, 5, 10, 25, 50, 100, 200, 500 years.
- **HEC-HMS ready output** — `HEC_HMS_hyetographs.xlsx`, one sheet per return period (Time HH:MM | Q50 mm | Q90 mm | Q50 cumulative mm).
- **Two technical reports** — stationarity analysis (TXT) and hyetograph / IDF metrics (TXT).
- **20 diagnostic figures** — 11 analytical plots + 9 individual hyetograph panels.
- **Self-installing packages** — automatically installs any missing CRAN packages on first run.

---

## Repository structure

```
MCM_Hyetographs/
├── MCM_Hyetographs_v1.0.0.R   ← main script (all modules)
├── Pmax_24h.xlsx               ← input: monthly maximum daily precipitation
├── ONI.xlsx                    ← input: NOAA CPC Oceanic Niño Index
├── README.md                   ← this file
├── LICENSE                     ← MIT license
└── results/                    ← created automatically on first run
    ├── plots/
    │   ├── 01_series_pmax.png
    │   ├── 02_oni_historical.png
    │   ├── 03_scatter_pmax_oni.png
    │   ├── 04_pettitt.png
    │   ├── 05_sneyers.png
    │   ├── 06_moving_window_mk.png
    │   ├── 07_white_heteroscedasticity.png
    │   ├── 08_frequency_curve.png
    │   ├── 09_IDF_curves.png
    │   ├── 10_hyetograph_T2.png  … T2_33, T5, T10, T25, T50, T100, T200, T500
    │   └── 11_comparative_panel.png
    ├── tables/
    │   ├── MCM_results.xlsx            (4 + N sheets: tests, quantiles, hyetographs)
    │   └── HEC_HMS_hyetographs.xlsx    (9 sheets, one per return period)
    ├── REPORT_MCM_Stationarity.txt
    └── REPORT_MCM_Hyetographs.txt
```

---

## Requirements

| Requirement | Details |
|---|---|
| R version | ≥ 4.0.0 (tested on 4.3.1) |
| Operating system | Windows, macOS, Linux |
| IDE | RStudio recommended; also runs via `Rscript` on the command line |
| Internet | Required on first run to install missing packages |

The script checks and installs all 25 required packages automatically.
For reference, the full list is:

```
here, readxl, dplyr, tidyr, lubridate, writexl,
ggplot2, gridExtra, grid, patchwork, ggrepel, scales, RColorBrewer,
Kendall, trend, strucchange, tseries, lmtest,
evd, fitdistrplus, gamlss, gamlss.dist, gamlssx,
knitr, kableExtra
```

---

## Input data format

### `Pmax_24h.xlsx` — Monthly maximum daily precipitation

Wide-format Excel file. **Column headers must be in English** (case-insensitive
after internal standardisation):

| YEAR | JAN | FEB | MAR | APR | MAY | JUN | JUL | AUG | SEP | OCT | NOV | DEC |
|------|-----|-----|-----|-----|-----|-----|-----|-----|-----|-----|-----|-----|
| 1954 | 45.2 | 38.0 | … | | | | | | | | | |
| …   | … | | | | | | | | | | | |

- Units: mm
- Missing values: leave blank or `NA`; the script tolerates up to 20 % missing data (alert above 10 %)
- Minimum recommended record: 20 years (GEV tail estimates are unreliable below this threshold)

### `ONI.xlsx` — Oceanic Niño Index (NOAA CPC)

Wide-format Excel file with the same column structure (YEAR + 12 months):

| YEAR | JAN | FEB | … | DEC |
|------|-----|-----|---|-----|
| 1950 | −1.4 | −1.2 | … | |

Download the latest ONI table from [NOAA CPC](https://origin.cpc.ncep.noaa.gov/products/analysis_monitoring/ensostuff/ONI_v5.php)
and save as `.xlsx`. If you have an `.xls` file, convert it once with:

```r
writexl::write_xlsx(readxl::read_xls("ONI.xls"), "ONI.xlsx")
```

Both input files must be placed **in the same folder as the script**.

---

## Quick start

### Option A — Interactive (RStudio)

1. Clone or download this repository.
2. Place `Pmax_24h.xlsx` and `ONI.xlsx` in the script folder.
3. Open `MCM_Hyetographs_v1.0.0.R` in RStudio.
4. Edit the `CONFIG` section at the top of the file for your study site (see [Configuration](#configuration)).
5. Click **Source** (or press `Ctrl+Shift+S`).
6. When prompted, enter the storm duration in hours and the time step in minutes.

### Option B — Batch / command line

Uncomment `storm_duration_h` and `time_step_min` inside `CONFIG`, then run:

```bash
Rscript MCM_Hyetographs_v1.0.0.R
```

No interactive prompts will appear.

---

## Configuration

All user-adjustable parameters are contained in the `CONFIG` list at the top
of the script. **Only modify this section** to adapt the script to a new site.

```r
CONFIG <- list(

  # ── File paths (relative to the script folder) ──────────────────────────
  file_pmax        = "Pmax_24h.xlsx",
  file_oni         = "ONI.xlsx",
  output_subfolder = "results",

  # ── Return periods (years) ───────────────────────────────────────────────
  return_periods   = c(2, 2.33, 5, 10, 25, 50, 100, 200, 500),

  # ── Batch mode (uncomment to skip interactive prompts) ───────────────────
  # storm_duration_h = 6,     # storm duration in hours
  # time_step_min    = 15,    # time step in minutes (must divide duration exactly)

  # ── Statistical significance level ───────────────────────────────────────
  alpha_stat       = 0.05,

  # ── IDF regional parameter ────────────────────────────────────────────────
  # P(1h)/P(24h) ratio — the ONLY regional parameter of the Bell (1969) method.
  # See regional guide below or compute from local overlapping records.
  ratio_60_1440    = 0.40,   # typical for inter-Andean valleys (Cali, Palmira)

  # ── MCM disaggregation ────────────────────────────────────────────────────
  n_sim            = 500,    # Monte Carlo realisations per return period
  mcm_seed         = NULL,   # integer seed for reproducibility (e.g. 42); NULL = random
  mcm_p0_factor    = 0.3,    # dry-fraction reduction factor (Molnar & Burlando 2005)
  mcm_max_frac     = 0.60,   # maximum fraction of Ptotal per time step; NULL = no limit

  # ── Huff quartile / coefficient of advance ───────────────────────────────
  # Option A (with rain gauge): provide observed median r → automatic quartile selection
  huff_r_median    = NULL,   # e.g. 0.167 for La Primavera station (Victoria Niño 2026)
  # Option B (no rain gauge): set quartile directly (1, 2, 3, or 4)
  huff_quartile    = 2,      # Q2 = Andean urban default

  # ── Plot export ───────────────────────────────────────────────────────────
  dpi              = 180,
  fig_width_in     = 10,
  fig_height_in    = 6
)
```

### Regional guide for `ratio_60_1440`

| Region / Zone (Colombia) | Range | Typical |
|---|---|---|
| Valle del Cauca, Cauca (Cali, Palmira) | 0.38–0.45 | **0.40** |
| Antioquia (Medellín, Aburrá Valley) | 0.36–0.44 | 0.40 |
| Cundinamarca (Bogotá and savanna) | 0.32–0.40 | 0.36 |
| Coffee Region (Manizales, Pereira) | 0.38–0.46 | 0.42 |
| Huila, Tolima (warm valleys) | 0.40–0.48 | 0.44 |
| Nariño (Pasto and highlands) | 0.30–0.38 | 0.34 |
| Dry Caribbean (Guajira, Valledupar) | 0.42–0.55 | 0.48 |
| Wet Caribbean (Barranquilla, Cartagena) | 0.40–0.50 | 0.45 |
| Pacific Region (Chocó, Buenaventura) | 0.35–0.48 | 0.42 |
| Orinoquía / Eastern Plains | 0.40–0.52 | 0.46 |
| Amazon region | 0.38–0.48 | 0.43 |

With overlapping hourly and daily records, compute directly:
`ratio_60_1440 = mean(Pmax_1h) / mean(Pmax_24h)` over common years.

### Huff quartile guide

| Quartile | Coefficient of advance *r* | Peak position | Typical regime (Colombia) |
|---|---|---|---|
| Q1 | 0.10–0.25 | 0–25 % of duration | High orography (> 1500 m), Caribbean coast |
| **Q2** *(default)* | 0.25–0.45 | 25–50 % | Andean urban: Cali, Medellín, Bogotá |
| Q3 | 0.45–0.60 | 50–75 % | Pacific region (Chocó, Buenaventura) |
| Q4 | 0.60–0.80 | 75–100 % | Eastern Plains (Apr–May transition) |

If sub-hourly rain-gauge records are available, compute the median *r* with
the methodology of [Victoria Niño (2026)](https://doi.org/10.31224/7062)
and set `huff_r_median`; the quartile is selected automatically.

---

## Workflow

The script executes 15 sequential modules:

```
[00] Package check & install
[01] Path detection & input validation
[02] Auxiliary functions
[03] Data loading (Pmax + ONI) ──────────────────────────────┐
[04] Exploratory visualisation                               │
[05] Stationarity tests (MK, Pettitt, Sneyers,               │
     MWMK, White, ADF)                                      │
[06] Majority-vote decision (≥ 3/6 → non-stationary)        ▼
[07] Frequency analysis ─── stationary GEV (evd::fgev)
                        └── non-stationary GEV-GAMLSS (ONI)
[08] Storm duration + time step (interactive or batch)
[09] IDF curves — Bell (1969) anchored on GEV quantiles
[10] MCM disaggregation — 500 Monte Carlo realisations
[11] Individual hyetograph plots (Q10 / Q50 / Q90)
[12] Comparative panel (all return periods)
[13] Excel export (MCM_results + HEC_HMS_hyetographs)
[14] Stationarity report (TXT)
[14b] Hyetograph + IDF metrics report (TXT)
[15] Execution summary
```

---

## Outputs

### Excel files

| File | Sheets | Contents |
|---|---|---|
| `MCM_results.xlsx` | 4 + 9 | Test summary, design precipitation, incremental Q10/Q50/Q90 hyetographs, cumulative Q50, one HEC-HMS sheet per return period |
| `HEC_HMS_hyetographs.xlsx` | 9 | One sheet per return period: Time HH:MM · Q50 (mm) · Q90 (mm) · Cumulative Q50 (mm) |

### Text reports

| File | Contents |
|---|---|
| `REPORT_MCM_Stationarity.txt` | Series statistics, six-test results table, final decision, design precipitation, MCM parameters, hyetograph summary (P50), methodological considerations, test power table |
| `REPORT_MCM_Hyetographs.txt` | GEV parameters and goodness of fit, design quantiles with ENSO scenarios, IDF table and equation, MCM parameters, ensemble metrics by return period, peak intensity by percentile, temporal concentration analysis, usage recommendations |

### Figures

| File | Description |
|---|---|
| `01_series_pmax.png` | Annual Pmax time series coloured by ONI, with linear trend |
| `02_oni_historical.png` | Annual mean ONI bar chart (1954–2024) |
| `03_scatter_pmax_oni.png` | Pmax vs. ONI scatter with linear regression and LOESS |
| `04_pettitt.png` | Pettitt U_t statistic with detected change-point |
| `05_sneyers.png` | Sequential Sneyers progressive/retrogressive statistics |
| `06_moving_window_mk.png` | Moving-Window Mann-Kendall (w = 21 yr), significant windows highlighted |
| `07_white_heteroscedasticity.png` | Residuals vs. fitted values (White test diagnostic) |
| `08_frequency_curve.png` | Empirical data vs. fitted GEV frequency curve |
| `09_IDF_curves.png` | IDF curves for all nine return periods (Bell 1969) |
| `10_hyetograph_T*.png` | Individual Q10/Q50/Q90 hyetograph + mass curve per return period (9 files) |
| `11_comparative_panel.png` | Multi-panel overview of all nine Q50 hyetographs |

---

## Reproducibility

To obtain exactly the results reported in the companion manuscript, set:

```r
mcm_seed = 42
```

in the `CONFIG` section before running. With `mcm_seed = NULL` (the default),
each run produces a stochastically equivalent but numerically distinct ensemble.

---

## Limitations and recommendations

1. **Scale extrapolation.** MCM Beta parameters (α, β) are calibrated from monthly
   Pmax data (~24 h scale) and applied at 15-min intervals — approximately 2 orders
   of magnitude of extrapolation. Ensemble CV of 72–77 % reflects this gap.
   For critical infrastructure, use the Q90 percentile hyetograph.
   Local recalibration from sub-hourly rain-gauge records eliminates the extrapolation.

2. **Non-stationary return periods.** Design quantiles are conditional on the neutral
   ENSO state (ONI = 0). In a non-stationary context the classical return period
   should be interpreted as an expected waiting time or replaced by the cumulative
   exceedance probability (CEP) over the design life (Salas & Obeysekera 2014).

3. **Minimum record length.** At least 20 years are required for GEV fitting;
   30+ years are recommended for T > 50-year quantiles.

4. **Single-station.** Results apply to the calibration station.
   Regionalisation across a network requires additional analysis.

5. **Huff quartile default.** Q2 is the default for Andean urban contexts.
   Verify against local rain-gauge records if available.

---

## Citation

If you use this script in your research or engineering work, please cite:

```bibtex
@misc{victoria2026mcm,
  author       = {Victoria Ni{\~n}o, Mauricio Javier},
  title        = {{MCM\_Hyetographs} v1.0.0: Synthetic Design Hyetographs
                  Under Non-Stationary Climate Conditions},
  year         = {2026},
  publisher    = {GitHub},
  url          = {https://github.com/YOUR_USERNAME/MCM_Hyetographs},
  note         = {Companion preprint: \doi{10.31224/7062}}
}
```

And the companion preprint:

```bibtex
@article{victoria2026preprint,
  author  = {Victoria Ni{\~n}o, Mauricio Javier},
  title   = {Synthetic Design Hyetographs Under Non-Stationary Climate
             Conditions: A Microcanonical Multiplicative Cascade Model
             with {GEV-GAMLSS} Frequency Analysis and {ENSO} Covariate
             for the Colombian Andes},
  journal = {EngrXiv (preprint)},
  year    = {2026},
  doi     = {10.31224/7062}
}
```

---

## References

Bell, F.C. (1969). Generalized rainfall-duration-frequency relationships.
*Journal of the Hydraulics Division, ASCE*, 95(1), 311–327.

Hingray, B. & Ben Haha, M. (2005). Statistical performances of various
deterministic and stochastic models for rainfall series disaggregation.
*Atmospheric Research*, 77(1–4), 169–185.
DOI: [10.1016/j.atmosres.2004.10.023](https://doi.org/10.1016/j.atmosres.2004.10.023)

Huff, F.A. (1967). Time distribution of rainfall in heavy storms.
*Water Resources Research*, 3(4), 1007–1019.
DOI: [10.1029/WR003i004p01007](https://doi.org/10.1029/WR003i004p01007)

Molnar, P. & Burlando, P. (2005). Preservation of rainfall properties in
stochastic disaggregation by a simple random cascade model.
*Atmospheric Research*, 77(1–4), 137–151.
DOI: [10.1016/j.atmosres.2004.10.024](https://doi.org/10.1016/j.atmosres.2004.10.024)

Rigby, R.A. & Stasinopoulos, D.M. (2005). Generalized additive models for
location, scale and shape. *Journal of the Royal Statistical Society: Series C*,
54(3), 507–554.
DOI: [10.1111/j.1467-9876.2005.00510.x](https://doi.org/10.1111/j.1467-9876.2005.00510.x)

Salas, J.D. & Obeysekera, J. (2014). Revisiting the concepts of return period
and risk for nonstationary hydrologic extreme events.
*Journal of Hydrologic Engineering*, 19(3), 554–568.
DOI: [10.1061/(ASCE)HE.1943-5584.0000820](https://doi.org/10.1061/(ASCE)HE.1943-5584.0000820)

Schertzer, D. & Lovejoy, S. (1987). Physical modeling and analysis of rain and
clouds by anisotropic scaling multiplicative processes.
*Journal of Geophysical Research*, 92(D8), 9693–9714.
DOI: [10.1029/JD092iD08p09693](https://doi.org/10.1029/JD092iD08p09693)

Victoria Niño, M.J. (2026). Design storm duration from hourly rainfall records
in a bimodal Andean climate. *EngrXiv* (preprint).
DOI: [10.31224/7062](https://doi.org/10.31224/7062)

Yue, S. & Wang, C. (2004). The Mann-Kendall test modified by effective sample
size to detect trend in serially correlated hydrological series.
*Water Resources Management*, 18(3), 201–218.
DOI: [10.1023/B:WARM.0000043140.61082.60](https://doi.org/10.1023/B:WARM.0000043140.61082.60)

---

## License

MIT License — Copyright © 2025 Mauricio Javier Victoria Niño.
See [`LICENSE`](LICENSE) for the full text.

---

## Author

**Mauricio Javier Victoria Niño**  
Independent Researcher · Cali, Colombia  
✉ [hidratecsa@gmail.com](mailto:hidratecsa@gmail.com)  
🔬 [ORCID 0009-0003-4328-5691](https://orcid.org/0009-0003-4328-5691)
