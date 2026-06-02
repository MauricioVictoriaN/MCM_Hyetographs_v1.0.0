# ==============================================================================
# MCM_Hyetographs_v1.0.0.R
# ==============================================================================
# ------------------------------------------------------------------------------
# Author  : Mauricio Javier Victoria Niño
#           Independent Researcher · Cali, Colombia
#           hidratecsa@gmail.com
#           ORCID: 0009-0003-4328-5691
# ------------------------------------------------------------------------------
# Description:
#   Generates design hyetograms using the Microcanonical Multiplicative Cascade
#   (MCM) model for sub-hourly temporal disaggregation of daily maximum
#   precipitation. Incorporates stationarity/non-stationarity analysis and
#   the El Niño–Southern Oscillation (ENSO) index (ONI) as a covariate for
#   non-stationary frequency analysis.
#
#   INCLUDES:
#   · Six stationarity tests (Mann-Kendall, Pettitt, Sequential Sneyers,
#     Moving-Window MK, White, Augmented Dickey-Fuller)
#   · Stationary GEV (MLE) and Non-stationary GEV (GAMLSS, ONI covariate)
#   · IDF curves derived from fitted GEV quantiles
#   · MCM disaggregation with Monte Carlo ensemble (n = 500)
#   · Design hyetograms for 9 return periods (T = 2, 2.33, 5, 10, 25,
#     50, 100, 200, 500 years), percentiles P10 / P50 / P90
#   · HEC-HMS ready Excel output (Specified Hyetograph format)
#   · Two detailed technical reports (stationarity + hyetograms + IDF metrics)
# ------------------------------------------------------------------------------
# Included modules:
#   00. Package installation and loading
#   01. User configuration (paths, design parameters)
#   02. Auxiliary functions and helpers
#   03. Data loading and preparation (Pmax 24h + ONI index)
#   04. Exploratory visualisation
#   05. Stationarity / non-stationarity analysis
#       05a. Mann-Kendall (classical)
#       05b. Pettitt change-point test
#       05c. Sequential Sneyers
#       05d. Moving-Window Mann-Kendall
#       05e. White heteroscedasticity test
#       05f. Augmented Dickey-Fuller (unit root)
#   06. Decision: stationary vs non-stationary (majority vote, 6 tests)
#   07. Frequency analysis — GEV stationary (evd::fgev) or
#       non-stationary GEV-GAMLSS with ONI covariate
#   08. Design storm parameters (interactive duration + interval input)
#   09. IDF curves (power-law, calibrated for Andean Colombia)
#   10. MCM disaggregation (binary multiplicative cascade, Beta(α, β))
#   11. Individual hyetogram plots (P10 / P50 / P90 bands)
#   12. WHO comparative panel (all return periods)
#   13. Results export — two Excel files (full analysis + HEC-HMS format)
#   14. Stationarity technical report (TXT)
#   14b. Hyetograms, IDF equation and metrics technical report (TXT)
#   15. Execution summary
# ------------------------------------------------------------------------------
# Input files (English headers required):
#   · Pmax_24h.xlsx  — Monthly maximum daily precipitation (mm), wide format:
#                      columns YEAR, JAN, FEB, ..., DEC
#   · ONI.xlsx       — Oceanic Niño Index (NOAA CPC), monthly, wide format:
#                      columns YEAR, JAN, FEB, MAR, APR, MAY, JUN, JUL, AUG, SEP, OCT, NOV, DEC
#
# Output folder: OUTPUT_DIR (set in CONFIG section below)
#   ├── plots/           PNG figures (01–11)
#   ├── tables/
#   │   ├── MCM_results.xlsx          Full results (4 + N sheets)
#   │   └── HEC_HMS_hyetographs.xlsx     HEC-HMS import-ready file
#   ├── REPORT_MCM_Stationarity.txt  Stationarity analysis report
#   └── REPORT_MCM_Hyetographs.txt      Hyetograms, IDF and metrics report
# ------------------------------------------------------------------------------
# References:
#   Cooley, D. (2013). Return periods and return levels under climate change.
#     In: AghaKouchak et al. (Eds.), Extremes in a Changing Climate.
#     Springer, Dordrecht. DOI: 10.1007/978-94-007-4479-0_4
#
#   Hingray, B. & Ben Haha, M. (2005). Statistical performances of various
#     deterministic and stochastic models for rainfall series disaggregation.
#     Atmospheric Research, 77(1–4), 169–185.
#     DOI: 10.1016/j.atmosres.2004.10.023
#
#   IDEAM (2014). Estudio Nacional del Agua 2014. Instituto de Hidrología,
#     Meteorología y Estudios Ambientales. Bogotá, Colombia.
#
#   Molnar, P. & Burlando, P. (2005). Preservation of rainfall properties in
#     stochastic disaggregation by a simple random cascade model.
#     Atmospheric Research, 77(1–4), 137–151.
#     DOI: 10.1016/j.atmosres.2004.10.024
#
#   Salas, J. D. & Obeysekera, J. (2014). Revisiting the concepts of return
#     period and risk for nonstationary hydrologic extreme events.
#     Journal of Hydrologic Engineering, 19(3), 554–568.
#     DOI: 10.1061/(ASCE)HE.1943-5584.0000820
#
#   Schertzer, D. & Lovejoy, S. (1987). Physical modeling and analysis of
#     rain and clouds by anisotropic scaling multiplicative processes.
#     Journal of Geophysical Research, 92(D8), 9693–9714.
#     DOI: 10.1029/JD092iD08p09693
#
#   Stedinger, J. R., Vogel, R. M. & Foufoula-Georgiou, E. (1993).
#     Frequency analysis of extreme events. In: Maidment, D. (Ed.),
#     Handbook of Hydrology. McGraw-Hill, New York. Chapter 18.
#
#   WMO (2008). Guide to Hydrological Practices, Volume I: Hydrology —
#     From Measurement to Hydrological Information. WMO-No. 168, 6th ed.
#     World Meteorological Organization, Geneva.
# ------------------------------------------------------------------------------
# License : MIT License
#           Copyright (c) 2025 Mauricio Javier Victoria Niño
#           Permission is hereby granted, free of charge, to any person
#           obtaining a copy of this software and associated documentation
#           files (the "Software"), to deal in the Software without
#           restriction, including without limitation the rights to use,
#           copy, modify, merge, publish, distribute, sublicense, and/or
#           sell copies of the Software, and to permit persons to whom the
#           Software is furnished to do so, subject to the following
#           conditions: The above copyright notice and this permission notice
#           shall be included in all copies or substantial portions of the
#           Software. THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF
#           ANY KIND, EXPRESS OR IMPLIED.
# ------------------------------------------------------------------------------
# Version : 1.0.0
# Date    : 2026-05-31
# Language: R (>= 4.0.0)
# ==============================================================================

# ==============================================================================
# MODULE 00: CENTRALISED CONFIGURATION
# ==============================================================================
# All user-adjustable parameters are concentrated here.
# Modify ONLY this section to adapt the script to a new study site.

CONFIG <- list(
  # ── Version ────────────────────────────────────────────────────────────────
  version           = "1.0.0",
  
  # ── Input file paths — relative to the script folder ───────────────────────
  # Place Pmax_24h.xlsx and ONI.xlsx in the same folder as this script,
  # or adjust the file names to match your directory structure.
  # ONI must be in .xlsx format (Excel 2007+). If you have ONI.xls,
  # convert it once with: writexl::write_xlsx(readxl::read_xls("ONI.xls"), "ONI.xlsx")
  file_pmax         = "Pmax_24h.xlsx",   # Pmax input file name (in script folder)
  file_oni          = "ONI.xlsx",         # ONI input file name  (in script folder)
  output_subfolder  = "results",          # output subfolder inside the script folder
  
  # ── Return periods (years) ─────────────────────────────────────────────────
  # 2.33 = mean annual flood (Langbein, 1949)
  return_periods    = c(2, 2.33, 5, 10, 25, 50, 100, 200, 500),
  
  # ── Design storm (set values to bypass interactive prompts) ───────────────
  # Uncomment to run in batch mode:
  # storm_duration_h  = 6,    # Storm duration (hours)
  # time_step_min     = 15,   # Time step (minutes); must divide duration exactly
  
  # ── Statistical significance level ────────────────────────────────────────
  alpha_stat        = 0.05,   # significance level for hypothesis tests
  
  # ── IDF — P(1h)/P(24h) ratio — regional parameter for Bell (1969) ──────────
  # The IDF curve is built with the Bell (1969) formula:
  #   P(t,T) = (0.54·t^0.25 - 0.50) · P(60min,T)   [t in minutes]
  #   P(60min,T) = P(24h,T) · ratio_60_1440
  #
  # ratio_60_1440 = P(1h) / P(24h) is the ratio of the 1-hour maximum
  # precipitation to the 24-hour maximum for the same return period.
  # It is the ONLY regional parameter of the method.
  #
  # SELECTION GUIDE BY COLOMBIAN REGION:
  # (Sources: IDEAM 2014, ENA; Vélez et al. 2002; Montoya et al. 2011;
  #  Chen 1983; IDF studies by UNAL, U. del Valle, U. de Antioquia)
  #
  # ┌──────────────────────────────────────────┬──────────┬─────────────────┐
  # │ Region / Zone (Colombia)                 │ Range    │ Typical value   │
  # ├──────────────────────────────────────────┼──────────┼─────────────────┤
  # │ Andean Region — inter-Andean valley      │          │                 │
  # │   Valle del Cauca, Cauca (Cali, Palmira) │ 0.38–0.45│ 0.40            │
  # │   Antioquia (Medellín, Aburrá Valley)    │ 0.36–0.44│ 0.40            │
  # │   Cundinamarca (Bogotá and savanna)      │ 0.32–0.40│ 0.36            │
  # │   Coffee Region (Manizales, Pereira)     │ 0.38–0.46│ 0.42            │
  # │   Huila, Tolima (warm valleys)           │ 0.40–0.48│ 0.44            │
  # │   Nariño (Pasto and highlands)           │ 0.30–0.38│ 0.34            │
  # ├──────────────────────────────────────────┼──────────┼─────────────────┤
  # │ Caribbean Coast                          │          │                 │
  # │   Dry Caribbean (Guajira, Valledupar)    │ 0.42–0.55│ 0.48            │
  # │   Wet Caribbean (Barranquilla, Cartagena)│ 0.40–0.50│ 0.45            │
  # │   Sierra Nevada de Santa Marta           │ 0.38–0.46│ 0.42            │
  # ├──────────────────────────────────────────┼──────────┼─────────────────┤
  # │ Pacific Region (Chocó, Buenaventura)     │ 0.35–0.48│ 0.42            │
  # ├──────────────────────────────────────────┼──────────┼─────────────────┤
  # │ Orinoquía / Eastern Plains               │ 0.40–0.52│ 0.46            │
  # ├──────────────────────────────────────────┼──────────┼─────────────────┤
  # │ Amazon region                            │ 0.38–0.48│ 0.43            │
  # └──────────────────────────────────────────┴──────────┴─────────────────┘
  #
  # USAGE NOTES:
  # · Higher values (>0.45): intense convective regimes, short storms with
  #   high peak intensity (Caribbean coast, Eastern Plains, warm valleys).
  # · Lower values (<0.38): stratiform or high-mountain regimes where
  #   rainfall is more uniform in time (Bogotá savanna, Nariño, páramos).
  # · With a rain gauge (continuous record): compute directly as
  #   ratio = mean(Pmax_1h) / mean(Pmax_24h) for the overlapping period
  #   and use that value instead of the tabular guidance.
  # · Without local data: use the typical value for the region.
  # · Plausible general range for Colombia: 0.28 – 0.58.
  #   Values outside this range should be justified.
  ratio_60_1440 = 0.40,   # P(1h)/P(24h) — see regional guide above
  
  # ── MCM disaggregation ────────────────────────────────────────────────────
  n_sim             = 500,    # Monte Carlo simulations per return period
  mcm_seed          = NULL,   # Set integer for reproducibility (e.g., 42)
  # Dry-fraction reduction factor for design storms.
  # Rationale: extreme events have fewer dry intervals than the calibrated
  # monthly mean. Value 0.3 from Molnar & Burlando (2005).
  mcm_p0_factor     = 0.3,    # P0 reduction: P0_design = P0_calibrated * mcm_p0_factor
  # Maximum fraction of Ptotal allowed in a single time step.
  # Prevents unrestricted peak concentration in one cascade block.
  # NULL = no limit (classic MCM behaviour). Suggested value: 0.60
  mcm_max_frac      = 0.60,   # Max fraction of Ptotal per interval (NULL to deactivate)
  
  # ── Hyetograph temporal structure — Huff quartile and advance coefficient r ─
  # The MCM generates rainfall MAGNITUDES; the TEMPORAL POSITION of the peak
  # is controlled by the Huff (1967) mass curve.
  #
  # QUARTILE SELECTION BASED ON THE COEFFICIENT OF ADVANCE (r):
  # ─────────────────────────────────────────────────────────────
  # With a local rain gauge, the quartile can be determined objectively
  # from the coefficient of advance r, defined as:
  #
  #   r = (t_peak − t_start) / D_physical    [dimensionless, 0 ≤ r ≤ 1]
  #
  # where t_peak is the timestamp of the maximum hourly intensity within the
  # event, t_start is the storm start time, and D_physical is the total
  # physical duration. The recommended design value is the MEDIAN of r
  # computed over all independent storms catalogued at the station (IETD method).
  #
  # Mapping between r (observed median) and Huff quartile:
  #
  #   r ∈ [0.10 – 0.25]  →  Q1: mass concentrated at start, very early peak
  #                          Impulsive convective storms, high orography,
  #                          mountain stations (e.g. La Primavera, Buga:
  #                          r_median = 0.167, Victoria Niño 2026)
  #
  #   r ∈ [0.25 – 0.45]  →  Q2: mass in first half, centred peak ← DEFAULT
  #                          Urban and peri-urban Andean convective storms.
  #                          Suitable for cities: Cali, Medellín, Bogotá,
  #                          Manizales, Pereira. Moderate orographic regime.
  #
  #   r ∈ [0.45 – 0.60]  →  Q3: mass in second half, late peak
  #                          Mixed convective-stratiform storms.
  #                          Pacific region (Chocó, Buenaventura), high
  #                          persistent humidity, long-duration events.
  #
  #   r ∈ [0.60 – 0.80]  →  Q4: mass at end, very late peak
  #                          Long-duration stratiform storms.
  #                          Eastern Plains (Apr–May transition). Rare.
  #
  # RECOMMENDED PROCEDURE:
  #   1. With rain gauge: compute the median r using the script by
  #      Victoria Niño (2026) — DOI: 10.31224/7062 — and set huff_r_median.
  #      The quartile is selected automatically from the mapping above.
  #   2. Without rain gauge: set huff_quartile directly from regional guide.
  #
  # QUICK REGIONAL GUIDE (no rain gauge, Colombia):
  #   Andean urban (Cali, Medellín, Bogotá, Coffee Region):  Q2
  #   Mountain stations / high orography (>1500 m, Valle):   Q1–Q2
  #   Caribbean coast:                                        Q1 or Q2
  #   Pacific region (Chocó, Buenaventura):                   Q3
  #   Eastern Plains / Orinoquía:                             Q2 or Q4
  #   Amazon region:                                          Q2 or Q3
  #
  # Ref: Huff, F.A. (1967). Water Resources Research, 3(4), 1007–1019.
  #      Victoria Niño, M.J. (2026). Design storm duration from hourly rainfall
  #      records in a bimodal Andean climate. EngrXiv.
  #      DOI: 10.31224/7062
  #
  # ── Parameters to set ───────────────────────────────────────────────────
  # Option A (WITH rain gauge): provide median r → automatic quartile.
  #   huff_r_median  = 0.167   # measured median r (e.g. La Primavera station)
  #   huff_quartile  = NULL    # NULL = automatic selection from r
  #
  # Option B (WITHOUT rain gauge): set quartile directly.
  #   huff_r_median  = NULL    # NULL = not available
  #   huff_quartile  = 2       # 1, 2, 3 or 4 from regional guide
  #
  # If both are provided, huff_r_median takes priority.
  huff_r_median     = NULL,   # observed median r (NULL = not available)
  huff_quartile     = 2,      # default quartile when huff_r_median = NULL
  
  # ── Plot export ───────────────────────────────────────────────────────────
  dpi               = 180,
  fig_width_in      = 10,
  fig_height_in     = 6
)


cat("╔══════════════════════════════════════════════════════════════════════╗\n")
cat("║   MCM — Microcanonical Multiplicative Cascade                       ║\n")
cat("║   Design Hyetographs with Non-Stationarity Analysis                 ║\n")
cat("╚══════════════════════════════════════════════════════════════════════╝\n\n")

cat("┌──────────────────────────────────────────────────────────────────────┐\n")
cat("│  USAGE NOTICE — READ BEFORE RUNNING                                 │\n")
cat("├──────────────────────────────────────────────────────────────────────┤\n")
cat("│  This script generates SYNTHETIC design hyetographs. Results are    │\n")
cat("│  ESTIMATES subject to the following limitations:                     │\n")
cat("│                                                                      │\n")
cat("│  1. MCM parameters (α, β) are calibrated from DAILY Pmax data,      │\n")
cat("│     not from sub-hourly rain gauges. Use as preliminary design;      │\n")
cat("│     calibrate locally for final engineering design.                  │\n")
cat("│                                                                      │\n")
cat("│  2. IDF uses Bell (1969) with regional P(1h)/P(24h) ratio.          │\n")
cat("│     Adjust CONFIG$ratio_60_1440 if local data are available.         │\n")
cat("│                                                                      │\n")
cat("│  3. Temporal structure uses Huff curves. The quartile can be         │\n")
cat("│     selected automatically from the coefficient of advance r         │\n")
cat("│     (CONFIG$huff_r_median, Victoria Niño 2026, DOI:10.31224/7062)    │\n")
cat("│     or set manually via CONFIG$huff_quartile. Default: Q2.           │\n")
cat("│                                                                      │\n")
cat("│  4. Minimum recommended record: 20 years. <20 years → high          │\n")
cat("│     GEV tail uncertainty for T > 50 years.                          │\n")
cat("│                                                                      │\n")
cat("│  Ref: Molnar & Burlando (2005), Bell (1969), Huff (1967),           │\n")
cat("│       IDEAM (2014), Salas & Obeysekera (2014).                       │\n")
cat("└──────────────────────────────────────────────────────────────────────┘\n\n")

# ── Package catalogue with descriptions ───────────────────────────────────────
packages_info <- list(
  # Portable paths
  here         = "Reproducible relative paths (here::here())",
  # Data handling
  readxl       = "Read Excel files (.xls / .xlsx)",
  dplyr        = "Data frame manipulation (filter, mutate, summarise...)",
  tidyr        = "Data reshaping (pivot_longer, pivot_wider...)",
  lubridate    = "Date handling and arithmetic",
  writexl      = "Write Excel files without Java dependencies",
  # Visualisation
  ggplot2      = "High-quality graphics system (Grammar of Graphics)",
  gridExtra    = "Multi-plot layout in a grid",
  grid         = "Base grid system for complex layouts",
  patchwork    = "Intuitive combination of ggplot objects",
  ggrepel      = "Non-overlapping labels in ggplot2",
  scales       = "Scales, formats and transformations in ggplot2",
  RColorBrewer = "Colour palettes for cartography and statistics",
  # Trend and stationarity tests
  Kendall      = "Mann-Kendall test and variants",
  trend        = "Pettitt, Sen slope, Sequential Sneyers",
  strucchange  = "Structural change detection (CUSUM, F-test)",
  # Time series
  tseries      = "ADF (Augmented Dickey-Fuller), KPSS, ARIMA",
  # Regression diagnostics
  lmtest       = "Breusch-Pagan / White (heteroscedasticity), DW",
  # Extreme value distributions
  evd          = "GEV, GPD: fitting and quantiles for extreme values",
  fitdistrplus = "Distribution fitting by MLE/MME with diagnostics",
  # Non-stationary models
  gamlss       = "GAMLSS: location-scale-shape regression models",
  gamlss.dist  = "Additional distributions for gamlss (GEV, etc.)",
  gamlssx      = "Standard GEV for maxima in GAMLSS (GEVfisher/GEVquasi)",
  # Reports
  knitr        = "Reproducible report generation (Rmarkdown)",
  kableExtra   = "Formatted tables in HTML/LaTeX/console"
)

required_packages <- names(packages_info)
n_total  <- length(required_packages)
n_ok     <- 0L
n_new    <- 0L
n_error  <- 0L
failures <- character(0)

cat(sprintf("  Checking %d required packages...\n\n", n_total))

for (pkg in required_packages) {
  desc <- packages_info[[pkg]]
  if (requireNamespace(pkg, quietly = TRUE)) {
    ver <- tryCatch(
      as.character(packageVersion(pkg)),
      error = function(e) "?"
    )
    cat(sprintf("  [OK]  %-15s v%-8s  %s\n", pkg, ver, desc))
    n_ok <- n_ok + 1L
  } else {
    cat(sprintf("  [--]  %-15s %-9s  %s  → installing...",
                pkg, "", desc))
    result <- tryCatch({
      install.packages(
        pkg,
        repos        = "https://cloud.r-project.org",
        quiet        = TRUE,
        dependencies = TRUE
      )
      requireNamespace(pkg, quietly = TRUE)
    }, error = function(e) FALSE)
    
    if (result) {
      ver <- tryCatch(as.character(packageVersion(pkg)), error = function(e) "?")
      cat(sprintf(" ✔ v%s\n", ver))
      n_new    <- n_new + 1L
    } else {
      cat(" ✗ FAILED\n")
      n_error  <- n_error + 1L
      failures <- c(failures, pkg)
    }
  }
}

# ── Package installation summary ─────────────────────────────────────────────
cat(sprintf(
  "\n  ┌─────────────────────────────────────────┐\n"
))
cat(sprintf(
  "  │ Total checked     : %-3d                  │\n", n_total
))
cat(sprintf(
  "  │ Already installed : %-3d                  │\n", n_ok
))
cat(sprintf(
  "  │ Newly installed   : %-3d                  │\n", n_new
))
cat(sprintf(
  "  │ With errors       : %-3d                  │\n", n_error
))
cat(sprintf(
  "  └─────────────────────────────────────────┘\n"
))

if (n_error > 0) {
  cat(sprintf(
    "\n  ⚠ The following packages could not be installed:\n    %s\n",
    paste(failures, collapse = ", ")
  ))
  cat("  Try installing them manually with:\n")
  cat(sprintf('    install.packages(c("%s"))\n',
              paste(failures, collapse = '", "')))
  cat("  The script will continue but may fail in modules that require them.\n\n")
} else {
  cat("\n  ✔ All packages available.\n\n")
}

# ── Load critical packages individually (errors visible) ──────────────────────
# Statistical analysis packages are loaded one by one so that any failure
# is immediately visible and not masked by suppressMessages.
library(here)
library(evd)
library(gamlss)
library(gamlss.dist)
library(gamlssx)        # Standard GEV for maxima in GAMLSS (GEVfisher)
library(fitdistrplus)
library(tseries)

# ── Load remaining libraries silently ───────────────────────────────────────
suppressMessages({
  library(dplyr);       library(tidyr);       library(ggplot2)
  library(gridExtra);   library(grid)
  library(Kendall);     library(trend)
  library(lmtest);      library(strucchange)
  library(scales);      library(RColorBrewer)
  library(lubridate);   library(patchwork)
  library(writexl)
})

cat("✔ Packages loaded successfully.\n\n")

# ==============================================================================
# MODULE 01: USER PARAMETERS
# ==============================================================================

# ── Robust script folder detection ─────────────────────────────────────────
# Three-level detection strategy, in order of reliability:
#   1. sys.frames(): works when script is run with source() or Rscript --file=
#   2. rstudioapi: works when running from RStudio (Source / Run button)
#   3. getwd(): last resort if the two above cannot resolve the path
#
.detect_script_folder <- function() {
  # Level 1: source() or Rscript --file=
  for (i in seq_along(sys.frames())) {
    env <- sys.frames()[[i]]
    if (exists("ofile", envir = env, inherits = FALSE)) {
      path <- get("ofile", envir = env)
      if (!is.null(path) && nchar(path) > 0)
        return(normalizePath(dirname(path), winslash = "/", mustWork = FALSE))      
    }
  }
  # Level 2: active RStudio with open document
  if (requireNamespace("rstudioapi", quietly = TRUE) &&
      rstudioapi::isAvailable()) {
    ctx <- tryCatch(rstudioapi::getActiveDocumentContext(), error = function(e) NULL)
    if (!is.null(ctx) && nchar(ctx$path) > 0)
      return(normalizePath(dirname(ctx$path), winslash = "/", mustWork = FALSE))
  }
  # Level 3: working directory (fallback — may not match the script folder)
  warning(paste0(
    "Could not automatically detect the script folder.\n",
    "Using working directory: ", getwd(), "\n",
    "If files are not found, set SCRIPT_DIR manually."))
  normalizePath(getwd(), winslash = "/", mustWork = FALSE)
}

SCRIPT_DIR <- .detect_script_folder()
cat(sprintf("  Script folder detected: %s\n", SCRIPT_DIR))

PATH_PMAX  <- file.path(SCRIPT_DIR, CONFIG$file_pmax)
PATH_ONI   <- file.path(SCRIPT_DIR, CONFIG$file_oni)
OUTPUT_DIR <- file.path(SCRIPT_DIR, CONFIG$output_subfolder)

# ── Verify input files exist ─────────────────────────────────────────────────
# Provides clear error messages before attempting to read missing files.
if (!file.exists(PATH_PMAX)) {
  stop(sprintf(paste0(
    "\n  ✗ Precipitation file not found:\n",
    "    %s\n\n",
    "  Check that CONFIG$file_pmax is correct and the file\n",
    "  is in the script folder:\n",
    "    %s\n"
  ), PATH_PMAX, SCRIPT_DIR))
}
if (!file.exists(PATH_ONI)) {
  stop(sprintf(paste0(
    "\n  ✗ ONI file not found:\n",
    "    %s\n\n",
    "  Check CONFIG$file_oni. If you have ONI.xls, convert it with:\n",
    "  writexl::write_xlsx(readxl::read_xls('ONI.xls'), 'ONI.xlsx')\n\n",
    "  Folder searched: %s\n"
  ), PATH_ONI, SCRIPT_DIR))
}
cat(sprintf("  ✔ Input files verified\n"))
cat(sprintf("  ✔ Pmax : %s\n", PATH_PMAX))
cat(sprintf("  ✔ ONI  : %s\n", PATH_ONI))

# ── Design parameters ────────────────────────────────────────────────────────
# Script prompts interactively; uncomment below to run in batch mode:
# STORM_DURATION_H <- CONFIG$storm_duration_h

RETURN_PERIODS <- CONFIG$return_periods  # Return periods (years)

# Alpha threshold for hypothesis tests
ALPHA <- CONFIG$alpha_stat

# ==============================================================================
# MODULE 02: AUXILIARY FUNCTIONS
# ==============================================================================

# ── Colour palette ──────────────────────────────────────────────────────
COL <- list(
  azul    = "#1B4F72",
  celeste = "#2E86C1",
  verde   = "#1E8449",
  naranja = "#E67E22",
  rojo    = "#C0392B",
  gris    = "#7F8C8D",
  amarillo= "#F4D03F"
)

theme_base <- theme_minimal(base_size = 11) +
  theme(
    plot.title    = element_text(face = "bold", size = 13, colour = COL$azul),
    plot.subtitle = element_text(size = 10, colour = COL$gris),
    axis.title    = element_text(size = 10),
    legend.position = "bottom",
    panel.grid.minor = element_blank(),
    plot.background = element_rect(fill = "white", colour = NA),
    panel.background = element_rect(fill = "white", colour = NA)
  )

save_plot <- function(p, name, width = 10, height = 6) {
  path <- file.path(OUTPUT_DIR, "plots", paste0(name, ".png"))
  ggsave(path, plot = p, width = width, height = height, dpi = 180,
         bg = "white")
  invisible(path)
}

cat_sep <- function(title) {
  cat(paste0("\n", strrep("─", 70), "\n"))
  cat(paste0("  ", toupper(title), "\n"))
  cat(paste0(strrep("─", 70), "\n"))
}

# Helper: clean label for column/sheet names (T=2.33 → "T2_33")
tr_label <- function(T_r) gsub("\\.", "_", sprintf("T%g", T_r))

# ==============================================================================
# MODULE 03: DATA LOADING AND PREPARATION
# ==============================================================================
cat_sep("3. Data loading and preparation")

dir.create(file.path(OUTPUT_DIR, "plots"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(OUTPUT_DIR, "tables"),   recursive = TRUE, showWarnings = FALSE)

# ── Maximum 24h Precipitation ─────────────────────────────────────────────────
pmax_wide <- readxl::read_excel(PATH_PMAX, sheet = 1)
# Standardise first column name (year) and 12 month columns
colnames(pmax_wide)[1] <- "YEAR"
months_en <- c("JAN","FEB","MAR","APR","MAY","JUN",
               "JUL","AUG","SEP","OCT","NOV","DEC")
colnames(pmax_wide)[2:13] <- months_en
pmax_wide$YEAR <- as.integer(pmax_wide$YEAR)
# Force all month columns to numeric
pmax_wide[months_en] <- lapply(pmax_wide[months_en], function(x) as.numeric(as.character(x)))

pmax_long <- pmax_wide %>%
  pivot_longer(cols = all_of(months_en), names_to = "MONTH", values_to = "PMAX") %>%
  mutate(
    MONTH_NUM = match(MONTH, months_en),
    DATE    = as.Date(paste(YEAR, MONTH_NUM, 15, sep = "-"))
  ) %>%
  dplyr::filter(!is.na(PMAX))

# Annual maxima series (for trend analysis)
pmax_annual <- pmax_long %>%
  group_by(YEAR) %>%
  summarise(PMAX_MAX = max(PMAX, na.rm = TRUE),
            PMAX_MEAN = mean(PMAX, na.rm = TRUE),
            .groups = "drop")

# ── ONI Index ────────────────────────────────────────────────────────────────
# read_excel() automatically detects the format (.xls / .xlsx).
# More robust than read_xls() which fails on Linux with certain .xls (libxls).
oni_wide <- readxl::read_excel(PATH_ONI, sheet = 1)
# Standard headers: YEAR + English month abbreviations.
# ONI.xlsx already uses these names from conversion; reassigned by
# position as a safeguard against typographic variations in the
# original NOAA file.
months_oni <- c("JAN","FEB","MAR","APR","MAY","JUN",
                "JUL","AUG","SEP","OCT","NOV","DEC")
colnames(oni_wide)[1]    <- "YEAR"
colnames(oni_wide)[2:13] <- months_oni
oni_wide$YEAR <- as.integer(oni_wide$YEAR)
# Remove rows with YEAR=NA (repeated headers inside .xls)
oni_wide <- oni_wide[!is.na(oni_wide$YEAR), ]
# Force all month columns to numeric (some may be read as
# character if the cell contains text, spaces, or a repeated header)
oni_wide[months_oni] <- lapply(oni_wide[months_oni], function(x) as.numeric(as.character(x)))

oni_long <- oni_wide %>%
  pivot_longer(cols = all_of(months_oni), names_to = "MONTH_ONI", values_to = "ONI") %>%
  mutate(MONTH_NUM = match(MONTH_ONI, months_oni))

# Annual ONI: mean and maximum
oni_annual <- oni_long %>%
  group_by(YEAR) %>%
  summarise(ONI_MEAN = mean(ONI, na.rm = TRUE),
            ONI_MAX  = max(ONI, na.rm = TRUE),
            ONI_MIN  = min(ONI, na.rm = TRUE),
            .groups = "drop")

# ── Join series ──────────────────────────────────────────────────────────
pmax_data <- pmax_annual %>%
  inner_join(oni_annual, by = "YEAR") %>%
  arrange(YEAR)

cat(sprintf("  Records loaded: %d common years (%d – %d)\n",
            nrow(pmax_data), min(pmax_data$YEAR), max(pmax_data$YEAR)))
cat(sprintf("  Annual Pmax range: %.1f – %.1f mm\n",
            min(pmax_data$PMAX_MAX), max(pmax_data$PMAX_MAX)))
cat(sprintf("  ONI range:         %.1f – %.1f\n",
            min(pmax_data$ONI_MEAN), max(pmax_data$ONI_MEAN)))

# ── Input data quality validation ────────────────────────────────────────────
cat("\n  Validating input data quality...\n")
data_alerts <- character(0)
data_warnings <- character(0)

# 1. Minimum series length
n_years <- nrow(pmax_data)
if (n_years < 20) {
  data_alerts <- c(data_alerts, sprintf(
    "ALERT: Only %d common years. Minimum recommended: 20 years for reliable GEV.", n_years))
} else if (n_years < 30) {
  data_warnings <- c(data_warnings, sprintf(
    "WARNING: %d common years. With <30 years GEV tail uncertainty is high.", n_years))
}

# 2. Missing Pmax data
n_total_meses <- nrow(pmax_long) + sum(is.na(pmax_wide[months_en]))
n_na_pmax <- sum(is.na(pmax_wide[months_en]))
pct_na <- n_na_pmax / (nrow(pmax_wide) * 12) * 100
if (pct_na > 20) {
  data_alerts <- c(data_alerts, sprintf(
    "ALERT: %.1f%% of monthly Pmax data are missing (>20%%).", pct_na))
} else if (pct_na > 10) {
  data_warnings <- c(data_warnings, sprintf(
    "WARNING: %.1f%% of monthly Pmax data are missing.", pct_na))
}

# 3. Hydrological plausibility of Pmax
pmax_max_abs <- max(pmax_data$PMAX_MAX, na.rm = TRUE)
pmax_min_abs <- min(pmax_data$PMAX_MAX, na.rm = TRUE)
if (pmax_max_abs > 400) {
  data_alerts <- c(data_alerts, sprintf(
    "ALERT: Maximum Pmax = %.1f mm/day exceeds 400 mm — check record.", pmax_max_abs))
}
if (pmax_min_abs < 5) {
  data_alerts <- c(data_alerts, sprintf(
    "ALERT: Minimum annual Pmax = %.1f mm — possible erroneous records.", pmax_min_abs))
}

# 4. Unusual coefficient of variation
cv_pmax <- sd(pmax_data$PMAX_MAX) / mean(pmax_data$PMAX_MAX)
if (cv_pmax > 0.8) {
  data_warnings <- c(data_warnings, sprintf(
    "WARNING: CV of Pmax = %.2f (>0.8) — high variability; check homogeneity.", cv_pmax))
}

# 5. Ratio P(1h)/P(24h) outside range
r <- CONFIG$ratio_60_1440
if (r < 0.28 || r > 0.58) {
  data_alerts <- c(data_alerts, sprintf(
    "ALERT: ratio_60_1440 = %.2f outside plausible range [0.28, 0.58]. Verify.", r))
} else if (r < 0.32 || r > 0.52) {
  data_warnings <- c(data_warnings, sprintf(
    "WARNING: ratio_60_1440 = %.2f unusual. Typical Colombia range: 0.32–0.52. Verify with local data.", r))
}

# 6. Huff quartile / valid r coefficient
huff_q <- CONFIG$huff_quartile
if (!is.null(CONFIG$huff_r_median)) {
  r_val <- CONFIG$huff_r_median
  if (!is.numeric(r_val) || r_val < 0 || r_val > 1) {
    data_alerts <- c(data_alerts,
                     "ALERT: huff_r_median must be numeric in [0,1]. huff_quartile will be used.")
  } else if (r_val < 0.10 || r_val > 0.80) {
    data_warnings <- c(data_warnings, sprintf(
      "WARNING: huff_r_median = %.3f outside usual range [0.10, 0.80].", r_val))
  } else {
    data_warnings <- c(data_warnings, sprintf(
      "INFO: huff_r_median = %.3f detected — Huff quartile will be selected automatically.", r_val))
  }
} else if (!huff_q %in% 1:4) {
  data_alerts <- c(data_alerts,
                   "ALERT: huff_quartile must be 1, 2, 3 or 4. Defaulting to 2.")
}

# ── Print validation result ─────────────────────────────────────────
if (length(data_alerts) == 0 && length(data_warnings) == 0) {
  cat("  ✔ Input data validated. No issues detected.\n\n")
} else {
  if (length(data_warnings) > 0) {
    cat("\n  ⚠ WARNINGS (execution continues):\n")
    for (av in data_warnings) cat(sprintf("    · %s\n", av))
  }
  if (length(data_alerts) > 0) {
    cat("\n  ✗ ALERTS (review before using results in design):\n")
    for (al in data_alerts) cat(sprintf("    · %s\n", al))
  }
  cat("\n")
}
# ==============================================================================
# MODULE 04: EXPLORATORY VISUALISATION
# ==============================================================================
cat_sep("4. Exploratory visualisation")

# ── 4.1 Pmax time series ────────────────────────────────────────────────────
p_series_pmax <- ggplot(pmax_data, aes(YEAR, PMAX_MAX)) +
  geom_line(colour = COL$celeste, linewidth = 0.8) +
  geom_point(aes(colour = ONI_MEAN), size = 2) +
  geom_smooth(method = "lm", se = TRUE, colour = COL$rojo,
              linetype = "dashed", linewidth = 0.9) +
  scale_colour_gradient2(low = COL$celeste, mid = "white",
                         high = COL$naranja, midpoint = 0,
                         name = "Mean ONI") +
  labs(title = "Historical 24-hour Maximum Daily Precipitation",
       subtitle = "Points coloured by ONI index; red line = linear trend",
       x = "Year", y = "Pmax 24h (mm)") +
  theme_base
save_plot(p_series_pmax, "01_series_pmax")

# ── 4.2 Historical ONI ─────────────────────────────────────────────────────────
p_oni <- ggplot(pmax_data, aes(YEAR, ONI_MEAN)) +
  geom_hline(yintercept = c(-0.5, 0.5), linetype = "dashed",
             colour = COL$gris, linewidth = 0.5) +
  geom_bar(stat = "identity",
           fill = ifelse(pmax_data$ONI_MEAN > 0, COL$naranja, COL$celeste),
           alpha = 0.8) +
  labs(title = "Oceanic Niño Index (ONI) — Annual Mean",
       subtitle = "Orange = El Niño  |  Blue = La Niña",
       x = "Year", y = "Mean ONI") +
  theme_base
save_plot(p_oni, "02_oni_historical")

# ── 4.3 Scatter Pmax vs ONI ────────────────────────────────────────────────
p_disp <- ggplot(pmax_data, aes(ONI_MEAN, PMAX_MAX)) +
  geom_point(colour = COL$celeste, size = 2, alpha = 0.7) +
  geom_smooth(method = "lm", se = TRUE, colour = COL$rojo,
              linewidth = 0.9) +
  geom_smooth(method = "loess", se = FALSE, colour = COL$verde,
              linewidth = 0.9, linetype = "dotted") +
  labs(title = "24-hour Maximum Precipitation vs. ONI",
       subtitle = "Red = linear regression  |  Green = LOESS",
       x = "Annual mean ONI", y = "Pmax 24h (mm)") +
  theme_base
save_plot(p_disp, "03_scatter_pmax_oni")

cat("  ✔ Exploratory plots generated\n")

# ==============================================================================
# MODULE 05: STATIONARITY / NON-STATIONARITY ANALYSIS
# ==============================================================================
cat_sep("5. Stationarity / non-stationarity analysis")

series <- pmax_data$PMAX_MAX
n <- length(series)
years <- pmax_data$YEAR

test_results <- list()

# ──────────────────────────────────────────────────────────────────────────────
# 5.1  Mann-Kendall (MK) — classical
# ──────────────────────────────────────────────────────────────────────────────
cat("\n  [5.1] Mann-Kendall — classical\n")
mk_res  <- MannKendall(series)
mk_tau  <- as.numeric(mk_res$tau)
mk_pval <- as.numeric(mk_res$sl)
mk_S    <- as.numeric(mk_res$S)
mk_varS <- as.numeric(mk_res$varS)
mk_z    <- mk_S / sqrt(mk_varS)

test_results$MK <- data.frame(
  Test = "Mann-Kendall",
  Statistic = round(mk_tau, 4),
  Statistic_name = "tau",
  p_value = round(mk_pval, 4),
  Conclusion = ifelse(mk_pval < ALPHA,
                      "SIGNIFICANT TREND (Non-stationary)",
                      "No significant trend (Stationary)")
)
cat(sprintf("    τ = %.4f  |  p-value = %.4f  →  %s\n",
            mk_tau, mk_pval, test_results$MK$Conclusion))

# ──────────────────────────────────────────────────────────────────────────────
# 5.2  Pettitt test (change-point detection)
# ──────────────────────────────────────────────────────────────────────────────
cat("\n  [5.2] Pettitt test\n")
pett_res  <- pettitt.test(series)
pett_pval <- pett_res$p.value
pett_U    <- pett_res$statistic
pett_t    <- pett_res$estimate   # change-point position

test_results$Pettitt <- data.frame(
  Test = "Pettitt",
  Statistic = round(as.numeric(pett_U), 4),
  Statistic_name = "U",
  p_value = round(pett_pval, 4),
  Conclusion = ifelse(pett_pval < ALPHA,
                      paste0("ABRUPT CHANGE detected in year ",
                             years[as.numeric(pett_t)], " (No stationary)"),
                      "No significant abrupt change (Stationary)")
)
cat(sprintf("    U = %.0f  |  p-value = %.4f  |  Change-point estimated: year %d\n    → %s\n",
            as.numeric(pett_U), pett_pval, years[as.numeric(pett_t)],
            test_results$Pettitt$Conclusion))

# ── Plot Pettitt ──────────────────────────────────────────────────────────
# Vectorised U_t: for each k, counts how many previous values are smaller.
uk_vals <- sapply(seq_len(n), function(k) sum(sign(series[k] - series[-k])))
uk_cum  <- cumsum(uk_vals)

df_pett <- data.frame(YEAR = years, U_t = uk_cum)
break_year <- years[which.max(abs(uk_cum))]

p_pettitt <- ggplot(df_pett, aes(YEAR, U_t)) +
  geom_line(colour = COL$azul, linewidth = 1) +
  geom_vline(xintercept = break_year,
             colour = COL$rojo, linetype = "dashed", linewidth = 1) +
  annotate("text", x = break_year + 1, y = max(uk_cum) * 0.9,
           label = paste0("Break ≈ ", break_year),
           colour = COL$rojo, hjust = 0, size = 3.5) +
  labs(title = "Pettitt Test — U_t statistic",
       subtitle = "Red line marks the point of maximum change",
       x = "Year", y = "Cumulative U_t") +
  theme_base
save_plot(p_pettitt, "04_pettitt")

# ──────────────────────────────────────────────────────────────────────────────
# 5.3  Sequential Sneyers (cumulative curvature analysis)
# ──────────────────────────────────────────────────────────────────────────────
cat("\n  [5.3] Sequential Sneyers\n")

# Progressive statistic u(t) and retrogressive u'(t)
prog_u <- function(x) {
  n <- length(x)
  S <- rep(0, n)
  for (k in 2:n) S[k] <- S[k-1] + sum(x[k] > x[1:(k-1)])
  # Variance of the Mann-Kendall S statistic (Kendall, 1975), used by
  # Sneyers (1990) in his sequential progressive/regressive test.
  # Ref: Sneyers, R. (1990). On the Statistical Analysis of Series of
  # Observations. WMO Technical Note No. 143. Geneva.
  u <- (S - (1:n) * ((1:n)-1) / 4) /
    sqrt((1:n) * ((1:n)-1) * (2*(1:n)+5) / 72)
  u[1] <- 0
  u
}

u_prog <- prog_u(series)
u_retr <- rev(prog_u(rev(series)))

df_sneyers <- data.frame(
  YEAR   = years,
  u_prog = u_prog,
  u_retr = -u_retr   # convention: retrogressive with opposite sign
)

# Detect crossings between u and u'
crossings <- which(diff(sign(u_prog - (-u_retr))) != 0)
crossing_years <- if (length(crossings) > 0) years[crossings] else NA

test_results$Sneyers <- data.frame(
  Test = "Sequential Sneyers",
  Statistic = round(tail(u_prog, 1), 4),
  Statistic_name = "u(t) final",
  p_value = NA,
  Conclusion = ifelse(any(abs(u_prog) > 1.96),
                      paste0("Possible change point. Crossings at: ",
                             paste(crossing_years, collapse = ", ")),
                      "No significant crossings — stationary behaviour")
)
cat(sprintf("    u(t) final = %.3f  |  Crossings detected: %s\n    → %s\n",
            tail(u_prog, 1),
            ifelse(all(is.na(crossing_years)), "none",
                   paste(crossing_years, collapse = ", ")),
            test_results$Sneyers$Conclusion))

p_sneyers <- ggplot(df_sneyers) +
  geom_hline(yintercept = c(-1.96, 1.96),
             linetype = "dashed", colour = COL$gris) +
  geom_hline(yintercept = 0, colour = "black", linewidth = 0.3) +
  geom_line(aes(YEAR, u_prog, colour = "Progressive"), linewidth = 1) +
  geom_line(aes(YEAR, u_retr, colour = "Retrogressive"),  linewidth = 1) +
  scale_colour_manual(values = c("Progressive" = COL$azul,
                                 "Retrogressive"  = COL$naranja),
                      name = "") +
  labs(title = "Sequential Sneyers — Progressive/retrogressive statistics",
       subtitle = "Dashed bands = 95% bounds (±1.96)",
       x = "Year", y = "u(t)") +
  theme_base
save_plot(p_sneyers, "05_sneyers")

# ──────────────────────────────────────────────────────────────────────────────
# 5.4  Moving-Window Mann-Kendall (sliding window)
# ──────────────────────────────────────────────────────────────────────────────
cat("\n  [5.4] Moving-Window Mann-Kendall\n")

window_size  <- max(10, round(n * 0.3))   # 30% of the series or minimum 10
tau_mw   <- numeric(n)
pval_mw  <- numeric(n)
years_mw <- numeric(n)

for (i in window_size:n) {
  sub_series <- series[(i - window_size + 1):i]
  res_mw    <- MannKendall(sub_series)
  tau_mw[i]  <- as.numeric(res_mw$tau)
  pval_mw[i] <- as.numeric(res_mw$sl)
  years_mw[i]<- years[i]
}

df_mw <- data.frame(
  YEAR  = years_mw[window_size:n],
  tau   = tau_mw[window_size:n],
  pval  = pval_mw[window_size:n],
  sig   = pval_mw[window_size:n] < ALPHA
)

n_sig_mw <- sum(df_mw$sig)
test_results$MWMK <- data.frame(
  Test = "Moving-Window MK",
  Statistic = round(mean(df_mw$tau), 4),
  Statistic_name = "mean tau",
  p_value = NA,
  Conclusion = ifelse(n_sig_mw > nrow(df_mw) * 0.3,
                      paste0("Variable trend detected (", n_sig_mw,
                             " significant windows out of ", nrow(df_mw), ")"),
                      "No persistent trend in sliding windows")
)
cat(sprintf("    Significant windows: %d / %d  →  %s\n",
            n_sig_mw, nrow(df_mw),
            test_results$MWMK$Conclusion))

p_mw <- ggplot(df_mw, aes(YEAR, tau, colour = sig)) +
  geom_hline(yintercept = 0, colour = "black", linewidth = 0.3) +
  geom_line(colour = COL$azul, linewidth = 0.7) +
  geom_point(size = 2.5) +
  scale_colour_manual(values = c("FALSE" = COL$gris, "TRUE" = COL$rojo),
                      labels = c("Not significant", "Significant"),
                      name = "Trend") +
  labs(title = paste0("Moving-Window Mann-Kendall (window = ", window_size, " years)"),
       subtitle = "Red dots: windows with statistically significant trend",
       x = "Year (right edge of window)", y = "Kendall τ") +
  theme_base
save_plot(p_mw, "06_moving_window_mk")

# ──────────────────────────────────────────────────────────────────────────────
# 5.5  White test (heteroscedasticity)
# ──────────────────────────────────────────────────────────────────────────────
cat("\n  [5.5] White test (heteroscedasticity)\n")

lm_model    <- lm(PMAX_MAX ~ YEAR, data = pmax_data)
white_res    <- bptest(lm_model, ~ fitted(lm_model) + I(fitted(lm_model)^2))
white_chi2   <- as.numeric(white_res$statistic)
white_pval   <- as.numeric(white_res$p.value)
white_df     <- as.numeric(white_res$parameter)

test_results$White <- data.frame(
  Test = "White (heteroscedasticity)",
  Statistic = round(white_chi2, 4),
  Statistic_name = "BP (chi2)",
  p_value = round(white_pval, 4),
  Conclusion = ifelse(white_pval < ALPHA,
                      "SIGNIFICANT HETEROSCEDASTICITY (non-constant variance)",
                      "No heteroscedasticity — homogeneous variance")
)
cat(sprintf("    BP = %.4f  (df=%d)  |  p-value = %.4f  →  %s\n",
            white_chi2, white_df, white_pval,
            test_results$White$Conclusion))

# Residuals plot
df_res <- data.frame(
  YEAR     = years,
  Residual = residuals(lm_model),
  Fitted   = fitted(lm_model)
)
p_white <- ggplot(df_res, aes(Fitted, Residual)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = COL$gris) +
  geom_point(colour = COL$celeste, size = 2, alpha = 0.7) +
  geom_smooth(method = "loess", se = FALSE,
              colour = COL$rojo, linewidth = 0.9) +
  labs(title = "Heteroscedasticity Diagnostic — White Test",
       subtitle = "Residuals vs. fitted values of the linear regression",
       x = "Fitted values (mm)", y = "Residuals (mm)") +
  theme_base
save_plot(p_white, "07_white_heteroscedasticity")

# ──────────────────────────────────────────────────────────────────────────────
# 5.6  Augmented Dickey-Fuller (ADF) – unit root
# ──────────────────────────────────────────────────────────────────────────────
cat("\n  [5.6] Augmented Dickey-Fuller\n")

adf_res    <- adf.test(series)
adf_stat   <- as.numeric(adf_res$statistic)
adf_pval   <- as.numeric(adf_res$p.value)
adf_lag    <- as.numeric(adf_res$parameter)

test_results$ADF <- data.frame(
  Test = "ADF (unit root)",
  Statistic = round(adf_stat, 4),
  Statistic_name = "Dickey-Fuller",
  p_value = round(adf_pval, 4),
  Conclusion = ifelse(adf_pval < ALPHA,
                      "STATIONARY (unit root rejected)",
                      "Possible unit root — non-stationary behaviour")
)
cat(sprintf("    ADF = %.4f  (lags=%d)  |  p-value = %.4f  →  %s\n",
            adf_stat, adf_lag, adf_pval,
            test_results$ADF$Conclusion))

# ==============================================================================
# MODULE 06: TEST SUMMARY AND STATIONARITY DECISION
# ==============================================================================
cat_sep("6. Test summary and decision")

summary_table <- bind_rows(test_results)
print(summary_table[, c("Test","Statistic","p_value","Conclusion")],
      row.names = FALSE)

# ── Automatic decision based on majority vote ───────────────────────────────────
non_stat_votes <- sum(c(
  mk_pval < ALPHA,                         # MK
  pett_pval < ALPHA,                        # Pettitt
  any(abs(u_prog) > 1.96),                  # Sneyers
  n_sig_mw > nrow(df_mw) * 0.3,            # MWMK
  white_pval < ALPHA,                       # White
  adf_pval >= ALPHA                         # ADF (does not reject unit root)
))
total_votes <- 6

cat(sprintf("\n  Votes for NON-STATIONARITY: %d / %d\n", non_stat_votes, total_votes))

IS_NON_STATIONARY <- non_stat_votes >= 3  # simple majority

if (IS_NON_STATIONARY) {
  cat("\n  ★ DECISION: NON-STATIONARY series")
  cat("\n    Non-stationary MCM model with ONI covariate will be used.\n")
} else {
  cat("\n  ★ DECISION: STATIONARY series")
  cat("\n    Stationary MCM model will be used.\n")
}

# ==============================================================================
# MODULE 07: FREQUENCY DISTRIBUTION FITTING
# ==============================================================================
cat_sep("7. Frequency distribution fitting")

if (!IS_NON_STATIONARY) {
  
  # ── Stationary model: classical GEV ────────────────────────────────────────
  cat("\n  Fitting stationary GEV (annual maxima)...\n")
  gev_fit      <- fgev(series)
  gev_fit_base <- gev_fit   # saved for reuse in Module 14b
  loc0    <- gev_fit$param["loc"]
  sc0     <- gev_fit$param["scale"]
  sh0     <- gev_fit$param["shape"]
  
  cat(sprintf("    μ = %.3f mm  |  σ = %.3f mm  |  ξ = %.4f\n",
              loc0, sc0, sh0))
  
  # Quantiles for return periods
  prob_exc <- 1 - 1 / RETURN_PERIODS
  q_pmax   <- qgev(prob_exc, loc = loc0, scale = sc0, shape = sh0)
  df_tr    <- data.frame(T_return = RETURN_PERIODS, Pmax_mm = round(q_pmax, 2))
  
} else {
  
  # ── Non-stationary model: GAMLSS with ONI covariate ────────────────────────
  cat("\n  Fitting non-stationary GEV (GAMLSS, ONI covariate)...\n")
  
  df_gamlss <- data.frame(
    y   = series,
    oni = pmax_data$ONI_MEAN
  )
  
  # NOTE on gamlss.dist vs gamlssx:
  # gamlss.dist v6.x does not export GEV() with standard maxima parameterisation.
  # We use gamlssx::GEVfisher() which implements GEV(μ, σ, ξ) for block maxima.
  # Ref: gamlssx package (CRAN), Rigby & Stasinopoulos (2005).
  # Model with location parameter as linear function of ONI
  fit_ns <- NULL   # declared before tryCatch to avoid superassignment (<<-)
  tryCatch({
    fit_ns <- gamlss(
      y ~ 1,
      sigma.formula = ~ 1,
      nu.formula    = ~ 1,
      mu.formula    = ~ oni,
      data    = df_gamlss,
      family  = gamlssx::GEVfisher(), # Standard GEV for maxima (gamlssx package)
      control = gamlss.control(n.cyc = 100, trace = FALSE)
    )
    cat("    ✔ Non-stationary GEV fitted with μ ~ ONI\n")
    cat(sprintf("    AIC = %.2f  |  BIC = %.2f\n", AIC(fit_ns), BIC(fit_ns)))
  }, error = function(e) {
    cat("    ⚠ Non-stationary GEV failed:", conditionMessage(e), "\n")
    cat("      Falling back to stationary GEV...\n")
    # fit_ns is already NULL (declared above), <<- is not used
  })
  
  # Predict quantiles for mean ONI (neutral climate scenario)
  oni_escenarios <- data.frame(
    scenario = c("La Nina (ONI=-1.0)", "Neutral (ONI=0.0)", "El Nino (ONI=+1.0)"),
    oni = c(-1.0, 0.0, 1.0)
  )
  
  if (!is.null(fit_ns)) {
    # Extract coefficients
    mu_int   <- fit_ns$mu.coefficients[1]
    mu_oni   <- if (length(fit_ns$mu.coefficients) > 1) fit_ns$mu.coefficients[2] else 0
    sigma_fit <- exp(fit_ns$sigma.coefficients[1])
    nu_fit    <- fit_ns$nu.coefficients[1]
    
    prob_exc <- 1 - 1 / RETURN_PERIODS
    
    df_tr_list <- lapply(seq_along(oni_escenarios$oni), function(s) {
      mu_s  <- mu_int + mu_oni * oni_escenarios$oni[s]
      q_s   <- qGEV(prob_exc, mu = mu_s, sigma = sigma_fit, nu = nu_fit)
      data.frame(
        Scenario = oni_escenarios$scenario[s],
        T_return = RETURN_PERIODS,
        Pmax_mm   = round(q_s, 2)
      )
    })
    df_tr_ns <- bind_rows(df_tr_list)
    
    # Use neutral scenario (ONI=0) for hyetograph generation
    q_pmax_neutral <- df_tr_ns %>%
      dplyr::filter(Scenario == "Neutral (ONI=0.0)") %>%
      dplyr::pull(Pmax_mm)
    df_tr <- data.frame(T_return = RETURN_PERIODS, Pmax_mm = q_pmax_neutral)
    
  } else {
    # Stationary fallback
    gev_fit_base <- fgev(series)
    loc0 <- gev_fit_base$param["loc"]
    sc0  <- gev_fit_base$param["scale"]
    sh0  <- gev_fit_base$param["shape"]
    prob_exc <- 1 - 1 / RETURN_PERIODS
    q_pmax   <- qgev(prob_exc, loc = loc0, scale = sc0, shape = sh0)
    df_tr    <- data.frame(T_return = RETURN_PERIODS, Pmax_mm = round(q_pmax, 2))
    df_tr_ns <- NULL
  }
}

# ── gev_fit_base: stationary reference for reports ──────────────────────────
# Always present regardless of path taken in Module 7.
# Does not affect design quantiles (which come from the NS model when applicable).
if (!exists("gev_fit_base")) {
  cat("  Fitting base GEV (stationary reference for report)...\n")
  gev_fit_base <- evd::fgev(series)
  cat(sprintf("    μ = %.3f | σ = %.3f | ξ = %.4f  (reference only; design uses NS model)\n",
              gev_fit_base$param["loc"],
              gev_fit_base$param["scale"],
              gev_fit_base$param["shape"]))
}

cat("\n  Design precipitation (neutral scenario or stationary):\n")
print(df_tr, row.names = FALSE)

# ── Distribution and quantile plot ──────────────────────────────────────────
df_empir <- data.frame(
  PMAX = sort(series),
  T_ret = 1 / (1 - (1:n) / (n + 1))
)

p_dist <- ggplot(df_empir, aes(T_ret, PMAX)) +
  geom_point(colour = COL$celeste, size = 2) +
  geom_line(data = data.frame(
    T_ret = RETURN_PERIODS,
    PMAX  = df_tr$Pmax_mm),
    aes(T_ret, PMAX), colour = COL$rojo, linewidth = 1) +
  geom_point(data = df_tr,
             aes(T_return, Pmax_mm),
             colour = COL$rojo, size = 3, shape = 18) +
  scale_x_log10(breaks = c(2, 5, 10, 25, 50, 100, 200, 500)) +
  labs(title = "Frequency Curve — 24-hour Maximum Precipitation",
       subtitle = "Blue dots = empirical data  |  Red line = fitted GEV",
       x = "Return period (years)", y = "Pmax 24h (mm)") +
  theme_base
save_plot(p_dist, "08_frequency_curve")
# ==============================================================================
# MODULE 08: DESIGN STORM PARAMETERS
# ==============================================================================
cat_sep("8. Design storm parameters")

# ── Interactive prompt: DURATION and INTERVAL ──────────────────────────────
if (exists("STORM_DURATION_H")) {
  Td <- STORM_DURATION_H
  cat(sprintf("  Design duration (batch mode): %g hours\n", Td))
} else {
  cat("\n  ENTER THE DESIGN STORM DURATION (in hours).\n")
  cat("  Typical values: 1, 2, 3, 6, 12, 24\n")
  cat("  Duration [hours]: ")
  input_raw <- readLines(con = stdin(), n = 1)
  Td <- as.numeric(trimws(input_raw))
  if (is.na(Td) || Td <= 0) {
    cat("  ⚠ Invalid input. Defaulting to 6 hours.\n")
    Td <- 6
  }
}

if (exists("TIME_STEP_MIN")) {
  dt_min <- TIME_STEP_MIN
  cat(sprintf("  Time step (batch mode): %g minutes\n", dt_min))
} else {
  cat("\n  ENTER THE TIME DISCRETISATION INTERVAL (in minutes).\n")
  cat("  Typical values: 5, 10, 15, 30, 60\n")
  cat("  Note: must divide the storm duration in minutes exactly.\n")
  cat("  Interval [minutes]: ")
  input_dt <- readLines(con = stdin(), n = 1)
  dt_min <- as.numeric(trimws(input_dt))
  if (is.na(dt_min) || dt_min <= 0) {
    cat("  ⚠ Invalid input. Defaulting to 15 minutes.\n")
    dt_min <- 15
  }
}

# Verify dt_min exactly divides the total duration in minutes
Td_min  <- Td * 60
if ((Td_min %% dt_min) != 0) {
  dt_min_orig <- dt_min
  # Adjust to the nearest integer divisor below
  divisors <- which(Td_min %% (1:Td_min) == 0)
  dt_min <- max(divisors[divisors <= dt_min_orig])
  cat(sprintf("  ⚠ %g min does not exactly divide %g h. Adjusted to %g min.\n",
              dt_min_orig, Td, dt_min))
}

N_steps <- as.integer(Td_min / dt_min)
if (N_steps < 2) {
  cat("  ⚠ Duration too short for the chosen interval. Adjusting to 4 minimum steps.\n")
  N_steps <- 4L
  dt_min  <- Td_min / N_steps
}

cat(sprintf("\n  Storm duration     : %g hours\n", Td))
cat(sprintf("  Time step         : %g minutes\n", dt_min))
cat(sprintf("  Number of steps   : %d\n", N_steps))

# ==============================================================================
# MODULE 09: IDF CURVES — BELL (1969) METHOD WITH STATION ANCHOR
# ==============================================================================
cat_sep("9. IDF curves — Bell (1969) anchored on station GEV quantiles")

# ── Rationale ──────────────────────────────────────────────────────────────────
# Bell (1969) expresses precipitation of duration t from P(60 min, T):
#
#   P(t, T) = (0.54 · t^0.25 - 0.50) · P(60min, T)     [t in minutes, P in mm]
#
# The ANCHOR P(60min, T) is calculated directly from station GEV quantiles:
#   P(60min, T) = P(24h, T) · ratio_60_1440
#
# where ratio_60_1440 = P(1h)/P(24h) is the only regional parameter needed.
# Everything else comes from station data.
# Intensity is simply I(t,T) = P(t,T) / (t/60).
#
# This formulation produces curves with REAL CURVATURE: function (0.54*t^0.25-0.50)
# is not a pure power law, so on log-log scale the curves appear
# as curves, not straight lines.
#
# Reference: Bell, F.C. (1969). Generalized rainfall-duration-frequency
#   relationships. J. Hydraul. Div. ASCE, 95(1), 311–327.

ratio_60_1440 <- CONFIG$ratio_60_1440   # P(1h)/P(24h) — only regional parameter

cat(sprintf("  Anchor: P(1h)/P(24h) = %.2f\n", ratio_60_1440))
cat(sprintf("  Formula: P(t,T) = (0.54·t^0.25 - 0.50) · P(60min,T)   [t in min]\n"))
cat(sprintf("  P(60min,T) = P(24h,T) · %.2f  [from station GEV quantiles]\n\n",
            ratio_60_1440))

# ── IDF function ───────────────────────────────────────────────────────────────
idf_func <- function(dur_min, Tr) {
  # Bell (1969): I(t,T) = P(60min,T) * (0.54*t^0.25 - 0.50) / f60 / (t/60)
  # Vectorized: dur_min and Tr can be vectors of the same length.
  P24_T <- df_tr$Pmax_mm[match(Tr, df_tr$T_return)]
  P60_T <- P24_T * ratio_60_1440
  f60   <- 0.54 * 60^0.25 - 0.50
  (P60_T * (0.54 * dur_min^0.25 - 0.50) / f60) / (dur_min / 60)
}

idf_P_mm <- function(dur_min, Tr) {
  idf_func(dur_min, Tr) * (dur_min / 60)
}

# ── Intensity table ─────────────────────────────────────────────────────
dur_table <- c(15, 30, 60, 120, 180, 360, 720, 1440)
cat("  Intensities I(mm/h):\n")
cat(sprintf("  %-6s", "T(a)"))
for (d in dur_table) cat(sprintf("  %5dmin", d))
cat("\n")
for (Tr in RETURN_PERIODS) {
  cat(sprintf("  %-6g", Tr))
  for (d in dur_table) cat(sprintf("  %8.1f", idf_func(d, Tr)))
  cat("\n")
}

# ── Plot IDF — linear scale, range 5-360 min (visible hyperbolic shape) ──
dur_idf_min <- seq(5, 360, by = 1)   # 1 min resolution, linear scale

df_idf <- expand.grid(T_return = RETURN_PERIODS, dur_min = dur_idf_min) %>%
  mutate(
    I_mm_h   = idf_func(dur_min, T_return),
    P_dur_mm = I_mm_h * (dur_min / 60)
  )

df_idf_pts <- expand.grid(T_return = RETURN_PERIODS, dur_min = dur_table) %>%
  mutate(I_mm_h = idf_func(dur_min, T_return)) %>%
  filter(dur_min <= 360)

p_idf <- ggplot(df_idf,
                aes(dur_min, I_mm_h,
                    colour = as.factor(T_return),
                    group  = T_return)) +
  geom_line(linewidth = 0.9) +
  geom_point(data = df_idf_pts,
             aes(dur_min, I_mm_h),
             size = 2.5, shape = 21, fill = "white") +
  scale_colour_brewer(palette = "Spectral", name = "Return period (years)") +
  scale_x_continuous(
    breaks = c(15, 30, 60, 120, 180, 360),
    labels = c("15min", "30min", "1h", "2h", "3h", "6h")
  ) +
  scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.05))) +
  labs(
    title = "Intensity-Duration-Frequency (IDF) Curves",
    subtitle = sprintf(
      "Bell (1969) | P(t,T) = (0.54·t^0.25-0.50)·P(60min,T) | P(1h)/P(24h) = %.2f (see CONFIG for regional guide)",
      ratio_60_1440),
    x = "Duration (minutes)", y = "Intensity (mm/h)"
  ) +
  theme_base +
  theme(panel.grid.major = element_line(colour = "grey92"))
save_plot(p_idf, "09_IDF_curves", width = 12, height = 7)
cat("  ✔ IDF curves generated\n")

# ==============================================================================
# MODULE 10: MCM MODEL — MICROCANONICAL MULTIPLICATIVE CASCADE
# ==============================================================================
cat_sep("10. MCM model — Hyetograph generation")

# ──────────────────────────────────────────────────────────────────────────────
# Theoretical foundation:
# The MCM disaggregates total precipitation through multipliers W at each
# level of the binary cascade. For each pair of sibling intervals (i, i+1):
#   P_i     = W_i * P_parent
#   P_{i+1} = (1-W_i) * P_parent
# where W ~ Beta(α, β) or mixed distribution with point mass at 0 and 1.
#
# Parameters (α, β) are estimated by MLE from the historical monthly
# precipitation series by sub-aggregating from 24h → 12h → ... → Δt
# ──────────────────────────────────────────────────────────────────────────────

# ── 10.1  MCM parameter estimation ────────────────────────────────────────
cat("\n  [10.1] MCM parameter estimation...\n")

# Build pairs (P_1h, P_2h, ...) from monthly pmax
# Use monthly values as base (24h resolution) and calculate
# Beta distribution parameters for disaggregation

monthly_series <- pmax_long$PMAX

# Beta parameter estimation function by maximum likelihood
# The W multipliers are generated from the observed asymmetry
# in sub-aggregations of synthetic hourly data (adapted for 24h data)

# A priori parameters calibrated for tropical-Andean regions
# Reference: Molnar & Burlando (2005), Hingray & Ben Haha (2005)
estimate_mcm_params <- function(series_pmax) {
  # Normalise: values strictly in (0, 1) for Beta fitting
  x <- series_pmax / max(series_pmax, na.rm = TRUE)
  x <- x[x > 0 & x < 1]   # exclude exactly 0 and 1 (outside Beta support)
  
  # Beta fitting by MLE with robust error handling
  fit_b <- tryCatch({
    if (length(x) < 10) stop("Fewer than 10 valid values for Beta fitting")
    suppressWarnings(fitdist(x, "beta", method = "mle"))
  }, error = function(e) {
    cat(sprintf(
      "    ⚠ fitdist Beta failed (%s).\n    → Using regional a priori parameters (α=1.5, β=1.5).\n",
      conditionMessage(e)))
    # A priori parameters for tropical-Andean regions
    # Equivalent to moderate symmetric Beta distribution (Molnar & Burlando, 2005)
    list(estimate = c(shape1 = 1.5, shape2 = 1.5))
  })
  
  list(
    alpha = unname(fit_b$estimate["shape1"]),
    beta  = unname(fit_b$estimate["shape2"]),
    P0    = mean(series_pmax == 0, na.rm = TRUE),  # dry fraction
    Pmax  = max(series_pmax, na.rm = TRUE)
  )
}

mcm_params <- estimate_mcm_params(monthly_series)
cat(sprintf("    α = %.4f  |  β = %.4f\n", mcm_params$alpha, mcm_params$beta))
cat(sprintf("    Dry fraction estimated P0 = %.3f\n", mcm_params$P0))
cat("\n  ⚠ METHODOLOGICAL NOTE (MCM):\n")
cat("    Parameters α and β calibrated from monthly Pmax data.\n")
cat("    Disaggregation to sub-hourly scale involves a jump of ~4 orders\n")
cat("    of magnitude. Results are plausible but parameters are\n")
cat("    regional a priori values, not locally calibrated at minute scale.\n")
cat("    For greater accuracy: recalibrate with hourly rain gauge data.\n\n")

# ── 10.2  MCM disaggregation function ────────────────────────────────────────
cat("\n  [10.2] Defining MCM disaggregation function...\n")
# ── Automatic quartile selection from median r (if available) ────────────────
huff_names <- c(
  "1" = "Q1 — Impulsive convective / high orography  (r: 0.10–0.25)",
  "2" = "Q2 — Andean convective / urban              (r: 0.25–0.45) [default]",
  "3" = "Q3 — Mixed / Pacific                        (r: 0.45–0.60)",
  "4" = "Q4 — Stratiform / Eastern Plains            (r: 0.60–0.80)"
)

r_med <- CONFIG$huff_r_median

if (!is.null(r_med) && is.numeric(r_med) && r_med >= 0 && r_med <= 1) {
  # Automatic selection from observed median r
  huff_q_auto <- if      (r_med < 0.25) 1L
  else if (r_med < 0.45) 2L
  else if (r_med < 0.60) 3L
  else                   4L
  cat(sprintf("  Input median r         : %.3f\n", r_med))
  cat(sprintf("  Huff quartile (auto)   : %s\n", huff_names[as.character(huff_q_auto)]))
  cat("  Ref: Victoria Nino (2026), DOI: 10.31224/7062\n")
  # Overwrite quartile in CONFIG for later use
  CONFIG$huff_quartile <- huff_q_auto
} else {
  # Fixed quartile from CONFIG
  huff_q_auto <- max(1L, min(4L, as.integer(round(CONFIG$huff_quartile))))
  CONFIG$huff_quartile <- huff_q_auto
  cat(sprintf("  Huff quartile (manual) : %s\n", huff_names[as.character(huff_q_auto)]))
  cat("  (For automatic selection: provide median r in CONFIG$huff_r_median)\n")
  cat("  (Ref for r methodology: Victoria Nino 2026, DOI: 10.31224/7062)\n")
}

mcm_disaggregate <- function(Ptotal, n_steps, alpha_b, beta_b, P0 = 0,
                             max_frac = NULL, seed = NULL) {
  # Disaggregates Ptotal into n_steps intervals using binary multiplicative cascade
  # with Huff Q2 temporal ordering (convective rainfall).
  #
  # TEMPORAL STRUCTURE:
  # The pure MCM generates magnitudes without temporal structure. For
  # tropical convective rainfall (Andean Colombia) the peak occurs in the
  # central third of the event. The Huff (1967) Q2 method is applied:
  # cascade-generated blocks are redistributed temporally
  # following the Huff Q2 mass curve, which concentrates the peak between
  # 30% and 60% of the total duration.
  # Ref: Huff, F.A. (1967). Time distribution of rainfall in heavy storms.
  #      Water Resources Research, 3(4), 1007-1019.
  
  if (!is.null(seed)) set.seed(seed)
  
  n_lev <- ceiling(log2(n_steps))
  
  cascade <- Ptotal
  
  for (lev in 1:n_lev) {
    n_parent <- length(cascade)
    children   <- numeric(n_parent * 2)
    for (k in 1:n_parent) {
      if (cascade[k] < 1e-10) {
        children[2*k-1] <- 0; children[2*k] <- 0
      } else if (runif(1) < P0) {
        if (runif(1) < 0.5) { children[2*k-1] <- cascade[k]; children[2*k] <- 0
        } else               { children[2*k-1] <- 0;          children[2*k] <- cascade[k] }
      } else {
        W <- rbeta(1, alpha_b, beta_b)
        # Random assignment of which child receives W (eliminates directional bias)
        if (runif(1) < 0.5) { children[2*k-1] <- W*cascade[k]; children[2*k] <- (1-W)*cascade[k]
        } else               { children[2*k-1] <- (1-W)*cascade[k]; children[2*k] <- W*cascade[k] }
      }
    }
    cascade <- children
  }
  
  magnitudes <- cascade[1:n_steps]
  
  # ── Temporal ordering — Huff curve per CONFIG$huff_quartile ─────────────────────
  # Huff (1967) curves describe the temporal rainfall distribution
  # according to the quartile in which the intensity peak occurs.
  # Values calibrated for Colombia from Huff (1967) and IDEAM (2014).
  t_huff <- c(0.00, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45,
              0.50, 0.55, 0.60, 0.65, 0.70, 0.75, 0.80, 0.85, 0.90, 0.95, 1.00)
  
  huff_tables <- list(
    # Q1: peak in first quarter — frontal storms / Caribbean coast
    "1" = c(0.00, 0.10, 0.22, 0.35, 0.48, 0.60, 0.70, 0.78, 0.84, 0.89,
            0.92, 0.95, 0.96, 0.97, 0.98, 0.985, 0.99, 0.993, 0.996, 0.998, 1.00),
    # Q2: peak in second quarter — convective Andean region (RECOMMENDED)
    "2" = c(0.00, 0.03, 0.07, 0.12, 0.18, 0.25, 0.34, 0.44, 0.55, 0.65,
            0.74, 0.81, 0.87, 0.91, 0.94, 0.96, 0.97, 0.98, 0.99, 0.995, 1.00),
    # Q3: peak in third quarter — mixed / Pacific region
    "3" = c(0.00, 0.02, 0.04, 0.07, 0.11, 0.15, 0.20, 0.26, 0.33, 0.41,
            0.51, 0.62, 0.72, 0.81, 0.88, 0.93, 0.96, 0.98, 0.99, 0.995, 1.00),
    # Q4: peak in fourth quarter — stratiform / Orinoquia transition
    "4" = c(0.00, 0.02, 0.03, 0.05, 0.08, 0.11, 0.14, 0.18, 0.23, 0.28,
            0.35, 0.43, 0.52, 0.62, 0.72, 0.82, 0.89, 0.94, 0.97, 0.99, 1.00)
  )
  
  huff_q_key <- as.character(max(1L, min(4L, as.integer(round(CONFIG$huff_quartile)))))
  p_huff <- huff_tables[[huff_q_key]]
  
  # Incremental Huff mass fraction for each time step
  t_norm   <- seq(1/n_steps, 1, length.out = n_steps)
  p_acum   <- approx(t_huff, p_huff, xout = t_norm)$y
  p_incr   <- diff(c(0, p_acum))           # Huff curve increments
  p_incr   <- pmax(p_incr, 0)
  p_incr   <- p_incr / sum(p_incr)         # normalise to 1
  
  # Sort magnitudes descending and assign to time positions
  # of highest Huff increment (peak position = maximum Huff increment)
  ord_huff <- order(p_incr, decreasing = TRUE)   # time positions sorted by weight
  ord_mag  <- order(magnitudes, decreasing = TRUE) # blocks sorted descending
  
  p_steps <- numeric(n_steps)
  p_steps[ord_huff] <- magnitudes[ord_mag]
  
  # Add stochastic perturbation: swap ~20% of adjacent positions
  # so the hyetograph is not deterministic for the same set of magnitudes
  n_swap <- max(1, round(n_steps * 0.20))
  for (i in seq_len(n_swap)) {
    pos <- sample(seq_len(n_steps - 1), 1)
    p_steps[c(pos, pos+1)] <- p_steps[c(pos+1, pos)]
  }
  
  # Exact mass correction
  if (sum(p_steps) > 1e-10) p_steps <- p_steps * (Ptotal / sum(p_steps))
  
  # Interval cap
  if (!is.null(max_frac) && Ptotal > 0) {
    cap          <- max_frac * Ptotal
    excess       <- pmax(p_steps - cap, 0)
    p_steps      <- pmin(p_steps, cap)
    excess_total <- sum(excess)
    below_cap    <- p_steps < cap
    if (any(below_cap) && excess_total > 0) {
      weight <- p_steps[below_cap] / sum(p_steps[below_cap])
      p_steps[below_cap] <- p_steps[below_cap] + excess_total * weight
    }
    if (sum(p_steps) > 1e-10) p_steps <- p_steps * (Ptotal / sum(p_steps))
  }
  
  p_steps
}

# ── 10.3  Generation of MCM hyetographs for each return period ─────────────
cat("\n  [10.3] Generating MCM hyetographs for all return periods...\n")

n_sim       <- CONFIG$n_sim
t_step_h  <- dt_min / 60
time_h    <- seq(t_step_h, Td, by = t_step_h)

hyet_results <- list()

for (T_r in RETURN_PERIODS) {
  Ptotal_T <- df_tr$Pmax_mm[df_tr$T_return == T_r]
  
  mat_sim <- matrix(0, nrow = n_sim, ncol = N_steps)
  for (s in 1:n_sim) {
    # Reproducible seed: varies by simulation and by return period T,
    # but is deterministic if CONFIG$mcm_seed is defined.
    seed_s <- if (!is.null(CONFIG$mcm_seed)) as.integer(CONFIG$mcm_seed) + s + as.integer(round(T_r * 100)) else NULL
    mat_sim[s, ] <- mcm_disaggregate(
      Ptotal  = Ptotal_T,
      n_steps = N_steps,
      alpha_b = mcm_params$alpha,
      beta_b  = mcm_params$beta,
      P0      = mcm_params$P0 * CONFIG$mcm_p0_factor,   # configurable factor from CONFIG
      max_frac= CONFIG$mcm_max_frac,                    # configurable interval cap from CONFIG
      seed    = seed_s
    )
  }
  
  # Ensemble quantiles: P50 (median), P90, P10
  hiet_p50  <- apply(mat_sim, 2, quantile, 0.50)
  hiet_p90  <- apply(mat_sim, 2, quantile, 0.90)
  hiet_p10  <- apply(mat_sim, 2, quantile, 0.10)
  hiet_mean <- colMeans(mat_sim)
  
  # Mass correction over quantiles:
  # Column-wise quantile does not conserve the sum exactly (known property
  # of the marginal quantile operator). Rescaled to guarantee that the
  # design hyetograph sums exactly to Ptotal, preserving the temporal shape.
  if (sum(hiet_p50) > 1e-10) hiet_p50 <- hiet_p50 * (Ptotal_T / sum(hiet_p50))
  if (sum(hiet_p90) > 1e-10) hiet_p90 <- hiet_p90 * (Ptotal_T / sum(hiet_p90))
  if (sum(hiet_p10) > 1e-10) hiet_p10 <- hiet_p10 * (Ptotal_T / sum(hiet_p10))
  
  hyet_results[[paste0("T", T_r)]] <- list(
    T_return = T_r,
    Ptotal    = Ptotal_T,
    time_h  = time_h,
    p50       = hiet_p50,
    p90       = hiet_p90,
    p10       = hiet_p10,
    mean_ens  = hiet_mean,
    ensemble  = mat_sim
  )
  
  cat(sprintf("    T = %6.2f yr  |  Ptotal = %6.2f mm  → OK\n", T_r, Ptotal_T))
}
# ==============================================================================
# MODULE 11: HYETOGRAPH PLOTS
# ==============================================================================
cat_sep("11. Hyetogram plots")

for (T_r in RETURN_PERIODS) {
  rh <- hyet_results[[paste0("T", T_r)]]
  
  df_hiet <- data.frame(
    t   = rh$time_h,
    p50 = rh$p50,
    p90 = rh$p90,
    p10 = rh$p10,
    med = rh$mean_ens
  ) %>%
    mutate(
      P_cum_50  = cumsum(p50),
      P_cum_90  = cumsum(p90),
      P_cum_10  = cumsum(p10)   # lower ensemble band
    )
  
  # ── Main hyetograph ──────────────────────────────────────────────────
  p_hiet <- ggplot(df_hiet, aes(x = t)) +
    geom_rect(aes(xmin = t - t_step_h, xmax = t,
                  ymin = 0, ymax = p90),
              fill = "#AED6F1", colour = NA, alpha = 0.6) +
    geom_rect(aes(xmin = t - t_step_h, xmax = t,
                  ymin = 0, ymax = p50),
              fill = COL$celeste, colour = NA, alpha = 0.9) +
    geom_rect(aes(xmin = t - t_step_h, xmax = t,
                  ymin = 0, ymax = p10),
              fill = COL$azul, colour = NA, alpha = 0.8) +
    geom_line(aes(y = p50), colour = COL$azul, linewidth = 0.7) +
    scale_x_continuous(breaks = seq(0, Td, by = max(1, Td/8))) +
    labs(
      title = sprintf("MCM Design Hyetograph  –  T = %g years", T_r),
      subtitle = sprintf(
        "Duration = %g h  |  Ptotal = %.2f mm  |  Interval = %g min\nDark blue = P10 | Mid blue = P50 (median) | Light blue = P90",
        Td, rh$Ptotal, dt_min),
      x = "Time (hours)",
      y = sprintf("Precipitation per interval (mm / %g min)", dt_min)
    ) +
    theme_base +
    theme(plot.subtitle = element_text(size = 8))
  
  # ── Accumulated mass curve ───────────────────────────────────────────────
  # P10–P90 band of accumulated mass using real ensemble quantiles.
  p_mass <- ggplot(df_hiet, aes(x = t)) +
    geom_ribbon(aes(ymin = P_cum_10, ymax = P_cum_90),
                fill = "#AED6F1", alpha = 0.5) +
    geom_line(aes(y = P_cum_50), colour = COL$azul, linewidth = 1) +
    geom_hline(yintercept = rh$Ptotal, linetype = "dashed",
               colour = COL$rojo, linewidth = 0.7) +
    annotate("text", x = Td * 0.05, y = rh$Ptotal * 1.02,
             label = paste0("Ptotal = ", round(rh$Ptotal, 1), " mm"),
             colour = COL$rojo, hjust = 0, size = 3.2) +
    scale_x_continuous(breaks = seq(0, Td, by = max(1, Td/8))) +
    labs(
      title = sprintf("Accumulated Mass Curve  –  T = %g years", T_r),
      subtitle = "Blue band = P10–P90 ensemble range",
      x = "Time (hours)",
      y = "Accumulated precipitation (mm)"
    ) +
    theme_base
  
  # ── Combine and save ────────────────────────────────────────────────────
  p_comb <- p_hiet / p_mass + plot_annotation(
    title = sprintf("MCM Hyetograph  |  T = %g yr  |  Duration = %g h",
                    T_r, Td),
    theme = theme(
      plot.title = element_text(face = "bold", size = 14, colour = COL$azul,
                                hjust = 0.5)
    )
  )
  
  save_plot(p_comb, paste0("10_hyetograph_", tr_label(T_r)), width = 10, height = 9)
  cat(sprintf("  ✔ Hyetograph T=%g yr saved\n", T_r))
}

# ==============================================================================
# MODULE 12: COMPARATIVE HYETOGRAPH PANEL (ALL RETURN PERIODS)
# ==============================================================================
cat_sep("12. Comparative panel — all return periods")

n_periods <- length(RETURN_PERIODS)
plot_list  <- list()

for (T_r in RETURN_PERIODS) {
  rh <- hyet_results[[paste0("T", T_r)]]
  df_h <- data.frame(t = rh$time_h, p50 = rh$p50, p90 = rh$p90,
                     p10 = rh$p10)
  plot_list[[length(plot_list) + 1]] <-
    ggplot(df_h, aes(x = t)) +
    geom_rect(aes(xmin = t - t_step_h, xmax = t, ymin = 0, ymax = p90),
              fill = "#AED6F1", alpha = 0.6) +
    geom_rect(aes(xmin = t - t_step_h, xmax = t, ymin = 0, ymax = p50),
              fill = COL$celeste, alpha = 0.9) +
    labs(title = paste0("T = ", T_r, " yr  |  ", round(rh$Ptotal, 1), " mm"),
         x = "h", y = "mm") +
    theme_base +
    theme(plot.title = element_text(size = 9, face = "bold"),
          axis.title = element_text(size = 8),
          axis.text  = element_text(size = 7),
          plot.margin = margin(4, 4, 4, 4))
}

p_who <- wrap_plots(plot_list, ncol = 4) +
  plot_annotation(
    title = sprintf("Comparative panel — MCM hyetographs for all return periods\nDuration = %g h | Interval = %g min", Td, dt_min),
    theme = theme(
      plot.title = element_text(face = "bold", size = 13, colour = COL$azul,
                                hjust = 0.5)
    )
  )

save_plot(p_who, "11_comparative_panel", width = 16, height = 12)
cat("  ✔ Comparative panel saved\n")

# ==============================================================================
# MODULE 13: RESULTS TABLES AND EXCEL EXPORT
# ==============================================================================
cat_sep("13. Results export")

# Helper: clean label for column/sheet names (T=2.33 → "T2_33")

# ── Table 1: Statistical test summary ───────────────────────────────────────
t1 <- dplyr::select(summary_table,
                    Test, Statistic, Statistic_name, p_value, Conclusion)

# ── Table 2: Design precipitation ───────────────────────────────────────────
t2 <- df_tr

# ── Table 3: P50 hyetographs per return period ───────────────────────────────
t3_list <- lapply(RETURN_PERIODS, function(T_r) {
  rh  <- hyet_results[[paste0("T", T_r)]]
  lbl <- tr_label(T_r)
  df  <- data.frame(
    Time_h    = round(rh$time_h, 4),
    Time_min  = round(rh$time_h * 60, 1)
  )
  df[[paste0("P_", lbl, "_P10_mm")]] <- round(rh$p10, 3)
  df[[paste0("P_", lbl, "_P50_mm")]] <- round(rh$p50, 3)
  df[[paste0("P_", lbl, "_P90_mm")]] <- round(rh$p90, 3)
  df
})
t3 <- t3_list[[1]]
for (i in 2:length(t3_list)) t3 <- cbind(t3, t3_list[[i]][, -(1:2)])

# ── Table 4: Cumulative P50 hyetographs ─────────────────────────────────────
t4_list <- lapply(RETURN_PERIODS, function(T_r) {
  rh  <- hyet_results[[paste0("T", T_r)]]
  lbl <- tr_label(T_r)
  data.frame(
    Time_h   = round(rh$time_h, 4),
    col_     = round(cumsum(rh$p50), 3)
  ) %>%
    dplyr::rename_with(~ paste0("Pcum_", lbl, "_mm"), "col_")
})
t4 <- t4_list[[1]]
for (i in 2:length(t4_list)) {
  lbl <- tr_label(RETURN_PERIODS[i])
  t4[[paste0("Pcum_", lbl, "_mm")]] <- t4_list[[i]][[2]]
}

# ── Table 5 (HEC-HMS): one sheet per return period ────────────────────────
# Required HEC-HMS format for precipitation hyetographs:
#   • Column 1 : Time  (hh:mm, starting at 00:00)
#   • Column 2 : Incremental Precipitation (mm) — MCM P50 median
#   • Column 3 : Incremental Precipitation P90 (mm) — conservative band
#   • Column 4 : Cumulative Precipitation (mm) — P50 accumulated
# Sheet name follows HEC-HMS convention: "Tr_XXX" (max 31 chars)

hec_sheets <- list()
for (T_r in RETURN_PERIODS) {
  rh  <- hyet_results[[paste0("T", T_r)]]
  lbl <- tr_label(T_r)
  
  # Time vector in HH:MM format (from 00:00 + dt)
  time_min_vec <- round(rh$time_h * 60)
  time_str <- sprintf("%02d:%02d",
                      time_min_vec %/% 60,
                      time_min_vec %%  60)
  
  hec_df <- data.frame(
    Time_hhmm              = time_str,
    Precip_incremental_P50_mm = round(rh$p50, 3),
    Precip_incremental_P90_mm = round(rh$p90, 3),
    Precip_accumulated_P50_mm = round(cumsum(rh$p50), 3),
    stringsAsFactors = FALSE
  )
  
  # Sheet name: max 31 characters, Excel compatible
  sheet_name <- paste0("HEC_Tr", gsub("\\.", "-", sprintf("%g", T_r)), "yr")
  hec_sheets[[sheet_name]] <- hec_df
}

# ── Save main Excel ───────────────────────────────────────────────────
excel_path <- file.path(OUTPUT_DIR, "tables", "MCM_results.xlsx")
write_xlsx(
  c(list(
    "Test_Summary"        = t1,
    "Design_Precipitation"     = t2,
    "Hyetographs_Incremental" = t3,
    "Hyetographs_Cumulative" = t4
  ), hec_sheets),
  path = excel_path
)
cat(sprintf("  ✔ Results exported to: %s\n", excel_path))
cat(sprintf("  ✔ HEC-HMS sheets generated: %d (one per return period)\n",
            length(hec_sheets)))

# ── Save independent HEC-HMS Excel ──────────────────────────────────────
# Clean file with only hyetograph sheets for direct import
# into HEC-HMS: Meteorologic Model → Specified Hyetograph → Import from Excel
hec_path <- file.path(OUTPUT_DIR, "tables", "HEC_HMS_hyetographs.xlsx")
write_xlsx(hec_sheets, path = hec_path)
cat(sprintf("  ✔ HEC-HMS file saved: %s\n", hec_path))
cat("\n  HEC-HMS file structure:\n")
cat("    Each sheet = one return period\n")
cat("    Col A: Time (HH:MM)\n")
cat("    Col B: Incremental precipitation P50 (mm)  ← use in HEC-HMS\n")
cat("    Col C: Incremental precipitation P90 (mm)  ← conservative scenario\n")
cat("    Col D: Cumulative precipitation P50 (mm)   ← verification\n")

# ==============================================================================
# REPORT AUXILIARY FUNCTIONS (shared by Modules 14 and 14b)
# ==============================================================================
# Extracted to avoid duplication between the two report blocks.

write_report_header <- function(title, subtitle) {
  cat("╔══════════════════════════════════════════════════════════════════════════╗\n")
  cat(sprintf("║   %-72s║\n", substr(title,   1, 72)))
  cat(sprintf("║   %-72s║\n", substr(subtitle, 1, 72)))
  cat("╚══════════════════════════════════════════════════════════════════════════╝\n\n")
  cat("REPORT GENERATED:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")
}

report_series_stats <- function() {
  cat(sprintf("  Record period         : %d – %d\n", min(pmax_data$YEAR), max(pmax_data$YEAR)))
  cat(sprintf("  Number of years       : %d\n", nrow(pmax_data)))
  cat(sprintf("  Mean annual Pmax 24h  : %.2f mm\n", mean(pmax_data$PMAX_MAX)))
  cat(sprintf("  Absolute max Pmax 24h : %.2f mm (year %d)\n",
              max(pmax_data$PMAX_MAX), pmax_data$YEAR[which.max(pmax_data$PMAX_MAX)]))
  cat(sprintf("  Absolute min Pmax 24h : %.2f mm (year %d)\n",
              min(pmax_data$PMAX_MAX), pmax_data$YEAR[which.min(pmax_data$PMAX_MAX)]))
  cat(sprintf("  Coeff. of variation   : %.3f\n", sd(pmax_data$PMAX_MAX)/mean(pmax_data$PMAX_MAX)))
  cat(sprintf("  Skewness              : %.3f\n",
              mean(((pmax_data$PMAX_MAX-mean(pmax_data$PMAX_MAX))/sd(pmax_data$PMAX_MAX))^3)))
  cat(sprintf("  Kurtosis              : %.3f\n",
              mean(((pmax_data$PMAX_MAX-mean(pmax_data$PMAX_MAX))/sd(pmax_data$PMAX_MAX))^4)))
}

report_gev_params <- function() {
  p <- gev_fit_base$param
  cat("  The GEV (Generalised Extreme Value) distribution has CDF:\n")
  cat("  F(x) = exp{ -[1 + ξ·(x-μ)/σ]^(-1/ξ) }   with 1 + ξ·(x-μ)/σ > 0\n\n")
  cat(sprintf("  Location parameter  μ = %.4f mm\n", p["loc"]))
  cat(sprintf("  Scale parameter     σ = %.4f mm\n", p["scale"]))
  cat(sprintf("  Shape parameter     ξ = %.4f\n\n",  p["shape"]))
}

# ==============================================================================
# MODULE 14: DETAILED STATIONARITY REPORT
# ==============================================================================
cat_sep("14. Stationarity technical report")

report_path <- file.path(OUTPUT_DIR, "REPORT_MCM_Stationarity.txt")
sink(report_path)
on.exit(sink(), add = TRUE)   # guarantees sink closure even on error

write_report_header(
  "TECHNICAL REPORT: MCM DESIGN HYETOGRAPHS",
  "Non-Stationarity Analysis and ENSO Influence — MCM Model"
)

cat("═══════════════════════════════════════════════════════════════════════════\n")
cat("1. SERIES INFORMATION\n")
cat("═══════════════════════════════════════════════════════════════════════════\n")
report_series_stats()

cat("\n═══════════════════════════════════════════════════════════════════════════\n")
cat("2. STATIONARITY ANALYSIS — DETAILED RESULTS\n")
cat("═══════════════════════════════════════════════════════════════════════════\n\n")

cat("┌─────────────────────────────┬──────────────┬──────────┬────────────────────────────────────────────────┐\n")
cat("│ Test                        │ Statistic    │ p-value  │ Conclusion                                     │\n")
cat("├─────────────────────────────┼──────────────┼──────────┼────────────────────────────────────────────────┤\n")
for (i in 1:nrow(summary_table)) {
  cat(sprintf("│ %-27s │ %-12.4f │ %-8s │ %-46s │\n",
              substr(summary_table$Test[i], 1, 27),
              as.numeric(summary_table$Statistic[i]),
              ifelse(is.na(summary_table$p_value[i]), "   N/A  ",
                     sprintf("%.4f", as.numeric(summary_table$p_value[i]))),
              substr(summary_table$Conclusion[i], 1, 46)))
}
cat("└─────────────────────────────┴──────────────┴──────────┴────────────────────────────────────────────────┘\n\n")

cat("TEST INTERPRETATION:\n\n")
cat("  • Mann-Kendall:\n")
cat("    Non-parametric test for monotonic trend based on the Kendall tau\n")
cat("    Kendall tau. H0: no trend. H0 is rejected if p < 0.05.\n")
cat(sprintf("    Result: τ = %.4f, p = %.4f. %s\n\n",
            mk_tau, mk_pval, test_results$MK$Conclusion))

cat("  • Pettitt:\n")
cat("    Non-parametric test for abrupt change (change point) in the mean.\n")
cat("    Detects whether a statistically significant break year exists.\n")
cat(sprintf("    Result: U = %.0f, p = %.4f. %s\n\n",
            as.numeric(pett_U), pett_pval, test_results$Pettitt$Conclusion))

cat("  • Sequential Sneyers:\n")
cat("    Analyses progressive u(t) and retrogressive u'(t) statistics.\n")
cat("    Crossings within the ±1.96 band suggest a regime change.\n")
cat(sprintf("    Result: %s\n\n", test_results$Sneyers$Conclusion))

cat("  • Moving-Window MK:\n")
cat("    Applies Mann-Kendall in sliding windows to detect local trends.\n")
cat(sprintf("    Window used: %d years.\n", window_size))
cat(sprintf("    Result: %s\n\n", test_results$MWMK$Conclusion))

cat("  • White (heteroscedasticity):\n")
cat("    Tests whether residual variance is constant. Heteroscedastic\n")
cat("    variance signals non-stationarity in dispersion.\n")
cat(sprintf("    Result: BP = %.4f, p = %.4f. %s\n\n",
            white_chi2, white_pval, test_results$White$Conclusion))

cat("  • Augmented Dickey-Fuller (ADF):\n")
cat("    Unit root test. H0: unit root exists (non-stationary).\n")
cat("    With p<0.05 H0 is rejected → stationary series.\n")
cat(sprintf("    Result: DF = %.4f (lags=%d), p = %.4f. %s\n\n",
            adf_stat, adf_lag, adf_pval, test_results$ADF$Conclusion))

cat("FINAL DECISION:\n")
cat(sprintf("  Non-stationarity votes: %d/%d\n", non_stat_votes, total_votes))
if (IS_NON_STATIONARY) {
  cat("  ★ NON-STATIONARY SERIES — Non-stationary MCM model adopted\n")
  cat("    with ONI covariate (Oceanic Niño Index).\n")
  cat("\n  RECOMMENDATIONS:\n")
  cat("  1. Prioritise non-stationary models for frequency analysis,\n")
  cat("     especially when climate change projections are available.\n")
  cat("  2. Update the model periodically as new historical data and\n")
  cat("     ONI projections become available.\n")
  cat("  3. Consider additional covariates such as sea surface temperature\n")
  cat("     (SST) or the AMO index to complement the analysis.\n")
  cat("  4. Update ONI with projections from global climate models (CMIP6)\n")
  cat("     to generate hyetographs under climate change scenarios.\n")
} else {
  cat("  ★ STATIONARY SERIES — Stationary MCM model adopted.\n")
  cat("\n  RECOMMENDATIONS:\n")
  cat("  1. Although the test indicates stationarity, it is recommended to\n")
  cat("     revisit the analysis periodically given the climate change context.\n")
  cat("  2. Keep the series updated to detect future changes.\n")
  cat("  3. Consider parameter uncertainty analysis (bootstrapping)\n")
  cat("     to quantify the robustness of design precipitation estimates.\n")
}

cat("\n═══════════════════════════════════════════════════════════════════════════\n")
cat("3. DESIGN PRECIPITATION\n")
cat("═══════════════════════════════════════════════════════════════════════════\n\n")

cat("┌─────────────────────┬──────────────────┐\n")
cat("│ T return (yr)       │ Pmax 24h (mm)    │\n")
cat("├─────────────────────┼──────────────────┤\n")
for (i in 1:nrow(df_tr)) {
  cat(sprintf("│ %-19g │ %-16.2f │\n", df_tr$T_return[i], df_tr$Pmax_mm[i]))
}
cat("└─────────────────────┴──────────────────┘\n\n")

cat("═══════════════════════════════════════════════════════════════════════════\n")
cat("4. MCM MODEL PARAMETERS\n")
cat("═══════════════════════════════════════════════════════════════════════════\n\n")
cat("  Foundation: Microcanonical Multiplicative Cascade model\n")
cat("  Multipliers W follow a Beta distribution with point mass\n")
cat("  mass at extremes 0 and 1 (dry/wet intervals).\n\n")
cat(sprintf("  α (shape1) = %.4f\n", mcm_params$alpha))
cat(sprintf("  β (shape2) = %.4f\n", mcm_params$beta))
cat(sprintf("  P0 (dry fraction) = %.4f\n", mcm_params$P0))
cat(sprintf("  Storm duration      = %g hours\n", Td))
cat(sprintf("  Time step           = %g minutes\n", dt_min))
cat(sprintf("  Number of steps     = %d\n", N_steps))
cat(sprintf("  Monte Carlo simulations = %d\n", n_sim))

cat("\n═══════════════════════════════════════════════════════════════════════════\n")
cat("5. HYETOGRAPH SUMMARY (P50 — MCM ENSEMBLE MEDIAN)\n")
cat("═══════════════════════════════════════════════════════════════════════════\n\n")

# Header
header_tr <- paste(sprintf("T=%gyr", RETURN_PERIODS), collapse = "  ")
cat(sprintf("  Time (h)     %s\n", header_tr))
cat(sprintf("  %s\n", strrep("-", 13 + 9 * length(RETURN_PERIODS))))

# Data
for (step in 1:N_steps) {
  vals <- sapply(RETURN_PERIODS, function(T_r) {
    hyet_results[[paste0("T", T_r)]]$p50[step]
  })
  cat(sprintf("  %8.3f h   %s\n",
              time_h[step],
              paste(sprintf("%6.2f", vals), collapse = "   ")))
}

cat("\n═══════════════════════════════════════════════════════════════════════════\n")
cat("6. METHODOLOGICAL CONSIDERATIONS AND LIMITATIONS\n")
cat("═══════════════════════════════════════════════════════════════════════════\n\n")
cat(sprintf(paste0(
  "  1. IDF CURVES: Bell (1969) anchored on station GEV quantiles.\n",
  "     P(t,T) = (0.54·t^0.25-0.50)·P(60min,T),  with P(1h)/P(24h)=%.2f.\n",
  "     If local rain gauge data are available, adjust ratio_60_1440 in CONFIG.\n\n"),
  CONFIG$ratio_60_1440))
cat("  2. MCM: Beta parameters were calibrated from monthly Pmax series.\n")
cat("     With hourly or sub-hourly data the estimation would be more\n")
cat("     precise. The model assumes statistical scale self-similarity.\n\n")
cat("  3. GEV DISTRIBUTION: GEV (Generalised Extreme Value) was adopted\n")
cat("     for frequency analysis, internationally recognised for\n")
cat("     hydrological maxima. Verify goodness of fit with KS/AD tests.\n\n")
cat("  4. NON-STATIONARITY: The non-stationary model uses ONI as a linear\n")
cat("     covariate in the GEV location parameter μ. More complex\n")
cat("     relationships (non-linear, interactions) may improve the fit.\n\n")
cat("  5. RETURN PERIODS: In a non-stationary context, the traditional\n")
cat("     return period loses its strict meaning. It is recommended to\n")
cat("     interpret it as 'expected waiting time' or use cumulative\n")
cat("     exceedance probabilities (Cooley 2013; Salas & Obeysekera 2014).\n\n")
cat("  6. UNCERTAINTY: The P10-P90 ensemble band represents the\n")
cat("     variability of the disaggregation process. For risk management\n")
cat("     it is recommended to use the P90 percentile as the conservative design.\n\n")

cat("  7. STATISTICAL POWER OF STATIONARITY TESTS\n")
cat("═══════════════════════════════════════════════════════════════════════════\n\n")
cat(sprintf("  n = %d years. α = %.2f. Approximate power for moderate trend (τ≈0.3):\n\n",
            n, ALPHA))
cat("  ┌───────────────────────────┬──────────────────────────┬──────────────┐\n")
cat("  │ Test                      │ Main strength             │ Power (n≈40) │\n")
cat("  ├───────────────────────────┼──────────────────────────┼──────────────┤\n")
cat("  │ Mann-Kendall              │ Monotonic trend           │ ~0.65–0.80   │\n")
cat("  │ Pettitt                   │ Abrupt change (point)     │ ~0.60–0.75   │\n")
cat("  │ Sequential Sneyers        │ Gradual change (qualit.)  │ Qualitative  │\n")
cat("  │ Moving-Window MK          │ Variable local trend      │ ~0.50–0.70   │\n")
cat("  │ White (heteroscedasticity)│ Variance change           │ ~0.55–0.70   │\n")
cat("  │ ADF (unit root)           │ Strong non-stationarity   │ ~0.40–0.60   │\n")
cat("  └───────────────────────────┴──────────────────────────┴──────────────┘\n\n")
cat("  Ref: Yue & Wang (2004), Hydrological Sciences Journal.\n")
cat("  Note: the majority vote system (3/6) is deliberately conservative;\n")
cat("  a consistent signal across tests of different nature is required.\n\n")

cat("═══════════════════════════════════════════════════════════════════════════\n")
cat("END OF REPORT\n")
cat("═══════════════════════════════════════════════════════════════════════════\n")

sink()
on.exit(NULL)   # reset on.exit after successful first sink closure
cat(sprintf("  ✔ Stationarity report saved: %s\n", report_path))
# ==============================================================================
# MODULE 14B: HYETOGRAPH, IDF CURVE AND METRICS REPORT
# ==============================================================================
cat_sep("14b. Hyetographs, IDF and metrics report")

report_path_hiet <- file.path(OUTPUT_DIR, "REPORT_MCM_Hyetographs.txt")
sink(report_path_hiet)
on.exit(sink(), add = TRUE)   # guarantees sink closes even on error

write_report_header(
  "TECHNICAL REPORT: MCM DESIGN HYETOGRAPHS",
  "IDF Curves · Hyetograph Metrics · Usage Recommendations"
)

# ── 1. Input data and configuration ─────────────────────────────────────────
cat("═══════════════════════════════════════════════════════════════════════════\n")
cat("1. INPUT DATA AND CONFIGURATION\n")
cat("═══════════════════════════════════════════════════════════════════════════\n\n")
report_series_stats()
cat(sprintf("\n  Storm duration      : %g hours\n", Td))
cat(sprintf("  Time step            : %g minutes\n", dt_min))
cat(sprintf("  Number of steps      : %d\n", N_steps))
cat(sprintf("  Monte Carlo sims    : %d per return period\n", n_sim))
cat(sprintf("  Frequency model     : %s\n\n",
            ifelse(IS_NON_STATIONARY,
                   "Non-stationary GEV (GAMLSS, ONI covariate) — fallback stationary GEV",
                   "Stationary GEV (MLE, evd package)")))

# ── 2. GEV parameters ──────────────────────────────────────────────────────────
cat("═══════════════════════════════════════════════════════════════════════════\n")
cat("2. FITTED GEV DISTRIBUTION PARAMETERS\n")
cat("═══════════════════════════════════════════════════════════════════════════\n\n")
report_gev_params()   # reuses gev_fit_base; does not re-fit the GEV

# Local variables for later use in this block (KS, tail type, bootstrap)
loc_r <- gev_fit_base$param["loc"]
sc_r  <- gev_fit_base$param["scale"]
sh_r  <- gev_fit_base$param["shape"]

gev_type <- ifelse(abs(sh_r) < 0.05, "Gumbel (ξ ≈ 0)",
                   ifelse(sh_r > 0, "Frechet (ξ > 0) — heavy right tail",
                          "Weibull (ξ < 0) — bounded upper tail"))
cat(sprintf("  Tail type: %s\n", gev_type))

# Bootstrap confidence interval for ξ (200 replicates)
set.seed(42)
sh_boot <- replicate(200, {
  s <- sample(series, length(series), replace = TRUE)
  tryCatch(fgev(s)$param["shape"], error = function(e) NA_real_)
})
sh_boot <- sh_boot[!is.na(sh_boot)]
ci_sh   <- quantile(sh_boot, c(0.025, 0.975))
cat(sprintf("  95%% bootstrap CI for ξ  : [%.4f, %.4f]\n\n", ci_sh[1], ci_sh[2]))

# KS goodness-of-fit
ks_res  <- ks.test(series, "pgev", loc = loc_r, scale = sc_r, shape = sh_r)
cat(sprintf("  Goodness of fit (KS)   : D = %.4f, p = %.4f  →  %s\n\n",
            ks_res$statistic, ks_res$p.value,
            ifelse(ks_res$p.value > 0.05,
                   "Adequate fit (H0 not rejected at 5%)",
                   "Questionable fit (H0 rejected at 5%)")))

# ── 3. Design precipitation ──────────────────────────────────────────────
cat("═══════════════════════════════════════════════════════════════════════════\n")
cat("3. DESIGN PRECIPITATION Pmax 24h (mm)\n")
cat("═══════════════════════════════════════════════════════════════════════════\n\n")
cat("  ┌──────────────┬──────────────┬──────────────┬──────────────┐\n")
cat("  │ T ret (yr)   │ Pmax 24h(mm) │ Exc. prob.(%)│ 50-yr risk(%)│\n")
cat("  ├──────────────┼──────────────┼──────────────┼──────────────┤\n")
for (i in seq_len(nrow(df_tr))) {
  Tr_i   <- df_tr$T_return[i]
  Pm_i   <- df_tr$Pmax_mm[i]
  prob_i <- round(100 / Tr_i, 2)
  riesg  <- round((1 - (1 - 1/Tr_i)^50) * 100, 1)
  cat(sprintf("  │ %-12g │ %-12.2f │ %-12.2f │ %-12.1f │\n",
              Tr_i, Pm_i, prob_i, riesg))
}
cat("  └──────────────┴──────────────┴──────────────┴──────────────┘\n")
cat("  Note: 50-yr risk = probability of exceedance at least once in 50 years.\n\n")

# ── 4. IDF equation and table ─────────────────────────────────────────
cat("═══════════════════════════════════════════════════════════════════════════\n")
cat("4. INTENSITY-DURATION-FREQUENCY (IDF) CURVES\n")
cat("═══════════════════════════════════════════════════════════════════════════\n\n")

cat("  IDF CURVES — BELL (1969) METHOD ANCHORED ON STATION DATA:\n\n")
cat("    P(t, T) = (0.54·t^0.25 - 0.50) · P(60min, T)   [t in min, P in mm]\n")
cat("    I(t, T) = P(t, T) / (t/60)                      [intensity in mm/h]\n\n")
cat(sprintf("    Anchor: P(60min,T) = P(24h,T) · %.2f   [ratio_60_1440 adopted]\n\n",
            CONFIG$ratio_60_1440))
cat("  REGIONAL PARAMETER ratio_60_1440 = P(1h)/P(24h):\n\n")
cat("  ┌──────────────────────────────────────────┬──────────┬──────────────┐\n")
cat("  │ Region / Zone (Colombia)                 │ Range    │ Typical val. │\n")
cat("  ├──────────────────────────────────────────┼──────────┼──────────────┤\n")
cat("  │ Valle del Cauca, Cauca (Cali, Palmira)   │ 0.38–0.45│ 0.40         │\n")
cat("  │ Antioquia (Medellin, Aburra Valley)      │ 0.36–0.44│ 0.40         │\n")
cat("  │ Cundinamarca (Bogota and savanna)        │ 0.32–0.40│ 0.36         │\n")
cat("  │ Coffee Region (Manizales, Pereira)       │ 0.38–0.46│ 0.42         │\n")
cat("  │ Huila, Tolima (warm valleys)             │ 0.40–0.48│ 0.44         │\n")
cat("  │ Narino (Pasto and highlands)            │ 0.30–0.38│ 0.34         │\n")
cat("  │ Dry Caribbean (Guajira, Valledupar)      │ 0.42–0.55│ 0.48         │\n")
cat("  │ Wet Caribbean (Barranquilla, Cartagena)  │ 0.40–0.50│ 0.45         │\n")
cat("  │ Pacific Region (Choco, Buenaventura)     │ 0.35–0.48│ 0.42         │\n")
cat("  │ Orinoquia / Eastern Plains               │ 0.40–0.52│ 0.46         │\n")
cat("  │ Amazon region                            │ 0.38–0.48│ 0.43         │\n")
cat("  └──────────────────────────────────────────┴──────────┴──────────────┘\n\n")
cat("  Higher values (>0.45): intense convective regime, short storms with\n")
cat("  high peak intensity (Caribbean coast, Eastern Plains, warm valleys).\n")
cat("  Lower values (<0.38): stratiform or high-mountain regime, rainfall\n")
cat("  more uniform in time (Bogota savanna, Narino, paramos).\n")
cat("  With local rain gauge: ratio = mean(Pmax_1h)/mean(Pmax_24h) measured.\n")
cat("  Ref: IDEAM (2014) ENA; Velez et al. (2002); Bell (1969).\n\n")
cat("  P(24h,T) quantiles come from the GEV fit to the historical series.\n")
cat("  The function (0.54·t^0.25-0.50) introduces real curvature in log-log.\n\n")

cat("  IDF TABLE — Intensities (mm/h) for selected durations:\n\n")
durs_tabla <- c(15, 30, 60, 120, 180, 360)
durs_label <- c("15min","30min","1h","2h","3h","6h")
encab <- paste(sprintf("%-8s", durs_label), collapse = " ")
cat(sprintf("  %-12s  %s\n", "T ret (yr)", encab))
cat(sprintf("  %s\n", strrep("-", 12 + 9*length(durs_tabla))))
for (i in seq_len(nrow(df_tr))) {
  Tr_i <- df_tr$T_return[i]
  Pm_i <- df_tr$Pmax_mm[i]
  vals_idf <- sapply(durs_tabla, function(d) {
    round(idf_func(d, Tr_i), 2)   # idf_func(dur_min, Tr) returns I in mm/h
  })
  cat(sprintf("  %-12g  %s\n", Tr_i,
              paste(sprintf("%-8.2f", vals_idf), collapse = " ")))
}
cat("\n")

# ── 5. MCM parameters ───────────────────────────────────────────────────────
cat("═══════════════════════════════════════════════════════════════════════════\n")
cat("5. MCM MODEL PARAMETERS\n")
cat("═══════════════════════════════════════════════════════════════════════════\n\n")
cat(sprintf("  Huff quartile adopted         : Q%d\n", CONFIG$huff_quartile))
cat(sprintf("  Stochastic perturbation       : %.0f%% of adjacent positions swapped\n",
            20))
cat("\n")
cat("  THEORETICAL FOUNDATION:\n")
cat("  The MCM (Schertzer & Lovejoy, 1987) disaggregates total precipitation\n")
cat("  through random multipliers W in a binary cascade:\n\n")
cat("    P_left  = W · P_parent\n")
cat("    P_right = (1-W) · P_parent,    W ~ Beta(α, β)\n\n")
cat("  The mass conservation property guarantees:\n")
cat("    Σ P_intervals = P_total     (exactly, by proportional correction)\n\n")
cat(sprintf("  α (shape1) of Beta = %.4f\n", mcm_params$alpha))
cat(sprintf("  β (shape2) of Beta = %.4f\n", mcm_params$beta))
cat(sprintf("  Multiplier mean E[W] = α/(α+β) = %.4f\n",
            mcm_params$alpha / (mcm_params$alpha + mcm_params$beta)))
cat(sprintf("  Multiplier variance Var[W] = αβ/[(α+β)²(α+β+1)] = %.6f\n",
            (mcm_params$alpha * mcm_params$beta) /
              ((mcm_params$alpha + mcm_params$beta)^2 *
                 (mcm_params$alpha + mcm_params$beta + 1))))
cat(sprintf("  Dry fraction P0 (calibrated)  = %.4f\n", mcm_params$P0))
cat(sprintf("  P0 reduction factor (design) = %.2f  → P0_design = %.4f\n",
            CONFIG$mcm_p0_factor, mcm_params$P0 * CONFIG$mcm_p0_factor))
cat(sprintf("  Rationale: extreme events have fewer dry intervals than the\n"))
cat(sprintf("  calibrated monthly mean (Molnar & Burlando, 2005).\n"))
if (!is.null(CONFIG$mcm_max_frac)) {
  cat(sprintf("  Interval cap (mcm_max_frac) = %.2f × Ptotal\n",
              CONFIG$mcm_max_frac))
} else {
  cat("  Interval cap: deactivated (mcm_max_frac = NULL)\n")
}
cat(sprintf("  Cascade levels = %d  (for %d steps over %g h)\n\n",
            ceiling(log2(N_steps)), N_steps, Td))

# ── 6. Statistical metrics by hyetograph ────────────────────────────────────
cat("═══════════════════════════════════════════════════════════════════════════\n")
cat("6. STATISTICAL METRICS BY HYETOGRAPH (MCM ENSEMBLE, n=500)\n")
cat("═══════════════════════════════════════════════════════════════════════════\n\n")

cat("  Definitions:\n")
cat("  · Peak P50   : maximum intensity of the ensemble median (mm/interval)\n")
cat("  · t_peak     : time at which the median peak occurs (hours)\n")
cat("  · Peak/Ptotal: rainfall fraction in the highest-intensity interval\n")
cat("  · CV_ens     : inter-simulation coefficient of variation (MCM uncertainty)\n")
cat("  · Skewness   : 3rd-order moment of the temporal distribution\n")
cat("  · P10-P90 band (mm) : mean width of the uncertainty band\n\n")

cat("  ┌──────────┬─────────┬────────┬───────────┬────────┬──────────┬──────────┐\n")
cat("  │T ret (yr)│Peak P50 │t_peak  │Peak/Ptotal│CV_ens  │Skewness  │P10-P90   │\n")
cat("  │          │  (mm)   │  (h)   │    (%)    │  (%)   │hyetograph│band (mm) │\n")
cat("  ├──────────┼─────────┼────────┼───────────┼────────┼──────────┼──────────┤\n")

for (T_r in RETURN_PERIODS) {
  rh   <- hyet_results[[paste0("T", T_r)]]
  p50  <- rh$p50
  mat  <- rh$ensemble
  
  pico_val   <- max(p50)
  t_pico     <- rh$time_h[which.max(p50)]   # peak time
  pico_frac  <- round(pico_val / rh$Ptotal * 100, 1)
  
  # CV inter-simulation: average CV per time step
  cv_paso    <- apply(mat, 2, function(x) sd(x)/mean(x) * 100)
  cv_ens     <- round(mean(cv_paso, na.rm = TRUE), 1)
  
  # Temporal skewness of P50 hyetograph
  t_norm <- rh$time_h / Td
  skew   <- round(sum(p50 * t_norm) / sum(p50), 3)  # normalised centre of mass
  
  # Mean width of P10-P90 band
  band   <- round(mean(rh$p90 - rh$p10), 3)
  
  cat(sprintf("  │ %-8g │ %7.3f │ %6.2f │ %9.1f │ %6.1f │ %8.3f │ %8.3f │\n",
              T_r, pico_val, t_pico, pico_frac, cv_ens, skew, band))
}
cat("  └──────────┴─────────┴────────┴───────────┴────────┴──────────┴──────────┘\n\n")

cat("  METRICS INTERPRETATION:\n\n")
cat("  · High CV_ens (>50%) indicates high MCM stochastic variability;\n")
cat("    for conservative design use the P90 percentile.\n")
cat("  · Temporal skewness: value <0.5 indicates rainfall concentrated in the first\n")
cat("    half of the event; >0.5 in the second half. The MCM generates distributions\n")
cat("    symmetric on average; variation between simulations is expected.\n")
cat("  · P10-P90 band: represents the intrinsic uncertainty of the disaggregation\n")
cat("    process. The band grows proportionally with return period.\n\n")

# ── 7. Peak intensity by percentile ─────────────────────────────────────────
cat("═══════════════════════════════════════════════════════════════════════════\n")
cat("7. PEAK INTENSITY ANALYSIS (mm/h) BY PERCENTILE\n")
cat("═══════════════════════════════════════════════════════════════════════════\n\n")
cat(sprintf("  Discretisation interval = %g min  →  conversion factor: ×%.4f\n\n",
            dt_min, 60/dt_min))

cat("  ┌──────────┬────────────────────────────────────────────────────────┐\n")
cat("  │T ret (yr)│   Peak intensity (mm/h) by ensemble percentile        │\n")
cat("  │          │   P10        P25        P50        P75        P90      │\n")
cat("  ├──────────┼────────────────────────────────────────────────────────┤\n")

factor_mmh <- 60 / dt_min  # conversion mm/interval → mm/h
for (T_r in RETURN_PERIODS) {
  rh  <- hyet_results[[paste0("T", T_r)]]
  mat <- rh$ensemble
  picos_sim <- apply(mat, 1, max) * factor_mmh
  q_picos   <- quantile(picos_sim, c(0.10, 0.25, 0.50, 0.75, 0.90))
  cat(sprintf("  │ %-8g │ %9.2f  %9.2f  %9.2f  %9.2f  %9.2f  │\n",
              T_r, q_picos[1], q_picos[2], q_picos[3], q_picos[4], q_picos[5]))
}
cat("  └──────────┴────────────────────────────────────────────────────────┘\n\n")

# ── 8. Temporal concentration (accumulated mass percentiles) ─────────────────
cat("═══════════════════════════════════════════════════════════════════════════\n")
cat("8. TEMPORAL CONCENTRATION — TIME TO REACH FRACTIONS OF PTOTAL\n")
cat("═══════════════════════════════════════════════════════════════════════════\n\n")
cat("  Time (h) at which cumulative rainfall reaches 25%, 50%, 75% and 90%\n")
cat("  of total precipitation (based on median P50 hyetograph):\n\n")
cat("  ┌──────────┬──────────┬──────────┬──────────┬──────────┐\n")
cat("  │T ret (yr)│  t(25%)  │  t(50%)  │  t(75%)  │  t(90%)  │\n")
cat("  ├──────────┼──────────┼──────────┼──────────┼──────────┤\n")
for (T_r in RETURN_PERIODS) {
  rh    <- hyet_results[[paste0("T", T_r)]]
  acum  <- cumsum(rh$p50) / rh$Ptotal
  t25   <- rh$time_h[which(acum >= 0.25)[1]]
  t50   <- rh$time_h[which(acum >= 0.50)[1]]
  t75   <- rh$time_h[which(acum >= 0.75)[1]]
  t90   <- rh$time_h[which(acum >= 0.90)[1]]
  cat(sprintf("  │ %-8g │ %8.3f │ %8.3f │ %8.3f │ %8.3f │\n",
              T_r, t25, t50, t75, t90))
}
cat("  └──────────┴──────────┴──────────┴──────────┴──────────┘\n\n")
cat("  Note: the MCM generates temporally symmetric structures on average.\n")
cat("  For hydraulic modelling it is recommended to compare with patterns\n")
cat("  observed locally (rain gauges) and adjust if necessary.\n\n")

# ── 9. Specific usage recommendations ───────────────────────────────────────
cat("═══════════════════════════════════════════════════════════════════════════\n")
cat("9. SPECIFIC USAGE RECOMMENDATIONS\n")
cat("═══════════════════════════════════════════════════════════════════════════\n\n")

cat("  A. RETURN PERIOD AND DESIGN PERCENTILE SELECTION:\n\n")
cat("  ┌────────────────────────────────────┬──────────┬──────────────────────┐\n")
cat("  │ Type of infrastructure             │ T ret.   │ MCM Percentile       │\n")
cat("  ├────────────────────────────────────┼──────────┼──────────────────────┤\n")
cat("  │ Urban storm drainage                │ 2–10yr   │ P50 (normal design)  │\n")
cat("  │ Secondary drainage channels         │ 10–25yr  │ P50–P75              │\n")
cat("  │ Main channels / collectors          │ 25–100yr │ P75–P90              │\n")
cat("  │ Flood control structures            │ 100–500yr│ P90 (conservative)  │\n")
cat("  │ Dams / spillways (NSF/USACE)        │ PMF/500yr│ P90                  │\n")
cat("  └────────────────────────────────────┴──────────┴──────────────────────┘\n\n")

cat("  B. USE IN HEC-HMS:\n\n")
cat("  1. Import the file HEC_HMS_hyetographs.xlsx.\n")
cat("  2. In Meteorologic Model, select 'Specified Hyetograph'.\n")
cat("  3. Use column B (P50) for normal design; column C (P90) for\n")
cat("     conservative verification.\n")
cat("  4. The time base must match the model computation interval.\n")
cat("  5. Verify that the model total depth matches\n")
cat("     with Pmax_mm from the design precipitation table.\n\n")

cat("  C. UNCERTAINTY AND SENSITIVITY ANALYSIS:\n\n")
cat("  1. Run the HEC-HMS model with P50 and P90 to quantify the range\n")
cat("     of expected peak flows.\n")
cat("  2. The P10-P90 band can be interpreted as a confidence interval\n")
cat("     for the design hyetograph.\n")
cat("  3. With longer series or hourly data, recalibrate the\n")
cat("     MCM parameters α and β for greater local accuracy.\n\n")

cat("  D. HUFF QUARTILE SELECTION — CLIMATE REGIME GUIDE:\n\n")
if (!is.null(CONFIG$huff_r_median)) {
  cat(sprintf("  Median r used           : %.3f (measured from rain gauge)\n",
              CONFIG$huff_r_median))
  cat(sprintf("  Quartile selected       : Q%d (automatic from r)\n\n",
              CONFIG$huff_quartile))
  cat("  Ref: Victoria Niño (2026). Design storm duration from hourly rainfall\n")
  cat("  records in a bimodal Andean climate. EngrXiv. DOI: 10.31224/7062\n\n")
} else {
  cat(sprintf("  Quartile adopted        : Q%d (manual, no r measured)\n\n",
              CONFIG$huff_quartile))
  cat("  For objective quartile determination from observed median r:\n")
  cat("  see: Victoria Nino (2026), DOI: 10.31224/7062\n\n")
}
cat("  ┌────────┬──────────────────┬─────────────┬──────────────────────────────────┐\n")
cat("  │Quartile│ Median r         │ Peak        │ Application (Colombia)           │\n")
cat("  │        │ (rain gauge)     │ position    │                                  │\n")
cat("  ├────────┼──────────────────┼─────────────┼──────────────────────────────────┤\n")
cat("  │  Q1    │ 0.10 – 0.25      │ 0–25% dur.  │ High orography, mountain         │\n")
cat("  │        │                  │ Very early  │ stations >1500 m (La Primavera   │\n")
cat("  │        │                  │             │ r=0.167), Caribbean coast        │\n")
cat("  ├────────┼──────────────────┼─────────────┼──────────────────────────────────┤\n")
cat("  │  Q2 ★  │ 0.25 – 0.45      │ 25–50% dur. │ ANDEAN URBAN: Cali,              │\n")
cat("  │        │                  │ 1st half    │ Medellin, Bogota, Manizales.     │\n")
cat("  │        │                  │ [DEFAULT]   │ Moderate orographic convective   │\n")
cat("  ├────────┼──────────────────┼─────────────┼──────────────────────────────────┤\n")
cat("  │  Q3    │ 0.45 – 0.60      │ 50–75% dur. │ Pacific (Choco, Buenaventura),   │\n")
cat("  │        │                  │ 2nd half    │ mixed convective-stratiform      │\n")
cat("  ├────────┼──────────────────┼─────────────┼──────────────────────────────────┤\n")
cat("  │  Q4    │ 0.60 – 0.80      │ 75–100% dur │ Eastern Plains (Apr–May),        │\n")
cat("  │        │                  │ Very late   │ long-duration stratiform events  │\n")
cat("  └────────┴──────────────────┴─────────────┴──────────────────────────────────┘\n\n")
cat("  r = (t_peak - t_start) / D_physical  [Victoria Nino 2026, DOI:10.31224/7062]\n")
cat("  Design value: median r over all independent IETD-catalogued storm events.\n")
cat("  Without rain gauge: use Q2 for Andean urban areas (conservative default).\n\n")
cat("  Ref: Huff (1967), Water Res. Res. 3(4); Victoria Niño (2026), EngrXiv,\n")
cat("       DOI: 10.31224/7062.\n\n")

cat("  E. LIMITATIONS AND PRECAUTIONS:\n\n")
cat(sprintf(paste0(
  "  1. IDF: Bell (1969) anchored on station GEV quantiles. ratio_60_1440=%.2f\n",
  "     (regional). With local rain gauge, calibrate this value in CONFIG.\n"),
  CONFIG$ratio_60_1440))
cat("  2. MCM: parameters α and β calibrated from DAILY Pmax data, not from\n")
cat("     sub-hourly rain gauges. Valid for preliminary design. For final\n")
cat("     engineering design, recalibrate with local hourly data.\n")
cat("  3. MCM: assumes statistical scale self-similarity. In catchments with\n")
cat("     organised convection (ITCZ) this assumption may not hold.\n")
cat("  4. In a non-stationary context, design quantiles correspond\n")
cat("     to the neutral ONI scenario. For El Nino/La Nina, use df_tr_ns.\n\n")

cat("  F. INPUT DATA VALIDATION:\n\n")
if (length(data_alerts) > 0 || length(data_warnings) > 0) {
  if (length(data_warnings) > 0) {
    cat("  Warnings detected during data loading:\n")
    for (av in data_warnings) cat(sprintf("    · %s\n", av))
    cat("\n")
  }
  if (length(data_alerts) > 0) {
    cat("  Alerts detected during data loading:\n")
    for (al in data_alerts) cat(sprintf("    · %s\n", al))
    cat("\n")
  }
} else {
  cat("  Data validated successfully. No issues detected.\n\n")
}

cat("  G. RECOMMENDATIONS FOR FUTURE IMPROVEMENT:\n\n")
cat("  1. Incorporate sub-hourly rain gauge data for direct calibration\n")
cat("     of MCM parameters directly at minute scale.\n")
cat("  2. Empirically verify the dominant Huff quartile with local rain\n")
cat("     gauges and adjust CONFIG$huff_quartile accordingly.\n")
cat("  3. Evaluate alternative disaggregation models: Bartlett-Lewis,\n")
cat("     Neyman-Scott, or higher-order Markov chains.\n")
cat("  4. Update ONI with CMIP6 projections to generate design hyetographs\n")
cat("     under climate change scenarios (SSP2-4.5, SSP5-8.5).\n")
cat("  5. Validate hyetographs against historical events recorded in\n")
cat("     rain gauges in the catchment of interest.\n\n")

# ── 10. Complete hyetograph table P50 and P90 ───────────────────────────────
cat("═══════════════════════════════════════════════════════════════════════════\n")
cat("10. HYETOGRAPH TABLE — INCREMENTAL PRECIPITATION P50 (mm)\n")
cat("═══════════════════════════════════════════════════════════════════════════\n\n")

header_tr2 <- paste(sprintf("T=%-5gyr", RETURN_PERIODS), collapse = "  ")
cat(sprintf("  %-8s   %s\n", "Time (h)", header_tr2))
cat(sprintf("  %s\n", strrep("-", 10 + 9 * length(RETURN_PERIODS))))
for (step in seq_len(N_steps)) {
  vals_p50 <- sapply(RETURN_PERIODS, function(T_r)
    hyet_results[[paste0("T", T_r)]]$p50[step])
  cat(sprintf("  %8.3f   %s\n", time_h[step],
              paste(sprintf("%7.3f", vals_p50), collapse = "  ")))
}

cat("\n═══════════════════════════════════════════════════════════════════════════\n")
cat("11. HYETOGRAPH TABLE — INCREMENTAL PRECIPITATION P90 (mm)\n")
cat("═══════════════════════════════════════════════════════════════════════════\n\n")

cat(sprintf("  %-8s   %s\n", "Time (h)", header_tr2))
cat(sprintf("  %s\n", strrep("-", 10 + 9 * length(RETURN_PERIODS))))
for (step in seq_len(N_steps)) {
  vals_p90 <- sapply(RETURN_PERIODS, function(T_r)
    hyet_results[[paste0("T", T_r)]]$p90[step])
  cat(sprintf("  %8.3f   %s\n", time_h[step],
              paste(sprintf("%7.3f", vals_p90), collapse = "  ")))
}

cat("\n═══════════════════════════════════════════════════════════════════════════\n")
cat("REFERENCES\n")
cat("═══════════════════════════════════════════════════════════════════════════\n\n")
cat("  Schertzer, D. & Lovejoy, S. (1987). Physical modeling and analysis of\n")
cat("  rain and clouds by anisotropic scaling multiplicative processes.\n")
cat("  J. Geophys. Res., 92(D8), 9693–9714.\n\n")
cat("  Molnar, P. & Burlando, P. (2005). Preservation of rainfall properties\n")
cat("  in stochastic disaggregation by a simple random cascade model.\n")
cat("  Atmospheric Research, 77(1–4), 137–151.\n\n")
cat("  Sneyers, R. (1990). On the Statistical Analysis of Series of Observations.\n")
cat("  WMO Technical Note No. 143. World Meteorological Organization, Geneva.\n\n")
cat("  Katz, R.W. (2013). Statistical methods for nonstationary extremes.\n")
cat("  In Extremes in a Changing Climate. Springer, Dordrecht.\n\n")
cat("  Salas, J.D. & Obeysekera, J. (2014). Revisiting the concepts of return\n")
cat("  period and risk for nonstationary hydrologic extreme events.\n")
cat("  J. Hydrol. Eng., 19(3), 554–568.\n\n")
cat("  IDEAM (2014). Estudio Nacional del Agua. Instituto de Hidrologia,\n")
cat("  Meteorologia y Estudios Ambientales. Bogota, Colombia.\n\n")
cat("  Yue, S. & Wang, C. (2004). The Mann-Kendall Test Modified by Effective\n")
cat("  Sample Size to Detect Trend in Serially Correlated Hydrological Series.\n")
cat("  Water Resources Management, 18(3), 201–218.\n\n")

cat("═══════════════════════════════════════════════════════════════════════════\n")
cat("VALIDATION NOTE\n")
cat("═══════════════════════════════════════════════════════════════════════════\n\n")
cat("  The hyetographs generated are SYNTHETIC, produced by the MCM model\n")
cat("  calibrated from daily precipitation data. Before using in final\n")
cat("  engineering design it is recommended to:\n\n")
cat("  1. Compare with historical events recorded in local rain gauges.\n")
cat("  2. Verify that the temporal shape (peak distribution) is consistent\n")
cat("     with patterns observed in the catchment of interest.\n")
cat("  3. If sub-hourly rain gauge data are available, recalibrate α and β\n")
cat("     of the MCM directly at minute scale (Module 10, function\n")
cat("     estimate_mcm_params). This removes the monthly-to-sub-hourly\n")
cat("     extrapolation assumption.\n\n")

cat("═══════════════════════════════════════════════════════════════════════════\n")
cat("END OF HYETOGRAPH REPORT\n")
cat("═══════════════════════════════════════════════════════════════════════════\n")

sink()
on.exit(NULL)   # reset on.exit after successful sink closure
cat(sprintf("  ✔ Hyetographs report saved: %s\n", report_path_hiet))

# ==============================================================================
# MODULE 15: FINAL EXECUTION SUMMARY
# ==============================================================================
cat_sep("FINAL SUMMARY")
cat(sprintf("\n  Results directory: %s/\n", OUTPUT_DIR))
cat("\n  Generated outputs:\n")
cat("  ├── plots/\n")
cat("  │   ├── 01_series_pmax.png\n")
cat("  │   ├── 02_oni_historical.png\n")
cat("  │   ├── 03_scatter_pmax_oni.png\n")
cat("  │   ├── 04_pettitt.png\n")
cat("  │   ├── 05_sneyers.png\n")
cat("  │   ├── 06_moving_window_mk.png\n")
cat("  │   ├── 07_white_heteroscedasticity.png\n")
cat("  │   ├── 08_frequency_curve.png\n")
cat("  │   ├── 09_IDF_curves.png\n")
for (T_r in RETURN_PERIODS) {
  cat(sprintf("  │   ├── 10_hyetograph_%s.png\n", tr_label(T_r)))
}
cat("  │   └── 11_comparative_panel.png\n")
cat("  ├── tables/\n")
cat("  │   ├── MCM_results.xlsx             (full analysis, 4+N sheets)\n")
cat("  │   └── HEC_HMS_hyetographs.xlsx     (one sheet per return period)\n")
cat("  ├── REPORT_MCM_Stationarity.txt\n")
cat("  └── REPORT_MCM_Hyetographs.txt\n")

cat(sprintf("\n  ✔ MCM script v%s completed successfully.\n", CONFIG$version))
cat(sprintf("  Design duration    : %g h\n", Td))
cat(sprintf("  Time step          : %g min\n", dt_min))
cat(sprintf("  Return periods     : %s years\n",
            paste(RETURN_PERIODS, collapse=", ")))
cat(sprintf("  Model adopted      : %s\n",
            ifelse(IS_NON_STATIONARY, "Non-stationary MCM (ONI)", "Stationary MCM")))
cat(sprintf("  MCM seed           : %s\n",
            ifelse(is.null(CONFIG$mcm_seed), "random (not reproducible)", as.character(CONFIG$mcm_seed))))
cat(sprintf("  P0 design factor   : %.2f × P0_calibrated\n", CONFIG$mcm_p0_factor))
cat(sprintf("  Interval cap       : %s\n",
            ifelse(is.null(CONFIG$mcm_max_frac), "no limit", sprintf("%.2f × Ptotal", CONFIG$mcm_max_frac))))
cat("\n  ⚠ VALIDATION REMINDER:\n")
cat("    Hyetographs are synthetic. Compare with local rain gauge records\n")
cat("    before use in final engineering design (see Hyetograph Report,\n")
cat("    Validation Note, and Section E Recommendations).\n\n")
