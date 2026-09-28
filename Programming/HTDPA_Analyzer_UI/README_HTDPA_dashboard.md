# HTD/PA Release Analysis Dashboard v5

Run with:

```bash
python -m pip install -r requirements_htdpa_dashboard.txt
streamlit run htdpa_release_dashboard.py
```

## Publication-matched model fitting

The dashboard now matches the publication scripts for transformed models:

- Korsmeyer-Peppas is fitted nonlinearly as `F = k*t^n` on the original release-fraction scale.
- KP and Peppas-Sahlin use recorded time since the experimental release origin.
- Residuals and R² are evaluated on the original release-fraction scale.
- An optional separate **KP / Peppas-Sahlin start-and-end window** can exclude
  near-background or out-of-domain observations without changing the common
  Zero/Higuchi/First-order interval.
- For the P23-MB water/PBS probe comparison, set the ordinary analysis window
  to 68–391 min and the KP/Peppas-Sahlin start to **85 min** to reproduce
  publication Figure F.
- The dashboard applies the entered LOQ before fitting and reports `n` as an apparent release exponent.
- Figure labels default directly to the uploaded filename without its extension and remain editable.

Use a common analysis window for Zero, Higuchi, and First-order fits so their
fit quality remains comparable. Use a different KP/Peppas-Sahlin window only
for a predefined model-domain or LOQ reason; both models should share that
same secondary window.

## Selectable analysis time window

Each uploaded condition now has an **Analysis time window** control with two modes:

- **Automatic (LOQ-valid start + terminal R² trimming)**
  - Starts at the first positive-release time point that passes the strict-majority LOQ rule (>50% of replicates at/above LOQ).
  - May remove terminal time points sequentially until the target linear R² is reached or the minimum number of points remains.
  - Internal points are never selectively deleted.

- **Manual start and end time**
  - Choose the exact starting time point and last time point from the uploaded dataset.
  - The selected boundaries are kept fixed; automatic terminal trimming is disabled.
  - Within the selected interval, the positive-release and >50% LOQ validity rules still apply.

The summary table records the window mode, requested start/end times, accepted fit start/end times, number of fit points, R², Kp, model R² values, CV statistics, and KP exponent.

## Core QC rules

1. Negative cumulative-release values are excluded.
2. A time point is quantitatively valid only if mean cumulative release is positive and a strict majority (>50%) of replicates are at/above assay LOQ.
3. CV is calculated using only quantifiable replicates at accepted time points.
4. LOQ is applied directly to measured receptor concentration when available.
5. If measured concentration is unavailable, LOQ is approximated as a release-% threshold and the dashboard issues a warning.
6. Kp and kinetic fits use the selected accepted analysis window.
7. KP and Peppas-Sahlin fits additionally use 0 < release <= 60% and the separately selected KP start time.

## Output

The dashboard generates the six-panel summary figure:

A. Mean cumulative release profiles
B. Median CV with IQR
C. Linear-fit R²
D. Kinetic-model R² heatmap
E. Apparent permeability Kp
F. Apparent Korsmeyer-Peppas exponent n

PNG export is 600 dpi, and the analysis summary can be downloaded as CSV.
