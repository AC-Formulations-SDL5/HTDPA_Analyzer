# SLF3X peak-to-peak flow versus RPM

Generate the figure with the standalone [Python script](SLF3X_peak_to_peak_vs_RPM.py).

![Peak-to-peak flow versus pump speed](SLF3X_peak_to_peak_vs_RPM.png)

## Installation

Requires Python 3.9 or newer, NumPy, pandas, and Matplotlib:

```bash
python3 -m pip install numpy pandas matplotlib
```

## Input files

Place these CSV files beside the script:

| Pump speed (RPM) | Source file |
| ---: | --- |
| 5 | `HTDPA_5RPM_30oC_1H.csv` |
| 10 | `HTDPA_10RPM_30oC_1H.csv` |
| 15 | `HTDPA_15RPM_30oC_1H.csv` |
| 20 | `HTDPA_20RPM_30oC_1H.csv` |
| 25 | `HTDPA_25rpm_30oC_1H.csv` |
| 30 | `HTDPA_30rpm_30oC_1.5H.csv` |

The script reads numeric flow columns `s1_flow_ml_min` through `s4_flow_ml_min`. At 30 RPM, sensor 2 is excluded because of the sensor issue documented in the combined-figure script; its column is unused. The source data contain one-second averaged readings. No Excel workbooks, backup files, existing images, or other project scripts are needed.

## Usage

From the script folder:

```bash
python3 SLF3X_peak_to_peak_vs_RPM.py
```

The output is `SLF3X_peak_to_peak_vs_RPM.png`, saved beside the script. Running the default command overwrites an existing image with that name. Source CSV files are unchanged.

To specify other locations:

```bash
python3 SLF3X_peak_to_peak_vs_RPM.py --data-dir "/path/to/data" --output "/path/to/figure.png"
```

The output directory must already exist. Show command-line help with:

```bash
python3 SLF3X_peak_to_peak_vs_RPM.py --help
```

## Calculations

### Flow recalibration

Each recorded flow value is corrected as follows:

```text
corrected flow = recorded flow / old embedded factor × gravimetric factor
```

| Sensor | Old embedded factor | Gravimetric factor |
| --- | ---: | ---: |
| S1 | 4.4 | 4.508343 |
| S2 | 4.2 | 4.474666 |
| S3 | 2.0 | 2.320488 |
| S4 | 3.0 | 3.340179 |

These constants assume the calibration basis of the supplied CSV files. Using already recalibrated values would apply the correction again.

### Plateau detection

For each sensor, the deadband is 15% of the 90th percentile of absolute finite flow values. Contiguous positive and negative runs outside that deadband must contain at least three samples.

A positive lobe plateau is the mean of values at or above that lobe's median. A negative lobe plateau is the mean of values at or below its median. Q+ and Q− are the respective means across detected lobe plateaus.

### Peak-to-peak flow and error bars

Each sensor's peak-to-peak flow is **Q+ − Q−**, in mL/min. It is a plateau-based measure, rather than the raw maximum minus minimum.

Each bar is the mean across four sensors at 5–25 RPM, and three sensors (S1, S3, S4) at 30 RPM. Error bars show ±1 sample standard deviation across sensors (`ddof=1`). Labels show the mean and coefficient of variation:

```text
CV (%) = 100 × sample SD / |mean|
```

## Output and validation

The output is a 19,200 × 10,800 pixel PNG rendered at 1200 DPI (16 × 9 inches). The script generates only the peak-to-peak figure.

The original 300 DPI version was checked against the existing figure using the supplied CSV files and reproduced it pixel for pixel. The script now renders the same plot at 1200 DPI. Different Matplotlib versions or fonts may affect rendering.
