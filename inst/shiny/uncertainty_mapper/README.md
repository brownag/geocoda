# Geocoda Uncertainty Mapper

Interactive web application for geostatistical uncertainty visualization.

## Quick Start

```r
library(geocoda)
run_uncertainty_mapper()
```

## Features

**Workflow**: Load CSV → Configure model → Run simulation → Explore results → Export

**Tabs**:
- **Data Upload**: Load CSV or example data. Auto-detects columns (composition, location, zone).
- **Model Setup**: Select backend (analytical, Stan, Nimble, 3D variants), configure parameters.
- **Simulation**: Run with specified nsim (100-10,000 realizations).
- **Visualization**: 
  - Probability maps: P(variable > threshold) by zone
  - Percentile maps: P10/P50/P90 uncertainty bands
  - Depth profiles: Uncertainty by depth level
  - 3D scatter: Interactive compositional space
  - Risk curves: Exceedance probability vs threshold
- **Download**: Export GeoTIFF, shapefile, JSON, CSV, HTML report

## Input Data Format

CSV with required columns: composition (sand, silt, clay %), spatial coordinates (x, y, optional z), zone assignment.

Data validation checks: compositions sum to 100% ± 1%, no negatives, no missing values, min 2 obs per zone.

## Backend Selection

| Backend | Speed | Use Case |
|---------|-------|----------|
| Analytical | Fastest | Screening |
| Stan | Moderate | Publication (full Bayesian) |
| Nimble | Moderate | Specialized models |
| Analytical 3D | Fast | Distance-weighted spatial |
| Stan 3D | Moderate | Rigorous 3D Bayesian |
| Nimble 3D | Moderate | Custom 3D MCMC |

## Common Parameters

- **Pooling Coefficient** (0-1): Shrinkage strength toward global estimates
- **nsim** (number of realizations): 100 = screening, 1,000 = typical, 5,000+ = publication

## Interpretation

- **Probability maps**: Red zones = high risk, blue zones = low risk
- **Percentile bands**: Narrow = low uncertainty, wide = high uncertainty  
- **Depth profiles**: Shaded region = 80% prediction interval (P10-P90)

## Input Data Example

```csv
zone,x,y,z,sand,silt,clay
A1,100,200,10,42,33,25
A1,105,205,15,45,30,25
B1,150,250,10,35,45,20
```

## Dependencies

Requires: shiny, shinydashboard, plotly, ggplot2, terra, sf, mvtnorm, geocoda

## Citation

If used in research, cite the geocoda package and this application.

---

**Version**: 0.2.0
**Status**: Production-ready


