
<!-- README.md is generated from README.Rmd. Please edit that file -->

# geocoda

<!-- badges: start -->
<!-- badges: end -->

## Overview

`geocoda` is an R package for geostatistical simulation of compositional
data. It implements a complete workflow for generating spatial
realizations of constrained multivariate compositions using Isometric
Log-Ratio (ILR) transformations and geostatistical kriging (Independent
Univariate or Linear Model of Coregionalization).

While designed to handle any compositional data (geochemistry,
vegetation proportions, etc.), the primary application is spatial
simulation of **soil texture separates** (sand, silt, clay) -
proportions that sum to a constant (typically 100%).

## Workflow

The package follows a five-step process:

1.  **Constrain** - Define validity bounds for each component (low,
    representative, high)
2.  **Transform** - Convert compositions to Isometric Log-Ratio (ILR)
    space to break the sum constraint
3.  **Model** - Fit a multivariate geostatistical model to ILR variables
    (univariate or LMC)
4.  **Simulate** - Generate spatial realizations in ILR space
5.  **Back-transform** - Return to original units, guaranteeing the sum
    constraint

## Installation

You can install the development version of geocoda from
[GitHub](https://github.com/brownag/geocoda) with:

``` r
# install.packages("remotes")
remotes::install_github("brownag/geocoda")
```

## Quick Start

``` r
library(geocoda)

# Define composition bounds (Sand, Silt, Clay percentages)
constraints <- list(
  SAND = list(min = 0, max = 40),
  SILT = list(min = 50, max = 80),
  CLAY = list(min = 10, max = 20)
)

# Expand into a valid composition grid
grid <- gc_expand_bounds(constraints, step = 0.5, target_sum = 100)

# Bootstrap samples from valid compositions
set.seed(123)
samples <- gc_resample_compositions(grid, n = 1000, method = "uniform")

# Estimate ILR parameters from samples
params <- gc_ilr_params(samples)

# Build a variogram model (user-specified)
library(gstat)
vgm_template <- vgm(psill = 1, model = "Exp", range = 30, nugget = 0.01)

# Construct the multivariate gstat model
model <- gc_ilr_model(params, variogram_model = vgm_template)

# Create a spatial grid
x.range <- seq(0, 100, by = 5)
y.range <- seq(0, 100, by = 5)
grid_sf <- sf::st_as_sf(
  expand.grid(x = x.range, y = y.range),
  coords = c("x", "y"),
  crs = "local"
)

# Simulate 10 realizations
sims <- gc_sim_composition(model, grid_sf, nsim = 10, target_names = c("sand", "silt", "clay"))

# Result is a SpatRaster with layers:
# sand.sim1, silt.sim1, clay.sim1, sand.sim2, silt.sim2, clay.sim2, ...
print(sims)
```

## Key Functions

- **`gc_expand_bounds()`** - Generate a grid of valid compositions from
  constraints
- **`gc_resample_compositions()`** - Draw samples from valid
  compositions (uniform or soil texture)
- **`gc_ilr_params()`** - Calculate ILR mean and covariance
- **`gc_ilr_model()`** - Construct multivariate gstat model with
  Independent Univariate kriging (default) or LMC; supports conditional
  simulation
- **`gc_vgm_defaults()`** - Suggest reasonable variogram parameters
  based on data
- **`gc_fit_vgm()`** - Fit empirical variograms per ILR dimension with
  optional aggregation
- **`gc_sim_composition()`** - Simulate and back-transform to original
  units; supports conditioning on observed data

## Interactive Web Application

**New in v0.2.0**: Geocoda Uncertainty Mapper provides a graphical interface for interactive uncertainty visualization without programming:

```r
library(geocoda)
run_uncertainty_mapper()  # Opens browser-based dashboard
```

**Features**:
- Upload your CSV data (or use built-in examples)
- Configure hierarchical model visually
- Run simulations with live progress tracking
- Five interactive visualizations:
  - Probability maps (threshold exceedance risk)
  - Percentile maps (P10/P50/P90 uncertainty bands)
  - Depth profiles (uncertainty by soil layer)
  - 3D interactive scatter (compositional relationships)
  - Risk curves (exceedance probabilities)
- Export results (GeoTIFF, shapefile, JSON, HTML report)

Perfect for decision-makers, educators, and rapid prototyping. See [inst/shiny/uncertainty_mapper/README.md](inst/shiny/uncertainty_mapper/README.md) for complete guide.

---

## Advanced Features

### Stratified Hierarchical Modeling (v3.0+)

The package now supports **domain-aware geostatistical simulation** via stratified 
hierarchical modeling. This extension handles non-stationarity common in real-world 
soil systems:

**Key capabilities**:
- **Automatic domain detection**: `gc_identify_strata()` detects zones of statistical similarity
- **Hierarchical estimation**: `gc_fit_hierarchical_model()` estimates zone-specific parameters with shrinkage pooling
- **Per-zone simulation**: `gc_sim_hierarchical()` respects domain boundaries in ensemble generation
- **Ensemble aggregation**: Per-zone statistics, percentile maps, probability thresholds
- **Advanced diagnostics**: Cross-validation, entropy analysis, bootstrap uncertainty per zone
- **Decision support**: Risk assessment under asymmetric costs; carbon stock accounting and compliance
- **3D modeling**: `gc_fit_vgm_3d_per_zone()` and `gc_sim_hierarchical_3d_per_zone()` for depth-stratified simulation

**Example**: Detect soil zones, model separately, then generate ensemble respecting boundaries:

``` r
# 1. Detect zones automatically
strata <- gc_identify_strata(
  data = soil_data,
  comp_cols = c("SAND", "SILT", "CLAY")
)

# 2. Define hierarchy
hierarchy <- gc_define_hierarchy(
  data = soil_data,
  comp_cols = c("SAND", "SILT", "CLAY"),
  group_col = "zone"
)

# 3. Fit zone-specific parameters
hz_fit <- gc_fit_hierarchical_model(
  data = soil_data,
  object = hierarchy
)

# 4. Ensemble simulation (zone-aware)
ensemble <- gc_sim_hierarchical(
  object = hz_fit,
  n = 500,
  zone_vector = strata$zone_assignments
)

# 5. Per-zone decisions (risk, carbon accounting)
risk <- gc_risk_assessment_per_zone(
  ensemble = ensemble,
  zone_col = "zone_id",
  thresholds = c(SAND = 50),
  cost_matrix = cost_df
)
```

See **Comprehensive Stratification Guide** in `.github/development/` for complete examples.

See **Comprehensive Stratification Guide** in `.github/development/` for complete examples.

---

## Learning Pathways

**Complete learning resources**: See [.github/development/VIGNETTE_INDEX.md](.github/development/VIGNETTE_INDEX.md) for 10 vignettes organized by skill level and topic.

### Quick Start (Beginners)
1. Run `vignette("Soil Texture Workflow")` - Core workflow concepts (20 min)
2. Explore `run_uncertainty_mapper()` - Interactive visualization without code (10 min)
3. Try `vignette("Ensemble Analysis & Risk Assessment")` - Decision-making (15 min)

### Regional Mapper (SSURGO Focus)
1. `vignette("SSURGO Integration")` - Query national soil database (20 min)
2. Example: `examples/01_ssurgo_sda_integration.Rmd` - Live SSURGO queries (30 min)
3. Example: `examples/03_stratified_hierarchical_modeling.Rmd` - Build regional model (45 min)

### Advanced Statistician
1. `vignette("Hierarchical Backends Deep Dive")` - Compare three backends (25 min)
2. `vignette("Advanced Diagnostics")` - Validation framework (25 min)
3. `vignette("3D Geostatistics with Stratigraphic Coordinates")` - Full 3D modeling (30 min)

### Interactive Learning
- **Shiny App**: `run_uncertainty_mapper()` - Visual parameter exploration
- **FAQ & Troubleshooting**: `vignette("FAQ & Troubleshooting")`
- **Parameter Guide**: `vignette("Parameter Selection")`

---

### Model Types: Univariate vs LMC

The `gc_ilr_model()` function supports two modeling approaches via the
`model_type` parameter:

**Independent Univariate Kriging** (`model_type = "univariate"`,
default): - Models each ILR dimension separately without
cross-covariance terms - Numerically stable and robust (avoids
positive-definite issues) - Standard practice in compositional
geostatistics - Efficient for large problems - Recommended for most
applications

**Linear Model of Coregionalization** (`model_type = "lmc"`): - Includes
cross-covariance terms between all pairs of ILR dimensions -
Theoretically more complete - More numerically complex - Useful when
cross-correlation structure is important

``` r
# Univariate (default, recommended)
model_univ <- gc_ilr_model(params, variogram_model = vgm_template)

# LMC (alternative with cross-covariance terms)
model_lmc <- gc_ilr_model(params, variogram_model = vgm_template, model_type = "lmc")

# Both produce equivalent results for most practical applications,
# but univariate is faster and more numerically stable
```

### Conditional Simulation

The `gc_ilr_model()` and `gc_sim_composition()` functions support
**conditioning on observed data**:

``` r
# Create conditioning data in ILR space
obs_ilr <- compositions::ilr(compositions::acomp(observed_samples))
colnames(obs_ilr) <- c("ilr1", "ilr2", "ilr3")
conditioning_data <- sf::st_as_sf(
  cbind(obs_locations, as.data.frame(obs_ilr)),
  coords = c("x", "y"),
  crs = "local"
)

# Build model with conditioning data
model_cond <- gc_ilr_model(
  params, 
  variogram_model = vgm_template,
  data = conditioning_data  # Pass observed data for kriging
)

# Simulate with conditioning
sims_cond <- gc_sim_composition(
  model_cond, 
  grid_sf, 
  nsim = 10,
  observed_data = conditioning_data,
  target_names = c("sand", "silt", "clay")
)
```

**Key properties:** - **Unconditional (default)**: Independent
realizations from spatial distribution - **Conditional**: Realizations
honor observed values exactly at sample locations - **Uncertainty
reduction**: Away from sample locations, uncertainty decreases with
distance - Maintains sum constraint in all cases

### Empirical Variogram Fitting

The `gc_fit_vgm()` function automatically fits empirical variograms to
each ILR dimension:

``` r
# Fit variograms to ILR values at observed locations
data_with_ilr <- data.frame(
  x = sample_locations$x,
  y = sample_locations$y,
  ilr1 = ilr_values[, 1],
  ilr2 = ilr_values[, 2],
  ilr3 = ilr_values[, 3]
)

# Per-dimension results
fitted_per_dim <- gc_fit_vgm(
  params, 
  data = data_with_ilr,
  aggregate = FALSE
)

# Or get aggregated template for model building
fitted_template <- gc_fit_vgm(
  params,
  data = data_with_ilr,
  aggregate = TRUE
)

model <- gc_ilr_model(params, variogram_model = fitted_template)
```

**Key advantages:** - Captures spatial structure unique to each
component - Covariance-weighted aggregation ensures high-variance
dimensions contribute appropriately - Returns per-dimension fitted
models for detailed inspection

### Vignette & Soil Texture Workflow

A comprehensive workflow vignette demonstrates a complete soil texture
simulation:

``` r
vignette("soil_texture_workflow", package = "geocoda")
```

The vignette covers: - Soil texture constraints (USDA triangle) -
Bootstrap sampling of realistic compositions - ILR parameter
estimation - Variogram fitting and optimization - Grid creation and
simulation - Validation (sum constraints, range checks) - Visualization
of results - Tips and troubleshooting

## Dependencies

**Imports:** - `compositions` - ILR transformations - `gstat` -
Geostatistical modeling and simulation - `terra` - Raster data
handling - `sf` - Spatial feature support

**Suggests:** - `aqp` - Soil texture bootstrapping utilities -
`testthat` - Unit testing - `knitr`, `rmarkdown` - Documentation

## Citation

If you use this package in your research, please cite it. Run
`citation("geocoda")` for citation information.

## License

MIT License. See [LICENSE.md](LICENSE.md) for details.

## Author

Andrew G. Brown

## References

- Aitchison, J. (1986). *The Statistical Analysis of Compositional
  Data*. Chapman and Hall, London. ISBN: 978-0-412-28060-3.

- Egozcue, J. J., Pawlowsky-Glahn, V., Mateu-Figueras, G., &
  Barceló-Vidal, C. (2003). Isometric Log-Ratio Transformations for
  Compositional Data Analysis. *Mathematical Geology*, 35(3), 279–300.
  <https://doi.org/10.1023/A:1023818214614>

- Chilès, J. P., & Delfiner, P. (2012). *Geostatistics: Modeling Spatial
  Uncertainty* (2nd ed.). Wiley, Hoboken, NJ.
  <https://doi.org/10.1002/9781118383087>

- Tolosana-Delgado, R., Mueller, U., & van den Boogaart, K. G. (2012).
  Geostatistics for Compositional Data: An Overview. *Mathematical
  Geosciences*, 44(4), 465–479.
  <https://doi.org/10.1007/s11004-011-9363-4>

- Pebesma, E. J. (2004). Multivariable geostatistics in S: the gstat
  package. *Computers & Geosciences*, 30, 683–691.
  <https://doi.org/10.1016/j.cageo.2004.03.012>
