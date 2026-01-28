# Simulate and Back-Transform Compositional Fields

Generate spatial realizations of compositional data by predicting in ILR
space using a gstat model, then back-transforming to the original units.
Supports both unconditional simulation (independent of data) and
conditional simulation (honoring observed values at sample locations).
Multiple realizations are stacked into a multi-layer
[`terra::SpatRaster`](https://rspatial.github.io/terra/reference/SpatRaster-class.html).

## Usage

``` r
gc_sim_composition(
  model,
  locations,
  nsim = 1,
  target_names = NULL,
  crs = "local",
  observed_data = NULL,
  nmax = NULL
)
```

## Arguments

- model:

  A `gstat` object (typically from [`gc_ilr_model()`](gc_ilr_model.md)).
  For conditional simulation, the model **must be built with the `data`
  parameter** so that conditioning data is embedded in the model.

- locations:

  An `sf` object or data frame defining the simulation grid. Must
  contain `x` and `y` columns. If a data frame is provided, it will be
  converted to an sf object with coordinates in the specified CRS.

- nsim:

  Number of realizations (default 1).

- target_names:

  Character vector of original component names. If `NULL`, inferred from
  the `gstat` model or defaults to `c("comp1", "comp2", ...)`.

- crs:

  Coordinate reference system (default `"local"`).

- observed_data:

  An optional `sf` object or data frame containing observed
  compositional samples for conditional simulation. Must have the same
  structure as the data used to build the model: columns `x`, `y` for
  coordinates, and columns for each ILR dimension named `ilr1`, `ilr2`,
  etc. (in ILR space). If `NULL` (default), unconditional simulation is
  performed.

- nmax:

  Maximum number of nearby observations to use for prediction at each
  location (default `NULL`, uses all observations). Setting a smaller
  `nmax` (e.g., 12-15) can:

  - Improve computational efficiency for large datasets

  - Reduce memory usage during prediction

  - Create more local uncertainty estimates (Screen Effect)

  - May introduce slight artifacts at spatial boundaries

  Recommended: Use `NULL` for small/medium datasets, `nmax=12-20` for
  large datasets.

## Value

A
[`terra::SpatRaster`](https://rspatial.github.io/terra/reference/SpatRaster-class.html)
with layers named according to the pattern: `<component>.sim<N>` where
`component` is from `target_names` and `sim<N>` indicates the
realization number. For example, with 3 components and 2 realizations:
`comp1.sim1, comp2.sim1, comp3.sim1, comp1.sim2, ...`

## Details

This function:

1.  Prepares spatial coordinates from `locations`

2.  Validates `observed_data` if provided (must have ilr columns and xy
    coords)

3.  Calls
    [`gstat::predict()`](https://r-spatial.github.io/gstat/reference/predict.gstat.html)
    to generate ILR values for `nsim` realizations. If the model was
    built with conditioning data, the predictions automatically honor
    that data at the sample locations. The `nmax` parameter limits the
    neighborhood used for kriging.

4.  Extracts simulated ILR columns (identified by pattern `.sim<N>`)

5.  Groups ILR values by realization

6.  Back-transforms each realization using
    [`compositions::ilrInv()`](https://rdrr.io/pkg/compositions/man/ilr.html)

7.  Rescales from 0-1 to 0-100

8.  Stacks results into a single `SpatRaster` with proper naming

**Conditional vs. Unconditional:** When the model is built with
`data = NULL`, all realizations are independent draws from the spatial
distribution (unconditional). When the model is built with
`data = <conditioning_data>`, all realizations are conditioned on those
values: predictions at sample locations exactly reproduce the observed
values, while predictions elsewhere reflect uncertainty updated by the
conditioning data.

**Important:** For conditional simulation to work, the model \*\*must
have been built with the conditioning data embedded via the `data`
parameter to [`gc_ilr_model()`](gc_ilr_model.md). The `observed_data`
parameter here is only for validation. Conditioning happens
automatically because the model was created with that data.

**Neighborhood Size Effects (Screen Effect):** When `nmax` is small
relative to data density, kriging typically uses only the nearest `nmax`
observations. This can create more localized uncertainty estimates and
may accelerate computation, but can also:

- Introduce discontinuities at boundaries between neighborhoods

- Miss spatial structure information from distant data

- Bias predictions if far-field data carries important variance
  information

The output strictly honors the sum constraint: all rows sum to
`target_sum` (typically 100 for percentages) within floating-point
precision.

## Examples

``` r
if (FALSE) { # \dontrun{
library(sf)
library(gstat)

# Assuming ilr_params and model are already defined
# (See workflow example in the package README)

# Create simulation grid
x.range <- seq(0, 100, by = 5)
y.range <- seq(0, 100, by = 5)
grid_df <- expand.grid(x = x.range, y = y.range)

# Unconditional simulation: 5 realizations
result <- gc_sim_composition(model, grid_df, nsim = 5,
                             target_names = c("sand", "silt", "clay"))
print(result)

# Conditional simulation: honor observed data
result_cond <- gc_sim_composition(model, grid_df, nsim = 5,
                                  target_names = c("sand", "silt", "clay"),
                                  observed_data = sample_df)

# Large dataset: use neighborhood limiting for efficiency
result_nmax <- gc_sim_composition(model, grid_df, nsim = 5,
                                  target_names = c("sand", "silt", "clay"),
                                  nmax = 15)

# Access results
terra::values(result)
} # }
```
