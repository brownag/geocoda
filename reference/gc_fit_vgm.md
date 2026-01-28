# Optimize Variogram Parameters

Fit empirical variograms to ILR values at observed locations and
estimate optimal variogram parameters for each ILR dimension. Results
can be returned per-dimension or aggregated into a single template
variogram for use in multivariate simulation. Includes diagnostics for
LMC admissibility.

## Usage

``` r
gc_fit_vgm(
  ilr_params,
  data,
  vgm_model_type = "Exp",
  maxdist = NULL,
  width = NULL,
  aggregate = FALSE,
  correct.diagonal = 1.01,
  fit.ranges = FALSE
)
```

## Arguments

- ilr_params:

  A list returned by [`gc_ilr_params()`](gc_ilr_params.md).

- data:

  A data frame or spatial object (sf) with columns `x`, `y` for spatial
  coordinates and columns named `ilr1`, `ilr2`, etc., containing ILR
  values at observed locations.

- vgm_model_type:

  Character string specifying variogram model type (e.g., `"Exp"`,
  `"Sph"`, `"Gau"`). Default `"Exp"`.

- maxdist:

  Maximum distance for empirical variogram calculation (default `NULL`,
  uses all pairs).

- width:

  Lag width for empirical variogram bins (default `NULL`, automatic).

- aggregate:

  Logical. If `TRUE`, aggregate fitted parameters across ILR dimensions
  using weighted average sills (weights = covariance diagonal). If
  `FALSE` (default), return per-dimension results.

- correct.diagonal:

  Numeric sill correction factor for LMC stability (default `1.01`).
  Multiplies each marginal (diagonal) sill by this factor to improve
  positive-definiteness of the LMC sill matrix and prevent numerical
  issues. Recommended range: 1.00-1.05. Use `1.01` for typical
  applications.

- fit.ranges:

  Logical. If `FALSE` (default), fixes ranges to fitted values during
  LMC construction to avoid over-parameterization. If `TRUE`, allows
  range re-optimization when building LMC models.

## Value

A list with:

- If `aggregate = FALSE`: A list of length `D-1` where each element
  contains `fitted_vgm` and `empirical_vgm` for that ILR dimension.

- If `aggregate = TRUE`: A single `vgm` object with aggregated
  parameters (mean range, weighted mean sill, weighted mean nugget).

- Attribute `ilr_dimension_names`: Names of ILR dimensions (ilr1, ilr2,
  ...)

- Attribute `fitted_params`: Data frame with per-dimension parameters

- Attribute `lmc_admissibility`: Logical indicating if sill matrix is
  positive-definite

## Details

This function:

1.  For each ILR dimension i:

    - Computes empirical variogram using
      [`gstat::variogram()`](https://r-spatial.github.io/gstat/reference/variogram.html)

    - Provides intelligent initial model (range ~ spatial extent / 3)

    - Fits parametric model using
      [`gstat::fit.variogram()`](https://r-spatial.github.io/gstat/reference/fit.variogram.html)

2.  Optionally aggregates results using covariance-weighted averaging

3.  Applies `correct.diagonal` to diagonal sills to ensure LMC
    positive-definiteness

Weighted aggregation combines results when building a single LMC model,
ensuring dimensions with higher variance contribute proportionally.

**LMC Admissibility:** The sill matrix (covariance at distance infinity)
must be positive-definite for valid spatial covariance. If eigenvalues
of the sill matrix have any zero or negative values, the LMC structure
is inadmissible. The `correct.diagonal` parameter helps stabilize the
sill matrix by inflating diagonal terms, which typically improves
admissibility without distorting the model significantly.

## Examples

``` r
if (FALSE) { # \dontrun{
# Example with observed ILR values on a spatial grid
library(gstat)

# Create sample data with spatial coordinates and ILR values
n_points <- 50
coords_df <- data.frame(
  x = runif(n_points, 0, 100),
  y = runif(n_points, 0, 100),
  ilr1 = rnorm(n_points, mean = 0.5, sd = 0.8),
  ilr2 = rnorm(n_points, mean = -0.2, sd = 0.6)
)

# Estimate parameters from bootstrap samples
samples <- data.frame(
  sand = c(20, 25, 30, 22),
  silt = c(60, 55, 50, 58),
  clay = c(20, 20, 20, 20)
)
params <- gc_ilr_params(samples)

# Fit empirical variograms
fitted <- gc_fit_vgm(
  params,
  coords_df,
  vgm_model_type = "Exp",
  aggregate = FALSE
)

# Or get single aggregated template
fitted_agg <- gc_fit_vgm(
  params,
  coords_df,
  vgm_model_type = "Exp",
  aggregate = TRUE
)
} # }
```
