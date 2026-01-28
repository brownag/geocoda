# Build Compositional Geostatistical Model

Constructs a gstat model for the ILR-transformed variables using either
Independent Univariate Kriging or Linear Model of Coregionalization
(LMC).

## Usage

``` r
gc_ilr_model(
  ilr_params,
  variogram_model,
  data = NULL,
  model_type = "univariate"
)
```

## Arguments

- ilr_params:

  A list from [`gc_ilr_params()`](gc_ilr_params.md) containing mean,
  cov, names.

- variogram_model:

  A [`vgm()`](https://r-spatial.github.io/gstat/reference/vgm.html)
  object defining the base variogram structure.

- data:

  Optional `sf` object with ILR columns for conditioning. If provided,
  Conditional Simulation is performed. If NULL (default), Unconditional
  Simulation is performed.

- model_type:

  Character string specifying the approach: `"univariate"` (default) or
  `"lmc"`. Univariate is numerically stable and standard practice. LMC
  includes cross-covariance terms between ILR dimensions.

## Value

A `gstat` object representing the model.

## Details

**Independent Univariate Kriging** (`model_type = "univariate"`): Models
each ILR dimension separately without cross-covariance terms. This is:

- Numerically stable (avoids positive-definite issues in LMC fitting)

- Standard practice in compositional geostatistics

- Robust across different datasets

- Efficient for large problems

The ILR transformation already decorrelates the data significantly, so
ignoring spatial cross-correlation between ILR coordinates has minimal
impact.

**Linear Model of Coregionalization** (`model_type = "lmc"`): Includes
cross-covariance terms between all pairs of ILR dimensions. This is:

- Theoretically more complete

- More numerically complex

- Useful when cross-correlation structure is important

For conditional simulation, the conditioning data must be passed to this
function (not to [`predict()`](https://rdrr.io/r/stats/predict.html)).
The model then automatically uses that data during prediction.
