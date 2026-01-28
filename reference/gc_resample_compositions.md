# Bootstrap Compositional Samples

Draw samples from a valid composition grid to create a "source
population" for estimating covariance structure in ILR space. Supports
both uniform random sampling and soil texture-aware bootstrapping (if
`aqp` is available).

## Usage

``` r
gc_resample_compositions(
  composition_grid,
  n = 1000,
  method = "uniform",
  seed = NULL
)
```

## Arguments

- composition_grid:

  A data frame of valid compositions (typically output from
  [`gc_expand_bounds()`](gc_expand_bounds.md)).

- n:

  Number of samples to draw (default 1000).

- method:

  Sampling method: `"uniform"` (simple random sampling) or
  `"soil_texture"` (uses
  [`aqp::bootstrapSoilTexture()`](https://ncss-tech.github.io/aqp/reference/bootstrapSoilTexture.html)
  if available, falls back to uniform if not).

- seed:

  Optional random seed for reproducibility.

## Value

A list with:

- `samples`: A data frame of sampled compositions.

- `method`: The method used.

- `n`: Number of samples requested.

## Details

- **Uniform sampling**: Performs simple random sampling of row indices
  without replacement. If `n > nrow(composition_grid)`, samples with
  replacement.

- **Soil texture sampling**: Attempts to delegate to
  [`aqp::bootstrapSoilTexture()`](https://ncss-tech.github.io/aqp/reference/bootstrapSoilTexture.html).
  If `aqp` is not installed, falls back to uniform sampling with a
  warning.

## Examples

``` r
# Create a simple composition grid
constraints <- list(
  SAND = list(min = 0, max = 40),
  SILT = list(min = 50, max = 80),
  CLAY = list(min = 10, max = 20)
)
grid <- gc_expand_bounds(constraints, step = 1.0, target_sum = 100)

# Uniform sampling
set.seed(42)
uniform_samps <- gc_resample_compositions(grid, n = 100, method = "uniform")
nrow(uniform_samps$samples)
#> [1] 100
```
