# Suggest Reasonable Variogram Default Parameters

Inspect data extent and ILR covariance structure to suggest ballpark
values for range, nugget, and sill of a variogram model. This is a
convenience function to help users get started before fine-tuning.

## Usage

``` r
gc_vgm_defaults(ilr_params, extent)
```

## Arguments

- ilr_params:

  A list returned by [`gc_ilr_params()`](gc_ilr_params.md).

- extent:

  A numeric vector of length 4: `c(xmin, ymin, xmax, ymax)`, or a
  spatial object (sf, terra) from which extent can be extracted.

## Value

A list containing:

- `range`: Suggested range parameter (approximately 1/3 of the diagonal
  extent)

- `nugget`: Suggested nugget ratio (typically 0.01 to 0.05 of average
  sill)

- `mean_sill`: Mean of the diagonal covariance terms

## Details

This function:

1.  Calculates the spatial extent diagonal (Euclidean distance from min
    to max)

2.  Suggests range as ~1/3 of the extent diagonal

3.  Calculates mean of diagonal covariance terms as representative sill

4.  Suggests nugget as 1% of the mean sill

These are heuristics and should be refined using empirical variography
(see [`gc_fit_vgm()`](gc_fit_vgm.md)).

## Examples

``` r
samples <- data.frame(
  sand = c(20, 25, 30, 22),
  silt = c(60, 55, 50, 58),
  clay = c(20, 20, 20, 20)
)

params <- gc_ilr_params(samples)
extent <- c(0, 0, 100, 100)

suggestions <- gc_vgm_defaults(params, extent)
print(suggestions)
#> $range
#> [1] 47.14045
#> 
#> $nugget
#> [1] 0.0001699233
#> 
#> $mean_sill
#> [1] 0.01699233
#> 
```
