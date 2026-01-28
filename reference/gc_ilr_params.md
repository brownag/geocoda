# Estimate ILR Parameters from Compositional Samples

Convert compositional samples to Isometric Log-Ratio (ILR) space and
compute the mean vector and covariance matrix. These statistics are used
to parameterize the multivariate geostatistical model.

## Usage

``` r
gc_ilr_params(samples)
```

## Arguments

- samples:

  A data frame of compositions (each row sums to constant, typically
  100% or 1).

## Value

A list containing:

- `mean`: Vector of ILR means (length D-1 where D is number of
  components)

- `cov`: Covariance matrix of ILR values (dimensions (D-1) x (D-1))

- `names`: Original component names (character vector)

- `base_class`: Class used for transformation ("acomp")

## Details

This function:

1.  Stores original column names for later use

2.  Converts to
    [compositions::acomp](https://rdrr.io/pkg/compositions/man/acomp.html)
    (absolute composition)

3.  Applies
    [compositions::ilr](https://rdrr.io/pkg/compositions/man/ilr.html)
    transformation

4.  Computes column means and covariance matrix of ILR values

The ILR transformation eliminates the sum constraint, allowing standard
multivariate geostatistics to be applied.

## Examples

``` r
# Simulate some simple compositions
samples <- data.frame(
  sand = c(20, 25, 30, 22),
  silt = c(60, 55, 50, 58),
  clay = c(20, 20, 20, 20)
)

params <- estimate_ilr_params(samples)
#> Error in estimate_ilr_params(samples): could not find function "estimate_ilr_params"
str(params)
#> Error: object 'params' not found
print(params$mean)
#> Error: object 'params' not found
print(params$cov)
#> Error: object 'params' not found
```
