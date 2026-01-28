# Expand Compositional Bounds into a Valid Grid

Generate a grid of valid compositions from component constraints. Uses
[`expand.grid()`](https://rdrr.io/r/base/expand.grid.html) to create
combinations of sequential component values within specified bounds,
then filters to rows where component sums equal the target sum (within
floating-point tolerance).

## Usage

``` r
gc_expand_bounds(constraints, step = 0.1, target_sum = 100, tol = 1e-06)
```

## Arguments

- constraints:

  A named list of lists or data frame containing `min` and `max` for
  each component. Example:
  `list(SAND = list(min = 0, max = 40), SILT = list(min = 50, max = 80), CLAY = list(min = 10, max = 20))`.

- step:

  Numeric resolution for the grid sequences (default 0.1).

- target_sum:

  Numeric constant sum for valid compositions (default 100).

- tol:

  Floating-point tolerance for sum constraint (default 1e-6).

## Value

A data frame where columns are components and rows are valid
compositions summing to `target_sum`.

## Details

This function performs:

1.  Validation that min values do not sum to more than `target_sum`

2.  Dynamic sequence generation for each component based on bounds

3.  Cartesian product expansion via
    [`expand.grid()`](https://rdrr.io/r/base/expand.grid.html)

4.  Filtering to compositions with sum approximately equal to
    `target_sum`

## Examples

``` r
# Define sand, silt, clay bounds
constraints <- list(
  SAND = list(min = 0, max = 40),
  SILT = list(min = 50, max = 80),
  CLAY = list(min = 10, max = 20)
)

# Expand with 1% resolution
grid <- gc_expand_bounds(constraints, step = 1.0, target_sum = 100)
nrow(grid)
#> [1] 341
head(grid)
#>   SAND SILT CLAY
#> 1   40   50   10
#> 2   39   51   10
#> 3   38   52   10
#> 4   37   53   10
#> 5   36   54   10
#> 6   35   55   10
```
