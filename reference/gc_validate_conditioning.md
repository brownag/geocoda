# Validate Conditional Simulation Data Honoring

Assess the quality and accuracy of conditional kriging by extracting
kriged predictions at the conditioning (observation) locations and
comparing them to the actual observed values. This provides a
quantitative measure of how well the model honors the conditioning data
in ILR space.

## Usage

``` r
gc_validate_conditioning(
  model,
  observed_data,
  metrics = c("rmse", "mae", "mean_error")
)
```

## Arguments

- model:

  A gstat kriging model (typically from
  [`gc_ilr_model()`](gc_ilr_model.md)) that was built with conditioning
  data via the `data` parameter. The model must contain observed ILR
  values at sample locations.

- observed_data:

  A data frame or sf object with the original conditioning data,
  containing columns `x`, `y` for spatial coordinates and columns
  `ilr1`, `ilr2`, etc., with ILR values in the same space as the model.

- metrics:

  Character vector specifying which error metrics to compute. Default
  `c("rmse", "mae", "mean_error")`. Options include:

  - `"rmse"`: Root mean squared error

  - `"mae"`: Mean absolute error

  - `"mean_error"`: Mean signed error (bias)

  - `"median_error"`: Median absolute error

  - `"sd_error"`: Standard deviation of errors

## Value

A list containing:

- `predictions_at_obs`: Data frame with kriged predictions at each
  observation location

- `observed_values`: Data frame with original observed ILR values

- `residuals`: Difference (observed - predicted) for each ILR dimension

- `error_metrics`: Data frame with error statistics (one row per ILR
  dimension)

- `overall_metrics`: Data frame with aggregated error statistics across
  all dimensions

## Details

**Conditional Honoring Assessment:**

Conditional kriging is intended to provide spatial predictions that
exactly (or very nearly) match observed values at the conditioning
locations. In practice, numerical precision limits and kriging
discretization may result in small residuals even for perfectly honored
data.

This function:

1.  Extracts kriged predictions at each observation location using the
    model

2.  Computes residuals (observed minus predicted) for each ILR dimension

3.  Calculates error metrics (RMSE, MAE, etc.) per dimension and overall

4.  Returns detailed diagnostics for reviewing data honoring quality

**Interpretation:**

- RMSE close to zero indicates excellent data honoring

- Mean error close to zero indicates unbiased predictions

- Non-zero residuals may indicate:

  - Numerical precision effects (typically \< 1e-10)

  - Model convergence issues (investigate with univariate kriging)

  - Grid resolution effects (predictions on coarse grids may miss exact
    locations)

For compositional data (original sand/silt/clay), back-transform RMSE in
ILR space to composition space using `ilrInv()` to assess error in
original units.

## Examples

``` r
if (FALSE) { # \dontrun{
# Assuming model and observed_data are available from gc_ilr_model()
# and original data frame with spatial coordinates

# Build model with conditioning data
model <- gc_ilr_model(ilr_params, data = conditioning_data)

# Validate conditioning
validation <- gc_validate_conditioning(model, conditioning_data)

# Review error metrics
print(validation$error_metrics)
print(validation$overall_metrics)

# Check for problematic observations
high_residuals <- which(abs(validation$residuals$ilr1) > 0.1)
if (length(high_residuals) > 0) {
  print("Observations with large residuals:")
  print(validation$residuals[high_residuals, ])
}

# Visualize residuals vs spatial location
plot(validation$predictions_at_obs$x,
     validation$residuals$ilr1,
     main = "ILR1 Residuals vs X Coordinate",
     xlab = "X", ylab = "Residual (Observed - Predicted)")
} # }
```
