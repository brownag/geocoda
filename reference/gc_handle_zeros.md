# Handle Zeros and Censored Data in Compositional Data

Process compositional data containing zeros or below-detection-limit
(censored) values using log-ratio imputation. Imputation is critical for
compositional data analysis because the geometric mean of zeros is
undefined, making standard covariance analysis unstable. This function
provides three strategies: multiplicative zero imputation, additive zero
imputation, or probabilistic imputation using the EM algorithm.

## Usage

``` r
gc_handle_zeros(comp_data, method = "mzero", dl = NULL, tolerance = 1e-10)
```

## Arguments

- comp_data:

  A data frame with columns representing compositional parts (e.g.,
  `sand`, `silt`, `clay`). Rows with at least one zero or negative value
  are candidates for imputation.

- method:

  Character string specifying imputation approach:

  - `"mzero"` (default): Multiplicative zero imputation using
    zCompositions::multLRM()

  - `"azero"`: Additive zero imputation using zCompositions::addLRM()

  - `"lrem"`: Log-ratio expectation-maximization using
    zCompositions::lrEM()

  Recommended: Use `"lrem"` for data with many zeros; use `"mzero"` for
  sparse zeros. Both are log-ratio based and preserve compositional
  geometry.

- dl:

  Numeric vector of detection limits (one per column) for handling
  censored data. Default `NULL` treats all positive values as observed.
  If provided, values below `dl[i]` are treated as left-censored in
  column `i`. Only used if `method = "lrem"`.

- tolerance:

  Numeric tolerance for identifying zeros (default 1e-10). Values with
  absolute value below this threshold are treated as zeros.

## Value

A list containing:

- `imputed_data`: Data frame with imputed compositional values

- `n_zeros_imputed`: Integer count of zero/censored values imputed

- `imputation_rate`: Proportion of values imputed (n_zeros /
  total_values)

- `method_used`: Character string with imputation method name

- `row_status`: Factor indicating which rows were modified:

  - `"observed"`: No imputation needed (no zeros/negative values)

  - `"imputed"`: At least one zero/censored value replaced

  - `"failed"`: Imputation failed (row excluded from analysis)

## Details

**Zero Imputation Strategies:**

**Multiplicative Zero Replacement (mzero)**: Replaces zeros with a small
multiple of the detection limit, then applies log-ratio closure. Fast,
appropriate for few isolated zeros.

**Additive Zero Replacement (azero)**: Adds a small constant to all
values before closure. Conservative and robust, but can distort
low-abundance components.

**Log-Ratio EM (lrem)**: Probabilistic imputation using
expectation-maximization on log-ratio transformed data. Respects
compositional geometry while honoring censoring patterns. Most
theoretically sound but slower than replacement methods.

**Detection Limits:** If `dl` is provided and method = `"lrem"`, values
below their detection limit are treated as left-censored (uncertainty in
exact value). EM iterates to estimate most likely imputed values
consistent with the censoring pattern and covariance structure.

## Examples

``` r
if (FALSE) { # \dontrun{
# Example: Soil texture data with some zeros or missing detections
soil_samples <- data.frame(
  sand = c(40, 35, 0, 42, 38),
  silt = c(35, 40, 45, 36, 40),
  clay = c(25, 25, 55, 22, 22)
)

# Check for zeros before imputation
has_zeros <- rowSums(soil_samples == 0) > 0
print(paste("Rows with zeros:", sum(has_zeros)))

# Impute using multiplicative zero replacement
result_mzero <- gc_handle_zeros(soil_samples, method = "mzero")
print(result_mzero$imputed_data)
print(paste("Imputation rate:", result_mzero$imputation_rate))

# Impute using log-ratio EM (more principled)
result_lrem <- gc_handle_zeros(soil_samples, method = "lrem")
print(result_lrem$imputed_data)

# With detection limits (censored measurements)
detection_limits <- c(sand = 1, silt = 1, clay = 1)
result_dl <- gc_handle_zeros(soil_samples, method = "lrem", dl = detection_limits)
} # }
```
