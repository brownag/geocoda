#' Expand Compositional Bounds into a Valid Grid
#'
#' Generate a grid of valid compositions from component constraints. Uses
#' `expand.grid()` to create combinations of sequential component values
#' within specified bounds, then filters to rows where component sums equal the
#' target sum (within floating-point tolerance).
#'
#' @param constraints A named list of lists or data frame containing `min` and
#'   `max` for each component. Example: `list(SAND = list(min = 0, max = 40),
#'   SILT = list(min = 50, max = 80), CLAY = list(min = 10, max = 20))`.
#' @param step Numeric resolution for the grid sequences (default 0.1).
#' @param target_sum Numeric constant sum for valid compositions (default 100).
#' @param tol Floating-point tolerance for sum constraint (default 1e-6).
#'
#' @return A data frame where columns are components and rows are valid
#'   compositions summing to `target_sum`.
#'
#' @details
#' This function performs:
#' 1. Validation that min values do not sum to more than `target_sum`
#' 2. Dynamic sequence generation for each component based on bounds
#' 3. Cartesian product expansion via `expand.grid()`
#' 4. Filtering to compositions with sum approximately equal to `target_sum`
#'
#' @references
#' Aitchison, J. (1986). The Statistical Analysis of Compositional Data.
#' Chapman and Hall, London.
#'
#' @examples
#' # Define sand, silt, clay bounds
#' constraints <- list(
#'   SAND = list(min = 0, max = 40),
#'   SILT = list(min = 50, max = 80),
#'   CLAY = list(min = 10, max = 20)
#' )
#'
#' # Expand with 1% resolution
#' grid <- gc_expand_bounds(constraints, step = 1.0, target_sum = 100)
#' nrow(grid)
#' head(grid)
#'
#' @importFrom stats setNames
#' @export
gc_expand_bounds <- function(constraints,
                             step = 0.1,
                             target_sum = 100,
                             tol = 1e-6) {
  # Extract component names and bounds
  comp_names <- names(constraints)
  n_comps <- length(comp_names)

  # Validate constraints
  min_sum <- sum(sapply(constraints, function(x) x$min))
  max_sum <- sum(sapply(constraints, function(x) x$max))

  if (min_sum > target_sum) {
    stop(
      "Minimum values sum to ",
      round(min_sum, 4),
      ", which exceeds target_sum = ",
      target_sum,
      ". Constraints are physically impossible."
    )
  }

  if (max_sum < target_sum) {
    stop(
      "Maximum values sum to ",
      round(max_sum, 4),
      ", which is less than target_sum = ",
      target_sum,
      ". Constraints are physically impossible."
    )
  }

  # Build sequences for each component
  sequences <- stats::setNames(vector("list", n_comps), comp_names)
  for (i in seq_along(comp_names)) {
    comp <- comp_names[i]
    min_val <- constraints[[comp]]$min
    max_val <- constraints[[comp]]$max
    sequences[[comp]] <- seq(from = min_val, to = max_val, by = step)
  }

  # Expand grid
  grid <- do.call(expand.grid, sequences)

  # Filter rows where sum =~ target_sum
  row_sums <- rowSums(grid)
  valid_idx <- abs(row_sums - target_sum) <= tol
  grid <- grid[valid_idx, ]

  # Reset row names
  rownames(grid) <- NULL

  return(grid)
}


#' Bootstrap Compositional Samples
#'
#' Draw samples from a valid composition grid to create a "source population"
#' for estimating covariance structure in ILR space. Supports both uniform
#' random sampling and soil texture-aware bootstrapping (if `aqp` is available).
#'
#' @param composition_grid A data frame of valid compositions (typically output
#'   from `gc_expand_bounds()`).
#' @param n Number of samples to draw (default 1000).
#' @param method Sampling method: `"uniform"` (simple random sampling) or
#'   `"soil_texture"` (uses [aqp::bootstrapSoilTexture()] if available,
#'   falls back to uniform if not).
#' @param seed Optional random seed for reproducibility.
#'
#' @return A list with:
#'   - `samples`: A data frame of sampled compositions.
#'   - `method`: The method used.
#'   - `n`: Number of samples requested.
#'
#' @details
#' - **Uniform sampling**: Performs simple random sampling of row indices without
#'   replacement. If `n > nrow(composition_grid)`, samples with replacement.
#' - **Soil texture sampling**: Attempts to delegate to
#'   `aqp::bootstrapSoilTexture()`. If `aqp` is not installed, falls back to
#'   uniform sampling with a warning.
#'
#' @examples
#' # Create a simple composition grid
#' constraints <- list(
#'   SAND = list(min = 0, max = 40),
#'   SILT = list(min = 50, max = 80),
#'   CLAY = list(min = 10, max = 20)
#' )
#' grid <- gc_expand_bounds(constraints, step = 1.0, target_sum = 100)
#'
#' # Uniform sampling
#' set.seed(42)
#' uniform_samps <- gc_resample_compositions(grid, n = 100, method = "uniform")
#' nrow(uniform_samps$samples)
#'
#' @importFrom compositions acomp ilr
#' @importFrom stats cov
#' @export
gc_resample_compositions <- function(composition_grid,
                                     n = 1000,
                                     method = "uniform",
                                     seed = NULL) {
  if (!is.null(seed)) {
    set.seed(seed)
  }

  method <- match.arg(method, c("uniform", "soil_texture"))

  if (method == "soil_texture") {
    # Try to use aqp::bootstrapSoilTexture
    if (!requireNamespace("aqp", quietly = TRUE)) {
      warning(
        "aqp package not installed. Falling back to uniform sampling.",
        call. = FALSE
      )
      method <- "uniform"
    }
  }

  if (method == "uniform") {
    # Simple random sampling
    n_grid <- nrow(composition_grid)
    replace_flag <- n > n_grid
    idx <- sample(n_grid, size = n, replace = replace_flag)
    samples_df <- composition_grid[idx, ]
    rownames(samples_df) <- NULL
  } else if (method == "soil_texture") {
    # Use aqp::bootstrapSoilTexture
    result <- aqp::bootstrapSoilTexture(composition_grid, n = n, method = "normal")
    samples_df <- result$samples
  }

  return(list(
    samples = samples_df,
    method = method,
    n = n
  ))
}


#' Handle Zeros and Censored Data in Compositional Data
#'
#' Process compositional data containing zeros or below-detection-limit (censored)
#' values using log-ratio imputation. Imputation is critical for compositional
#' data analysis because the geometric mean of zeros is undefined, making standard
#' covariance analysis unstable. This function provides three strategies:
#' multiplicative zero imputation, additive zero imputation, or probabilistic
#' imputation using the EM algorithm.
#'
#' @param comp_data A data frame with columns representing compositional parts
#'   (e.g., `sand`, `silt`, `clay`). Rows with at least one zero or negative
#'   value are candidates for imputation.
#' @param method Character string specifying imputation approach:
#'   - `"mzero"` (default): Multiplicative zero imputation using zCompositions::multLRM()
#'   - `"azero"`: Additive zero imputation using zCompositions::addLRM()
#'   - `"lrem"`: Log-ratio expectation-maximization using zCompositions::lrEM()
#'
#'   Recommended: Use `"lrem"` for data with many zeros; use `"mzero"` for sparse
#'   zeros. Both are log-ratio based and preserve compositional geometry.
#'
#' @param dl Numeric vector of detection limits (one per column) for handling
#'   censored data. Default `NULL` treats all positive values as observed.
#'   If provided, values below `dl[i]` are treated as left-censored in column `i`.
#'   Only used if `method = "lrem"`.
#'
#' @param tolerance Numeric tolerance for identifying zeros (default 1e-10).
#'   Values with absolute value below this threshold are treated as zeros.
#'
#' @return A list containing:
#'   - `imputed_data`: Data frame with imputed compositional values
#'   - `n_zeros_imputed`: Integer count of zero/censored values imputed
#'   - `imputation_rate`: Proportion of values imputed (n_zeros / total_values)
#'   - `method_used`: Character string with imputation method name
#'   - `row_status`: Factor indicating which rows were modified:
#'     - `"observed"`: No imputation needed (no zeros/negative values)
#'     - `"imputed"`: At least one zero/censored value replaced
#'     - `"failed"`: Imputation failed (row excluded from analysis)
#'
#' @details
#' **Zero Imputation Strategies:**
#'
#' **Multiplicative Zero Replacement (mzero)**:
#' Replaces zeros with a small multiple of the detection limit, then applies
#' log-ratio closure. Fast, appropriate for few isolated zeros.
#'
#' **Additive Zero Replacement (azero)**:
#' Adds a small constant to all values before closure. Conservative and robust,
#' but can distort low-abundance components.
#'
#' **Log-Ratio EM (lrem)**:
#' Probabilistic imputation using expectation-maximization on log-ratio transformed
#' data. Respects compositional geometry while honoring censoring patterns.
#' Most theoretically sound but slower than replacement methods.
#'
#' **Detection Limits:**
#' If `dl` is provided and method = `"lrem"`, values below their detection limit
#' are treated as left-censored (uncertainty in exact value). EM iterates to
#' estimate most likely imputed values consistent with the censoring pattern
#' and covariance structure.
#'
#' @examples
#' \dontrun{
#' # Example: Soil texture data with some zeros or missing detections
#' soil_samples <- data.frame(
#'   sand = c(40, 35, 0, 42, 38),
#'   silt = c(35, 40, 45, 36, 40),
#'   clay = c(25, 25, 55, 22, 22)
#' )
#'
#' # Check for zeros before imputation
#' has_zeros <- rowSums(soil_samples == 0) > 0
#' print(paste("Rows with zeros:", sum(has_zeros)))
#'
#' # Impute using multiplicative zero replacement
#' result_mzero <- gc_handle_zeros(soil_samples, method = "mzero")
#' print(result_mzero$imputed_data)
#' print(paste("Imputation rate:", result_mzero$imputation_rate))
#'
#' # Impute using log-ratio EM (more principled)
#' result_lrem <- gc_handle_zeros(soil_samples, method = "lrem")
#' print(result_lrem$imputed_data)
#'
#' # With detection limits (censored measurements)
#' detection_limits <- c(sand = 1, silt = 1, clay = 1)
#' result_dl <- gc_handle_zeros(soil_samples, method = "lrem", dl = detection_limits)
#' }
#'
#' @references
#' Aitchison, J. (1986). The Statistical Analysis of Compositional Data.
#' Chapman and Hall, London.
#'
#' Thió-Henestrosa, S., & Martín-Fernández, J. A. (2003).
#' Dealing with Compositional Data: The Fcompositions Package.
#' *Computational Statistics & Data Analysis*, 43, 523-536.
#'
#' @importFrom methods is
#' @export
gc_handle_zeros <- function(comp_data,
                            method = "mzero",
                            dl = NULL,
                            tolerance = 1e-10) {
  # Validate input
  if (!is.data.frame(comp_data)) {
    stop("comp_data must be a data frame")
  }

  if (nrow(comp_data) == 0) {
    stop("comp_data must have at least one row")
  }

  if (ncol(comp_data) < 2) {
    stop("comp_data must have at least 2 columns (compositional parts)")
  }

  comp_data <- as.data.frame(lapply(comp_data, as.numeric))

  has_zero <- rowSums(comp_data < tolerance | is.na(comp_data)) > 0
  has_negative <- rowSums(comp_data < 0, na.rm = TRUE) > 0
  problematic_rows <- has_zero | has_negative

  n_zeros_total <- sum(comp_data < tolerance | is.na(comp_data))
  n_cols <- ncol(comp_data)
  n_rows <- nrow(comp_data)
  total_values <- n_cols * n_rows
  imputation_rate <- n_zeros_total / total_values

  row_status <- factor(
    rep("observed", nrow(comp_data)),
    levels = c("observed", "imputed", "failed")
  )
  row_status[problematic_rows] <- "imputed"

  if (!requireNamespace("zCompositions", quietly = TRUE)) {
    stop(
      "zCompositions package is required for zero imputation. ",
      "Install with: install.packages('zCompositions')"
    )
  }

  method <- match.arg(method, c("mzero", "azero", "lrem"))

  tryCatch({
    if (method == "mzero") {
      # Multiplicative zero replacement
      imputed_data <- zCompositions::multLRM(comp_data, label = NA)
    } else if (method == "azero") {
      # Additive zero replacement
      imputed_data <- zCompositions::addLRM(comp_data, label = NA)
    } else if (method == "lrem") {
      # Log-ratio EM imputation
      # Build detection limit matrix if provided
      if (!is.null(dl)) {
        if (length(dl) != n_cols) {
          stop("dl must have same length as number of columns in comp_data")
        }
        # Create matrix of detection limits (one per column, repeated for all rows)
        dl_matrix <- matrix(
          rep(dl, n_rows),
          nrow = n_rows,
          byrow = TRUE
        )
        colnames(dl_matrix) <- colnames(comp_data)
        imputed_data <- zCompositions::lrEM(comp_data,
          dl = dl_matrix,
          label = NA
        )
      } else {
        # No detection limits provided, treat zeros as missing only
        imputed_data <- zCompositions::lrEM(comp_data, label = NA)
      }
    }

    imputed_data <- as.data.frame(imputed_data)
    colnames(imputed_data) <- colnames(comp_data)

    if (any(is.na(imputed_data)) || any(imputed_data < 0)) {
      row_status[problematic_rows] <- "failed"
      warning(
        "Imputation produced invalid values (NAs or negative). ",
        "Review imputed_data and consider using a different method."
      )
    }

    return(list(
      imputed_data = imputed_data,
      n_zeros_imputed = n_zeros_total,
      imputation_rate = imputation_rate,
      method_used = method,
      row_status = row_status
    ))
  }, error = function(e) {
    warning(
      "Imputation failed with method '", method, "': ",
      e$message, ". Returning original data."
    )
    row_status[problematic_rows] <- "failed"

    return(list(
      imputed_data = comp_data,
      n_zeros_imputed = n_zeros_total,
      imputation_rate = imputation_rate,
      method_used = paste0(method, " (failed)"),
      row_status = row_status
    ))
  })
}
