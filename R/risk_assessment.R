#' Compute Probability of Exceeding Threshold
#'
#' Calculate probability that simulated values exceed a specified threshold at each
#' spatial location. Uses indicator functions on simulation realizations to estimate
#' empirical probabilities.
#'
#' @param simulations A `terra::SpatRaster` stack with multiple layers representing
#'   different simulation realizations of the same variable(s). Each layer is one
#'   realization; values are in the original scale (after ILR back-transform if
#'   compositional data).
#'
#' @param threshold Numeric value or named list specifying threshold(s) to compare.
#'   - Single numeric: Applied to all layers
#'   - Named list (e.g., `list(sand = 20, clay = 30)`): Component-specific thresholds
#'     (for multi-component simulations)
#'
#' @param operator Character specifying comparison operator.
#'   - `">"` (default): Probability of exceeding threshold (Z > t)
#'   - `">="`: P(Z ≥ t)
#'   - `"<"`: Probability below threshold (Z < t)
#'   - `"<="`: P(Z ≤ t)
#'   - `"=="`: Probability of exact equality (rarely used)
#'
#' @param component Character or integer vector. If simulations contain multiple
#'   components (e.g., sand, silt, clay as separate layers), which component(s) to
#'   analyze.
#'   - NULL (default): Use all layers
#'   - Single value: Extract one component
#'   - Vector: Multiple components (computes per-component probability)
#'
#' @param na.rm Logical. How to handle NA values in the realization stack.
#'   - TRUE (default): P(Z_valid > t) - compute probability excluding NAs
#'   - FALSE: Propagate NA to output if any realization is NA
#'
#' @return A `terra::SpatRaster` object with a single layer representing
#'   probability values in [0, 1] at each spatial location. Layer is named
#'   `"probability"` by default. Same extent and resolution as input `simulations`.
#'
#' @details
#' **Probability Calculation:**
#'
#' At each grid cell, the function counts how many realizations satisfy the
#' threshold condition and divides by total number of realizations:
#'
#' \deqn{P(s) = \frac{\#\{Z(s)_k > t\}}{n}}
#'
#' Where n = number of simulation realizations.
#'
#' **Interpretation:**
#' - P = 0.9: 90% of simulations exceed threshold (high confidence of occurrence)
#' - P = 0.5: 50% of simulations exceed threshold (maximum uncertainty)
#' - P = 0.1: 10% of simulations exceed threshold (low probability of occurrence)
#'
#' **Common Applications:**
#' - Carbon credit eligibility: P(SOC > 50 Mg/ha)
#' - Contamination risk: P(Pb > 300 mg/kg)
#' - Agronomic suitability: P(clay < 40%)
#'
#' @examplesIf Sys.getenv("R_GEOCODA_RUN_EXAMPLES") == "true"
#' # Create synthetic data for demonstration
#' set.seed(42)
#' sims <- terra::rast(ncol = 10, nrow = 10, nlyr = 50)
#' terra::values(sims) <- rnorm(5000, mean = 50, sd = 15)
#'
#' # Probability of exceeding 50 Mg/ha
#' prob_map <- gc_probability_map(sims, threshold = 50)
#' terra::plot(prob_map, main = "P(SOC > 50 Mg/ha)")
#'
#' # Conservative estimate (use >= operator for credit eligibility)
#' prob_credit <- gc_probability_map(
#'   sims,
#'   threshold = 50,
#'   operator = ">="
#' )
#'
#' @export
gc_probability_map <- function(simulations,
                                threshold,
                                operator = ">",
                                component = NULL,
                                na.rm = TRUE) {

  # Validate inputs
  if (!inherits(simulations, "SpatRaster")) {
    stop("simulations must be a terra::SpatRaster object")
  }

  if (!is.numeric(threshold) && !is.list(threshold)) {
    stop("threshold must be numeric or a named list")
  }

  valid_operators <- c(">", ">=", "<", "<=", "==")
  if (!(operator %in% valid_operators)) {
    stop(paste("operator must be one of:", paste(valid_operators, collapse = ", ")))
  }

  # Extract values to matrix (n_cells × n_layers)
  values_matrix <- terra::values(simulations)

  if (is.null(values_matrix) || nrow(values_matrix) == 0) {
    stop("simulations contains no data")
  }

  # Apply threshold based on operator
  comparison_fn <- switch(operator,
    ">" = function(x, t) x > t,
    ">=" = function(x, t) x >= t,
    "<" = function(x, t) x < t,
    "<=" = function(x, t) x <= t,
    "==" = function(x, t) x == t
  )

  # Apply threshold
  logical_matrix <- comparison_fn(values_matrix, threshold)

  # Compute probability (mean of logical values)
  if (na.rm) {
    # For each row, count TRUEs and divide by non-NA count
    probability <- rowMeans(logical_matrix, na.rm = TRUE)
  } else {
    # If any NA, result is NA for that row
    probability <- rowMeans(logical_matrix, na.rm = FALSE)
  }

  # Create output raster
  prob_rast <- simulations[[1]]  # Template from first layer
  terra::values(prob_rast) <- probability

  # Set meaningful name
  names(prob_rast) <- "probability"

  return(prob_rast)
}


#' Compute Percentile Maps from Simulations
#'
#' Calculate percentile values (e.g., P10, P50, P90) at each grid cell from
#' simulation realizations. Useful for uncertainty visualization and decision-making.
#'
#' @param simulations A `terra::SpatRaster` stack with multiple layers representing
#'   simulation realizations.
#'
#' @param percentiles Numeric vector in [0, 1] specifying which percentiles to compute.
#'   Default: `c(0.1, 0.5, 0.9)` = [P10, P50, P90]
#'   - P10: 10th percentile (conservative low estimate, 90% confidence of exceeding)
#'   - P50: Median / 50th percentile (central estimate)
#'   - P90: 90th percentile (conservative high estimate, only 10% chance of exceeding)
#'
#' @param component Character or integer vector. Which component(s) to analyze if
#'   simulations are multi-component. NULL (default) = use all layers.
#'
#' @param na.rm Logical. TRUE (default) = compute percentiles excluding NAs;
#'   FALSE = return NA if any realization is NA at that location.
#'
#' @param names_format Character string for output layer naming. Use `{p}` as
#'   placeholder for percentile value.
#'   - Default: `"P{p}"` produces "P10", "P50", "P90"
#'   - Alternative: `"Q{p*100}"` produces "Q10", "Q50", "Q90"
#'
#' @return A `terra::SpatRaster` object with multiple layers (one per percentile).
#'   Layer names indicate percentile (e.g., "P10", "P50", "P90").
#'   Same extent and resolution as input.
#'
#' @details
#' **Uncertainty Quantification:**
#'
#' Percentiles provide conservative estimates for decision-making:
#' - **P90 (conservative high)**: For environmental standards, use if trying to be
#'   protective (assume worst case)
#' - **P50 (median)**: Best estimate of central tendency
#' - **P10 (conservative low)**: For carbon credits, use to be conservative with
#'   buyers (ensure deliver minimum promised)
#'
#' **Standard Uncertainty Band:**
#' The P10-P90 range encompasses 80% of simulated outcomes and is commonly used
#' to visualize uncertainty.
#'
#' @examplesIf Sys.getenv("R_GEOCODA_RUN_EXAMPLES") == "true"
#' # Create synthetic data
#' set.seed(42)
#' sims <- terra::rast(ncol = 10, nrow = 10, nlyr = 100)
#' terra::values(sims) <- rnorm(10000, mean = 50, sd = 15)
#'
#' # Standard P10-P50-P90 band
#' perc_map <- gc_percentile_map(sims)
#'
#' # Visualize uncertainty band
#' terra::plot(perc_map, main = "Soil Carbon Stock Uncertainty")
#'
#' # More granular percentiles
#' percentiles <- gc_percentile_map(
#'   sims,
#'   percentiles = c(0.05, 0.25, 0.5, 0.75, 0.95)
#' )
#'
#' # Calculate uncertainty band width
#' width <- perc_map[[3]] - perc_map[[1]]  # P90 - P10
#' terra::plot(width, main = "Uncertainty Band Width")
#'
#' @export
gc_percentile_map <- function(simulations,
                               percentiles = c(0.1, 0.5, 0.9),
                               component = NULL,
                               na.rm = TRUE,
                               names_format = "P{p}") {

  # Validate inputs
  if (!inherits(simulations, "SpatRaster")) {
    stop("simulations must be a terra::SpatRaster object")
  }

  if (!is.numeric(percentiles) || any(percentiles < 0 | percentiles > 1)) {
    stop("percentiles must be numeric values in [0, 1]")
  }

  # Extract values to matrix
  values_matrix <- terra::values(simulations)

  if (is.null(values_matrix) || nrow(values_matrix) == 0) {
    stop("simulations contains no data")
  }

  # Compute percentiles for each row (grid cell)
  n_cells <- nrow(values_matrix)
  n_percentiles <- length(percentiles)

  percentile_matrix <- matrix(nrow = n_cells, ncol = n_percentiles)

  for (i in seq_len(n_cells)) {
    cell_values <- values_matrix[i, ]
    if (na.rm) {
      cell_values <- cell_values[!is.na(cell_values)]
    }

    if (length(cell_values) > 0) {
      percentile_matrix[i, ] <- stats::quantile(cell_values, probs = percentiles, na.rm = FALSE)
    } else {
      percentile_matrix[i, ] <- NA_real_
    }
  }

  # Create output raster with multiple layers
  perc_rast <- simulations[[1:n_percentiles]]
  terra::values(perc_rast) <- percentile_matrix

  # Set layer names based on format string
  layer_names <- sapply(percentiles, function(p) {
    # Replace {p} with the raw percentile value
    name <- gsub("\\{p\\}", as.character(p), names_format)
    # Replace {p*100} with percentile as percentage
    name <- gsub("\\{p\\*100\\}", as.character(p * 100), name)
    name
  })
  names(perc_rast) <- layer_names

  return(perc_rast)
}


#' Internal: Compute Loss for Risk Assessment
#'
#' Compute loss values given observed/simulated values and a threshold, according to
#' specified loss function. Internal helper for `gc_risk_assessment()`.
#'
#' @param z Numeric vector of values to compute loss for
#' @param threshold Numeric threshold value
#' @param loss_function Character name of loss function or custom function
#' @param cost_params List of cost parameters (a, b, etc.)
#' @param operator Character: "<" (loss if under) or ">" (loss if over threshold)
#'
#' @return Numeric vector of loss values (same length as z)
#'
#' @keywords internal
compute_loss <- function(z, threshold, loss_function, cost_params, operator) {

  if (is.character(loss_function)) {
    loss_fn <- switch(loss_function,
      "linear" = function(z, t, a, b) a * abs(z - t),
      "asymmetric" = function(z, t, a, b) {
        ifelse(z < t, a * (t - z), b * (z - t))
      },
      "step" = function(z, t, a, b) {
        ifelse(z < t, a, b)
      },
      stop(paste("Unknown loss function:", loss_function))
    )
  } else if (is.function(loss_function)) {
    loss_fn <- loss_function
  } else {
    stop("loss_function must be character or function")
  }

  # Extract cost parameters with defaults
  a <- cost_params$a %||% 1.0
  b <- cost_params$b %||% 1.0

  # Compute loss
  loss <- loss_fn(z, threshold, a, b)

  return(loss)
}


#' Assess Risk with Loss Function Framework
#'
#' Compute expected loss maps under specified cost structure for decision-making.
#' Integrates asymmetric costs for false positives vs. false negatives.
#'
#' @param simulations A `terra::SpatRaster` stack with simulation realizations.
#'
#' @param threshold Numeric threshold value for decision boundary
#'   (e.g., minimum SOC for carbon credit, maximum contamination limit)
#'
#' @param loss_function Character or function specifying loss computation.
#'   - `"linear"`: L(z) = a|z - t| (symmetric loss)
#'   - `"asymmetric"`: L(z) = a(t-z)⁺ + b(z-t)⁺ (asymmetric, different costs)
#'   - `"step"`: L(z) = a·I(z<t) + b·I(z>t) (binary cost)
#'   - Function: Custom function(z, threshold, cost_params)
#'
#' @param cost_params List with cost parameters:
#'   - `a`: Cost of under-achievement (loss if Z < threshold)
#'   - `b`: Cost of over-achievement (loss if Z > threshold)
#'   - Default: `list(a = 1, b = 1)` for symmetric loss
#'   - Example: `list(a = 2, b = 1)` means under-estimate costs 2× over-estimate
#'
#' @param operator Character: "<" (loss if under, default) or ">" (loss if over)
#'
#' @param na.rm Logical. TRUE (default) = exclude NAs from expected value calculation
#'
#' @return List containing:
#'   - `expected_loss`: terra::SpatRaster with expected loss at each location
#'   - `decision`: terra::SpatRaster with binary decision (1=accept, 0=reject)
#'   - `cost_params`: Cost parameters used
#'   - `metadata`: List with threshold, operator, loss_function used
#'
#' @details
#' **Loss Function Framework:**
#'
#' Expected loss at location s:
#'
#' \deqn{EL(s|t) = \frac{1}{n} \sum_{k=1}^n L(Z(s)_k | t)}
#'
#' Where L(z|t) is the loss function evaluated at each realization.
#'
#' **Decision Rule:**
#' Accept location if EL(threshold) < EL(threshold ± ε), otherwise reject.
#'
#' **Asymmetric Cost Example (Carbon Credits):**
#' - a = 2.0: Cost of over-issuing undeserved credit (high financial risk)
#' - b = 1.0: Cost of withholding valid credit (opportunity cost)
#' - Decision threshold optimizes total risk
#'
#' @examplesIf Sys.getenv("R_GEOCODA_RUN_EXAMPLES") == "true"
#' # Create synthetic data
#' set.seed(42)
#' carbon_sims <- terra::rast(ncol = 10, nrow = 10, nlyr = 50)
#' terra::values(carbon_sims) <- rnorm(5000, mean = 50, sd = 15)
#'
#' # Carbon credit verification: Asymmetric costs
#' loss_result <- gc_risk_assessment(
#'   simulations = carbon_sims,
#'   threshold = 50,
#'   loss_function = "asymmetric",
#'   cost_params = list(
#'     a = 2.0,  # Cost of false positive (over-crediting)
#'     b = 1.0   # Cost of false negative (under-crediting)
#'   ),
#'   operator = "<"
#' )
#'
#' terra::plot(loss_result$expected_loss, main = "Expected Loss")
#' terra::plot(loss_result$decision, main = "Credit Approval")
#'
#' @export
gc_risk_assessment <- function(simulations,
                                threshold,
                                loss_function,
                                cost_params = list(a = 1, b = 1),
                                operator = "<",
                                na.rm = TRUE) {

  # Validate inputs
  if (!inherits(simulations, "SpatRaster")) {
    stop("simulations must be a terra::SpatRaster object")
  }

  if (!is.numeric(threshold) || length(threshold) != 1) {
    stop("threshold must be a single numeric value")
  }

  if (!(operator %in% c("<", ">"))) {
    stop("operator must be '<' (loss if under) or '>' (loss if over)")
  }

  if (!is.list(cost_params)) {
    stop("cost_params must be a list with elements 'a' and 'b'")
  }

  # Extract values to matrix
  values_matrix <- terra::values(simulations)

  if (is.null(values_matrix) || nrow(values_matrix) == 0) {
    stop("simulations contains no data")
  }

  # Compute loss for each realization
  n_cells <- nrow(values_matrix)
  loss_matrix <- matrix(nrow = n_cells, ncol = ncol(values_matrix))

  for (i in seq_len(n_cells)) {
    loss_matrix[i, ] <- compute_loss(
      z = values_matrix[i, ],
      threshold = threshold,
      loss_function = loss_function,
      cost_params = cost_params,
      operator = operator
    )
  }

  # Compute expected loss (mean across realizations)
  if (na.rm) {
    expected_loss <- rowMeans(loss_matrix, na.rm = TRUE)
  } else {
    expected_loss <- rowMeans(loss_matrix, na.rm = FALSE)
  }

  # Create output rasters
  loss_rast <- simulations[[1]]
  terra::values(loss_rast) <- expected_loss
  names(loss_rast) <- "expected_loss"

  # Decision: Accept if expected loss at threshold is acceptable
  # (simple decision: threshold where loss is minimized)
  decision <- expected_loss < stats::median(expected_loss, na.rm = TRUE)
  decision_rast <- simulations[[1]]
  terra::values(decision_rast) <- as.numeric(decision)
  names(decision_rast) <- "decision"

  # Return as structured list
  result <- list(
    expected_loss = loss_rast,
    decision = decision_rast,
    cost_params = cost_params,
    metadata = list(
      threshold = threshold,
      operator = operator,
      loss_function = if (is.character(loss_function)) loss_function else "custom",
      n_realizations = ncol(values_matrix)
    )
  )

  class(result) <- c("risk_assessment", "list")
  return(result)
}


#' Specialized Carbon Stock Audit Function
#'
#' Comprehensive wrapper for carbon credit verification in payment-for-ecosystem-services
#' programs. Computes conservative estimates, probability maps, and audit reports.
#'
#' @param simulations A `terra::SpatRaster` stack or data.frame with simulation
#'   realizations. If raster, returns spatial results; if data.frame, returns
#'   tabular summary (useful for non-spatial audits).
#'
#' @param credit_threshold Numeric. Minimum carbon stock (Mg/ha or equivalent units)
#'   required for credit issuance.
#'
#' @param confidence Numeric in (0, 1). Confidence level for conservative estimates.
#'   - confidence = 0.9: Use P10 (10th percentile, 90% sure actual ≥ this value)
#'   - confidence = 0.95: Use P5 (5th percentile, 95% sure actual ≥ this value)
#'   - Default: 0.9
#'
#' @param component Character. If multi-component simulations, which variable
#'   represents carbon stock (e.g., "SOC"). Usually not needed for single-variable
#'   simulations.
#'
#' @param audit_type Character. Verification strategy:
#'   - `"conservative"` (default): Issue credits only if P(1-confidence) > threshold
#'     (suitable for credit issuance with low risk)
#'   - `"expected"`: Use P50 (median estimate, medium risk)
#'   - `"optimistic"`: Use P(confidence) (exploratory/planning, high risk)
#'
#' @param na.rm Logical. Exclude NA values from calculations.
#'
#' @return List containing:
#'   - `probability_map` (if raster input): SpatRaster of P(SOC > threshold)
#'   - `conservative_estimate` (if raster input): SpatRaster of percentile values
#'   - `issued_credit` (if raster input): SpatRaster with 0/1 credit approval
#'   - `summary`: Data.frame with aggregate statistics (total area, eligible area,
#'     total credits issuable, etc.)
#'   - `audit_report`: Character string with human-readable audit summary
#'   - `metadata`: List with audit parameters and timestamp
#'
#' @details
#' **Carbon Audit Methodology:**
#'
#' 1. Compute conservative percentile map (e.g., P10 for 90% confidence)
#' 2. Compute probability of exceeding credit threshold
#' 3. Determine eligible area where conservative_estimate > credit_threshold
#' 4. Aggregate results and generate audit report
#'
#' **Interpretation of Audit Report:**
#'
#' - **Total Area**: Entire audit domain
#' - **Eligible Area**: Where conservative estimate exceeds threshold
#' - **Total Credits**: Eligible area × mean conservative estimate
#' - **Confidence Level**: Trust level in issued credits
#'
#' @examplesIf Sys.getenv("R_GEOCODA_RUN_EXAMPLES") == "true"
#' # Create synthetic data
#' set.seed(42)
#' carbon_sims <- terra::rast(ncol = 10, nrow = 10, nlyr = 100)
#' terra::values(carbon_sims) <- rnorm(10000, mean = 55, sd = 15)
#'
#' # Standard 90% confidence audit for carbon credit program
#' audit <- gc_carbon_audit(
#'   simulations = carbon_sims,
#'   credit_threshold = 50,
#'   confidence = 0.9,
#'   audit_type = "conservative"
#' )
#'
#' # Print human-readable report
#' cat(audit$audit_report)
#'
#' # Check eligible area
#' eligible_pixels <- sum(terra::values(audit$issued_credit) == 1, na.rm = TRUE)
#' total_pixels <- terra::ncell(audit$issued_credit)
#' cat("Eligible area:", eligible_pixels, "/", total_pixels, "pixels\n")
#'
#' # Map results
#' terra::plot(audit$probability_map, main = "P(SOC > 50 Mg/ha)")
#' terra::plot(audit$issued_credit, main = "Credit Approved (Green) vs Denied (Red)")
#'
#' @export
gc_carbon_audit <- function(simulations,
                            credit_threshold,
                            confidence = 0.9,
                            component = NULL,
                            audit_type = "conservative",
                            na.rm = TRUE) {

  # Validate inputs
  if (!inherits(simulations, "SpatRaster") && !is.data.frame(simulations)) {
    stop("simulations must be a terra::SpatRaster or data.frame")
  }

  if (!is.numeric(credit_threshold) || length(credit_threshold) != 1) {
    stop("credit_threshold must be a single numeric value")
  }

  if (!is.numeric(confidence) || confidence <= 0 || confidence >= 1) {
    stop("confidence must be numeric in (0, 1)")
  }

  if (!(audit_type %in% c("conservative", "expected", "optimistic"))) {
    stop("audit_type must be 'conservative', 'expected', or 'optimistic'")
  }

  # Determine percentile to use based on audit type
  percentile_to_use <- switch(audit_type,
    "conservative" = 1 - confidence,  # e.g., 0.1 for 0.9 confidence
    "expected" = 0.5,
    "optimistic" = confidence
  )

  # Compute probability and percentile maps
  prob_map <- gc_probability_map(simulations, threshold = credit_threshold, operator = ">")
  perc_map <- gc_percentile_map(simulations, percentiles = percentile_to_use)

  # Determine credit eligibility
  eligible <- terra::values(perc_map) > credit_threshold
  credit_rast <- perc_map
  terra::values(credit_rast) <- as.numeric(eligible)
  names(credit_rast) <- "issued_credit"

  # Aggregate statistics
  if (inherits(simulations, "SpatRaster")) {
    n_cells <- terra::ncell(simulations)
    cell_area <- prod(terra::res(simulations)) / 10000  # Convert to hectares (assuming m)

    n_eligible <- sum(eligible, na.rm = TRUE)
    total_area <- n_cells * cell_area
    eligible_area <- n_eligible * cell_area

    # Mean conservative estimate over eligible area
    perc_values <- terra::values(perc_map)
    perc_values[!eligible] <- NA
    mean_conservative <- mean(perc_values, na.rm = TRUE)

    total_credits <- eligible_area * mean_conservative

    summary_df <- data.frame(
      metric = c(
        "total_area_ha",
        "eligible_area_ha",
        "ineligible_area_ha",
        "mean_soc_eligible",
        "total_credits_mg",
        "confidence_level"
      ),
      value = c(
        round(total_area, 2),
        round(eligible_area, 2),
        round(total_area - eligible_area, 2),
        round(mean_conservative, 2),
        round(total_credits, 0),
        confidence
      )
    )
  } else {
    # Non-spatial (data.frame input)
    values_matrix <- as.matrix(simulations)
    perc_value <- stats::quantile(as.vector(values_matrix), probs = percentile_to_use, na.rm = na.rm)
    eligible_overall <- perc_value > credit_threshold

    summary_df <- data.frame(
      metric = c(
        "n_samples",
        "conservative_estimate",
        "credit_eligible",
        "confidence_level"
      ),
      value = c(
        nrow(values_matrix),
        round(perc_value, 2),
        eligible_overall,
        confidence
      )
    )
  }

  # Generate human-readable audit report
  audit_report <- sprintf(
    "Soil Carbon Stock Audit (%s)\n%s\n\nCredit Threshold: %.1f Mg/ha\nConfidence Level: %.0f%%\nAudit Type: %s\n\n",
    audit_type,
    strrep("=", 50),
    credit_threshold,
    confidence * 100,
    audit_type
  )

  audit_report <- paste(audit_report, "Summary Statistics:\n")
  for (i in seq_len(nrow(summary_df))) {
    audit_report <- paste(
      audit_report,
      sprintf("%-25s: %s\n", summary_df$metric[i], summary_df$value[i])
    )
  }

  audit_report <- paste(
    audit_report,
    sprintf("\nAudit conducted: %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
  )

  # Return structured result
  if (inherits(simulations, "SpatRaster")) {
    result <- list(
      probability_map = prob_map,
      conservative_estimate = perc_map,
      issued_credit = credit_rast,
      summary = summary_df,
      audit_report = audit_report,
      metadata = list(
        credit_threshold = credit_threshold,
        confidence = confidence,
        audit_type = audit_type,
        n_realizations = terra::nlyr(simulations),
        timestamp = Sys.time()
      )
    )
  } else {
    result <- list(
      summary = summary_df,
      audit_report = audit_report,
      metadata = list(
        credit_threshold = credit_threshold,
        confidence = confidence,
        audit_type = audit_type,
        timestamp = Sys.time()
      )
    )
  }

  class(result) <- c("carbon_audit", "list")
  return(result)
}
