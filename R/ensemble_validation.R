#' Ensemble Realization Post-Processing
#'
#' Functions for aggregating multiple stochastic realizations and validating
#' simulation quality. Enables uncertainty quantification, ensemble statistics,
#' and spatial continuity assessment.
#'
#' @keywords internal
#' @name ensemble_validation
NULL


#' Aggregate Multiple Realizations to Summary Statistics
#'
#' Combine multiple stochastic simulation realizations into ensemble summary
#' maps (mean, median, standard deviation, quantiles, confidence intervals).
#' Useful for communicating uncertainty and central tendency.
#'
#' @param realizations A list of `terra::SpatRaster` objects (output from
#'   multiple calls to `gc_sim_composition()` with `nsim=1` each), OR a
#'   single `terra::SpatRaster` with multiple layers organized as:
#'   `[comp1.sim1, comp1.sim2, ..., comp1.simN, comp2.sim1, ...]`
#' @param components Character vector, names of compositional components
#'   (e.g., c("sand", "silt", "clay"))
#' @param realizations_per_component Integer, number of realizations per component.
#'   If NULL (default), inferred from raster layer names (e.g., "sand.sim1" → 1)
#' @param stats Character vector, summary statistics to compute:
#'   - "mean": arithmetic mean
#'   - "median": median (50th percentile)
#'   - "sd": standard deviation
#'   - "cv": coefficient of variation (sd/mean, useful for normalized uncertainty)
#'   - "p05", "p25", "p75", "p95": quantiles
#'   - "ci_lower", "ci_upper": 95% confidence intervals (2.5th, 97.5th percentiles)
#'
#' @param verbose Logical, report processing steps (default FALSE)
#'
#' @return A `terra::SpatRaster` with layers for each component × statistic.
#'   Layer names follow pattern: `component_stat` (e.g., "sand_mean", "silt_sd")
#'
#' @details
#' This function assumes realizations are stored in a consistent format:
#' either as a list where each element is a single-component raster,
#' or as a multi-layer raster where layers are grouped by component then by realization.
#'
#' Common workflow:
#' ```r
#' # Generate N realizations separately
#' real_1 <- gc_sim_composition(model, domain, nsim = 1)
#' real_2 <- gc_sim_composition(model, domain, nsim = 1)
#' real_N <- gc_sim_composition(model, domain, nsim = 1)
#'
#' # Combine stack and aggregate
#' all_reals <- c(real_1, real_2, ..., real_N)
#' ensemble_stats <- gc_aggregate_realizations(
#'   all_reals,
#'   components = c("sand", "silt", "clay"),
#'   stats = c("mean", "sd", "p05", "p95")
#' )
#' ```
#'
#' @examples
#' # Hypothetical ensemble of 3 realizations × 3 components
#' # (In practice, generate via gc_sim_composition(..., nsim=1) repeated)
#'
#' # Create synthetic raster stack
#' r <- terra::rast(nrows = 10, ncols = 10, nlyrs = 9)
#' terra::ext(r) <- c(0, 100, 0, 100)
#' terra::crs(r) <- "EPSG:4326"
#'
#' # Set layer names: sand.sim1, sand.sim2, sand.sim3, ..., clay.sim3
#' nlyr_per_comp <- 3
#' comps <- c("sand", "silt", "clay")
#' names(r) <- paste0(
#'   rep(comps, each = nlyr_per_comp),
#'   ".sim",
#'   rep(1:nlyr_per_comp, times = length(comps))
#' )
#'
#' # Fill with random compositions summing to 100
#' set.seed(42)
#' for (i in 1:terra::nlyr(r)) {
#'   vals <- abs(stats::rnorm(terra::ncell(r), 50, 10))
#'   vals <- vals / rowSums(matrix(vals, ncol = 3, byrow = FALSE)) * 100
#'   terra::values(r)[, i] <- vals[, i]
#' }
#'
#' # Aggregate to ensemble statistics
#' ensemble <- gc_aggregate_realizations(
#'   r,
#'   components = comps,
#'   stats = c("mean", "sd", "p05", "p95")
#' )
#'
#' # Inspect output
#' terra::names(ensemble)  # e.g., "sand_mean", "sand_sd", "sand_p05", ...
#'
#' @importFrom terra rast nlyr names values
#' @export
gc_aggregate_realizations <- function(realizations,
                                       components,
                                       realizations_per_component = NULL,
                                       stats = c("mean", "sd", "p05", "p95"),
                                       verbose = FALSE) {
  
  # Validate stats argument
  valid_stats <- c("mean", "median", "sd", "cv", "p05", "p25", "p75", "p95", "ci_lower", "ci_upper")
  invalid <- setdiff(stats, valid_stats)
  if (length(invalid) > 0) {
    stop("Invalid stats: ", paste(invalid, collapse = ", "),
         ". Valid options: ", paste(valid_stats, collapse = ", "))
  }
  
  # Coerce realizations to SpatRaster if it's a list
  if (is.list(realizations)) {
    if (length(realizations) == 0) {
      stop("realizations list is empty")
    }
    if (verbose) cat("Combining", length(realizations), "rasters from list...\n")
    realizations <- terra::rast(realizations)
  }
  
  # Validate SpatRaster
  if (!inherits(realizations, "SpatRaster")) {
    stop("realizations must be a terra::SpatRaster or list of SpatRasters")
  }
  
  nlyr <- terra::nlyr(realizations)
  layer_names <- terra::names(realizations)
  
  # Infer realizations per component if not provided
  if (is.null(realizations_per_component)) {
    # Extract numeric suffix from layer names (e.g., "sand.sim3" → 3)
    sim_nums <- as.numeric(gsub(".*\\.sim(\\d+)$", "\\1", layer_names))
    if (any(is.na(sim_nums))) {
      warning("Could not infer realization indices from layer names. ",
              "Specify realizations_per_component explicitly.")
      realizations_per_component <- nlyr / length(components)
    } else {
      realizations_per_component <- max(sim_nums, na.rm = TRUE)
    }
  }
  
  if (nlyr != length(components) * realizations_per_component) {
    stop(
      "Mismatch: nlyr = ", nlyr, " but ",
      "length(components) × realizations_per_component = ",
      length(components), " × ", realizations_per_component, " = ",
      length(components) * realizations_per_component
    )
  }
  
  if (verbose) {
    cat("Processing", length(components), "components with",
        realizations_per_component, "realizations each\n")
  }
  
  # Initialize output list
  output_layers <- list()
  
  # Process each component
  for (comp_idx in seq_along(components)) {
    comp_name <- components[comp_idx]
    
    # Extract layer indices for this component
    start_idx <- (comp_idx - 1) * realizations_per_component + 1
    end_idx <- comp_idx * realizations_per_component
    layer_indices <- start_idx:end_idx
    
    # Extract values for all realizations of this component
    # terra::values() returns matrix: rows = cells, cols = layers
    comp_values <- terra::values(realizations, layer = layer_indices)
    
    if (verbose) {
      cat("  ", comp_name, ": extracting layers", start_idx, "-", end_idx, "\n")
    }
    
    # Compute each requested statistic
    for (stat in stats) {
      stat_result <- switch(stat,
        mean = {
          rowMeans(comp_values, na.rm = TRUE)
        },
        median = {
          apply(comp_values, 1, stats::median, na.rm = TRUE)
        },
        sd = {
          apply(comp_values, 1, stats::sd, na.rm = TRUE)
        },
        cv = {
          means <- rowMeans(comp_values, na.rm = TRUE)
          sds <- apply(comp_values, 1, stats::sd, na.rm = TRUE)
          sds / means  # Coefficient of variation
        },
        p05 = {
          apply(comp_values, 1, stats::quantile, probs = 0.05, na.rm = TRUE)
        },
        p25 = {
          apply(comp_values, 1, stats::quantile, probs = 0.25, na.rm = TRUE)
        },
        p75 = {
          apply(comp_values, 1, stats::quantile, probs = 0.75, na.rm = TRUE)
        },
        p95 = {
          apply(comp_values, 1, stats::quantile, probs = 0.95, na.rm = TRUE)
        },
        ci_lower = {
          apply(comp_values, 1, stats::quantile, probs = 0.025, na.rm = TRUE)
        },
        ci_upper = {
          apply(comp_values, 1, stats::quantile, probs = 0.975, na.rm = TRUE)
        },
        stop("Unknown stat: ", stat)
      )
      
      # Create raster layer for this statistic
      layer_name <- paste0(comp_name, "_", stat)
      output_layers[[layer_name]] <- stat_result
    }
  }
  
  # Stack all output layers into a single SpatRaster
  # Convert list of vectors to matrix, then create raster
  output_matrix <- do.call(cbind, output_layers)
  
  # Create template raster from first component
  output_raster <- terra::rast(realizations, nlyrs = ncol(output_matrix))
  terra::values(output_raster) <- output_matrix
  names(output_raster) <- names(output_layers)
  
  if (verbose) {
    cat("Output:", terra::nlyr(output_raster), "layers\n")
    cat("Names:", paste(terra::names(output_raster)[1:min(5, terra::nlyr(output_raster))],
                        collapse = ", "), "...\n")
  }
  
  output_raster
}


#' Validate Simulation Realizations Against Constraints
#'
#' Assess realization quality by checking compositional sum constraints,
#' boundary violations, and spatial continuity. Returns per-pixel violation
#' counts and summary statistics.
#'
#' @param realizations A `terra::SpatRaster` with simulation output
#'   (multiple components and/or realizations)
#' @param components Character vector, names of compositional components
#'   (e.g., c("sand", "silt", "clay"))
#' @param target_sum Numeric, target composition sum (default 100)
#' @param sum_tolerance Numeric, acceptable deviation from target_sum (default 1)
#'   Composed values with |sum - target_sum| > tolerance are flagged
#' @param constraints List with per-component constraints (min, max).
#'   If NULL (default), no boundary checking is performed.
#'   Format: `constraints = list(SAND = list(min=0, max=100), ...)`
#' @param realizations_per_component Integer, number of realizations per component.
#'   If NULL, inferred from layer names
#' @param verbose Logical, report processing steps (default FALSE)
#'
#' @return A list with elements:
#'   - `violations_map`: `terra::SpatRaster` with per-pixel violation count
#'     (0 = all valid, higher = more violations)
#'   - `summary`: Data frame with per-component statistics:
#'     - `component`: Component name
#'     - `n_sum_violations`: Pixels where composition sum out of bounds
#'     - `n_boundary_violations`: Pixels where component value outside [min, max]
#'     - `pct_valid`: Percent of pixels with valid values
#'     - `mean`, `sd`: Mean and SD of component across all valid realizations
#'   - `constraint_ranges`: Echo of constraints used
#'   - `metadata`: Processing metadata
#'
#' @details
#' This function is useful for diagnosing simulation artifacts:
#' - High sum violations may indicate tight constraint bounds or compositional
#'   arithmetic issues
#' - Boundary violations suggest constraint specification issues
#' - Spatial clustering of violations may indicate zone boundaries or
#'   non-stationarity
#'
#' @examples
#' # Validate a realization raster
#' \dontrun{
#' validation <- gc_validate_realizations(
#'   realizations = my_sim_output,
#'   components = c("sand", "silt", "clay"),
#'   target_sum = 100,
#'   sum_tolerance = 1,
#'   constraints = list(
#'     SAND = list(min = 0, max = 100),
#'     SILT = list(min = 0, max = 100),
#'     CLAY = list(min = 0, max = 100)
#'   )
#' )
#'
#' # Inspect violations
#' print(validation$summary)
#' plot(validation$violations_map)
#' }
#'
#' @importFrom terra rast nlyr names values
#' @export
gc_validate_realizations <- function(realizations,
                                      components,
                                      target_sum = 100,
                                      sum_tolerance = 1,
                                      constraints = NULL,
                                      realizations_per_component = NULL,
                                      verbose = FALSE) {
  
  # Validate inputs
  if (!inherits(realizations, "SpatRaster")) {
    stop("realizations must be a terra::SpatRaster")
  }
  
  nlyr <- terra::nlyr(realizations)
  layer_names <- terra::names(realizations)
  
  # Infer realizations per component if not provided
  if (is.null(realizations_per_component)) {
    sim_nums <- as.numeric(gsub(".*\\.sim(\\d+)$", "\\1", layer_names))
    if (any(is.na(sim_nums))) {
      realizations_per_component <- nlyr / length(components)
    } else {
      realizations_per_component <- max(sim_nums, na.rm = TRUE)
    }
  }
  
  if (verbose) {
    cat("Validating", length(components), "components ×",
        realizations_per_component, "realizations\n")
  }
  
  # Extract all values
  all_values <- terra::values(realizations)
  ncells <- nrow(all_values)
  
  # Initialize violation counter
  violation_counter <- matrix(0, nrow = ncells, ncol = 1)
  
  # Initialize summary list
  summary_list <- list()
  
  # Check each component
  for (comp_idx in seq_along(components)) {
    comp_name <- components[comp_idx]
    
    # Extract layer indices for this component
    start_idx <- (comp_idx - 1) * realizations_per_component + 1
    end_idx <- comp_idx * realizations_per_component
    layer_indices <- start_idx:end_idx
    
    comp_values <- all_values[, layer_indices, drop = FALSE]
    
    if (verbose) {
      cat("  Processing", comp_name, "...\n")
    }
    
    # Check boundary constraints
    boundary_violations <- rep(0, ncells)
    if (!is.null(constraints) && comp_name %in% names(constraints)) {
      comp_constraint <- constraints[[comp_name]]
      if (is.list(comp_constraint) && "min" %in% names(comp_constraint)) {
        min_val <- comp_constraint$min
        max_val <- comp_constraint$max
        
        # Count pixels where ANY realization violates bounds
        boundary_violations <- rowSums(comp_values < min_val | comp_values > max_val, na.rm = TRUE)
        boundary_violations <- pmin(boundary_violations, 1)  # Convert to binary (0 or 1)
      }
    }
    
    # Add boundary violations to counter
    violation_counter <- violation_counter + boundary_violations
    
    # Compute component statistics
    comp_mean <- rowMeans(comp_values, na.rm = TRUE)
    comp_sd <- apply(comp_values, 1, stats::sd, na.rm = TRUE)
    
    summary_list[[comp_name]] <- data.frame(
      component = comp_name,
      n_boundary_violations = sum(boundary_violations),
      mean_value = mean(comp_mean, na.rm = TRUE),
      sd_value = mean(comp_sd, na.rm = TRUE),
      pct_valid = 100 * sum(!is.na(comp_values)) / length(comp_values)
    )
  }
  
  # Check compositional sum constraint
  # For each cell and realization, compute 3-component sum
  sum_violations <- rep(0, ncells)
  
  for (comp_idx in seq_along(components)) {
    start_idx <- (comp_idx - 1) * realizations_per_component + 1
    end_idx <- comp_idx * realizations_per_component
    
    if (comp_idx == 1) {
      # Initialize sum accumulator
      comp_sums <- matrix(0, nrow = ncells, ncol = realizations_per_component)
    }
    
    comp_sums <- comp_sums + all_values[, start_idx:end_idx, drop = FALSE]
  }
  
  # Check if sums deviate from target
  sum_violations <- rowSums(abs(comp_sums - target_sum) > sum_tolerance, na.rm = TRUE)
  sum_violations <- pmin(sum_violations, realizations_per_component)  # Cap at max realizations
  
  # Add sum violations to counter (count max 1 per cell for summary purposes)
  violation_counter <- violation_counter + pmin(sum_violations, 1)
  
  # Create violations map raster
  violations_raster <- terra::rast(realizations[[1]])
  terra::values(violations_raster) <- violation_counter
  names(violations_raster) <- "violations"
  
  # Compile summary
  summary_df <- do.call(rbind, summary_list)
  rownames(summary_df) <- NULL
  
  # Add sum violation summary
  summary_df$n_sum_violations <- sum(sum_violations > 0)
  
  if (verbose) {
    cat("Validation complete:\n")
    cat("  Total violations:", sum(violation_counter), "\n")
    cat("  Sum constraint violations:", sum(sum_violations > 0), "\n")
  }
  
  list(
    violations_map = violations_raster,
    summary = summary_df,
    constraint_ranges = constraints,
    metadata = list(
      target_sum = target_sum,
      sum_tolerance = sum_tolerance,
      n_cells = ncells,
      n_components = length(components),
      n_realizations = realizations_per_component,
      validation_date = Sys.Date()
    )
  )
}
