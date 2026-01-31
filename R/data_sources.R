#' Data Source Integration
#'
#' Functions for integrating external data sources into the geostatistical
#' simulation workflow. Supports SSURGO soil survey data, depth-stratified
#' parametric specifications, and component weighting schemes.
#'
#' @details
#' These functions accept soil data in a tabular format (data frame). The specific
#' source of that data (Soil Data Access API, NASIS, CSV, etc.) is flexible—what
#' matters is that columns match the expected schema (e.g., mukey, cokey, comppct,
#' sandtotal, silttotal, claytotal, hzdept_r, hzdepb_r).
#'
#' @keywords internal
#' @name data_sources
NULL


#' Convert Soil Component Data to Composition Parameters
#'
#' Transform soil component data (sand, silt, clay percentages) from any source
#' into constraint and samples for simulation. Aggregates multiple components per
#' map unit using specified weighting scheme.
#'
#' @param x Data frame with soil component data. Must contain:
#'   - `mukey` (character/numeric): map unit key, identifies soil mapping unit
#'   - `cokey` (character/numeric): component key, identifies soil component
#'   - `comppct` (numeric): percent of map unit occupied by component (0-100)
#'   - `sandtotal` (numeric): sand content (%)
#'   - `silttotal` (numeric): silt content (%)
#'   - `claytotal` (numeric): clay content (%)
#'   - `hzdept_r` (numeric, optional): depth to top of horizon (cm)
#'   - `hzdepb_r` (numeric, optional): depth to bottom of horizon (cm)
#'
#' @param weight_method Character, aggregation weight method for components:
#'   - "comppct" (default): weight by component percent (comppct)
#'   - "area": equal weight to all components in map unit
#'   - "custom": use external weight column
#' @param weight_col Character, name of weight column (only if weight_method = "custom")
#' @param depth_range Numeric vector c(min_depth, max_depth) in cm to filter horizons.
#'   If NULL (default), uses all horizons.
#' @param target_sum Numeric, target sum for compositions (default 100)
#' @param agg_method Character, aggregation method by map unit:
#'   - "weighted_mean" (default): weighted mean across components
#'   - "percentile": use specified percentile (e.g., "p50" for median)
#' @param percentile Numeric (0-100), percentile to use if agg_method = "percentile"
#' @param verbose Logical, report processing steps (default FALSE)
#'
#' @return List with elements:
#'   - `constraints`: Map unit-level constraints (min, representative, max)
#'   - `samples`: Point-level component samples (unweighted)
#'   - `samples_weighted`: Component samples with weights applied
#'   - `mu_stats`: Summary statistics by map unit
#'   - `metadata`: Processing metadata (source, method, depth_range, n_records)
#'   - `input_data`: Original input data for reference
#'
#' @details
#' The constraint bounds are derived from:
#' - Minimum: 10th percentile across all components
#' - Representative: Weighted mean of component percentages
#' - Maximum: 90th percentile across all components
#'
#' This approach preserves spatial and component heterogeneity while providing
#' realistic bounds for the simulation domain.
#'
#' @examples
#' # Soil component data from any source (SSURGO, SDA, NASIS, CSV, etc.)
#' soil_data <- data.frame(
#'   mukey = c(rep(1001, 3), rep(1002, 2)),
#'   cokey = c(1001.1, 1001.2, 1001.3, 1002.1, 1002.2),
#'   comppct = c(60, 25, 15, 50, 50),
#'   sandtotal = c(35, 15, 50, 40, 20),
#'   silttotal = c(50, 60, 30, 45, 65),
#'   claytotal = c(15, 25, 20, 15, 15),
#'   hzdept_r = c(0, 0, 0, 0, 0),
#'   hzdepb_r = c(25, 25, 25, 25, 25)
#' )
#'
#' # Convert to parameters
#' params <- gc_ssurgo_to_params(
#'   soil_data,
#'   weight_method = "comppct",
#'   depth_range = c(0, 25)
#' )
#'
#' # Extract constraints for simulation
#' constraints <- params$constraints
#'
#' @importFrom stats quantile weighted.mean aggregate
#' @export
gc_ssurgo_to_params <- function(x,
                                weight_method = "comppct",
                                weight_col = NULL,
                                depth_range = NULL,
                                target_sum = 100,
                                agg_method = "weighted_mean",
                                percentile = 50,
                                verbose = FALSE) {
  
  # Validate required columns
  required_cols <- c("mukey", "cokey", "comppct", "sandtotal", "silttotal", "claytotal")
  missing_cols <- setdiff(required_cols, names(x))
  
  if (length(missing_cols) > 0) {
    stop(
      "Missing required columns: ",
      paste(missing_cols, collapse = ", "),
      ". Required: mukey, cokey, comppct, sandtotal, silttotal, claytotal"
    )
  }
  
  data <- x
  
  # Filter by depth range if provided
  if (!is.null(depth_range)) {
    if (!"hzdept_r" %in% names(data) || !"hzdepb_r" %in% names(data)) {
      warning("depth_range specified but horizon depth columns not found. Skipping depth filtering.")
    } else {
      if (verbose) cat("Filtering to depth range:", paste(depth_range, collapse = "-"), "cm\n")
      
      # Keep horizons that overlap with depth_range
      # Overlap occurs when: horizon_top < range_bottom AND horizon_bottom > range_top
      data <- data[
        data$hzdept_r < depth_range[2] & data$hzdepb_r > depth_range[1],
        ,
        drop = FALSE
      ]
      
      if (nrow(data) == 0) {
        stop("No data matches the specified depth_range.")
      }
    }
  }
  
  # Validate compositional columns sum to approximately target_sum
  comp_cols <- c("sandtotal", "silttotal", "claytotal")
  comp_sums <- rowSums(data[, comp_cols])
  
  invalid_sums <- abs(comp_sums - target_sum) > 0.5
  if (any(invalid_sums)) {
    warning(
      "Found ",
      sum(invalid_sums),
      " records with composition sums outside [",
      target_sum - 0.5,
      ", ",
      target_sum + 0.5,
      "]. These will be rescaled."
    )
    
    # Rescale to target_sum
    data[, comp_cols] <- sweep(data[, comp_cols], 1, comp_sums, "/") * target_sum
  }
  
  if (verbose) {
    cat("Processing", nrow(data), "soil component records\n")
    cat("Number of map units:", n_distinct(data$mukey), "\n")
  }
  
  # Create weight column if needed
  if (weight_method == "comppct") {
    data$weight <- data$comppct / 100  # Normalize to 0-1
  } else if (weight_method == "area") {
    # Count components per map unit and equal-weight them
    data$weight <- 1 / as.numeric(table(data$mukey)[match(data$mukey, names(table(data$mukey)))])
  } else if (weight_method == "custom") {
    if (is.null(weight_col) || !weight_col %in% names(data)) {
      stop("weight_col must be specified and exist in data when weight_method = 'custom'")
    }
    data$weight <- data[[weight_col]]
  } else {
    stop("weight_method must be one of: 'comppct', 'area', 'custom'")
  }
  
  # Normalize weights within each map unit
  data <- data[order(data$mukey), ]
  mu_group <- data$mukey
  for (mu in unique(mu_group)) {
    mu_idx <- mu_group == mu
    data$weight[mu_idx] <- data$weight[mu_idx] / sum(data$weight[mu_idx])
  }
  
  # Aggregate to map unit level
  mu_stats <- data.frame()
  constraints <- list()
  
  for (component in comp_cols) {
    if (verbose && component == comp_cols[1]) {
      cat("Computing aggregate statistics...\n")
    }
    
    # Percentile bounds across all components (global)
    all_vals <- data[[component]]
    p10 <- quantile(all_vals, probs = 0.1, na.rm = TRUE)
    p90 <- quantile(all_vals, probs = 0.9, na.rm = TRUE)
    
    comp_name <- gsub("total", "", component)
    comp_name <- toupper(comp_name)
    
    constraints[[comp_name]] <- list(min = as.numeric(p10), max = as.numeric(p90))
  }
  
  # Compute per-map-unit statistics
  mu_list <- split(data, data$mukey)
  
  mu_stats_list <- lapply(mu_list, function(mu_data) {
    list(
      mukey = unique(mu_data$mukey),
      n_components = nrow(mu_data),
      sand_mean = stats::weighted.mean(mu_data$sandtotal, w = mu_data$weight, na.rm = TRUE),
      silt_mean = stats::weighted.mean(mu_data$silttotal, w = mu_data$weight, na.rm = TRUE),
      clay_mean = stats::weighted.mean(mu_data$claytotal, w = mu_data$weight, na.rm = TRUE),
      sand_sd = sd(mu_data$sandtotal),
      silt_sd = sd(mu_data$silttotal),
      clay_sd = sd(mu_data$claytotal)
    )
  })
  
  mu_stats <- do.call(rbind, lapply(mu_stats_list, as.data.frame))
  
  if (verbose) {
    cat("\nConstraints computed:\n")
    for (comp in names(constraints)) {
      cat(
        "  ",
        comp,
        ": [",
        round(constraints[[comp]]$min, 2),
        ", ",
        round(constraints[[comp]]$max, 2),
        "]\n",
        sep = ""
      )
    }
  }
  
  list(
    constraints = constraints,
    samples = data[, c("mukey", "cokey", "sandtotal", "silttotal", "claytotal")],
    samples_weighted = cbind(
      data[, c("mukey", "cokey", "sandtotal", "silttotal", "claytotal")],
      weight = data$weight
    ),
    mu_stats = mu_stats,
    metadata = list(
      source = "SSURGO",
      weight_method = weight_method,
      depth_range = depth_range,
      n_records = nrow(data),
      n_map_units = length(unique(data$mukey)),
      processing_date = Sys.Date()
    ),
    input_data = x  # Keep original for reference/debugging
  )
}


#' Stratify Compositional Parameters by Soil Horizon Depth
#'
#' Create depth-stratified constraint and sample sets from soil composition data
#' with multiple horizons. Enables depth-specific simulation workflows.
#'
#' @param x Data frame with soil composition by horizon. Must contain:
#'   - Compositional columns (e.g., sand, silt, clay)
#'   - `hzdept_r` (numeric): depth to top of horizon (cm)
#'   - `hzdepb_r` (numeric): depth to bottom of horizon (cm)
#'
#' @param comp_cols Character vector, names of compositional columns
#'   (default: c("sandtotal", "silttotal", "claytotal"))
#' @param horizon_breaks Numeric vector of depth breakpoints (cm) to define strata
#'   (default: c(0, 30, 60, 100) for surface, subsurface-1, subsurface-2, deep)
#' @param horizon_labels Character vector, names for strata
#'   (default: c("surface", "subsurface_1", "subsurface_2", "deep"))
#' @param agg_method Character, aggregation method:
#'   - "weighted" (default): weight by horizon thickness
#'   - "mean": simple mean
#'   - "median": median values
#'
#' @return List with one element per stratum. Each contains:
#'   - `constraints`: Constraint list (min, representative, max)
#'   - `samples`: Component samples from that stratum
#'   - `depth_range`: Depth bounds (cm) for the stratum
#'   - `n_samples`: Number of samples in stratum
#'
#' @details
#' Horizons spanning multiple depth strata are split proportionally. For example,
#' a horizon from 20-50 cm split between surface (0-30) and subsurface_1 (30-60)
#' would include:
#'   - 10 cm (20-30) in surface stratum
#'   - 20 cm (30-50) in subsurface_1 stratum
#'
#' @examples
#' # Composition data with multiple horizons
#' horizon_data <- data.frame(
#'   sandtotal = c(35, 30, 25),
#'   silttotal = c(50, 55, 60),
#'   claytotal = c(15, 15, 15),
#'   hzdept_r = c(0, 25, 60),
#'   hzdepb_r = c(25, 60, 100)
#' )
#'
#' strata <- gc_horizon_stratify(
#'   horizon_data,
#'   horizon_breaks = c(0, 30, 60, 100)
#' )
#'
#' # Use surface stratum for shallow simulation
#' surface_params <- strata$surface
#'
#' @importFrom stats quantile weighted.mean
#' @export
gc_horizon_stratify <- function(x,
                                comp_cols = c("sandtotal", "silttotal", "claytotal"),
                                horizon_breaks = c(0, 30, 60, 100),
                                horizon_labels = NULL,
                                agg_method = "weighted") {
  
  # Validate inputs
  if (!all(c("hzdept_r", "hzdepb_r") %in% names(x))) {
    stop("x must contain 'hzdept_r' and 'hzdepb_r' columns")
  }
  
  if (!all(comp_cols %in% names(x))) {
    missing <- setdiff(comp_cols, names(x))
    stop("x missing columns: ", paste(missing, collapse = ", "))
  }
  
  if (is.null(horizon_labels)) {
    n_breaks <- length(horizon_breaks) - 1
    horizon_labels <- paste0("stratum_", seq_len(n_breaks))
  }
  
  if (length(horizon_labels) != length(horizon_breaks) - 1) {
    stop(
      "horizon_labels length (",
      length(horizon_labels),
      ") must equal number of strata (",
      length(horizon_breaks) - 1,
      ")"
    )
  }
  
  # Initialize output list
  strata <- stats::setNames(vector("list", length(horizon_labels)), horizon_labels)
  
  # Process each stratum
  for (stratum_idx in seq_along(horizon_labels)) {
    depth_min <- horizon_breaks[stratum_idx]
    depth_max <- horizon_breaks[stratum_idx + 1]
    label <- horizon_labels[stratum_idx]
    
    # Find horizons that overlap this stratum
    overlaps <- x$hzdept_r < depth_max & x$hzdepb_r > depth_min
    stratum_data <- x[overlaps, , drop = FALSE]
    
    if (nrow(stratum_data) == 0) {
      # Empty stratum
      strata[[label]] <- list(
        constraints = setNames(
          vector("list", length(comp_cols)),
          comp_cols
        ),
        samples = data.frame(),
        depth_range = c(depth_min, depth_max),
        n_samples = 0
      )
      next
    }
    
    # Compute horizon thickness within stratum
    stratum_data$thickness <- pmin(stratum_data$hzdepb_r, depth_max) -
      pmax(stratum_data$hzdept_r, depth_min)
    
    # Aggregate compositional values
    if (agg_method == "weighted") {
      contrib_weights <- stratum_data$thickness / sum(stratum_data$thickness)
      
      agg_comps <- colSums(
        sweep(stratum_data[, comp_cols], 1, contrib_weights, "*")
      )
    } else if (agg_method == "mean") {
      agg_comps <- colMeans(stratum_data[, comp_cols], na.rm = TRUE)
    } else if (agg_method == "median") {
      agg_comps <- apply(stratum_data[, comp_cols], 2, median, na.rm = TRUE)
    } else {
      stop("agg_method must be one of: 'weighted', 'mean', 'median'")
    }
    
    # Build constraints from stratum data
    constraints <- list()
    for (comp in comp_cols) {
      comp_vals <- stratum_data[[comp]]
      constraints[[comp]] <- list(
        min = quantile(comp_vals, 0.1, na.rm = TRUE),
        representative = as.numeric(agg_comps[comp]),
        max = quantile(comp_vals, 0.9, na.rm = TRUE)
      )
    }
    
    strata[[label]] <- list(
      constraints = constraints,
      samples = stratum_data[, comp_cols],
      depth_range = c(depth_min, depth_max),
      n_samples = nrow(stratum_data)
    )
  }
  
  strata
}


#' Apply Component Weighting to Soil Samples
#'
#' Weight compositional samples (typically soil components) based on their
#' relative importance (e.g., component percent, expert opinion, etc.).
#' Useful for aggregating multi-component map units to representative
#' single-component simulations.
#'
#' @param x Data frame with compositional columns and a weight column
#' @param comp_cols Character vector, names of compositional columns
#' @param weight_col Character, name of weight column (must exist in samples)
#' @param normalize Logical, normalize weights to sum to 1 (default TRUE)
#' @param method Character, weighting application method:
#'   - "sample_weights" (default): return samples with weights column
#'   - "weighted_mean": aggregate to single weighted composition
#'   - "quantile_weighted": return weighted quantiles for constraints
#'
#' @return
#' - If method = "sample_weights": data frame with added/modified weight column
#' - If method = "weighted_mean": single-row data frame of weighted composition
#' - If method = "quantile_weighted": list with quantile-based constraints
#'
#' @examples
#' soil_samples <- data.frame(
#'   sand = c(35, 25, 45),
#'   silt = c(50, 60, 40),
#'   clay = c(15, 15, 15),
#'   comppct = c(60, 25, 15)  # Component percent as weight
#' )
#'
#' # Apply weighting
#' weighted_samples <- gc_weight_components(
#'   soil_samples,
#'   comp_cols = c("sand", "silt", "clay"),
#'   weight_col = "comppct",
#'   method = "sample_weights"
#' )
#'
#' @importFrom stats weighted.mean
#' @export
gc_weight_components <- function(x,
                                  comp_cols,
                                  weight_col,
                                  normalize = TRUE,
                                  method = "sample_weights") {
  
  # Validate inputs
  if (!weight_col %in% names(x)) {
    stop("weight_col '", weight_col, "' not found in x")
  }
  
  if (!all(comp_cols %in% names(x))) {
    missing <- setdiff(comp_cols, names(x))
    stop("comp_cols not found: ", paste(missing, collapse = ", "))
  }
  
  # Extract and validate weights
  weights <- x[[weight_col]]
  
  if (any(is.na(weights))) {
    warning("Found NA weights; removing affected rows")
    valid_idx <- !is.na(weights)
    x <- x[valid_idx, ]
    weights <- weights[valid_idx]
  }
  
  if (any(weights < 0)) {
    stop("Weights must be non-negative")
  }
  
  # Normalize weights if requested
  if (normalize && sum(weights) > 0) {
    weights <- weights / sum(weights)
  }
  
  # Apply method
  if (method == "sample_weights") {
    # Return samples with normalized weight column
    result <- x
    result$weight <- weights
    return(result)
    
  } else if (method == "weighted_mean") {
    # Compute single weighted composition
    weighted_comps <- sapply(comp_cols, function(col) {
      stats::weighted.mean(x[[col]], w = weights, na.rm = TRUE)
    })
    
    return(as.data.frame(as.list(weighted_comps)))
    
  } else if (method == "quantile_weighted") {
    # Compute weighted quantiles for each component
    constraints <- list()
    
    for (col in comp_cols) {
      comp_vals <- x[[col]]
      
      # Sort by value
      sorted_idx <- order(comp_vals)
      sorted_vals <- comp_vals[sorted_idx]
      sorted_weights <- weights[sorted_idx]
      
      # Cumulative weights
      cum_weights <- cumsum(sorted_weights)
      
      # Weighted quantiles
      p10_idx <- which(cum_weights >= 0.1)[1]
      p50_idx <- which(cum_weights >= 0.5)[1]
      p90_idx <- which(cum_weights >= 0.9)[1]
      
      constraints[[col]] <- list(
        min = sorted_vals[p10_idx],
        representative = sorted_vals[p50_idx],
        max = sorted_vals[p90_idx]
      )
    }
    
    return(constraints)
    
  } else {
    stop("method must be one of: 'sample_weights', 'weighted_mean', 'quantile_weighted'")
  }
}


#' Define Compositional Constraints from Parametric Specifications
#'
#' Create constraint bounds for geostatistical simulation based on expert
#' knowledge, literature values, or parametric distributions rather than
#' empirical data. Enables "prior-driven" simulation modes.
#'
#' @param components Character vector, names of compositional components
#'   (default: c("SAND", "SILT", "CLAY"))
#' @param means Numeric vector, central values (representative composition) for components
#' @param sds Numeric vector, standard deviations for components
#' @param mins Numeric vector, hard minimum bounds (default: means - 2*sds)
#' @param maxs Numeric vector, hard maximum bounds (default: means + 2*sds)
#' @param target_sum Numeric, target composition sum (default 100)
#' @param dist Character, distribution assumption ("normal", "uniform", or "custom")
#' @param correlations Numeric matrix, correlation structure between components
#'   (optional, for documentation/future use)
#'
#' @return List with elements:
#'   - `constraints`: Constraint list suitable for gc_expand_bounds()
#'   - `parameters`: Vector of means and SDs
#'   - `metadata`: Distribution and correlation info
#'
#' @details
#' This function constructs compositional constraints from parametric specifications,
#' useful when:
#'   - Simulation domain has no empirical data
#'   - Expert judgment is the primary knowledge source
#'   - Performing sensitivity analyses around nominal values
#'
#' The function does NOT enforce the compositional sum constraint strictly;
#' users should verify that specified means and bounds are coherent before
#' using the constraints in gc_expand_bounds().
#'
#' @examples
#' # Loamy texture: loamy-sand typical composition
#' constraints <- gc_parametric_constraints(
#'   components = c("SAND", "SILT", "CLAY"),
#'   means = c(35, 50, 15),
#'   sds = c(5, 5, 3),
#'   dist = "normal"
#' )
#'
#' # Verify constraint bounds
#' str(constraints$constraints)
#'
#' @export
gc_parametric_constraints <- function(components = c("SAND", "SILT", "CLAY"),
                                       means,
                                       sds,
                                       mins = NULL,
                                       maxs = NULL,
                                       target_sum = 100,
                                       dist = "normal",
                                       correlations = NULL) {
  
  # Validate inputs
  if (length(components) != length(means)) {
    stop("components and means must have same length")
  }
  
  if (length(components) != length(sds)) {
    stop("components and sds must have same length")
  }
  
  if (any(sds < 0)) {
    stop("sds must be non-negative")
  }
  
  # Set default bounds if not provided
  if (is.null(mins)) {
    mins <- pmax(means - 2 * sds, 0)  # No negative percentages
  }
  
  if (is.null(maxs)) {
    maxs <- pmin(means + 2 * sds, target_sum)  # Can't exceed 100%
  }
  
  # Validate min/max
  if (any(mins < 0) || any(maxs > target_sum)) {
    warning("mins/maxs adjusted to valid range [0, target_sum]")
    mins <- pmax(mins, 0)
    maxs <- pmin(maxs, target_sum)
  }
  
  if (any(mins > maxs)) {
    stop("mins > maxs for some components; unable to build valid constraints")
  }
  
  # Build constraint list
  constraints <- setNames(vector("list", length(components)), components)
  for (i in seq_along(components)) {
    constraints[[components[i]]] <- list(
      min = mins[i],
      max = maxs[i]
    )
  }
  
  # Check if mean falls within bounds
  means_valid <- means >= mins & means <= maxs
  if (!all(means_valid)) {
    warning(
      "Some means fall outside [min, max] bounds; ",
      "recommend adjusting means or expanding bounds"
    )
  }
  
  # Check compositional consistency
  mean_sum <- sum(means)
  if (abs(mean_sum - target_sum) > 1) {
    warning(
      "Mean composition sums to ",
      round(mean_sum, 2),
      " instead of target_sum = ",
      target_sum,
      "; ",
      "consider rescaling means"
    )
  }
  
  list(
    constraints = constraints,
    parameters = data.frame(
      component = components,
      mean = means,
      sd = sds,
      min = mins,
      max = maxs
    ),
    metadata = list(
      distribution = dist,
      target_sum = target_sum,
      mean_sum = mean_sum,
      correlations = correlations,
      creation_date = Sys.Date()
    )
  )
}
