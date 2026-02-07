#' Zone-Stratified Ensemble Aggregation
#'
#' Compute zone-specific ensemble statistics (mean, SD, confidence intervals) from
#' hierarchical model realizations. Essential for understanding zone-scale uncertainty
#' and decision support in spatially heterogeneous domains.
#'
#' @param ensemble List returned by [gc_sim_hierarchical()] containing ensemble realizations
#' @param zone_col Character name of zone column in ensemble$data, or NULL (see Details)
#' @param zone_definition sf polygon, SpatialPolygons, or factor vector defining zones.
#'   If ensemble data includes zone_id column, this may be omitted.
#' @param conf_level Numeric confidence level (default 0.95 for 95% confidence intervals)
#' @param na.rm Logical; remove NA values before computing statistics
#' @param verbose Logical; print per-zone summary to console
#'
#' @details
#' ## Zone Definition Flexibility
#'
#' Zones can be specified three ways:
#'
#' 1. **Factor in ensemble data** (fastest): Set `zone_col` to column name in ensemble$data
#'    - Example: `gc_ensemble_per_zone(ensemble, zone_col = "zone_id")`
#'
#' 2. **sf polygon geometry**: Provide SpatialPolygons object; points assigned by intersection
#'    - Requires ensemble$data to have coordinates (x, y columns)
#'    - Example: `gc_ensemble_per_zone(ensemble, zone_definition = zones_sf)`
#'
#' 3. **Factor vector**: Length must match ensemble location count
#'    - Example: `gc_ensemble_per_zone(ensemble, zone_definition = zone_factor)`
#'
#' ## Computation Details
#'
#' Per-zone statistics computed for all compositional columns in ensemble:
#' - **Mean**: `mean(X_zone)`
#' - **SD**: `sd(X_zone)`
#' - **95% CI**: quantile(X_zone, c(0.025, 0.975))
#'
#' When ensemble has multiple realizations per location, all points (location × realization)
#' are included in zone aggregation. This preserves spatial uncertainty structure.
#'
#' ## Shrinkage Intensity Assessment
#'
#' For hierarchical models, zone SD reflects both data uncertainty and shrinkage pooling:
#' - Large zones (n >> 5): SD ≈ true parameter uncertainty
#' - Small zones (n << 5): SD reduced by shrinkage; compare with global SD for pooling intensity
#'
#' @return List with elements:
#'   - `zone_stats`: Data frame with columns:
#'     - `zone`: Zone name/ID
#'     - `n`: Count of ensemble points in zone
#'     - `COMPONENT_mean`, `COMPONENT_sd`: Mean and SD for each component
#'     - `COMPONENT_lower`, `COMPONENT_upper`: Confidence interval bounds
#'   - `zone_names`: Character vector of unique zone identifiers
#'   - `n_realizations`: Number of ensemble realizations (if detectable)
#'   - `components`: Names of compositional columns analyzed
#'
#' @seealso [gc_fit_hierarchical_model()], [gc_sim_hierarchical()],
#'   [gc_aggregate_realizations()] (global ensemble aggregation),
#'   [gc_percentile_map_per_zone()], [gc_probability_map_per_zone()]
#'
#' @examples
#' \dontrun{
#' # Generate ensemble from stratified hierarchical model
#' data(soil_hierarchy_example)
#'
#' hz_fit <- gc_fit_hierarchical_model(
#'   data = soil_hierarchy_example$comp,
#'   hierarchy = soil_hierarchy_example$hierarchy,
#'   backend = "analytical"
#' )
#'
#' ensemble <- gc_sim_hierarchical(
#'   object = hz_fit,
#'   n = 100,
#'   zone_vector = soil_hierarchy_example$zone,
#'   coords = soil_hierarchy_example$coords
#' )
#'
#' # Using zone column already in ensemble data
#' per_zone_stats <- gc_ensemble_per_zone(
#'   ensemble = ensemble,
#'   zone_col = "zone_id",
#'   conf_level = 0.95,
#'   verbose = TRUE
#' )
#'
#' # View zone-specific means and uncertainty
#' print(per_zone_stats$zone_stats)
#' }
#'
#' @export
gc_ensemble_per_zone <- function(ensemble,
                                  zone_col = NULL,
                                  zone_definition = NULL,
                                  conf_level = 0.95,
                                  na.rm = TRUE,
                                  verbose = FALSE) {

  if (!is.list(ensemble) || !"data" %in% names(ensemble)) {
    stop("ensemble must be a list with 'data' element (from gc_sim_hierarchical)")
  }

  ensemble_data <- ensemble$data

  # Determine zone assignments
  zone_assignment <- NULL

  if (!is.null(zone_col)) {
    if (!(zone_col %in% names(ensemble_data))) {
      stop(sprintf("Zone column '%s' not found in ensemble$data", zone_col))
    }
    zone_assignment <- ensemble_data[[zone_col]]
  } else if (!is.null(zone_definition)) {
    if (is.factor(zone_definition) || is.character(zone_definition)) {
      if (length(zone_definition) != nrow(ensemble_data)) {
        stop("zone_definition length must equal ensemble point count")
      }
      zone_assignment <- zone_definition
    } else if (inherits(zone_definition, "sf")) {
      # Spatial assignment via sf
      message("Performing spatial point-in-polygon overlay...")
      coords_cols <- c("x", "y")
      if (!all(coords_cols %in% names(ensemble_data))) {
        stop("ensemble$data must have 'x', 'y' columns for spatial zone assignment")
      }

      coords_sf <- sf::st_as_sf(
        ensemble_data,
        coords = coords_cols,
        crs = sf::st_crs(zone_definition)
      )

      # Assign zones via intersection
      overlay <- sf::st_join(coords_sf, zone_definition, join = sf::st_within)

      # Use first geometry column name as zone identifier
      geom_col <- attr(zone_definition, "sf_column")
      zone_cols <- names(zone_definition)[names(zone_definition) != geom_col]

      if (length(zone_cols) > 0) {
        zone_assignment <- overlay[[zone_cols[1]]]
      } else {
        zone_assignment <- seq_len(nrow(overlay))
      }
    } else {
      stop("zone_definition must be factor, character, or sf object")
    }
  } else {
    stop("Either zone_col or zone_definition must be specified")
  }

  if (is.null(zone_assignment)) {
    stop("Could not assign zones to ensemble points")
  }

  # Identify compositional columns
  non_comp_cols <- c("x", "y", "zone_id", zone_col)
  comp_cols <- !(names(ensemble_data) %in% non_comp_cols)
  comp_data <- ensemble_data[, comp_cols]

  if (ncol(comp_data) == 0) {
    stop("No compositional columns found in ensemble$data")
  }

  zone_names <- unique(zone_assignment)
  n_zones <- length(zone_names)

  # Pre-allocate results data frame
  alpha <- 1 - conf_level
  result_cols <- c("zone", "n")

  for (comp in names(comp_data)) {
    result_cols <- c(
      result_cols,
      paste0(comp, "_mean"),
      paste0(comp, "_sd"),
      paste0(comp, "_lower"),
      paste0(comp, "_upper")
    )
  }

  zone_stats <- data.frame(matrix(
    nrow = n_zones,
    ncol = length(result_cols)
  ))
  names(zone_stats) <- result_cols

  # Compute per-zone statistics
  for (i in seq_along(zone_names)) {
    zone_name <- zone_names[i]
    zone_mask <- zone_assignment == zone_name

    zone_stats$zone[i] <- as.character(zone_name)
    zone_stats$n[i] <- sum(zone_mask)

    for (comp in names(comp_data)) {
      zone_vals <- comp_data[zone_mask, comp]

      if (na.rm) {
        zone_vals <- zone_vals[!is.na(zone_vals)]
      }

      if (length(zone_vals) > 0) {
        zone_stats[[paste0(comp, "_mean")]][i] <- mean(zone_vals)
        zone_stats[[paste0(comp, "_sd")]][i] <- sd(zone_vals)

        ci <- quantile(zone_vals, c(alpha / 2, 1 - alpha / 2), na.rm = na.rm)
        zone_stats[[paste0(comp, "_lower")]][i] <- ci[1]
        zone_stats[[paste0(comp, "_upper")]][i] <- ci[2]
      }
    }
  }

  if (verbose) {
    cat("Zone-wise Ensemble Aggregation Summary:\n")
    cat(sprintf("  Total zones: %d\n", n_zones))
    cat(sprintf("  Total ensemble points: %d\n", nrow(ensemble_data)))
    cat(sprintf("  Components analyzed: %s\n", paste(names(comp_data), collapse = ", ")))
    cat(sprintf("  Confidence level: %.1f%%\n", conf_level * 100))
    cat("\nPer-zone sample sizes:\n")
    for (i in seq_along(zone_names)) {
      cat(sprintf("  %s: %d points\n", zone_stats$zone[i], zone_stats$n[i]))
    }
  }

  structure(
    list(
      zone_stats = zone_stats,
      zone_names = zone_names,
      n_zones = n_zones,
      conf_level = conf_level,
      components = names(comp_data)
    ),
    class = c("ensemble_per_zone", "list")
  )
}


#' Per-Zone Percentile Maps
#'
#' Extract percentiles (e.g., P10, P50, P90) of ensemble realizations stratified by zone.
#' Useful for visualizing uncertainty bounds and median predictions within each zone.
#'
#' @param ensemble List from [gc_sim_hierarchical()]
#' @param zone_col Character zone column name in ensemble$data, or NULL
#' @param zone_definition sf polygon, factor, or character zone definitions (see [gc_ensemble_per_zone()])
#' @param probs Numeric vector of percentiles (default c(0.10, 0.50, 0.90) for P10, P50, P90)
#' @param components Character vector of components to compute (default: all compositional columns)
#'
#' @details
#' This function computes quantiles separately for each zone. Useful for:
#' - Visualizing zone-specific prediction intervals (P10-P90)
#' - Risk mapping with zone-informed bounds
#' - Communicating uncertainty variability across zones
#'
#' ## Output Structure
#'
#' Returns a list with nested structure:
#' ```
#' $zone_percentiles
#'   $zone_name_1
#'     $SAND: data frame with columns p10, p50, p90, ...
#'     $SILT: ...
#'   $zone_name_2
#'     ...
#' ```
#'
#' @return List with elements:
#'   - `zone_percentiles`: Named list of percentile data frames by zone
#'   - `zones`: Character vector of zone names
#'   - `percentiles`: Numeric vector of computed percentiles
#'   - `components`: Character vector of components
#'
#' @seealso [gc_ensemble_per_zone()], [gc_percentile_map()] (global percentiles)
#'
#' @examples
#' \dontrun{
#' # Get P10, P50, P90 percentiles by zone
#' zone_percentiles <- gc_percentile_map_per_zone(
#'   ensemble = ensemble,
#'   zone_col = "zone_id",
#'   probs = c(0.10, 0.50, 0.90)
#' )
#' }
#'
#' @export
gc_percentile_map_per_zone <- function(ensemble,
                                        zone_col = NULL,
                                        zone_definition = NULL,
                                        probs = c(0.10, 0.50, 0.90),
                                        components = NULL) {

  if (!is.list(ensemble) || !"data" %in% names(ensemble)) {
    stop("ensemble must be from gc_sim_hierarchical()")
  }

  ensemble_data <- ensemble$data

  # Get zone assignment (reuse logic from gc_ensemble_per_zone)
  if (!is.null(zone_col)) {
    if (!(zone_col %in% names(ensemble_data))) {
      stop(sprintf("Zone column '%s' not found", zone_col))
    }
    zone_assignment <- ensemble_data[[zone_col]]
  } else if (!is.null(zone_definition)) {
    if (is.factor(zone_definition) || is.character(zone_definition)) {
      zone_assignment <- zone_definition
    } else {
      stop("zone_definition must be factor or character for percentile extraction")
    }
  } else {
    stop("Either zone_col or zone_definition required")
  }

  zone_names <- unique(zone_assignment)

  # Identify compositional columns
  non_comp_cols <- c("x", "y", "zone_id", zone_col)
  comp_cols <- !(names(ensemble_data) %in% non_comp_cols)
  comp_names <- names(ensemble_data)[comp_cols]

  if (!is.null(components)) {
    comp_names <- intersect(components, comp_names)
  }

  # Compute percentiles per zone
  zone_percentiles <- list()

  for (zone_name in zone_names) {
    zone_mask <- zone_assignment == zone_name
    zone_data <- ensemble_data[zone_mask, comp_names, drop = FALSE]

    zone_percentiles[[as.character(zone_name)]] <- list()

    for (comp in comp_names) {
      vals <- zone_data[[comp]]
      zone_percentiles[[as.character(zone_name)]][[comp]] <-
        setNames(as.list(quantile(vals, probs, na.rm = TRUE)),
                 paste0("p", round(probs * 100)))
    }
  }

  list(
    zone_percentiles = zone_percentiles,
    zones = zone_names,
    percentiles = probs,
    components = comp_names
  )
}


#' Per-Zone Probability Maps (Risk Assessment)
#'
#' Compute probability that ensemble realizations exceed specified thresholds,
#' stratified by zone. Essential for risk-based decision making in heterogeneous domains.
#'
#' @param ensemble List from [gc_sim_hierarchical()]
#' @param zone_col Character zone column name, or NULL
#' @param zone_definition sf polygon, factor, or character zones
#' @param thresholds Named numeric vector of thresholds (e.g., c(SAND = 50, CLAY = 30))
#'   indicating "risk" levels per component
#' @param operators Character vector of comparison operators (">" [default], "<", ">=", "<=")
#'   (recycled if length 1, else must match threshold length)
#'
#' @details
#' ## Risk Definition
#'
#' Probability computed as proportion of ensemble realizations meeting thresholds:
#'
#' $$P(\text{risk}) = \frac{\text{count}(X > \theta)}{\text{total count}}$$
#'
#' For multiple thresholds, risk is computed independently per component
#' (not as joint probability).
#'
#' ## Management Interpretation
#'
#' - P = 0.95: Very likely to exceed threshold (95% of ensemble realizes risk)
#' - P = 0.50: Uncertain; median risk
#' - P = 0.05: Unlikely to exceed (only 5% of ensemble)
#'
#' Use zone-specific probabilities to tailor risk communication:
#' - Upland (P = 0.05): "Low risk in upland"
#' - Lowland (P = 0.80): "High risk in lowland, recommend mitigation"
#'
#' @return Data frame with columns:
#'   - `zone`: Zone identifier
#'   - `component`: Compositional component
#'   - `threshold`: Threshold value
#'   - `n_exceed`: Count of realizations exceeding threshold
#'   - `n_total`: Total realizations in zone
#'   - `probability`: Proportion exceeding (0-1)
#'   - `probability_pct`: Percentage (0-100)
#'
#' @seealso [gc_ensemble_per_zone()], [gc_probability_map()] (global risk)
#'
#' @examples
#' \dontrun{
#' # Define risk thresholds: high sand (>60%) or high clay (>35%)
#' risk <- gc_probability_map_per_zone(
#'   ensemble = ensemble,
#'   zone_col = "zone_id",
#'   thresholds = c(SAND = 60, CLAY = 35),
#'   operators = ">"
#' )
#'
#' # Display by zone
#' print(risk)
#' }
#'
#' @export
gc_probability_map_per_zone <- function(ensemble,
                                         zone_col = NULL,
                                         zone_definition = NULL,
                                         thresholds = NULL,
                                         operators = ">") {

  if (!is.list(ensemble) || !"data" %in% names(ensemble)) {
    stop("ensemble must be from gc_sim_hierarchical()")
  }

  if (is.null(thresholds)) {
    stop("thresholds must be specified (e.g., c(SAND = 50, CLAY = 30))")
  }

  if (!is.numeric(thresholds) || is.null(names(thresholds))) {
    stop("thresholds must be a named numeric vector")
  }

  ensemble_data <- ensemble$data

  # Get zone assignment
  if (!is.null(zone_col)) {
    zone_assignment <- ensemble_data[[zone_col]]
  } else if (!is.null(zone_definition)) {
    zone_assignment <- zone_definition
  } else {
    stop("Either zone_col or zone_definition required")
  }

  zone_names <- unique(zone_assignment)

  # Recycle operators if needed
  if (length(operators) == 1) {
    operators <- rep(operators, length(thresholds))
  } else if (length(operators) != length(thresholds)) {
    stop("operators length must be 1 or match thresholds length")
  }

  # Compute per-zone probabilities
  results <- list()

  for (zone_name in zone_names) {
    zone_mask <- zone_assignment == zone_name
    zone_data <- ensemble_data[zone_mask, names(thresholds), drop = FALSE]
    n_total <- sum(zone_mask)

    for (j in seq_along(thresholds)) {
      comp <- names(thresholds)[j]
      thresh <- thresholds[j]
      op <- operators[j]

      if (!(comp %in% names(zone_data))) {
        warning(sprintf("Component '%s' not found in ensemble", comp))
        next
      }

      vals <- zone_data[[comp]]

      # Apply comparison operator
      exceed <- switch(op,
        ">" = vals > thresh,
        "<" = vals < thresh,
        ">=" = vals >= thresh,
        "<=" = vals <= thresh,
        stop(sprintf("Unknown operator '%s'", op))
      )

      n_exceed <- sum(exceed, na.rm = TRUE)
      prob <- n_exceed / n_total

      results[[length(results) + 1]] <- list(
        zone = as.character(zone_name),
        component = comp,
        threshold = thresh,
        operator = op,
        n_exceed = n_exceed,
        n_total = n_total,
        probability = prob,
        probability_pct = round(prob * 100, 1)
      )
    }
  }

  as.data.frame(do.call(rbind, lapply(results, as.data.frame)))
}


#' Per-Zone Ensemble Quality Report
#'
#' Generate diagnostic summaries of ensemble quality and shrinkage pooling intensity
#' for each zone. Assesses convergence, coverage, and whether hierarchical borrowing
#' is effective.
#'
#' @param ensemble List from [gc_sim_hierarchical()]
#' @param zone_col Character zone column name, or NULL
#' @param zone_definition sf polygon, factor, or character zones
#' @param global_stats Optional list from [gc_aggregate_realizations()] for comparison
#' @param verbose Logical; print detailed report to console
#'
#' @details
#' ## Quality Metrics
#'
#' - **Coverage**: Proportion of realizations producing valid compositions (0-1)
#' - **Effective Sample Size**: Zone n × coverage (accounts for invalid realizations)
#' - **Shrinkage Intensity**: Compare zone SD to global SD; ratio < 1 indicates pooling
#' - **Range**: Max - Min; wide ranges suggest high ensemble diversity
#'
#' ## Interpretation
#'
#' | Metric | Good | Poor | Action |
#' |--------|------|------|--------|
#' | Coverage | > 0.99 | < 0.95 | Check constraint satisfaction |
#' | Effective n | > 10 | < 5 | Consider combining zones |
#' | Shrinkage ratio | 0.5-0.8 | > 1.0 | Increase hierarchy prior strength |
#' | Range | Narrow | Extreme | Check ensemble diversity |
#'
#' @return List with elements:
#'   - `quality_summary`: Data frame of per-zone diagnostics
#'   - `zones`: Character vector of zone names
#'   - `convergence_notes`: Character vector of interpretation strings
#'
#' @seealso [gc_ensemble_per_zone()], [gc_ensemble_quality_report()] (global report)
#'
#' @examples
#' \dontrun{
#' global_stats <- gc_aggregate_realizations(ensemble)
#'
#' quality <- gc_ensemble_quality_report_per_zone(
#'   ensemble = ensemble,
#'   zone_col = "zone_id",
#'   global_stats = global_stats,
#'   verbose = TRUE
#' )
#'
#' print(quality$quality_summary)
#' }
#'
#' @export
gc_ensemble_quality_report_per_zone <- function(ensemble,
                                                 zone_col = NULL,
                                                 zone_definition = NULL,
                                                 global_stats = NULL,
                                                 verbose = FALSE) {

  if (!is.list(ensemble) || !"data" %in% names(ensemble)) {
    stop("ensemble must be from gc_sim_hierarchical()")
  }

  ensemble_data <- ensemble$data

  # Get zone assignment
  if (!is.null(zone_col)) {
    zone_assignment <- ensemble_data[[zone_col]]
  } else if (!is.null(zone_definition)) {
    zone_assignment <- zone_definition
  } else {
    stop("Either zone_col or zone_definition required")
  }

  zone_names <- unique(zone_assignment)

  # Identify compositional columns
  non_comp_cols <- c("x", "y", "zone_id", zone_col)
  comp_cols <- !(names(ensemble_data) %in% non_comp_cols)
  comp_data <- ensemble_data[, comp_cols]

  # Compute global statistics for shrinkage comparison
  if (is.null(global_stats)) {
    global_sd <- sapply(comp_data, sd, na.rm = TRUE)
  } else if (is.list(global_stats)) {
    # Extract from ensemble_quality_report output if available
    global_sd <- sapply(names(comp_data), function(comp) {
      if (!is.null(global_stats[[paste0(comp, "_sd")]])) {
        global_stats[[paste0(comp, "_sd")]]
      } else {
        sd(comp_data[[comp]], na.rm = TRUE)
      }
    })
  } else {
    global_sd <- sapply(comp_data, sd, na.rm = TRUE)
  }

  # Build quality report
  quality_list <- list()

  for (zone_name in zone_names) {
    zone_mask <- zone_assignment == zone_name
    zone_data <- comp_data[zone_mask, ]
    n_zone <- sum(zone_mask)

    # Coverage: proportion of non-NA values
    coverage <- mean(!is.na(rowSums(zone_data)))

    # Effective sample size
    eff_n <- n_zone * coverage

    # Compute statistics per component
    zone_sd <- sapply(zone_data, sd, na.rm = TRUE)
    zone_range <- sapply(zone_data, function(x) {
      max(x, na.rm = TRUE) - min(x, na.rm = TRUE)
    })

    # Shrinkage intensity (zone SD / global SD)
    shrinkage_ratio <- zone_sd / global_sd

    quality_list[[as.character(zone_name)]] <- list(
      zone = as.character(zone_name),
      n = n_zone,
      coverage = round(coverage, 4),
      eff_n = round(eff_n, 1),
      mean_sd = round(mean(zone_sd, na.rm = TRUE), 3),
      mean_global_sd = round(mean(global_sd, na.rm = TRUE), 3),
      mean_shrinkage_ratio = round(mean(shrinkage_ratio, na.rm = TRUE), 3),
      min_zone_sd = round(min(zone_sd, na.rm = TRUE), 3),
      max_zone_sd = round(max(zone_sd, na.rm = TRUE), 3)
    )
  }

  quality_df <- as.data.frame(do.call(rbind, lapply(quality_list, as.data.frame)))

  # Generate convergence notes
  convergence_notes <- character()

  for (i in seq_len(nrow(quality_df))) {
    zone_name <- quality_df$zone[i]
    eff_n <- quality_df$eff_n[i]
    coverage <- quality_df$coverage[i]
    shrink <- quality_df$mean_shrinkage_ratio[i]

    note <- sprintf("%s:", zone_name)

    if (coverage < 0.99) {
      note <- paste(note, "low coverage")
    }
    if (eff_n < 5) {
      note <- paste(note, "small effective n (consider combining zones)")
    }
    if (shrink > 1.0) {
      note <- paste(note, "zone SD > global (weak hierarchical pooling)")
    }
    if (shrink < 0.5) {
      note <- paste(note, "strong shrinkage (zone estimates dominated by global prior)")
    }

    convergence_notes <- c(convergence_notes, note)
  }

  if (verbose) {
    cat("Per-Zone Ensemble Quality Report\n")
    cat("==================================\n\n")
    cat("Quality Summary:\n")
    print(quality_df)
    cat("\nConvergence Assessment:\n")
    for (note in convergence_notes) {
      cat(sprintf("  %s\n", note))
    }
  }

  list(
    quality_summary = quality_df,
    zones = zone_names,
    convergence_notes = convergence_notes
  )
}
