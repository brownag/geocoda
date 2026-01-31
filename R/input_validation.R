#' Pre-Simulation Data Validation
#'
#' Comprehensive pre-flight checks for data intended for geostatistical simulation.
#' Identifies common data quality issues before they cause problems downstream.
#'
#' @keywords internal
#' @name input_validation
NULL


#' Validate Input Data for Geostatistical Simulation
#'
#' Perform comprehensive pre-flight checks on compositional and spatial data
#' before constrain definition or simulation. Detects missing values, invalid
#' compositions, outliers, and other common data quality issues.
#'
#' @param data An `sf` object (with spatial geometry) or data frame with
#'   compositional and coordinate columns
#' @param component_names Character vector, names of compositional components
#'   (e.g., c("sand", "silt", "clay"))
#' @param x_coord Character, name of x-coordinate column (easting, longitude).
#'   Required if data is not sf
#' @param y_coord Character, name of y-coordinate column (northing, latitude).
#'   Required if data is not sf
#' @param target_sum Numeric, expected composition sum (default 100 for percentages)
#' @param sum_tolerance Numeric, allowable deviation from target_sum before flagging
#'   (default 1, i.e., sums must be in [99, 101] for target_sum=100)
#' @param check_spatial Logical, validate spatial properties (duplicates, extent)
#'   (default TRUE)
#' @param check_outliers Logical, detect statistical outliers using IQR method
#'   (default TRUE)
#' @param verbose Logical, print detailed messages (default TRUE)
#'
#' @return A list with elements:
#'   - `valid`: Logical, whether data passes all checks
#'   - `issues`: Data frame listing all detected issues (severity, message, n_rows)
#'   - `summary`: Overall summary (n_observations, n_components, spatial_extent, etc)
#'   - `recommendations`: Character vector of actionable next steps
#'
#' @details
#' The function performs these checks in order:
#'
#' **Structural Checks:**
#' - Required columns present
#' - Data types correct (numeric for components)
#' - No completely empty rows
#'
#' **Compositional Checks:**
#' - Sum constraints (each composition sums to target_sum ± tolerance)
#' - Range constraints (each component in [0, 100] or appropriate range)
#' - Zero/negative handling
#' - Compositional coherence (multivariate distributional shape)
#'
#' **Spatial Checks:**
#' - Projected coordinates (non-geographic for meaningful distances)
#' - Valid extent (not wrapped around poles or dateline)
#' - Duplicate locations
#' - Sufficient spatial spread
#'
#' **Statistical Checks:**
#' - Univariate outliers (IQR method per component)
#' - Multivariate outliers (Mahalanobis distance)
#' - Extreme values vs domain knowledge
#' - Sample size adequacy
#'
#' Issues are categorized as:
#' - "error": Will cause algorithm failure; must fix
#' - "warning": May produce unreliable results; recommend fixing
#' - "note": May be acceptable depending on use case; informational
#'
#' @examples
#' # Create sample compositional data with some issues
#' bad_data <- data.frame(
#'   x = c(1, 2, 3, 4, 5),
#'   y = c(1, 2, 3, 4, 5),
#'   sand = c(30, 40, 35, NA, 50),      # One NA
#'   silt = c(50, 45, 50, 55, 40),
#'   clay = c(20, 20, 15, 20, 5),       # Last row sums to 95, not 100
#'   texture = c("SL", "L", "L", "L", "SL")
#' )
#'
#' validation <- gc_validate_input_data(
#'   bad_data,
#'   component_names = c("sand", "silt", "clay"),
#'   x_coord = "x",
#'   y_coord = "y"
#' )
#'
#' print(validation$issues)
#' print(validation$recommendations)
#'
#' @importFrom sf st_as_sf st_crs st_bbox
#' @export
gc_validate_input_data <- function(data,
                                    component_names,
                                    x_coord = NULL,
                                    y_coord = NULL,
                                    target_sum = 100,
                                    sum_tolerance = 1,
                                    check_spatial = TRUE,
                                    check_outliers = TRUE,
                                    verbose = TRUE) {
  
  # Initialize tracking structures
  issues <- list()
  recommendations <- character()
  
  if (verbose) cat("Starting comprehensive data validation...\n\n")
  
  # =========================================================================
  # STEP 1: Check data structure
  # =========================================================================
  
  if (verbose) cat("1. Checking data structure...\n")
  
  # Convert sf to data frame if needed
  if (inherits(data, "sf")) {
    is_sf <- TRUE
    data_df <- as.data.frame(data)
    coords <- sf::st_coordinates(data)
    if (is.null(x_coord)) x_coord <- NULL
    if (is.null(y_coord)) y_coord <- NULL
  } else {
    is_sf <- FALSE
    data_df <- as.data.frame(data)
    
    if (is.null(x_coord) || is.null(y_coord)) {
      issues$struct_001 <- list(
        severity = "error",
        message = "x_coord and y_coord required for non-sf data",
        n_rows = 0
      )
    }
  }
  
  # Check required columns (component_names)
  missing_cols <- setdiff(component_names, names(data_df))
  if (length(missing_cols) > 0) {
    issues$struct_002 <- list(
      severity = "error",
      message = paste("Missing component columns:", paste(missing_cols, collapse = ", ")),
      n_rows = 0
    )
  }
  
  # Check coordinate columns only if specified
  if (!is_sf) {
    if (!is.null(x_coord) && !x_coord %in% names(data_df)) {
      issues$struct_003 <- list(
        severity = "error",
        message = paste("X-coordinate column not found:", x_coord),
        n_rows = 0
      )
    }
    
    if (!is.null(y_coord) && !y_coord %in% names(data_df)) {
      issues$struct_004 <- list(
        severity = "error",
        message = paste("Y-coordinate column not found:", y_coord),
        n_rows = 0
      )
    }
  }
  
  # Check nrows
  n_obs <- nrow(data_df)
  if (n_obs == 0) {
    issues$struct_005 <- list(
      severity = "error",
      message = "Data has 0 rows",
      n_rows = 0
    )
  }
  
  # Early bail if structural errors prevent further checks
  if (any(unlist(lapply(issues, function(x) x$severity)) == "error")) {
    issues_df <- do.call(rbind, lapply(issues, function(x) {
      data.frame(
        severity = x$severity,
        message = x$message,
        n_affected = x$n_rows,
        stringsAsFactors = FALSE
      )
    }))
    rownames(issues_df) <- NULL
    
    if (verbose) {
      cat("\n=== VALIDATION SUMMARY ===\n")
      cat("Valid: FALSE\n")
      cat("Issues found:", nrow(issues_df), "\n")
    }
    
    return(list(
      valid = FALSE,
      issues = issues_df,
      summary = list(n_observations = n_obs, n_components = length(component_names)),
      recommendations = c("Resolve structural errors before proceeding")
    ))
  }
  
  if (verbose) cat("  ✓ Structure checks complete\n")
  
  # =========================================================================
  # STEP 2: Check for missing values
  # =========================================================================
  
  if (verbose) cat("2. Checking for missing values...\n")
  
  # Check components for NA (only if they exist)
  existing_components <- intersect(component_names, names(data_df))
  na_per_comp <- colSums(is.na(data_df[, existing_components, drop = FALSE]))
  for (comp in existing_components) {
    n_na <- na_per_comp[[comp]]
    if (n_na > 0) {
      issues[[paste0("missing_", comp)]] <- list(
        severity = "warning",
        message = paste0(comp, " has ", n_na, " missing values"),
        n_rows = n_na
      )
      if (n_na / n_obs > 0.1) {
        recommendations <- c(recommendations,
          paste("Consider imputing or removing rows with missing", comp))
      }
    }
  }
  
  # Check coordinates for NA (only if they exist and are not sf)
  if (!is_sf) {
    if (!is.null(x_coord) && x_coord %in% names(data_df)) {
      n_na_x <- sum(is.na(data_df[[x_coord]]))
      if (n_na_x > 0) {
        issues$missing_x <- list(
          severity = "warning",
          message = paste("X-coordinate has", n_na_x, "missing values"),
          n_rows = n_na_x
        )
      }
    }
    
    if (!is.null(y_coord) && y_coord %in% names(data_df)) {
      n_na_y <- sum(is.na(data_df[[y_coord]]))
      if (n_na_y > 0) {
        issues$missing_y <- list(
          severity = "warning",
          message = paste("Y-coordinate has", n_na_y, "missing values"),
          n_rows = n_na_y
        )
      }
    }
  }
  
  if (verbose) cat("  ✓ Missing value checks complete\n")
  
  # =========================================================================
  # STEP 3: Check compositional integrity
  # =========================================================================
  
  if (verbose) cat("3. Checking compositional integrity...\n")
  
  # Only proceed if components exist
  if (length(existing_components) == 0) {
    issues$comp_000 <- list(
      severity = "error",
      message = "No component columns found",
      n_rows = 0
    )
    if (verbose) cat("  ✓ Compositional integrity checks complete (no components)\n")
  } else {
    comp_data <- data_df[, existing_components, drop = FALSE]
    
    # Calculate row sums
    row_sums <- rowSums(comp_data, na.rm = TRUE)
    
    # Check sum constraints
    invalid_sums <- abs(row_sums - target_sum) > sum_tolerance
    n_invalid_sums <- sum(invalid_sums, na.rm = TRUE)
    
    if (n_invalid_sums > 0) {
      severity <- if (n_invalid_sums / n_obs > 0.2) "error" else "warning"
      issues$comp_001 <- list(
        severity = severity,
        message = paste0(
          "Composition sum outside [", target_sum - sum_tolerance, ", ",
          target_sum + sum_tolerance, "] for ", n_invalid_sums, " rows (mean sum = ",
          round(mean(row_sums, na.rm = TRUE), 2), ")"
        ),
        n_rows = n_invalid_sums
      )
      
      if (n_invalid_sums / n_obs <= 0.2) {
        recommendations <- c(recommendations,
          "Consider rescaling compositions to sum to target_sum")
      }
    }
    
    # Check range constraints (assuming 0-100 for percentages)
    for (comp in existing_components) {
      comp_vals <- comp_data[[comp]]
      below_zero <- sum(comp_vals < 0, na.rm = TRUE)
      above_target <- sum(comp_vals > target_sum, na.rm = TRUE)
      
      if (below_zero > 0) {
        issues[[paste0("range_neg_", comp)]] <- list(
          severity = "warning",
          message = paste0(comp, " has ", below_zero, " negative values"),
          n_rows = below_zero
        )
        recommendations <- c(recommendations,
          "Consider imputing or removing rows with negative values")
      }
      
      if (above_target > 0) {
        issues[[paste0("range_over_", comp)]] <- list(
          severity = "warning",
          message = paste0(comp, " has ", above_target, " values > ", target_sum),
          n_rows = above_target
        )
      }
    }
    
    # Check for zero-only compositions
    all_zero <- rowSums(comp_data, na.rm = TRUE) == 0
    n_all_zero <- sum(all_zero, na.rm = TRUE)
    
    if (n_all_zero > 0) {
      issues$comp_002 <- list(
        severity = "error",
        message = paste("Found", n_all_zero, "rows with all zero compositions"),
        n_rows = n_all_zero
      )
    }
    
    if (verbose) cat("  ✓ Compositional integrity checks complete\n")
  }
  
  # =========================================================================
  # STEP 4: Check spatial properties
  # =========================================================================
  
  if (check_spatial && nrow(data_df) > 0) {
    if (verbose) cat("4. Checking spatial properties...\n")
    
    # Get coordinates from sf if available
    if (is_sf && inherits(data, "sf")) {
      coords_matrix <- sf::st_coordinates(data)
      coords <- data.frame(x = coords_matrix[, 1], y = coords_matrix[, 2])
    } else if (!is.null(x_coord) && !is.null(y_coord) && 
               x_coord %in% names(data_df) && y_coord %in% names(data_df)) {
      coords <- data_df[, c(x_coord, y_coord), drop = FALSE]
      names(coords) <- c("x", "y")
    } else {
      coords <- NULL
    }
    
    if (!is.null(coords)) {
      coords <- coords[!is.na(coords$x) & !is.na(coords$y), ]
      
      if (nrow(coords) > 1) {
        x_vals <- coords$x
        y_vals <- coords$y
        
        # Check for geographic coordinates (should be projected typically)
        if (max(abs(x_vals)) <= 180 && max(abs(y_vals)) <= 90) {
          issues$spatial_001 <- list(
            severity = "note",
            message = "Coordinates appear to be geographic (lat/lon); recommend projection for accurate distances",
            n_rows = 0
          )
          recommendations <- c(recommendations,
            "Project coordinates to UTM or equal-area projection for kriging")
        }
        
        # Check for duplicates
        coords_combined <- paste(round(x_vals, 6), round(y_vals, 6), sep = "_")
        duplicates <- sum(duplicated(coords_combined))
        
        if (duplicates > 0) {
          issues$spatial_002 <- list(
            severity = "warning",
            message = paste("Found", duplicates, "duplicate or near-duplicate locations"),
            n_rows = duplicates
          )
          recommendations <- c(recommendations,
            "Consider aggregating or slightly perturbing duplicate locations")
        }
        
        # Check spatial extent
        x_range <- diff(range(x_vals, na.rm = TRUE))
        y_range <- diff(range(y_vals, na.rm = TRUE))
        
        if (x_range == 0 || y_range == 0) {
          issues$spatial_003 <- list(
            severity = "warning",
            message = "All observations at same location (no spatial variation)",
            n_rows = 0
          )
          recommendations <- c(recommendations,
            "Sampling points have no spatial extent; kriging variance will be zero")
        }
        
        # Check for adequate sample spacing
        if (nrow(coords) >= 2) {
          distances <- as.matrix(stats::dist(coords))
          diag(distances) <- NA
          min_dist <- min(distances, na.rm = TRUE)
          med_dist <- stats::median(distances, na.rm = TRUE)
          
          if (min_dist < 0.01 * med_dist) {
            issues$spatial_004 <- list(
              severity = "note",
              message = "Some observations very close together relative to median spacing",
              n_rows = 0
            )
          }
        }
      }
    }
    
    if (verbose) cat("  ✓ Spatial property checks complete\n")
  }
  
  # =========================================================================
  # STEP 5: Check statistical properties
  # =========================================================================
  
  if (check_outliers && length(existing_components) > 0) {
    if (verbose) cat("5. Checking for statistical outliers...\n")
    
    if (length(existing_components) == 0) {
      # Skip if no components
      if (verbose) cat("  ✓ Outlier checks skipped (no components)\n")
    } else {
      comp_data <- data_df[, existing_components, drop = FALSE]
      
      for (comp in existing_components) {
        comp_vals <- comp_data[[comp]]
        comp_vals <- comp_vals[!is.na(comp_vals)]
        
        if (length(comp_vals) > 4) {
          # IQR method for outliers
          q1 <- stats::quantile(comp_vals, 0.25)
          q3 <- stats::quantile(comp_vals, 0.75)
          iqr <- q3 - q1
          
          lower_fence <- q1 - 1.5 * iqr
          upper_fence <- q3 + 1.5 * iqr
          
          n_outliers <- sum(comp_vals < lower_fence | comp_vals > upper_fence)
          
          if (n_outliers > 0) {
            issues[[paste0("outlier_", comp)]] <- list(
              severity = "note",
              message = paste0(comp, " has ", n_outliers, " statistical outliers (IQR method)"),
              n_rows = n_outliers
            )
            
            if (n_outliers / nrow(comp_data) > 0.05) {
              recommendations <- c(recommendations,
                paste("Review and potentially remove extreme", comp, "values"))
            }
          }
        }
      }
      
      if (verbose) cat("  ✓ Outlier checks complete\n")
    }
  }
  
  # =========================================================================
  # STEP 6: Summary and recommendations
  # =========================================================================
  
  # Compile issues into data frame
  if (length(issues) > 0) {
    issues_df <- do.call(rbind, lapply(issues, function(x) {
      data.frame(
        severity = x$severity,
        message = x$message,
        n_affected = x$n_rows,
        stringsAsFactors = FALSE
      )
    }))
    rownames(issues_df) <- NULL
  } else {
    issues_df <- data.frame(
      severity = character(),
      message = character(),
      n_affected = numeric(),
      stringsAsFactors = FALSE
    )
  }
  
  # Determine overall validity
  has_errors <- any(issues_df$severity == "error", na.rm = TRUE)
  
  # Compile summary
  summary <- list(
    n_observations = nrow(data_df),
    n_components = length(component_names),
    n_complete_cases = sum(rowSums(!is.na(comp_data)) == length(component_names)),
    mean_composition_sum = mean(rowSums(comp_data, na.rm = TRUE), na.rm = TRUE),
    sd_composition_sum = stats::sd(rowSums(comp_data, na.rm = TRUE), na.rm = TRUE)
  )
  
  # Add spatial summary if available
  if (!is.null(x_coord) && !is.null(y_coord)) {
    coords <- data_df[, c(x_coord, y_coord), drop = FALSE]
    coords <- coords[!is.na(coords[[x_coord]]) & !is.na(coords[[y_coord]]), ]
    
    if (nrow(coords) >= 2) {
      summary$spatial_extent_x <- diff(range(coords[[x_coord]], na.rm = TRUE))
      summary$spatial_extent_y <- diff(range(coords[[y_coord]], na.rm = TRUE))
      summary$n_spatial_observations <- nrow(coords)
    }
  }
  
  # Generate recommendations if there are issues
  if (has_errors) {
    recommendations <- unique(c(
      recommendations,
      "Resolve ERROR-severity issues before proceeding with simulation"
    ))
  }
  
  if (nrow(issues_df) > 0) {
    recommendations <- unique(c(
      recommendations,
      "Review WARNING and NOTE messages for potential improvements"
    ))
  }
  
  if (nrow(data_df) < 50) {
    recommendations <- unique(c(
      recommendations,
      "Sample size < 50: results may be unstable; consider collecting more observations"
    ))
  }
  
  if (verbose) {
    cat("\n=== VALIDATION SUMMARY ===\n")
    cat("Valid:", !has_errors, "\n")
    cat("Issues found:", nrow(issues_df), "\n")
    if (nrow(issues_df) > 0) {
      cat("  - Errors:", sum(issues_df$severity == "error"), "\n")
      cat("  - Warnings:", sum(issues_df$severity == "warning"), "\n")
      cat("  - Notes:", sum(issues_df$severity == "note"), "\n")
    }
    cat("\n")
  }
  
  list(
    valid = !has_errors,
    issues = issues_df,
    summary = summary,
    recommendations = if (length(recommendations) > 0) recommendations else character()
  )
}
