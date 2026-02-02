#' Standardize Soil Profiles to Regular Depth Intervals
#'
#' Convert irregular soil horizons to regular depth intervals using mass-preserving
#' splines (quadratic spline interpolation via the aqp package).
#'
#' @param profiles Data frame with soil profile data including:
#'   - hzdept_r: Horizon depth top (cm)
#'   - hzdepb_r: Horizon depth bottom (cm)
#'   - id_col column: Profile identifier
#'   - Component columns with compositional or numeric values
#'
#' @param method Character specifying interpolation method:
#'   - "spline": Mass-preserving quadratic splines (requires aqp package)
#'   - "linear": Simple linear interpolation
#'   Default: "spline"
#'
#' @param intervals Numeric vector of standard depth slices (cm)
#'   Default: c(0, 5, 15, 30, 60, 100) (GlobalSoilMap standard)
#'   These are depth interval BOTTOMS: 0-5, 5-15, 15-30, 30-60, 60-100, 100+ cm
#'
#' @param comp_cols Character vector of column names with values to standardize
#'   (e.g., c("sandtotal", "silttotal", "claytotal") for texture)
#'
#' @param id_col Character: Column name for profile identifier
#'   (e.g., "mukey", "profile_id", "pedon_id")
#'
#' @param verbose Logical: If TRUE, print progress messages
#'
#' @return Data frame with standardized depth intervals containing:
#'   - depth_top: Top of interval (cm)
#'   - depth_bot: Bottom of interval (cm)
#'   - id_col column: Original profile ID
#'   - Standardized values for each comp_cols at each interval
#'
#'   Attributes:
#'   - method: Interpolation method used
#'   - intervals: Depth intervals used
#'
#' @details
#' **Mass-Preserving Splines (aqp::spc2mpspline):**
#' - Ensures total mass (sum) is conserved: ∫ f(z) dz = measured_value
#' - Quadratic interpolation for smooth depth functions
#' - Ideal for soil properties where conservation of mass matters
#' - Example: Total carbon (mg/cm²) in profile must equal sum of horizon carbon
#'
#' **GlobalSoilMap Standard Intervals:**
#' - 0-5 cm (surface layer, influenced by land use, management)
#' - 5-15 cm (topsoil, root zone)
#' - 15-30 cm (subsoil, transitional)
#' - 30-60 cm (subsurface, reduced activity)
#' - 60-100 cm (deep subsoil)
#' - 100+ cm (parent material, if available)
#'
#' **Alternatives to Specify Intervals:**
#' - Custom intervals: c(10, 20, 30, 40, 50) for 10 cm slices
#' - Two-depth system: c(30, 100) for 0-30, 30-100 cm
#' - Single depth: c(100) for 0-100 cm aggregate
#'
#' **Mass Conservation Validation:**
#' Use [gc_validate_depth_standardization()] to check that:
#' - Original horizon values sum to standardized interval values
#' - No mass was gained or lost in interpolation
#'
#' @examples
#' \dontrun{
#' # Load SSURGO component data
#' comp_data <- ssurgo_components  # Your SSURGO data frame
#'
#' # Standardize to GlobalSoilMap depths
#' standardized <- gc_standardize_depths(
#'   profiles = comp_data,
#'   method = "spline",
#'   intervals = c(5, 15, 30, 60, 100),
#'   comp_cols = c("sandtotal", "silttotal", "claytotal"),
#'   id_col = "mukey",
#'   verbose = TRUE
#' )
#'
#' # Check mass conservation
#' validation <- gc_validate_depth_standardization(comp_data, standardized)
#' print(validation)
#' }
#'
#' @export
gc_standardize_depths <- function(profiles, method = "spline",
                                   intervals = c(5, 15, 30, 60, 100),
                                   comp_cols = NULL, id_col = "mukey",
                                   verbose = FALSE) {
  # Input validation
  stopifnot(is.data.frame(profiles))
  stopifnot(all(c("hzdept_r", "hzdepb_r") %in% names(profiles)))
  stopifnot(id_col %in% names(profiles))
  stopifnot(method %in% c("spline", "linear"))
  stopifnot(is.numeric(intervals), length(intervals) > 0, all(intervals > 0))

  # Determine columns to standardize
  if (is.null(comp_cols)) {
    # Auto-detect: exclude coordinate/ID columns
    exclude_cols <- c("hzdept_r", "hzdepb_r", id_col, "x", "y", "gid", "lon", "lat")
    comp_cols <- setdiff(names(profiles), exclude_cols)
    if (length(comp_cols) == 0) {
      stop("No suitable columns to standardize. Specify comp_cols explicitly.")
    }
  }

  stopifnot(all(comp_cols %in% names(profiles)))

  # Check method availability
  if (method == "spline") {
    if (!requireNamespace("aqp", quietly = TRUE)) {
      stop("Method 'spline' requires aqp package. Install with: install.packages('aqp')",
           "Or use method='linear' as fallback.")
    }
  }

  # Split by profile
  profiles_list <- split(profiles, profiles[[id_col]])
  n_profiles <- length(profiles_list)

  if (verbose) {
    message("Standardizing ", n_profiles, " profiles to intervals: ",
            paste(c(0, intervals), collapse = "-"), " cm")
  }

  # Standardize each profile
  standardized_list <- lapply(seq_along(profiles_list), function(i) {
    prof <- profiles_list[[i]]
    prof_id <- unique(prof[[id_col]])

    if (method == "spline") {
      gc_standardize_depths_spline(prof, intervals, comp_cols, id_col, prof_id)
    } else {
      gc_standardize_depths_linear(prof, intervals, comp_cols, id_col, prof_id)
    }
  })

  # Combine results
  standardized <- do.call(rbind, standardized_list)
  rownames(standardized) <- NULL

  # Set attributes
  attr(standardized, "method") <- method
  attr(standardized, "intervals") <- intervals
  attr(standardized, "comp_cols") <- comp_cols

  if (verbose) {
    message("Standardization complete. Output: ", nrow(standardized), " depth slices")
  }

  return(standardized)
}


#' Internal: Spline-based depth standardization
#'
#' @keywords internal
gc_standardize_depths_spline <- function(prof, intervals, comp_cols, id_col, prof_id) {
  # For each component column, fit spline and evaluate at intervals
  standardized_values <- lapply(comp_cols, function(col) {
    tryCatch({
      spline_fit <- aqp::spc2mpspline(
        obj = prof,
        var_name = col,
        d = intervals,
        vhigh = 100,  # Maximum possible value
        vlow = 0      # Minimum possible value
      )
      # Extract fitted values at interval bottoms
      values <- spline_fit$fitted_values
      return(values)
    },
    error = function(e) {
      warning("Spline fitting failed for column '", col, "' in profile ",
              prof_id, ". Returning NA.")
      return(rep(NA_real_, length(intervals)))
    })
  })

  # Build output data frame
  result <- data.frame(
    depth_top = c(0, intervals[-length(intervals)]),
    depth_bot = intervals
  )
  result[[id_col]] <- prof_id

  # Add standardized values
  for (j in seq_along(comp_cols)) {
    result[[comp_cols[j]]] <- standardized_values[[j]]
  }

  return(result)
}


#' Internal: Linear interpolation depth standardization
#'
#' @keywords internal
gc_standardize_depths_linear <- function(prof, intervals, comp_cols, id_col, prof_id) {
  standardized_values <- lapply(comp_cols, function(col) {
    # For each interval, average values from overlapping horizons
    interval_values <- numeric(length(intervals))

    for (i in seq_along(intervals)) {
      interval_top <- if (i == 1) 0 else intervals[i - 1]
      interval_bot <- intervals[i]

      # Find overlapping horizons
      overlap <- prof$hzdept_r < interval_bot & prof$hzdepb_r > interval_top
      overlap_indices <- which(overlap)

      if (length(overlap_indices) > 0) {
        # Weighted average by overlap fraction
        weights <- numeric(length(overlap_indices))
        for (j in seq_along(overlap_indices)) {
          idx <- overlap_indices[j]
          overlap_top <- max(prof$hzdept_r[idx], interval_top)
          overlap_bot <- min(prof$hzdepb_r[idx], interval_bot)
          weights[j] <- overlap_bot - overlap_top
        }
        weights <- weights / sum(weights)
        interval_values[i] <- sum(prof[[col]][overlap_indices] * weights)
      } else {
        interval_values[i] <- NA_real_
      }
    }

    return(interval_values)
  })

  # Build output data frame
  result <- data.frame(
    depth_top = c(0, intervals[-length(intervals)]),
    depth_bot = intervals
  )
  result[[id_col]] <- prof_id

  # Add standardized values
  for (j in seq_along(comp_cols)) {
    result[[comp_cols[j]]] <- standardized_values[[j]]
  }

  return(result)
}


#' Validate Mass Conservation in Depth Standardization
#'
#' Check that depth standardization conserved total mass (sum) of properties.
#' Important for quantities like total carbon, clay mass, etc.
#'
#' @param profiles_original Data frame with original horizon data
#' @param profiles_standardized Data frame from [gc_standardize_depths()]
#' @param comp_cols Character vector of columns to validate
#'   (should match columns used in standardization)
#' @param id_col Character: Profile identifier column
#' @param tolerance Numeric: Relative tolerance for mass conservation check
#'   Default: 0.01 (1% error allowed)
#'
#' @return Data frame with validation results:
#'   - Column: Component column name
#'   - n_profiles: Number of profiles checked
#'   - n_conserved: Number of profiles meeting tolerance
#'   - pct_conserved: Percentage of profiles conserving mass
#'   - max_error_pct: Maximum error observed (%)
#'   - mean_error_pct: Mean error across profiles (%)
#'
#' @details
#' **Mass conservation check:**
#' For each profile and component:
#' - Original mass = sum of (horizon_value × horizon_thickness)
#' - Standardized mass = sum of (interval_value × interval_thickness)
#' - Error% = |Original - Standardized| / Original × 100
#'
#' Good mass conservation: error% < 1%
#'
#' @examples
#' \dontrun{
#' # Standardize and validate
#' original <- ssurgo_components
#' standardized <- gc_standardize_depths(original, comp_cols = c("sandtotal", "claytotal"))
#' validation <- gc_validate_depth_standardization(original, standardized)
#' print(validation)
#' }
#'
#' @export
gc_validate_depth_standardization <- function(profiles_original,
                                              profiles_standardized,
                                              comp_cols = NULL, id_col = "mukey",
                                              tolerance = 0.01) {
  stopifnot(is.data.frame(profiles_original), is.data.frame(profiles_standardized))

  # Get comp_cols from standardized if not provided
  if (is.null(comp_cols)) {
    comp_cols <- attr(profiles_standardized, "comp_cols")
    if (is.null(comp_cols)) {
      stop("comp_cols not found in standardized data attributes. Specify explicitly.")
    }
  }

  # Validation results
  results <- list()

  for (col in comp_cols) {
    if (!col %in% names(profiles_original) || !col %in% names(profiles_standardized)) {
      warning("Column '", col, "' not found in both data frames. Skipping.")
      next
    }

    # Get unique profiles
    profiles_ids <- unique(profiles_original[[id_col]])
    n_profiles <- length(profiles_ids)

    errors <- numeric(n_profiles)

    for (i in seq_len(n_profiles)) {
      prof_id <- profiles_ids[i]

      # Original mass (sum of value × thickness)
      orig_prof <- profiles_original[profiles_original[[id_col]] == prof_id, ]
      if (nrow(orig_prof) > 0) {
        thickness <- orig_prof$hzdepb_r - orig_prof$hzdept_r
        orig_mass <- sum(orig_prof[[col]] * thickness)
      } else {
        orig_mass <- NA_real_
      }

      # Standardized mass
      std_prof <- profiles_standardized[profiles_standardized[[id_col]] == prof_id, ]
      if (nrow(std_prof) > 0) {
        thickness <- std_prof$depth_bot - std_prof$depth_top
        std_mass <- sum(std_prof[[col]] * thickness)
      } else {
        std_mass <- NA_real_
      }

      # Error percentage
      if (!is.na(orig_mass) && !is.na(std_mass) && orig_mass != 0) {
        errors[i] <- abs(orig_mass - std_mass) / orig_mass
      } else {
        errors[i] <- NA_real_
      }
    }

    # Summarize
    valid_errors <- errors[!is.na(errors)]
    n_conserved <- sum(valid_errors <= tolerance)

    results[[col]] <- data.frame(
      Column = col,
      n_profiles = n_profiles,
      n_conserved = n_conserved,
      pct_conserved = ifelse(n_profiles > 0, n_conserved / n_profiles * 100, NA),
      max_error_pct = max(valid_errors) * 100,
      mean_error_pct = mean(valid_errors) * 100
    )
  }

  # Combine results
  if (length(results) == 0) {
    return(data.frame())
  }

  validation <- do.call(rbind, results)
  rownames(validation) <- NULL

  return(validation)
}
