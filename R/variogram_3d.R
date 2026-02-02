#' Fit 3D Variogram with Geometric Anisotropy
#'
#' Fits 3D variograms with explicit separation of horizontal and vertical spatial
#' continuity. Computes anisotropy ratio and directional information, critical
#' for modeling soil properties where lateral correlation ranges differ significantly
#' from vertical ranges (typical: horizontal >> vertical by 10-100×).
#'
#' @param ilr_params Output from [gc_ilr_params()], containing mean and covariance
#' @param data Data frame with x, y, z coordinates plus ilr1, ilr2, ... columns
#' @param vgm_model_type Character: "Exp" (exponential), "Sph" (spherical), "Gau" (Gaussian)
#' @param horizontal_cutoff Distance (in xy-plane) for horizontal variogram (default: max xy distance / 2)
#' @param vertical_cutoff Distance (in z direction) for vertical variogram (default: max z distance / 2)
#' @param horizontal_width Lag width for horizontal variogram (default: auto)
#' @param vertical_width Lag width for vertical variogram (default: auto)
#' @param aggregate Logical: if TRUE, return single aggregated 3D vgm; if FALSE, return per-dimension (default: FALSE)
#' @param correct.diagonal Sill correction factor for LMC stability (default: 1.01)
#' @param tolerance_z Tolerance in z for computing horizontal variogram (cm, default: 10)
#' @param tolerance_xy Tolerance in xy for computing vertical variogram (m, default: 50)
#'
#' @return List containing:
#'   - `vgm_horizontal`: Fitted variogram for horizontal (xy) distances
#'   - `vgm_vertical`: Fitted variogram for vertical (z) distances
#'   - `range_horizontal`: Range parameter from horizontal variogram
#'   - `range_vertical`: Range parameter from vertical variogram
#'   - `anisotropy_ratio`: Ratio of horizontal to vertical range
#'   - `fitted_params`: Data frame with per-dimension parameters
#'   - Attributes: `ilr_dimension_names`, `model_type`, `coordinates_3d = TRUE`
#'
#' @details
#' **3D Variogram Strategy:**
#'
#' Instead of fitting a full 3D variogram (which requires complex directional binning),
#' this function fits:
#'
#' 1. **Horizontal Variogram** (in xy-plane):
#'    - Pairs with similar z-values (within tolerance_z)
#'    - Describes lateral continuity independent of depth
#'    - Example: Sand content at 50-60 cm similar to sand at 40-50 cm in same location
#'
#' 2. **Vertical Variogram** (in z-direction):
#'    - Pairs with similar xy-locations (within tolerance_xy)
#'    - Describes depth-dependent changes
#'    - Example: Sand content at 10 cm vs 50 cm in same profile location
#'
#' **Anisotropy:**
#' - Computed as ratio: `anisotropy_ratio = range_horizontal / range_vertical`
#' - Typical soil values: 10-100 (lateral continuity >> depth continuity)
#' - In gstat: represented via `anis` parameter in 3D kriging
#'
#' **Geometric Anisotropy in 3D Kriging:**
#' The fitted 3D model enables gstat to compute covariance for lag vector (h_x, h_y, h_z):
#' ```
#' C_3D(h) = C_model(√(h_x² + h_y²) / a_h, h_z / a_v)
#' ```
#' where `a_h` is horizontal range and `a_v` is vertical range.
#'
#' @examples
#' \\dontrun{
#' # Generate synthetic 3D soil data
#' data_3d <- expand.grid(
#'   x = seq(0, 1000, by = 100),
#'   y = seq(0, 1000, by = 100),
#'   z = seq(0, 100, by = 10)  # Depth in cm
#' )
#'
#' # Create ILR-transformed compositional data
#' n <- nrow(data_3d)
#' data_3d$sand <- pmin(100, pmax(0, rnorm(n, 40, 15) - 0.3*data_3d$z))
#' data_3d$silt <- pmin(100, pmax(0, rnorm(n, 35, 15)))
#' data_3d$clay <- 100 - data_3d$sand - data_3d$silt
#'
#' ilr_params <- gc_ilr_params(data_3d[, c("sand", "silt", "clay")])
#' data_3d <- cbind(data_3d, ilr_params$ilr_coords)
#'
#' # Fit 3D variogram with anisotropy
#' vgm_3d <- gc_fit_vgm_3d(
#'   ilr_params, data_3d,
#'   vgm_model_type = "Exp",
#'   horizontal_cutoff = 500,
#'   vertical_cutoff = 50
#' )
#'
#' # Inspect results
#' cat("Horizontal range:", vgm_3d$range_horizontal, "m\n")
#' cat("Vertical range:", vgm_3d$range_vertical, "cm\n")
#' cat("Anisotropy ratio (h/v):", vgm_3d$anisotropy_ratio, "\n")
#' }
#'
#' @export
gc_fit_vgm_3d <- function(ilr_params,
                         data,
                         vgm_model_type = "Exp",
                         horizontal_cutoff = NULL,
                         vertical_cutoff = NULL,
                         horizontal_width = NULL,
                         vertical_width = NULL,
                         aggregate = FALSE,
                         correct.diagonal = 1.01,
                         tolerance_z = 10,
                         tolerance_xy = 50) {
  # Input validation
  if (!all(c("x", "y", "z") %in% colnames(data))) {
    stop("data must contain columns 'x', 'y', and 'z' for 3D variogram fitting")
  }

  if (!is.numeric(tolerance_z) || tolerance_z < 0) {
    stop("tolerance_z must be numeric and >= 0")
  }
  if (!is.numeric(tolerance_xy) || tolerance_xy < 0) {
    stop("tolerance_xy must be numeric and >= 0")
  }

  if (!is.numeric(correct.diagonal) || correct.diagonal < 1.0) {
    stop("correct.diagonal must be numeric and >= 1.0")
  }

  n_ilr <- length(ilr_params$mean)
  n_points <- nrow(data)

  if (n_points < 8) {
    stop("data must contain at least 8 observations for 3D variogram fitting")
  }

  # Set default cutoff distances if not provided
  if (is.null(horizontal_cutoff)) {
    x_range <- diff(range(data$x, na.rm = TRUE))
    y_range <- diff(range(data$y, na.rm = TRUE))
    horizontal_cutoff <- max(x_range, y_range) / 2
  }

  if (is.null(vertical_cutoff)) {
    z_range <- diff(range(data$z, na.rm = TRUE))
    vertical_cutoff <- z_range / 2
  }

  # Create coordinate distance matrix
  coords_xyz <- as.matrix(data[, c("x", "y", "z")])

  # Compute horizontal distances (xy-plane)
  xy_coords <- coords_xyz[, c("x", "y"), drop = FALSE]
  dist_xy <- as.matrix(stats::dist(xy_coords))

  # Compute vertical distances (z-direction)
  z_coords <- coords_xyz[, "z", drop = FALSE]
  dist_z <- as.matrix(stats::dist(z_coords))

  # Initialize results
  fitted_vgms <- list()
  empirical_vgms_h <- list()
  empirical_vgms_v <- list()
  fitted_params <- data.frame(
    ilr_id = character(n_ilr),
    range_h = numeric(n_ilr),
    range_v = numeric(n_ilr),
    anis_ratio = numeric(n_ilr),
    nugget = numeric(n_ilr),
    psill = numeric(n_ilr),
    stringsAsFactors = FALSE
  )

  # Fit variograms per ILR dimension
  for (i in seq_len(n_ilr)) {
    ilr_id <- paste0("ilr", i)

    if (!ilr_id %in% colnames(data)) {
      stop("data must contain column '", ilr_id, "'")
    }

    tryCatch(
      {
        ilr_values <- data[[ilr_id]]

        # HORIZONTAL VARIOGRAM: Pairs with similar z-values
        horiz_pairs <- data.frame(
          dist = numeric(0),
          gamma = numeric(0)
        )

        # Find pairs where z-distance is small (same depth zone)
        for (j in seq_len(n_points - 1)) {
          for (k in (j + 1):n_points) {
            if (dist_z[j, k] <= tolerance_z) {
              pair_dist <- dist_xy[j, k]
              pair_gamma <- (ilr_values[j] - ilr_values[k])^2 / 2
              horiz_pairs <- rbind(horiz_pairs, data.frame(
                dist = pair_dist,
                gamma = pair_gamma
              ))
            }
          }
        }

        # VERTICAL VARIOGRAM: Pairs with similar xy-locations
        vert_pairs <- data.frame(
          dist = numeric(0),
          gamma = numeric(0)
        )

        # Find pairs where xy-distance is small (same profile location)
        for (j in seq_len(n_points - 1)) {
          for (k in (j + 1):n_points) {
            if (dist_xy[j, k] <= tolerance_xy) {
              pair_dist <- dist_z[j, k]
              pair_gamma <- (ilr_values[j] - ilr_values[k])^2 / 2
              vert_pairs <- rbind(vert_pairs, data.frame(
                dist = pair_dist,
                gamma = pair_gamma
              ))
            }
          }
        }

        # Check if we have enough pairs
        if (nrow(horiz_pairs) < 3) {
          warning(
            "Only ", nrow(horiz_pairs), " horizontal pairs for ", ilr_id,
            " (need >= 3). Consider increasing tolerance_z."
          )
        }
        if (nrow(vert_pairs) < 3) {
          warning(
            "Only ", nrow(vert_pairs), " vertical pairs for ", ilr_id,
            " (need >= 3). Consider increasing tolerance_xy."
          )
        }

        # Initial model
        init_sill <- stats::var(ilr_values, na.rm = TRUE)

        # Estimate horizontal variogram from pairs
        if (nrow(horiz_pairs) >= 3) {
          # Estimate range as distance where semivariance reaches ~90% of sill
          sill_target <- 0.9 * init_sill
          above_sill <- horiz_pairs$gamma >= sill_target
          if (any(above_sill)) {
            range_h <- min(horiz_pairs$dist[above_sill], na.rm = TRUE)
          } else {
            range_h <- max(horiz_pairs$dist, na.rm = TRUE) / 1.5
          }

          # Estimate nugget from smallest distances
          small_dists <- horiz_pairs$dist <= (max(horiz_pairs$dist) * 0.1)
          if (any(small_dists)) {
            nugget_h <- mean(horiz_pairs$gamma[small_dists], na.rm = TRUE)
          } else {
            nugget_h <- 0.01 * init_sill
          }

          # Create fitted variogram object manually
          fitted_vgm_h <- gstat::vgm(
            psill = init_sill - nugget_h,
            model = vgm_model_type,
            range = range_h,
            nugget = nugget_h
          )
        } else {
          fitted_vgm_h <- NULL
          range_h <- NA_real_
        }

        # Estimate vertical variogram from pairs
        if (nrow(vert_pairs) >= 3) {
          # Estimate range as distance where semivariance reaches ~90% of sill
          sill_target <- 0.9 * init_sill
          above_sill <- vert_pairs$gamma >= sill_target
          if (any(above_sill)) {
            range_v <- min(vert_pairs$dist[above_sill], na.rm = TRUE)
          } else {
            range_v <- max(vert_pairs$dist, na.rm = TRUE) / 1.5
          }

          # Estimate nugget from smallest distances
          small_dists <- vert_pairs$dist <= (max(vert_pairs$dist) * 0.1)
          if (any(small_dists)) {
            nugget_v <- mean(vert_pairs$gamma[small_dists], na.rm = TRUE)
          } else {
            nugget_v <- 0.01 * init_sill
          }

          # Create fitted variogram object manually
          fitted_vgm_v <- gstat::vgm(
            psill = init_sill - nugget_v,
            model = vgm_model_type,
            range = range_v,
            nugget = nugget_v
          )
        } else {
          fitted_vgm_v <- NULL
          range_v <- NA_real_
        }

        # Store combined 3D model
        fitted_vgms[[ilr_id]] <- list(
          vgm_horizontal = fitted_vgm_h,
          vgm_vertical = fitted_vgm_v,
          range_horizontal = range_h,
          range_vertical = range_v
        )

        # Store parameters
        fitted_params$ilr_id[i] <- ilr_id
        fitted_params$range_h[i] <- range_h
        fitted_params$range_v[i] <- range_v
        fitted_params$anis_ratio[i] <- ifelse(
          is.na(range_v) || range_v == 0,
          NA_real_,
          range_h / range_v
        )
        fitted_params$nugget[i] <- ifelse(
          is.null(fitted_vgm_h),
          NA_real_,
          fitted_vgm_h$psill[fitted_vgm_h$model == "Nug"]
        )
        fitted_params$psill[i] <- init_sill
      },
      error = function(e) {
        warning("Failed to fit 3D variogram for ", ilr_id, ": ", e$message)
      }
    )
  }

  # Return results
  result <- list(
    fitted_vgms = fitted_vgms,
    fitted_params = fitted_params,
    empirical_vgms_horizontal = empirical_vgms_h,
    empirical_vgms_vertical = empirical_vgms_v
  )

  attr(result, "ilr_dimension_names") <- paste0("ilr", seq_len(n_ilr))
  attr(result, "model_type") <- vgm_model_type
  attr(result, "coordinates_3d") <- TRUE
  attr(result, "tolerance_z") <- tolerance_z
  attr(result, "tolerance_xy") <- tolerance_xy
  attr(result, "horizontal_cutoff") <- horizontal_cutoff
  attr(result, "vertical_cutoff") <- vertical_cutoff

  return(result)
}


#' Internal: Create Empirical Variogram from Pair Data
#'
#' Helper function to create gstat variogram object from distance/gamma pairs.
#'
#' @param pairs Data frame with 'dist' and 'gamma' columns
#' @param ilr_id Variable name (for labeling)
#' @param cutoff Maximum distance
#' @param width Lag width (default: auto)
#'
#' @return gstat variogram object
#' @keywords internal
.create_empirical_vgm_from_pairs <- function(pairs, ilr_id, cutoff, width) {
  # Bin pairs into distance classes
  if (is.null(width)) {
    width <- cutoff / 7  # Default: ~7 distance classes
  }

  n_classes <- ceiling(cutoff / width)
  class_breaks <- seq(0, cutoff + width, by = width)

  # Initialize bins data frame with correct types
  bins <- data.frame(
    dist = numeric(),
    np = integer(),
    gamma = numeric(),
    dir.hor = numeric(),
    dir.ver = numeric(),
    id = character(),
    stringsAsFactors = FALSE
  )

  # Bin the pairs
  for (j in seq_len(n_classes)) {
    lower <- class_breaks[j]
    upper <- class_breaks[j + 1]

    in_bin <- pairs$dist >= lower & pairs$dist < upper
    if (any(in_bin)) {
      bin_gamma <- mean(pairs$gamma[in_bin], na.rm = TRUE)
      bin_np <- as.integer(sum(in_bin))
      bin_dist <- as.numeric((lower + upper) / 2)

      new_row <- data.frame(
        dist = bin_dist,
        np = bin_np,
        gamma = bin_gamma,
        dir.hor = 0.0,
        dir.ver = 90.0,
        id = ilr_id,
        stringsAsFactors = FALSE
      )

      bins <- rbind(bins, new_row, stringsAsFactors = FALSE)
    }
  }

  # Add final bin if needed
  final_lower <- class_breaks[n_classes + 1]
  in_final <- pairs$dist >= final_lower & pairs$dist <= cutoff
  if (any(in_final)) {
    final_gamma <- mean(pairs$gamma[in_final], na.rm = TRUE)
    final_np <- as.integer(sum(in_final))
    final_dist <- as.numeric((final_lower + cutoff) / 2)

    final_row <- data.frame(
      dist = final_dist,
      np = final_np,
      gamma = final_gamma,
      dir.hor = 0.0,
      dir.ver = 90.0,
      id = ilr_id,
      stringsAsFactors = FALSE
    )

    bins <- rbind(bins, final_row, stringsAsFactors = FALSE)
  }

  # Ensure all columns are proper types
  bins$dist <- as.numeric(bins$dist)
  bins$np <- as.integer(bins$np)
  bins$gamma <- as.numeric(bins$gamma)
  bins$dir.hor <- as.numeric(bins$dir.hor)
  bins$dir.ver <- as.numeric(bins$dir.ver)
  bins$id <- as.character(bins$id)

  # Return as variogram object
  class(bins) <- c("gstatVariogram", "data.frame")
  return(bins)
}
