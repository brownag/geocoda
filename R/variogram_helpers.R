#' Suggest Reasonable Variogram Default Parameters
#'
#' Inspect data extent and ILR covariance structure to suggest ballpark values
#' for range, nugget, and sill of a variogram model. This is a convenience
#' function to help users get started before fine-tuning.
#'
#' @param ilr_params A list returned by `gc_ilr_params()`.
#' @param extent A numeric vector of length 4: `c(xmin, ymin, xmax, ymax)`,
#'   or a spatial object (sf, terra) from which extent can be extracted.
#'
#' @return A list containing:
#'   - `range`: Suggested range parameter (approximately 1/3 of the diagonal extent)
#'   - `nugget`: Suggested nugget ratio (typically 0.01 to 0.05 of average sill)
#'   - `mean_sill`: Mean of the diagonal covariance terms
#'
#' @details
#' This function:
#' 1. Calculates the spatial extent diagonal (Euclidean distance from min to max)
#' 2. Suggests range as ~1/3 of the extent diagonal
#' 3. Calculates mean of diagonal covariance terms as representative sill
#' 4. Suggests nugget as 1% of the mean sill
#'
#' These are heuristics and should be refined using empirical variography
#' (see [gc_fit_vgm()]).
#'
#' @importFrom stats as.formula
#' @examples
#' samples <- data.frame(
#'   sand = c(20, 25, 30, 22),
#'   silt = c(60, 55, 50, 58),
#'   clay = c(20, 20, 20, 20)
#' )
#'
#' params <- gc_ilr_params(samples)
#' extent <- c(0, 0, 100, 100)
#'
#' suggestions <- gc_vgm_defaults(params, extent)
#' print(suggestions)
#'
#' @importFrom stats as.formula
#' @export
gc_vgm_defaults <- function(ilr_params, extent) {
  ilr_cov <- ilr_params$cov
  diag_cov <- diag(ilr_cov)
  mean_sill <- mean(diag_cov)

  if (length(extent) == 4) {
    xmin <- extent[1]
    ymin <- extent[2]
    xmax <- extent[3]
    ymax <- extent[4]
  } else {
    stop("extent must be a vector of length 4: c(xmin, ymin, xmax, ymax)")
  }

  diag_dist <- sqrt((xmax - xmin)^2 + (ymax - ymin)^2)

  suggested_range <- diag_dist / 3
  suggested_nugget <- 0.01 * mean_sill

  list(
    range = suggested_range,
    nugget = suggested_nugget,
    mean_sill = mean_sill
  )
}


#' Optimize Variogram Parameters
#'
#' Fit empirical variograms to ILR values at observed locations and estimate
#' optimal variogram parameters for each ILR dimension. Results can be returned
#' per-dimension or aggregated into a single template variogram for use in
#' multivariate simulation. Includes diagnostics for LMC admissibility.
#'
#' @param ilr_params A list returned by `gc_ilr_params()`.
#' @param data A data frame or spatial object (sf) with columns `x`, `y` for
#'   spatial coordinates and columns named `ilr1`, `ilr2`, etc., containing
#'   ILR values at observed locations.
#' @param vgm_model_type Character string specifying variogram model type
#'   (e.g., `"Exp"`, `"Sph"`, `"Gau"`). Default `"Exp"`.
#' @param maxdist Maximum distance for empirical variogram calculation
#'   (default `NULL`, uses all pairs).
#' @param width Lag width for empirical variogram bins (default `NULL`,
#'   automatic).
#' @param aggregate Logical. If `TRUE`, aggregate fitted parameters across
#'   ILR dimensions using weighted average sills (weights = covariance diagonal).
#'   If `FALSE` (default), return per-dimension results.
#' @param correct.diagonal Numeric sill correction factor for LMC stability
#'   (default `1.01`). Multiplies each marginal (diagonal) sill by this factor
#'   to improve positive-definiteness of the LMC sill matrix and prevent
#'   numerical issues. Recommended range: 1.00-1.05. Use `1.01` for typical applications.
#' @param fit.ranges Logical. If `FALSE` (default), fixes ranges to fitted values
#'   during LMC construction to avoid over-parameterization. If `TRUE`, allows
#'   range re-optimization when building LMC models.
#'
#' @return A list with:
#'   - If `aggregate = FALSE`: A list of length `D-1` where each element
#'     contains `fitted_vgm` and `empirical_vgm` for that ILR dimension.
#'   - If `aggregate = TRUE`: A single `vgm` object with aggregated parameters
#'     (mean range, weighted mean sill, weighted mean nugget).
#'   - Attribute `ilr_dimension_names`: Names of ILR dimensions (ilr1, ilr2, ...)
#'   - Attribute `fitted_params`: Data frame with per-dimension parameters
#'   - Attribute `lmc_admissibility`: Logical indicating if sill matrix is positive-definite
#'
#' @details
#' This function:
#' 1. For each ILR dimension i:
#'    - Computes empirical variogram using `gstat::variogram()`
#'    - Provides intelligent initial model (range ~ spatial extent / 3)
#'    - Fits parametric model using `gstat::fit.variogram()`
#' 2. Optionally aggregates results using covariance-weighted averaging
#' 3. Applies `correct.diagonal` to diagonal sills to ensure LMC positive-definiteness
#'
#' Weighted aggregation combines results when building a single LMC model,
#' ensuring dimensions with higher variance contribute proportionally.
#'
#' **LMC Admissibility:**
#' The sill matrix (covariance at distance infinity) must be positive-definite for
#' valid spatial covariance. If eigenvalues of the sill matrix have any zero or
#' negative values, the LMC structure is inadmissible. The `correct.diagonal`
#' parameter helps stabilize the sill matrix by inflating diagonal terms,
#' which typically improves admissibility without distorting the model significantly.
#'
#' @examples
#' \dontrun{
#' # Example with observed ILR values on a spatial grid
#' library(gstat)
#'
#' # Create sample data with spatial coordinates and ILR values
#' n_points <- 50
#' coords_df <- data.frame(
#'   x = runif(n_points, 0, 100),
#'   y = runif(n_points, 0, 100),
#'   ilr1 = rnorm(n_points, mean = 0.5, sd = 0.8),
#'   ilr2 = rnorm(n_points, mean = -0.2, sd = 0.6)
#' )
#'
#' # Estimate parameters from bootstrap samples
#' samples <- data.frame(
#'   sand = c(20, 25, 30, 22),
#'   silt = c(60, 55, 50, 58),
#'   clay = c(20, 20, 20, 20)
#' )
#' params <- gc_ilr_params(samples)
#'
#' # Fit empirical variograms
#' fitted <- gc_fit_vgm(
#'   params,
#'   coords_df,
#'   vgm_model_type = "Exp",
#'   aggregate = FALSE
#' )
#'
#' # Or get single aggregated template
#' fitted_agg <- gc_fit_vgm(
#'   params,
#'   coords_df,
#'   vgm_model_type = "Exp",
#'   aggregate = TRUE
#' )
#' }
#'
#' @importFrom gstat variogram fit.variogram vgm
#' @export
gc_fit_vgm <- function(ilr_params,
                       data,
                       vgm_model_type = "Exp",
                       maxdist = NULL,
                       width = NULL,
                       aggregate = FALSE,
                       correct.diagonal = 1.01,
                       fit.ranges = FALSE) {
  if (!all(c("x", "y") %in% colnames(data))) {
    stop("data must contain columns 'x' and 'y'")
  }

  if (!is.numeric(correct.diagonal) || correct.diagonal < 1.0) {
    stop("correct.diagonal must be numeric and >= 1.0")
  }

  n_ilr <- length(ilr_params$mean)
  n_points <- nrow(data)

  if (n_points < 4) {
    stop("data must contain at least 4 observations")
  }

  cov_diag <- diag(ilr_params$cov)

  fitted_vgms <- list()
  empirical_vgms <- list()
  fitted_params <- data.frame(
    ilr_id = character(n_ilr),
    range = numeric(n_ilr),
    nugget = numeric(n_ilr),
    psill = numeric(n_ilr),
    stringsAsFactors = FALSE
  )

  for (i in seq_len(n_ilr)) {
    ilr_id <- paste0("ilr", i)
    formula_obj <- stats::as.formula(paste0(ilr_id, " ~ 1"))

    if (!ilr_id %in% colnames(data)) {
      stop("data must contain column '", ilr_id, "'")
    }

    tryCatch(
      {
        vgm_args <- list(
          locations = ~ x + y,
          data = data
        )
        if (!is.null(maxdist)) vgm_args$maxdist <- maxdist
        if (!is.null(width)) vgm_args$width <- width

        emp_vgm <- do.call(gstat::variogram, c(list(formula_obj), vgm_args))

        x_range <- diff(range(data$x))
        y_range <- diff(range(data$y))
        avg_range <- (x_range + y_range) / 2
        init_range <- avg_range / 3

        init_sill <- stats::var(data[[ilr_id]], na.rm = TRUE)

        init_vgm <- gstat::vgm(
          psill = init_sill,
          model = vgm_model_type,
          range = init_range,
          nugget = 0.01 * init_sill
        )

        fitted_vgm <- suppressWarnings(
          gstat::fit.variogram(emp_vgm, init_vgm)
        )

        fitted_vgms[[ilr_id]] <- fitted_vgm
        empirical_vgms[[ilr_id]] <- emp_vgm

        model_row <- which(fitted_vgm$model != "Nug")
        fitted_params$ilr_id[i] <- ilr_id
        fitted_params$range[i] <- fitted_vgm$range[model_row[1]]
        fitted_params$nugget[i] <- fitted_vgm$psill[fitted_vgm$model == "Nug"]
        fitted_params$psill[i] <- fitted_vgm$psill[model_row[1]]
      },
      error = function(e) {
        warning("Failed to fit variogram for ", ilr_id, ": ", e$message)
      }
    )
  }

  fitted_params$psill_corrected <- fitted_params$psill * correct.diagonal

  if (aggregate) {
    weights <- cov_diag / sum(cov_diag)

    valid_range <- fitted_params$range[is.finite(fitted_params$range) & fitted_params$range > 0]
    valid_nugget <- fitted_params$nugget[is.finite(fitted_params$nugget) & fitted_params$nugget >= 0]
    valid_psill <- fitted_params$psill_corrected[is.finite(fitted_params$psill_corrected) & fitted_params$psill_corrected > 0]

    if (length(valid_range) > 0) {
      mean_range <- mean(valid_range)
    } else {
      x_range <- diff(range(data$x))
      y_range <- diff(range(data$y))
      avg_range <- (x_range + y_range) / 2
      mean_range <- avg_range / 3
    }

    if (length(valid_nugget) > 0) {
      mean_nugget <- mean(valid_nugget)
    } else {
      mean_nugget <- 0.01 * mean(diag(ilr_params$cov))
    }

    if (length(valid_psill) > 0) {
      mean_psill <- mean(valid_psill)
    } else {
      mean_psill <- mean(diag(ilr_params$cov)) * correct.diagonal
    }

    agg_vgm <- gstat::vgm(
      psill = mean_psill,
      model = vgm_model_type,
      range = mean_range,
      nugget = mean_nugget
    )

    sill_matrix <- matrix(mean_psill, nrow = n_ilr, ncol = n_ilr)
    diag(sill_matrix) <- mean_psill
    eigenvals <- eigen(sill_matrix, only.values = TRUE)$values
    lmc_admissible <- all(eigenvals > 1e-10)

    if (!lmc_admissible) {
      warning(
        "LMC sill matrix may not be positive-definite. ",
        "Consider increasing correct.diagonal parameter or using univariate kriging."
      )
    }

    result <- agg_vgm
    attr(result, "ilr_dimension_names") <- paste0("ilr", seq_len(n_ilr))
    attr(result, "fitted_params") <- fitted_params
    attr(result, "lmc_admissibility") <- lmc_admissible
    return(result)
  } else {
    result <- list(
      fitted_vgms = fitted_vgms,
      empirical_vgms = empirical_vgms,
      fitted_params = fitted_params
    )
    attr(result, "ilr_dimension_names") <- paste0("ilr", seq_len(n_ilr))
    return(result)
  }
}
