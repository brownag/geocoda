#' Convert Variogram Model to Covariance at Specific Distance
#'
#' Evaluate covariance C(h) from a fitted variogram model at specified distances.
#' Uses the relationship: C(h) = C(0) - γ(h), where C(0) is the sill and γ(h) is the variogram.
#'
#' @param vgm_object A gstat vgm object (fitted variogram model)
#' @param distance Scalar distance at which to evaluate covariance
#'
#' @return Numeric covariance value C(h)
#'
#' @importFrom gstat variogramLine
#'
#' @details
#' The total sill C(0) = sum of psill values + nugget effect.
#' For a variogram model with nugget and structure psills, this recovers the
#' point-lag covariance that is needed for MAF computations.
#'
#' @examples
#' \dontrun{
#' # Fit a variogram
#' vgm_fit <- gstat::fit.variogram(
#'   gstat::variogram(log(zinc) ~ 1, meuse),
#'   gstat::vgm(1, "Exp", 300, 0.1)
#' )
#'
#' # Evaluate covariance at lag 100
#' C_100 <- gc_vgm_to_cov(vgm_fit, distance = 100)
#' }
#'
#' @export
gc_vgm_to_cov <- function(vgm_object, distance) {
  stopifnot(is.numeric(distance), length(distance) == 1, distance >= 0)

  # Total sill: sum of all psill values (includes nugget as first row's psill)
  C_0 <- sum(vgm_object$psill)

  # Evaluate variogram at this distance using gstat::variogramLine
  vgm_line <- gstat::variogramLine(vgm_object, dist_vector = distance)

  # Extract gamma at this distance (first element of returned data frame)
  gamma_h <- vgm_line$gamma[1]

  # Covariance: C(h) = C(0) - γ(h)
  C_h <- C_0 - gamma_h

  return(as.numeric(C_h))
}

#' Build Multivariate Covariance Matrix at Specified Lag
#'
#' Construct the spatial covariance matrix between all ILR dimensions at distance lag_h.
#' Uses univariate covariances on diagonal and simplified approximation for off-diagonal terms.
#'
#' @param ilr_params List from [gc_ilr_params()] containing mean, cov, names
#' @param variogram_models List of fitted gstat vgm objects (one per ILR dimension)
#' @param lag_h Distance at which to evaluate covariances
#'
#' @return D×D covariance matrix at lag h, where D is the number of ILR dimensions
#'
#' @details
#' **Diagonal terms** C_ii(h) are computed directly from the univariate variogram
#' using [gc_vgm_to_cov()].
#'
#' **Off-diagonal terms** C_ij(h) are approximated using:
#' ```
#' C_ij(h) ≈ C_ij(0) * sqrt(C_ii(h)/C_ii(0)) * sqrt(C_jj(h)/C_jj(0))
#' ```
#'
#' This assumes that the correlation structure between dimensions is relatively stable
#' across lags (i.e., dimensions that are correlated at lag 0 remain similarly correlated
#' at lag h). This is a simplification that works well when all dimensions have similar
#' spatial structure, but may underestimate cross-covariance if dimensions have very
#' different spatial ranges.
#'
#' @examples
#' \dontrun{
#' # ILR parameters from soil texture data
#' ilr_params <- gc_ilr_params(texture_data[, c("sand", "silt", "clay")])
#'
#' # Fit variograms per dimension
#' vgm_list <- list(
#'   gc_fit_vgm(ilr_params, texture_data, aggregate = TRUE),
#'   gc_fit_vgm(ilr_params, texture_data, aggregate = TRUE)
#' )
#'
#' # Covariance matrix at lag 50 (meters)
#' C_50 <- gc_covariance_matrix(ilr_params, vgm_list, lag_h = 50)
#' }
#'
#' @export
gc_covariance_matrix <- function(ilr_params, variogram_models, lag_h) {
  stopifnot(is.list(ilr_params), !is.null(ilr_params$cov))
  stopifnot(is.list(variogram_models), length(variogram_models) > 0)
  stopifnot(is.numeric(lag_h), length(lag_h) == 1, lag_h >= 0)

  n_ilr <- ncol(ilr_params$cov)

  # Initialize covariance matrix
  C_h <- matrix(0, n_ilr, n_ilr)
  colnames(C_h) <- rownames(C_h) <- paste0("ilr", seq_len(n_ilr))

  # Fill in diagonal: univariate covariances at lag h
  for (i in seq_len(n_ilr)) {
    C_h[i, i] <- gc_vgm_to_cov(variogram_models[[i]], distance = lag_h)
  }

  # Fill in off-diagonal: approximated cross-covariances
  for (i in seq_len(n_ilr - 1)) {
    for (j in (i + 1):n_ilr) {
      # Collocated (lag 0) covariance
      C_ij_0 <- ilr_params$cov[i, j]

      # Total sills (covariance at lag 0)
      C_ii_0 <- ilr_params$cov[i, i]
      C_jj_0 <- ilr_params$cov[j, j]

      # Covariance at lag h
      C_ii_h <- C_h[i, i]
      C_jj_h <- C_h[j, j]

      # Approximation: preserve correlation structure
      # C_ij(h) ≈ C_ij(0) * sqrt(C_ii(h)/C_ii(0)) * sqrt(C_jj(h)/C_jj(0))
      correlation_factor <- sqrt(C_ii_h / C_ii_0) * sqrt(C_jj_h / C_jj_0)
      C_h[i, j] <- C_h[j, i] <- C_ij_0 * correlation_factor
    }
  }

  return(C_h)
}
