#' Validate Conditional Simulation Data Honoring
#'
#' Assess the quality and accuracy of conditional kriging by extracting
#' kriged predictions at the conditioning (observation) locations and comparing
#' them to the actual observed values. This provides a quantitative measure of
#' how well the model honors the conditioning data in ILR space.
#'
#' @param model A gstat kriging model (typically from `gc_ilr_model()`) that was
#'   built with conditioning data via the `data` parameter. The model must contain
#'   observed ILR values at sample locations.
#' @param observed_data A data frame or sf object with the original conditioning
#'   data, containing columns `x`, `y` for spatial coordinates and columns `ilr1`,
#'   `ilr2`, etc., with ILR values in the same space as the model.
#' @param metrics Character vector specifying which error metrics to compute.
#'   Default `c("rmse", "mae", "mean_error")`. Options include:
#'   - `"rmse"`: Root mean squared error
#'   - `"mae"`: Mean absolute error
#'   - `"mean_error"`: Mean signed error (bias)
#'   - `"median_error"`: Median absolute error
#'   - `"sd_error"`: Standard deviation of errors
#'
#' @return A list containing:
#'   - `predictions_at_obs`: Data frame with kriged predictions at each observation location
#'   - `observed_values`: Data frame with original observed ILR values
#'   - `residuals`: Difference (observed - predicted) for each ILR dimension
#'   - `error_metrics`: Data frame with error statistics (one row per ILR dimension)
#'   - `overall_metrics`: Data frame with aggregated error statistics across all dimensions
#'
#' @details
#' **Conditional Honoring Assessment:**
#'
#' Conditional kriging is intended to provide spatial predictions that exactly
#' (or very nearly) match observed values at the conditioning locations.
#' In practice, numerical precision limits and kriging discretization may result
#' in small residuals even for perfectly honored data.
#'
#' This function:
#' 1. Extracts kriged predictions at each observation location using the model
#' 2. Computes residuals (observed minus predicted) for each ILR dimension
#' 3. Calculates error metrics (RMSE, MAE, etc.) per dimension and overall
#' 4. Returns detailed diagnostics for reviewing data honoring quality
#'
#' **Interpretation:**
#' - RMSE close to zero indicates excellent data honoring
#' - Mean error close to zero indicates unbiased predictions
#' - Non-zero residuals may indicate:
#'   - Numerical precision effects (typically < 1e-10)
#'   - Model convergence issues (investigate with univariate kriging)
#'   - Grid resolution effects (predictions on coarse grids may miss exact locations)
#'
#' For compositional data (original sand/silt/clay), back-transform RMSE in ILR
#' space to composition space using `ilrInv()` to assess error in original units.
#'
#' @examples
#' \dontrun{
#' # Assuming model and observed_data are available from gc_ilr_model()
#' # and original data frame with spatial coordinates
#'
#' # Build model with conditioning data
#' model <- gc_ilr_model(ilr_params, data = conditioning_data)
#'
#' # Validate conditioning
#' validation <- gc_validate_conditioning(model, conditioning_data)
#'
#' # Review error metrics
#' print(validation$error_metrics)
#' print(validation$overall_metrics)
#'
#' # Check for problematic observations
#' high_residuals <- which(abs(validation$residuals$ilr1) > 0.1)
#' if (length(high_residuals) > 0) {
#'   print("Observations with large residuals:")
#'   print(validation$residuals[high_residuals, ])
#' }
#'
#' # Visualize residuals vs spatial location
#' plot(validation$predictions_at_obs$x,
#'      validation$residuals$ilr1,
#'      main = "ILR1 Residuals vs X Coordinate",
#'      xlab = "X", ylab = "Residual (Observed - Predicted)")
#' }
#'
#' @importFrom sf st_as_sf st_coordinates st_drop_geometry
#' @importFrom methods is
#' @export
gc_validate_conditioning <- function(model,
                                     observed_data,
                                     metrics = c("rmse", "mae", "mean_error")) {
  requireNamespace("gstat", quietly = TRUE)

  if (methods::is(observed_data, "data.frame")) {
    if (!all(c("x", "y") %in% colnames(observed_data))) {
      stop("observed_data must contain columns 'x' and 'y'")
    }
    obs_sf <- sf::st_as_sf(observed_data, coords = c("x", "y"), crs = "local")
  } else if (methods::is(observed_data, "sf")) {
    obs_sf <- observed_data
  } else {
    stop("observed_data must be a data frame or sf object")
  }

  ilr_cols <- grep("^ilr\\d+$", colnames(obs_sf), value = TRUE)
  if (length(ilr_cols) == 0) {
    stop("observed_data must contain ILR columns (ilr1, ilr2, ...)")
  }

  obs_coords <- sf::st_coordinates(obs_sf)
  obs_data_no_geom <- sf::st_drop_geometry(obs_sf)

  pred_result <- predict(
    model,
    newdata = obs_sf,
    nsim = 0,
    debug.level = 0
  )

  pred_data_no_geom <- sf::st_drop_geometry(pred_result)

  pred_list <- list()
  for (ilr_col in ilr_cols) {
    if (ilr_col %in% colnames(pred_data_no_geom)) {
      pred_list[[ilr_col]] <- pred_data_no_geom[[ilr_col]]
    } else {
      pred_list[[ilr_col]] <- rep(NA_real_, nrow(pred_data_no_geom))
    }
  }
  
  if (length(pred_list) > 0) {
    pred_values <- as.data.frame(pred_list)
  } else {
    pred_values <- data.frame(matrix(NA_real_, nrow = nrow(pred_data_no_geom), ncol = length(ilr_cols)))
    colnames(pred_values) <- ilr_cols
  }
  
  if (nrow(pred_values) != nrow(obs_sf)) {
    warning(
      "Prediction result has ", nrow(pred_values), " rows but observed data has ",
      nrow(obs_sf), " rows. Predictions may be misaligned."
    )
  }

  obs_list <- list()
  for (ilr_col in ilr_cols) {
    if (ilr_col %in% colnames(obs_data_no_geom)) {
      obs_list[[ilr_col]] <- obs_data_no_geom[[ilr_col]]
    } else {
      obs_list[[ilr_col]] <- rep(NA_real_, nrow(obs_data_no_geom))
    }
  }
  obs_values <- as.data.frame(obs_list)

  residuals <- obs_values - pred_values
  residuals$x <- obs_coords[, 1]
  residuals$y <- obs_coords[, 2]

  error_metrics <- data.frame(
    ilr_dimension = ilr_cols,
    stringsAsFactors = FALSE
  )

  for (metric in metrics) {
    metric_lower <- tolower(metric)
    if (metric_lower == "rmse") {
      error_metrics$RMSE <- sqrt(colMeans(residuals[, ilr_cols]^2, na.rm = TRUE))
    } else if (metric_lower == "mae") {
      error_metrics$MAE <- colMeans(abs(residuals[, ilr_cols]), na.rm = TRUE)
    } else if (metric_lower == "mean_error") {
      error_metrics$Mean_Error <- colMeans(residuals[, ilr_cols], na.rm = TRUE)
    } else if (metric_lower == "median_error") {
      error_metrics$Median_Error <- apply(
        residuals[, ilr_cols],
        2,
        function(x) median(abs(x), na.rm = TRUE)
      )
    } else if (metric_lower == "sd_error") {
      error_metrics$SD_Error <- apply(
        residuals[, ilr_cols],
        2,
        function(x) sd(x, na.rm = TRUE)
      )
    } else {
      warning("Unknown metric: ", metric)
    }
  }

  overall_resid <- unlist(residuals[, ilr_cols])
  overall_metrics <- data.frame(
    metric = character(0),
    value = numeric(0),
    stringsAsFactors = FALSE
  )

  if ("rmse" %in% tolower(metrics)) {
    overall_metrics <- rbind(
      overall_metrics,
      data.frame(metric = "RMSE", value = sqrt(mean(overall_resid^2, na.rm = TRUE)))
    )
  }
  if ("mae" %in% tolower(metrics)) {
    overall_metrics <- rbind(
      overall_metrics,
      data.frame(metric = "MAE", value = mean(abs(overall_resid), na.rm = TRUE))
    )
  }
  if ("mean_error" %in% tolower(metrics)) {
    overall_metrics <- rbind(
      overall_metrics,
      data.frame(metric = "Mean_Error", value = mean(overall_resid, na.rm = TRUE))
    )
  }
  if ("median_error" %in% tolower(metrics)) {
    overall_metrics <- rbind(
      overall_metrics,
      data.frame(metric = "Median_Error", value = median(abs(overall_resid), na.rm = TRUE))
    )
  }
  if ("sd_error" %in% tolower(metrics)) {
    overall_metrics <- rbind(
      overall_metrics,
      data.frame(metric = "SD_Error", value = sd(overall_resid, na.rm = TRUE))
    )
  }

  pred_at_obs <- cbind(
    data.frame(x = obs_coords[, 1], y = obs_coords[, 2]),
    pred_values
  )

  return(list(
    predictions_at_obs = pred_at_obs,
    observed_values = cbind(
      data.frame(x = obs_coords[, 1], y = obs_coords[, 2]),
      obs_values
    ),
    residuals = residuals,
    error_metrics = error_metrics,
    overall_metrics = overall_metrics,
    n_observations = nrow(obs_sf),
    n_ilr_dimensions = length(ilr_cols)
  ))
}
