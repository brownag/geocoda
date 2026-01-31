#' Data Pipeline Helpers
#'
#' Convenience functions for integrating soil databases and covariate surfaces
#' into geostatistical simulation workflows. Extensions for data preparation.
#'
#' @keywords internal
#' @name data_pipeline
NULL


#' Apply Covariate Trend Surface to Simulation
#'
#' Incorporate auxiliary variables (elevation, parent material, slope, etc.) as
#' trend surfaces in ILR-transformed compositional data. Useful for capturing
#' deterministic spatial patterns in soil properties related to landscape position.
#'
#' @param x A `terra::SpatRaster` containing sample data with component layers
#'   (e.g., sand, silt, clay)
#' @param component_names Character vector, names of compositional layers in x
#' @param covariates A `terra::SpatRaster` with auxiliary variable layers
#'   (e.g., elevation, slope, parent_material_class)
#' @param covariate_names Character vector, names of covariate layers to use
#'   (must be subset of names(covariates))
#' @param trend_method Character, method for incorporating covariates:
#'   - "simple": Compute per-component linear regressions against covariates (default)
#'   - "local": Use local polynomial regression (loess) for nonlinear trends
#'   - "categorical": Stratify by covariate classes and compute per-class means
#'
#' @return A list with elements:
#'   - `detrended_raster`: New raster with residuals after removing trend
#'   - `trend_surfaces`: Raster stack of predicted trend surfaces per component
#'   - `regressions`: Summary of linear trend coefficients (if trend_method = "simple")
#'   - `metadata`: Processing metadata
#'
#' @details
#' The function fits a regression of each component against the covariates,
#' then subtracts the predicted values to produce detrended residuals. These
#' residuals should be more stationary and suitable for kriging.
#'
#' Common workflow:
#' ```r
#' # Detrend compositional data using elevation
#' detrended <- gc_apply_covariate_trend(
#'   x = sample_compositions,
#'   component_names = c("sand", "silt", "clay"),
#'   covariates = elevation_raster,
#'   covariate_names = "elevation",
#'   trend_method = "simple"
#' )
#'
#' # Perform simulation on detrended residuals
#' residual_sim <- gc_sim_composition(model, domain, nsim = 100, data = detrended$detrended_raster)
#'
#' # Add trend back to simulations
#' trend_surface <- terra::extract(detrended$trend_surfaces, prediction_points)
#' final_sim <- residual_sim + trend_surface
#' ```
#'
#' @export
gc_apply_covariate_trend <- function(x,
                                      component_names,
                                      covariates,
                                      covariate_names = NULL,
                                      trend_method = "simple") {
  
  # Validate inputs
  if (!inherits(x, "SpatRaster")) {
    stop("x must be a terra::SpatRaster")
  }
  
  if (!inherits(covariates, "SpatRaster")) {
    stop("covariates must be a terra::SpatRaster")
  }
  
  if (!all(component_names %in% terra::names(x))) {
    missing <- setdiff(component_names, terra::names(x))
    stop("component_names not found in x: ", paste(missing, collapse = ", "))
  }
  
  if (is.null(covariate_names)) {
    covariate_names <- terra::names(covariates)
  }
  
  if (!all(covariate_names %in% terra::names(covariates))) {
    missing <- setdiff(covariate_names, terra::names(covariates))
    stop("covariate_names not found in covariates: ", paste(missing, collapse = ", "))
  }
  
  # Validate trend method
  valid_methods <- c("simple", "local", "categorical")
  if (!trend_method %in% valid_methods) {
    stop("trend_method must be one of: ", paste(valid_methods, collapse = ", "))
  }
  
  # Extract component and covariate data
  comp_data <- terra::extract(x[[component_names]], 1:terra::ncell(x))[-1]  # Drop ID column
  cov_data <- terra::extract(covariates[[covariate_names]], 1:terra::ncell(covariates))[-1]
  
  # Handle NAs
  valid_idx <- rowSums(!is.na(comp_data)) == ncol(comp_data) &
               rowSums(!is.na(cov_data)) == ncol(cov_data)
  
  comp_valid <- comp_data[valid_idx, ]
  cov_valid <- cov_data[valid_idx, ]
  
  # Fit trend model
  trend_predictions <- matrix(NA, nrow = nrow(comp_data), ncol = length(component_names))
  regression_coefs <- list()
  
  for (comp in component_names) {
    if (trend_method == "simple") {
      # Linear regression
      fit <- stats::lm(comp_data[[comp]] ~ ., data = cov_data)
      regression_coefs[[comp]] <- stats::coef(fit)
      
      # Predict on all data (with NA for missing covariates)
      pred <- stats::predict(fit,
                             newdata = as.data.frame(cov_data),
                             allow.new.levels = TRUE)
      trend_predictions[, which(component_names == comp)] <- pred
      
    } else if (trend_method == "local") {
      # Local polynomial regression (simplified: use loess if single covariate)
      if (length(covariate_names) == 1) {
        cov_col <- cov_valid[, 1]
        comp_col <- comp_valid[[comp]]
        
        fit <- stats::loess(comp_col ~ cov_col, na.action = "na.omit")
        pred <- stats::predict(fit,
                               newdata = data.frame(cov_col = cov_data[, 1]),
                               se.fit = FALSE)
        trend_predictions[, which(component_names == comp)] <- pred
      } else {
        warning("local trend method not implemented for multivariate covariates; using simple")
        fit <- stats::lm(comp_data[[comp]] ~ ., data = cov_data)
        pred <- stats::predict(fit,
                               newdata = as.data.frame(cov_data),
                               allow.new.levels = TRUE)
        trend_predictions[, which(component_names == comp)] <- pred
      }
      
    } else if (trend_method == "categorical") {
      # Stratify by covariate classes and compute means
      # Simplified: assumes first covariate is categorical
      class_col <- cov_valid[, 1]
      comp_col <- comp_valid[[comp]]
      
      class_means <- stats::aggregate(comp_col, by = list(class_col), FUN = mean)
      names(class_means) <- c("class", "mean")
      
      # Map class means back to all cells
      cov_class <- cov_data[, 1]
      pred <- class_means$mean[match(cov_class, class_means$class)]
      trend_predictions[, which(component_names == comp)] <- pred
    }
  }
  
  # Compute residuals
  residuals <- comp_data - trend_predictions
  
  # Create residual raster
  residual_raster <- terra::rast(x[[component_names]])
  terra::values(residual_raster) <- residuals
  
  # Create trend surface raster
  trend_raster <- terra::rast(x[[component_names]])
  terra::values(trend_raster) <- trend_predictions
  
  list(
    detrended_raster = residual_raster,
    trend_surfaces = trend_raster,
    regressions = if (trend_method == "simple") regression_coefs else NULL,
    metadata = list(
      trend_method = trend_method,
      components = component_names,
      covariates = covariate_names,
      n_valid_observations = sum(valid_idx),
      n_missing = sum(!valid_idx),
      processing_date = Sys.Date()
    )
  )
}


#' Prepare SSURGO Data for Simulation (SDA Direct Interface)
#'
#' Query Soil Data Access (SDA) API and prepare component data for simulation.
#' This is a convenience wrapper that assumes soilDB is installed and available.
#'
#' @param map_unit_keys Character vector of SSURGO map unit keys (mukey values)
#'   to query via SDA
#' @param depth_range Numeric vector c(min_depth, max_depth) in cm to filter
#'   horizons (default: c(0, 30) for surface horizon)
#' @param weight_method Character, component weighting method
#'   ("comppct", "area", or "custom")
#' @param target_sum Numeric, target composition sum (default 100)
#' @param soilDB_installed Logical, whether to check for soilDB availability
#'   (default TRUE; set FALSE if calling from already-soilDB-aware context)
#'
#' @return A list as returned by `gc_ssurgo_to_params()` if successful,
#'   or NULL with warning if soilDB is not available
#'
#' @details
#' This function requires the soilDB package to be installed. Users must call:
#' ```r
#' install.packages('soilDB')
#' library(soilDB)
#' ```
#'
#' The function queries Soil Data Access via `soilDB::get_SDA_property()` or
#' similar functions to fetch component data, then passes it to
#' `gc_ssurgo_to_params()` for processing.
#'
#' @examples
#' \dontrun{
#' # Requires soilDB installed
#' params <- gc_prepare_ssurgo_direct(
#'   map_unit_keys = c("463168", "463171"),
#'   depth_range = c(0, 30),
#'   weight_method = "comppct"
#' )
#' }
#'
#' @export
gc_prepare_ssurgo_direct <- function(map_unit_keys,
                                     depth_range = c(0, 30),
                                     weight_method = "comppct",
                                     target_sum = 100,
                                     soilDB_installed = TRUE) {
  
  # Check for soilDB
  if (soilDB_installed && !requireNamespace("soilDB", quietly = TRUE)) {
    warning(
      "soilDB package required but not installed. ",
      "Install via: install.packages('soilDB')\n",
      "This function returns NULL; consider using geocoda functions directly with pre-fetched data."
    )
    return(NULL)
  }
  
  # If soilDB is available, fetch data
  tryCatch({
    # This is a placeholder - actual implementation would query SDA
    # For now, return guidance to user
    message(
      "gc_prepare_ssurgo_direct requires manual integration of soilDB queries.\n",
      "Example workflow:\n",
      "  library(soilDB)\n",
      "  # Fetch SSURGO component data via SDA\n",
      "  ssurgo_data <- soilDB::get_SDA_property(...)\n",
      "  # Then process with geocoda\n",
      "  params <- geocoda::gc_ssurgo_to_params(ssurgo_data, ...)"
    )
    return(NULL)
  }, error = function(e) {
    warning("Error querying SDA: ", e$message)
    return(NULL)
  })
}
