#' Estimate ILR Parameters from Compositional Samples
#'
#' Convert compositional samples to Isometric Log-Ratio (ILR) space and compute
#' the mean vector and covariance matrix. These statistics are used to parameterize
#' the multivariate geostatistical model.
#'
#' @param samples A data frame of compositions (each row sums to constant,
#'   typically 100% or 1).
#'
#' @return A list containing:
#'   - `mean`: Vector of ILR means (length D-1 where D is number of components)
#'   - `cov`: Covariance matrix of ILR values (dimensions (D-1) x (D-1))
#'   - `names`: Original component names (character vector)
#'   - `base_class`: Class used for transformation ("acomp")
#'
#' @details
#' This function:
#' 1. Stores original column names for later use
#' 2. Converts to [compositions::acomp] (absolute composition)
#' 3. Applies [compositions::ilr] transformation
#' 4. Computes column means and covariance matrix of ILR values
#'
#' The ILR transformation eliminates the sum constraint, allowing standard
#' multivariate geostatistics to be applied.
#'
#' @examples
#' # Simulate some simple compositions
#' samples <- data.frame(
#'   sand = c(20, 25, 30, 22),
#'   silt = c(60, 55, 50, 58),
#'   clay = c(20, 20, 20, 20)
#' )
#'
#' params <- estimate_ilr_params(samples)
#' str(params)
#' print(params$mean)
#' print(params$cov)
#'
#' @importFrom compositions acomp ilr
#' @importFrom stats cov
#' @export
gc_ilr_params <- function(samples) {
  comp_names <- colnames(samples)

  comp <- compositions::acomp(samples)
  ilr_vals <- compositions::ilr(comp)

  ilr_mean <- colMeans(ilr_vals)
  ilr_cov <- stats::cov(ilr_vals)

  list(
    mean = ilr_mean,
    cov = ilr_cov,
    names = comp_names,
    base_class = "acomp"
  )
}


#' Build Compositional Geostatistical Model
#'
#' Constructs a gstat model for the ILR-transformed variables using either
#' Independent Univariate Kriging or Linear Model of Coregionalization (LMC).
#'
#' @param ilr_params A list from `gc_ilr_params()` containing mean, cov, names.
#' @param variogram_model A `vgm()` object defining the base variogram structure.
#' @param data Optional `sf` object with ILR columns for conditioning.
#'   If provided, Conditional Simulation is performed.
#'   If NULL (default), Unconditional Simulation is performed.
#' @param model_type Character string specifying the approach: `"univariate"` 
#'   (default) or `"lmc"`. Univariate is numerically stable and standard practice.
#'   LMC includes cross-covariance terms between ILR dimensions.
#'
#' @return A `gstat` object representing the model.
#'
#' @details
#' **Independent Univariate Kriging** (`model_type = "univariate"`):
#' Models each ILR dimension separately without cross-covariance terms. This is:
#' - Numerically stable (avoids positive-definite issues in LMC fitting)
#' - Standard practice in compositional geostatistics
#' - Robust across different datasets
#' - Efficient for large problems
#'
#' The ILR transformation already decorrelates the data significantly, so ignoring
#' spatial cross-correlation between ILR coordinates has minimal impact.
#'
#' **Linear Model of Coregionalization** (`model_type = "lmc"`):
#' Includes cross-covariance terms between all pairs of ILR dimensions. This is:
#' - Theoretically more complete
#' - More numerically complex
#' - Useful when cross-correlation structure is important
#'
#' For conditional simulation, the conditioning data must be passed to this function
#' (not to `predict()`). The model then automatically uses that data during prediction.
#'
#' @importFrom gstat gstat
#' @export
gc_ilr_model <- function(ilr_params, variogram_model, data = NULL, model_type = "univariate") {
  model_type <- match.arg(model_type, c("univariate", "lmc"))
  
  g <- NULL
  n_ilr <- nrow(ilr_params$cov)

  if (model_type == "univariate") {
    for (i in seq_len(n_ilr)) {
      var_name <- paste0("ilr", i)
      
      local_model <- variogram_model
      local_model$psill <- local_model$psill * ilr_params$cov[i, i]
      
      f <- as.formula(paste(var_name, "~ 1"))
      
      if (is.null(data)) {
        g <- gstat::gstat(g, id = var_name, formula = f,
                          dummy = TRUE, beta = ilr_params$mean[i],
                          model = local_model)
      } else {
        g <- gstat::gstat(g, id = var_name, formula = f,
                          data = data, beta = ilr_params$mean[i],
                          model = local_model)
      }
    }
  } else {
    for (i in seq_len(n_ilr)) {
      var_name <- paste0("ilr", i)
      
      local_model <- variogram_model
      local_model$psill <- local_model$psill * ilr_params$cov[i, i]
      
      f <- as.formula(paste(var_name, "~ 1"))
      
      if (is.null(data)) {
        g <- gstat::gstat(g, id = var_name, formula = f,
                          dummy = TRUE, beta = ilr_params$mean[i],
                          model = local_model)
      } else {
        g <- gstat::gstat(g, id = var_name, formula = f,
                          data = data, beta = ilr_params$mean[i],
                          model = local_model)
      }
    }
    
    for (i in seq_len(n_ilr - 1)) {
      for (j in (i + 1):n_ilr) {
        var_i <- paste0("ilr", i)
        var_j <- paste0("ilr", j)
        
        cross_model <- variogram_model
        cross_model$psill <- variogram_model$psill * ilr_params$cov[i, j]
        
        g <- gstat::gstat(g, id = c(var_i, var_j), model = cross_model)
      }
    }
  }
  
  return(g)
}
