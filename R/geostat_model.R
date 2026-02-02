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
#' @references
#' Egozcue, J. J., Pawlowsky-Glahn, V., Mateu-Figueras, G., & Barceló-Vidal, C. (2003).
#' Isometric Log-Ratio Transformations for Compositional Data Analysis.
#' *Mathematical Geology*, 35(3), 279-300. https://doi.org/10.1023/A:1023818214614
#'
#' Aitchison, J. (1986). The Statistical Analysis of Compositional Data.
#' Chapman and Hall, London.
#'
#' @examples
#' # Simulate some simple compositions
#' samples <- data.frame(
#'   sand = c(20, 25, 30, 22),
#'   silt = c(60, 55, 50, 58),
#'   clay = c(20, 20, 20, 20)
#' )
#'
#' params <- gc_ilr_params(samples)
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
#' Optionally applies Minimum/Maximum Autocorrelation Factors (MAF) for improved
#' spatial decorrelation.
#'
#' @param ilr_params A list from `gc_ilr_params()` containing mean, cov, names.
#' @param variogram_model A `vgm()` object defining the base variogram structure.
#' @param data Optional `sf` object with ILR columns for conditioning.
#'   If provided, Conditional Simulation is performed.
#'   If NULL (default), Unconditional Simulation is performed.
#' @param model_type Character string specifying the approach: `"univariate"`
#'   (default) or `"lmc"`. Univariate is numerically stable and standard practice.
#'   LMC includes cross-covariance terms between ILR dimensions.
#' @param use_maf Logical. If TRUE, apply MAF decorrelation (default FALSE).
#' @param maf_lag_h Numeric. Distance for MAF computation. If NULL, uses mean_range/3.
#'
#' @return A `gstat` object representing the model with attributes:
#'   - `ilr_params`: The input ILR parameters
#'   - `maf_object`: MAF transformation object (if use_maf = TRUE)
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
#' @references
#' Wackernagel, H. (2003). Multivariate Geostatistics: An Introduction with Applications (3rd ed.).
#' Springer-Verlag, Berlin.
#'
#' Chilès, J. P., & Delfiner, P. (2012). Geostatistics: Modeling Spatial Uncertainty (2nd ed.).
#' Wiley, Hoboken, NJ.
#'
#' Tolosana-Delgado, R., Mueller, U., & van den Boogaart, K. G. (2012).
#' Geostatistics for Compositional Data: An Overview.
#' *Mathematical Geosciences*, 44(4), 465-479. https://doi.org/10.1007/s11004-011-9363-4
#'
#' @examples
#' \dontrun{
#' ilr_prm <- gc_ilr_params(compositions::acomp(data.frame(
#'   sand = c(50, 60, 45), silt = c(30, 25, 35), clay = c(20, 15, 20)
#' )))
#' vgm <- gstat::vgm(psill = 0.5, model = "Exp", range = 50, nugget = 0.1)
#' model <- gc_ilr_model(ilr_prm, vgm, model_type = "univariate")
#' }
#'
#' @importFrom gstat gstat
#' @export
gc_ilr_model <- function(ilr_params, variogram_model, data = NULL,
                         model_type = "univariate", use_maf = FALSE,
                         maf_lag_h = NULL) {
  model_type <- match.arg(model_type, c("univariate", "lmc"))

  # MAF decorrelation (optional)
  maf_object <- NULL
  ilr_mean_for_inv <- ilr_params$mean  # Store for back-transform
  if (use_maf) {
    # Prepare variogram models (one per ILR dimension)
    n_ilr_temp <- nrow(ilr_params$cov)
    vgm_list <- replicate(n_ilr_temp, variogram_model, simplify = FALSE)

    # Compute MAF rotation
    maf_object <- gc_compute_maf(ilr_params, vgm_list, lag_h = maf_lag_h)

    # Transform data to MAF space if conditioning data provided
    if (!is.null(data)) {
      data <- gc_transform_maf(data, maf_object, ilr_mean = ilr_params$mean)
    }

    # Update mean vector to MAF means (all zero after centering and rotation)
    ilr_params$mean <- rep(0, nrow(ilr_params$cov))
  }

  g <- NULL
  n_ilr <- nrow(ilr_params$cov)

  if (model_type == "univariate") {
    for (i in seq_len(n_ilr)) {
      var_name <- if (use_maf) paste0("maf", i) else paste0("ilr", i)
      
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
      var_name <- if (use_maf) paste0("maf", i) else paste0("ilr", i)
      
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
        var_i <- if (use_maf) paste0("maf", i) else paste0("ilr", i)
        var_j <- if (use_maf) paste0("maf", j) else paste0("ilr", j)
        
        cross_model <- variogram_model
        cross_model$psill <- variogram_model$psill * ilr_params$cov[i, j]
        
        g <- gstat::gstat(g, id = c(var_i, var_j), model = cross_model)
      }
    }
  }

  # Store parameters in attributes for later use in simulation
  attr(g, "ilr_params") <- ilr_params
  if (!is.null(maf_object)) {
    attr(g, "maf_object") <- maf_object
    attr(g, "ilr_mean_original") <- ilr_mean_for_inv  # For back-transform
  }

  return(g)
}
