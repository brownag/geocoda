#' Zone-Stratified Cross-Validation
#'
#' Compute k-fold cross-validation metrics stratified by zone.
#' Identifies which zones have well-constrained or data-poor model fits,
#' enabling targeted sampling and highlighting zones needing more data.
#'
#' @param data Data frame with compositional columns and zone assignment
#' @param zone_vector Character or factor zone assignments (length = nrow(data))
#' @param comp_cols Character vector of compositional column names
#' @param k Numeric number of cross-validation folds (default 5)
#' @param method Character cross-validation method: "kfold" or "loo" (leave-one-out)
#' @param model_func Function that fits model to training data and predicts on test data
#'   (default: uses gc_ilr_model for compositional predictions)
#' @param verbose Logical; print per-zone CV results
#'
#' @details
#' ## Cross-Validation Strategy
#'
#' For each zone, separately stratify k-fold CV to ensure zone representation
#' in both training and test sets. Prevents data-poor zones from being excluded
#' by random fold assignment.
#'
#' **CV Metrics computed per zone**:
#' - RMSE: Root mean square error on held-out test set
#' - MAE: Mean absolute error (robust to outliers)
#' - R²: Coefficient of determination
#' - n_test: Count of test observations per zone
#'
#' ## Interpretation
#'
#' - **Good CV (low RMSE, high R²)**: Zone is well-fitted; sufficient data
#' - **Poor CV (high RMSE, low R²)**: Zone is underfitted; needs more sampling
#' - **Asymmetric CV**: Some folds much worse (suggests outliers or inhomogeneity)
#'
#' @return Data frame with columns:
#'   - `zone`: Zone identifier
#'   - `fold`: Fold number (1 to k)
#'   - `n_test`: Test set size in fold
#'   - `rmse`: Root mean square error
#'   - `mae`: Mean absolute error
#'   - `r_squared`: Coefficient of determination (per component)
#'   - Zone-level summary statistics (mean RMSE, SD, consistency)
#'
#' @seealso [gc_ensemble_per_zone()], [gc_compute_entropy_per_zone()],
#'   [gc_cross_validate()] (global cross-validation)
#'
#' @examples
#' \dontrun{
#' data(soil_hierarchy_example)
#'
#' cv_results <- gc_cross_validate_per_zone(
#'   data = soil_hierarchy_example$comp,
#'   zone_vector = soil_hierarchy_example$zone,
#'   comp_cols = c("SAND", "SILT", "CLAY"),
#'   k = 5,
#'   verbose = TRUE
#' )
#'
#' # Identify poorly-validated zones
#' poor_zones <- cv_results[cv_results$mean_rmse > quantile(cv_results$mean_rmse, 0.75), ]
#' cat("Zones needing more data:", unique(poor_zones$zone), "\n")
#' }
#'
#' @export
gc_cross_validate_per_zone <- function(data,
                                        zone_vector,
                                        comp_cols = NULL,
                                        k = 5,
                                        method = "kfold",
                                        model_func = NULL,
                                        verbose = FALSE) {

  if (!is.data.frame(data)) {
    stop("data must be a data frame")
  }

  if (length(zone_vector) != nrow(data)) {
    stop("zone_vector length must equal nrow(data)")
  }

  if (is.null(comp_cols)) {
    # Guess compositional columns (common soil names)
    possible_comp <- c("SAND", "SILT", "CLAY", "OM", "OC")
    comp_cols <- intersect(possible_comp, names(data))
    if (length(comp_cols) == 0) {
      stop("Could not identify compositional columns; specify comp_cols")
    }
  }

  if (is.null(model_func)) {
    # Default: ILR-based model
    model_func <- function(train_data, test_data, comps) {
      tryCatch(
        {
          # Fit ILR model on training data
          ilr_fit <- gc_ilr_model(train_data[, comps])
          # Predict on test data
          pred <- predict(ilr_fit, test_data[, comps])
          list(predictions = pred, obs = test_data[, comps])
        },
        error = function(e) list(predictions = NA, obs = NA)
      )
    }
  }

  zone_names <- unique(zone_vector)
  results <- list()

  for (zone_name in zone_names) {
    zone_mask <- zone_vector == zone_name
    zone_data <- data[zone_mask, comp_cols]
    zone_n <- nrow(zone_data)

    # Create folds within zone (or LOO if method = "loo")
    if (method == "loo") {
      folds <- seq_len(zone_n)
    } else {
      fold_size <- max(1, floor(zone_n / k))
      folds <- rep(seq_len(k), length.out = zone_n)
      folds <- folds[sample.int(zone_n)]
    }

    fold_results <- list()

    for (fold_id in unique(folds)) {
      test_mask <- folds == fold_id
      train_data <- zone_data[!test_mask, ]
      test_data <- zone_data[test_mask, ]

      if (nrow(train_data) == 0 || nrow(test_data) == 0) {
        next
      }

      # Fit and predict
      fit_result <- model_func(train_data, test_data, comp_cols)

      if (all(is.na(fit_result$predictions))) {
        next
      }

      # Compute metrics
      obs <- fit_result$obs
      pred <- fit_result$predictions

      rmse <- sqrt(mean((obs - pred)^2, na.rm = TRUE))
      mae <- mean(abs(obs - pred), na.rm = TRUE)

      # R² per component
      r_squared <- sapply(comp_cols, function(comp) {
        ss_res <- sum((obs[[comp]] - pred[[comp]])^2, na.rm = TRUE)
        ss_tot <- sum((obs[[comp]] - mean(obs[[comp]], na.rm = TRUE))^2, na.rm = TRUE)
        if (ss_tot == 0) NA else 1 - (ss_res / ss_tot)
      })

      fold_results[[length(fold_results) + 1]] <- list(
        zone = as.character(zone_name),
        fold = fold_id,
        n_test = nrow(test_data),
        rmse = round(rmse, 4),
        mae = round(mae, 4),
        r_squared_mean = round(mean(r_squared, na.rm = TRUE), 4),
        r_squared_sd = round(sd(r_squared, na.rm = TRUE), 4)
      )
    }

    # Aggregate per zone
    if (length(fold_results) > 0) {
      fold_df <- as.data.frame(do.call(rbind, lapply(fold_results, as.data.frame)))

      zone_summary <- list(
        zone = as.character(zone_name),
        n_samples = zone_n,
        n_folds = nrow(fold_df),
        mean_rmse = round(mean(fold_df$rmse), 4),
        sd_rmse = round(sd(fold_df$rmse), 4),
        mean_r_squared = round(mean(fold_df$r_squared_mean), 4),
        sd_r_squared = round(sd(fold_df$r_squared_mean), 4)
      )

      results[[length(results) + 1]] <- zone_summary
    }
  }

  result_df <- as.data.frame(do.call(rbind, lapply(results, as.data.frame)))

  if (verbose) {
    cat("Zone-Stratified Cross-Validation Summary\n")
    cat("=========================================\n")
    print(result_df)
    cat("\nInterpretation:\n")
    cat("- Zones with high mean_rmse: underfitted (need more data)\n")
    cat("- Zones with high sd_rmse: unstable (heterogeneous conditions)\n")
    cat("- Zones with low mean_r_squared: poor predictability\n")
  }

  result_df
}


#' Zone-Stratified Compositional Entropy
#'
#' Compute Shannon entropy of compositional uncertainty by zone.
#' Entropy quantifies compositional diversity and uncertainty; high entropy
#' zones indicate greater sampling difficulty or inherent variability.
#'
#' @param ensemble List from [gc_sim_hierarchical()] or [gc_aggregate_realizations()]
#' @param zone_col Character zone column name, or NULL
#' @param zone_definition sf polygon, factor, or character zones
#' @param method Character; entropy computation method: "shannon" (default), "simpson"
#' @param normalize Logical; normalize entropy to [0, 1] range
#'
#' @details
#' ## Shannon Entropy
#'
#' For compositional data on the simplex, entropy measures uncertainty in
#' the distribution of proportions:
#'
#' $$H = -\sum_i p_i \log(p_i)$$
#'
#' where p_i is the mean proportion of component i in the zone.
#'
#' **Interpretation**:
#' - High entropy: Components nearly equally abundant (uncertain composition)
#' - Low entropy: One component dominates (certain/homogeneous composition)
#'
#' High-entropy zones are "hot spots" of compositional uncertainty and may
#' benefit from targeted additional sampling.
#'
#' @return Data frame with columns:
#'   - `zone`: Zone identifier
#'   - `n_ensemble`: Ensemble points in zone
#'   - `shannon_entropy`: Shannon entropy value
#'   - `entropy_normalized`: Entropy scaled to [0, 1]
#'   - `dominant_component`: Component with highest mean proportion
#'   - `dominant_proportion`: Mean proportion of dominant component
#'   - `uncertainty_class`: "Low" / "Medium" / "High" based on entropy quantiles
#'
#' @seealso [gc_ensemble_per_zone()], [gc_cross_validate_per_zone()],
#'   [gc_compute_entropy()] (global entropy)
#'
#' @examples
#' \dontrun{
#' ensemble <- gc_sim_hierarchical(
#'   object = hz_fit,
#'   n = 100,
#'   zone_vector = zone_vector,
#'   coords = coords
#' )
#'
#' entropy_results <- gc_compute_entropy_per_zone(
#'   ensemble = ensemble,
#'   zone_col = "zone_id",
#'   method = "shannon",
#'   verbose = TRUE
#' )
#'
#' # Identify high-entropy zones (uncertainty hotspots)
#' hotspots <- entropy_results[entropy_results$uncertainty_class == "High", ]
#' }
#'
#' @export
gc_compute_entropy_per_zone <- function(ensemble,
                                         zone_col = NULL,
                                         zone_definition = NULL,
                                         method = "shannon",
                                         normalize = TRUE,
                                         verbose = FALSE) {

  if (!is.list(ensemble) || !"data" %in% names(ensemble)) {
    stop("ensemble must be from gc_sim_hierarchical() or gc_aggregate_realizations()")
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
  comp_names <- names(comp_data)

  # Compute entropy per zone
  results <- list()

  for (zone_name in zone_names) {
    zone_mask <- zone_assignment == zone_name
    zone_ensemble <- comp_data[zone_mask, ]
    n_ensemble <- sum(zone_mask)

    # Compute mean proportions
    mean_props <- colMeans(zone_ensemble, na.rm = TRUE)

    # Normalize to sum to 1 (in case of rounding)
    mean_props <- mean_props / sum(mean_props)

    # Compute entropy
    if (method == "shannon") {
      # Shannon entropy: -sum(p * log(p))
      entropy <- -sum(mean_props * log(mean_props + 1e-10), na.rm = TRUE)
    } else if (method == "simpson") {
      # Simpson diversity: 1 - sum(p^2)
      entropy <- 1 - sum(mean_props^2)
    } else {
      entropy <- NA
    }

    # Normalize to [0, 1]
    n_components <- length(comp_names)
    if (normalize && method == "shannon") {
      max_entropy <- log(n_components)
      entropy_normalized <- entropy / max_entropy
    } else {
      entropy_normalized <- entropy
    }

    # Identify dominant component
    dominant_idx <- which.max(mean_props)
    dominant_component <- comp_names[dominant_idx]
    dominant_proportion <- mean_props[dominant_idx]

    results[[length(results) + 1]] <- list(
      zone = as.character(zone_name),
      n_ensemble = n_ensemble,
      shannon_entropy = round(entropy, 4),
      entropy_normalized = round(entropy_normalized, 4),
      dominant_component = dominant_component,
      dominant_proportion = round(dominant_proportion, 4)
    )
  }

  result_df <- as.data.frame(do.call(rbind, lapply(results, as.data.frame)))

  # Classify uncertainty
  entropy_quantiles <- quantile(result_df$shannon_entropy)
  result_df$uncertainty_class <- cut(
    result_df$shannon_entropy,
    breaks = c(0, entropy_quantiles[2], entropy_quantiles[4], Inf),
    labels = c("Low", "Medium", "High"),
    include.lowest = TRUE
  )

  if (verbose) {
    cat("Zone-Stratified Compositional Entropy\n")
    cat("====================================\n")
    print(result_df)
    cat("\nUncertainty Hotspots (High Entropy):\n")
    hotspots <- result_df[result_df$uncertainty_class == "High", ]
    if (nrow(hotspots) > 0) {
      for (i in seq_len(nrow(hotspots))) {
        cat(sprintf("%s: entropy = %.3f (recommend additional sampling)\n",
                    hotspots$zone[i], hotspots$shannon_entropy[i]))
      }
    }
  }

  result_df
}


#' Zone-Stratified Bootstrap Parameter Uncertainty
#'
#' Compute parameter confidence intervals via bootstrap resampling, stratified by zone.
#' Assesses robustness of zone-specific parameter estimates under data uncertainty.
#'
#' @param data Data frame with compositional and zone columns
#' @param zone_vector Character or factor zone assignments
#' @param param_func Function that computes parameters from data
#'   (Default: computes ILR means and variances)
#' @param n_bootstrap Numeric number of bootstrap samples (default 1000)
#' @param conf_level Numeric confidence level (default 0.95)
#' @param verbose Logical; print per-zone bootstrap results
#'
#' @details
#' ## Bootstrap Strategy
#'
#' For each zone, independently resample with replacement n_bootstrap times.
#' Compute parameters each time; confidence intervals from bootstrap quantiles.
#'
#' **Parameters estimated** (default):
#' - Mean ILR coordinates per component
#' - Variance of ILR coordinates
#' - Closure-checked composition (inverse ILR transform)
#'
#' ## Interpretation
#'
#' - **Narrow CI**: Parameter well-estimated (stable)
#' - **Wide CI**: Parameter poorly-estimated (high data uncertainty)
#' - **Asymmetric CI**: Skewed distribution (possible outliers)
#'
#' @return List with per-zone bootstrap results:
#'   - `zone_estimates`: Data frame with central estimate, lower, upper CI
#'   - `zones`: Character vector of zone names
#'   - `n_bootstrap`: Number of bootstrap samples
#'   - `conf_level`: Confidence level used
#'
#' @seealso [gc_cross_validate_per_zone()], [gc_compute_entropy_per_zone()],
#'   [gc_bootstrap_uncertainty()] (global bootstrap)
#'
#' @export
gc_bootstrap_uncertainty_per_zone <- function(data,
                                               zone_vector,
                                               param_func = NULL,
                                               n_bootstrap = 1000,
                                               conf_level = 0.95,
                                               verbose = FALSE) {

  if (!is.data.frame(data)) {
    stop("data must be a data frame")
  }

  if (length(zone_vector) != nrow(data)) {
    stop("zone_vector length must equal nrow(data)")
  }

  if (is.null(param_func)) {
    # Default: compute ILR means and variance
    param_func <- function(zone_data) {
      comp_cols <- c("SAND", "SILT", "CLAY")
      comp_cols <- intersect(comp_cols, names(zone_data))

      if (length(comp_cols) < 2) {
        return(rep(NA, 3))
      }

      ilr_params <- gc_ilr_params(zone_data[, comp_cols])
      c(
        mean_ilr_1 = ilr_params$mean[1],
        mean_ilr_2 = if (length(ilr_params$mean) > 1) ilr_params$mean[2] else NA,
        var_ilr = mean(diag(ilr_params$cov), na.rm = TRUE)
      )
    }
  }

  zone_names <- unique(zone_vector)
  alpha <- 1 - conf_level

  results <- list()

  for (zone_name in zone_names) {
    zone_mask <- zone_vector == zone_name
    zone_data <- data[zone_mask, ]
    zone_n <- nrow(zone_data)

    # Bootstrap
    boot_params <- matrix(nrow = n_bootstrap, ncol = 3)

    for (b in seq_len(n_bootstrap)) {
      boot_idx <- sample.int(zone_n, replace = TRUE)
      boot_data <- zone_data[boot_idx, ]
      boot_params[b, ] <- param_func(boot_data)
    }

    # Compute quantiles
    ci_lower <- apply(boot_params, 2, quantile, alpha / 2, na.rm = TRUE)
    ci_upper <- apply(boot_params, 2, quantile, 1 - alpha / 2, na.rm = TRUE)
    central <- apply(boot_params, 2, mean, na.rm = TRUE)

    results[[length(results) + 1]] <- list(
      zone = as.character(zone_name),
      n_obs = zone_n,
      param_1_est = round(central[1], 4),
      param_1_lower = round(ci_lower[1], 4),
      param_1_upper = round(ci_upper[1], 4),
      param_2_est = round(central[2], 4),
      param_2_lower = round(ci_lower[2], 4),
      param_2_upper = round(ci_upper[2], 4),
      variance_est = round(central[3], 4),
      variance_lower = round(ci_lower[3], 4),
      variance_upper = round(ci_upper[3], 4)
    )
  }

  result_df <- as.data.frame(do.call(rbind, lapply(results, as.data.frame)))

  if (verbose) {
    cat("Zone-Stratified Bootstrap Parameter Uncertainty\n")
    cat("===============================================\n")
    print(result_df)
    cat("\nInterpretation:\n")
    cat("- Narrow CI: parameter is well-estimated\n")
    cat("- Wide CI: parameter is poorly-estimated (needs more data)\n")
    cat("- Asymmetric CI: possible outliers or skewed distribution\n")
  }

  list(
    zone_estimates = result_df,
    zones = zone_names,
    n_bootstrap = n_bootstrap,
    conf_level = conf_level
  )
}
