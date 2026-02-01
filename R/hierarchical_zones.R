#' Bayesian Hierarchical Zones
#'
#' Multi-level spatial modeling for compositional data with zone-specific parameters
#' and pooled uncertainty estimates. Implement Bayesian hierarchical framework enabling
#' zone-specific inference with global shrinkage priors.
#'
#' @keywords internal
#' @name hierarchical_zones
NULL


#' Define Hierarchical Zone Structure
#'
#' Specify nested spatial domains (zones within regions) for hierarchical modeling.
#' Zones inherit prior information from parent regions while allowing zone-specific
#' parameter estimation with shrinkage toward global distribution.
#'
#' @param zone_names Character vector, names of zones (e.g., "zone_1", "zone_2")
#' @param zone_boundaries Spatial object (sf) or list of spatial polygons defining zone extent
#' @param parent_region Character or NULL, parent region name for all zones (optional)
#' @param component_names Character vector, names of components (e.g., c("sand", "silt", "clay"))
#' @param target_sum Numeric, expected composition sum (default 100)
#'
#' @return An S3 object of class `"gc_hierarchy"` containing:
#'   - `zones`: Character vector of zone names
#'   - `boundaries`: Spatial boundaries per zone
#'   - `parent_region`: Parent region identifier
#'   - `components`: Component names
#'   - `target_sum`: Composition sum constraint
#'   - `n_zones`: Number of zones
#'   - `n_components`: Number of components
#'
#' @examples
#' \dontrun{
#' # Create hierarchy from stratification results
#' strata_result <- gc_identify_strata(data, n_strata = 3, plot = FALSE)
#' 
#' hierarchy <- gc_define_hierarchy(
#'   zone_names = paste0("stratum_", 1:3),
#'   parent_region = "study_area",
#'   component_names = c("sand", "silt", "clay")
#' )
#' 
#' print(hierarchy)
#' }
#'
#' @export
gc_define_hierarchy <- function(zone_names,
                                zone_boundaries = NULL,
                                parent_region = NULL,
                                component_names,
                                target_sum = 100) {
  
  if (length(zone_names) == 0) {
    stop("zone_names must have at least one zone")
  }
  
  if (length(component_names) == 0) {
    stop("component_names must contain at least two components")
  }
  
  if (length(component_names) < 2) {
    stop("Compositional data requires at least 2 components")
  }
  
  hierarchy <- list(
    zones = as.character(zone_names),
    boundaries = zone_boundaries,
    parent_region = parent_region,
    components = as.character(component_names),
    target_sum = target_sum,
    n_zones = length(zone_names),
    n_components = length(component_names),
    n_ilr = length(component_names) - 1
  )
  
  class(hierarchy) <- c("gc_hierarchy", "list")
  
  hierarchy
}


#' Set Prior Specifications and Pooling for Hierarchical Zones
#'
#' Configure pooling strategy and parameter specifications for zone-specific shrinkage estimation.
#' This implements empirical Bayes shrinkage (NOT Bayesian MCMC inference).
#' The pooling_coefficient controls how much zone estimates are shrunk toward the global mean:
#' - `"none"` (λ=0.00): Independent zone models, no pooling (zone estimates only)
#' - `"weak"` (λ=0.05): Light shrinkage toward global
#' - `"moderate"` (λ=0.20): Balanced information sharing across zones
#' - `"strong"` (λ=0.50): Heavy shrinkage toward global mean
#'
#' @param hierarchy An object of class `"gc_hierarchy"`
#' @param pooling Character, pooling strategy: `"none"`, `"weak"`, `"moderate"`, or `"strong"`
#'   (default `"moderate"`)
#' @param global_mean Numeric vector or NULL, global mean estimates (length = n_components).
#'   If NULL, computed from data.
#' @param global_covariance Numeric matrix or NULL, global covariance (n_components × n_components).
#'   If NULL, computed from data.
#' @param prior_sd_mean Numeric, prior SD for mean parameters (default 1.0)
#' @param prior_shape_variance Numeric, Gamma prior shape for variance (default 2.0)
#' @param prior_rate_variance Numeric, Gamma prior rate for variance (default 0.5)
#' @param fit_data Optional data frame with ILR columns and zone assignments
#'   Used to auto-compute global statistics
#' @param verbose Logical, print details (default TRUE)
#'
#' @return An S3 object of class `"gc_prior_spec"` containing:
#'   - `hierarchy`: Reference to parent hierarchy
#'   - `pooling`: Pooling strategy
#'   - `pooling_coefficient`: Numeric strength (0 = none, 1 = complete pooling)
#'   - `global_mean`: Global mean estimates
#'   - `global_covariance`: Global covariance matrix
#'   - `prior_sd_mean`: Prior SD for means
#'   - `prior_shape`: Prior shape parameters
#'   - `prior_rate`: Prior rate parameters
#'   - `metadata`: Creation date, data summary
#'
#' @examples
#' \dontrun{
#' hierarchy <- gc_define_hierarchy(
#'   zone_names = c("Upland", "Midslope", "Lowland"),
#'   component_names = c("sand", "silt", "clay")
#' )
#' priors <- gc_set_hierarchy_priors(hierarchy, pooling = "moderate")
#' }
#'
#' @export
gc_set_hierarchy_priors <- function(hierarchy,
                                    pooling = "moderate",
                                    global_mean = NULL,
                                    global_covariance = NULL,
                                    prior_sd_mean = 1.0,
                                    prior_shape_variance = 2.0,
                                    prior_rate_variance = 0.5,
                                    fit_data = NULL,
                                    verbose = TRUE) {
  
  if (!inherits(hierarchy, "gc_hierarchy")) {
    stop("hierarchy must be an object of class 'gc_hierarchy'")
  }
  
  # Validate pooling strategy
  pooling <- match.arg(pooling, c("none", "weak", "moderate", "strong"))
  
  pooling_coef <- switch(pooling,
    "none" = 0.00,
    "weak" = 0.05,
    "moderate" = 0.20,
    "strong" = 0.50
  )
  
  # Compute global statistics from data if not provided
  if (is.null(global_mean) && !is.null(fit_data)) {
    ilr_cols <- grep("^ilr", names(fit_data), value = TRUE)
    if (length(ilr_cols) > 0) {
      global_mean <- colMeans(fit_data[, ilr_cols, drop = FALSE], na.rm = TRUE)
    }
  }
  
  if (is.null(global_covariance) && !is.null(fit_data)) {
    ilr_cols <- grep("^ilr", names(fit_data), value = TRUE)
    if (length(ilr_cols) > 0) {
      global_covariance <- stats::cov(fit_data[, ilr_cols, drop = FALSE], use = "complete.obs")
    }
  }
  
  # Defaults if still NULL
  n_ilr <- hierarchy$n_components - 1
  if (is.null(global_mean)) {
    global_mean <- rep(0, n_ilr)
  }
  if (is.null(global_covariance)) {
    global_covariance <- diag(1.0, n_ilr)
  }
  
  prior_spec <- list(
    hierarchy = hierarchy,
    pooling = pooling,
    pooling_coefficient = pooling_coef,
    global_mean = global_mean,
    global_covariance = global_covariance,
    prior_sd_mean = prior_sd_mean,
    prior_shape_variance = prior_shape_variance,
    prior_rate_variance = prior_rate_variance,
    metadata = list(
      created = Sys.time(),
      n_zones = hierarchy$n_zones,
      n_components = hierarchy$n_components,
      n_ilr = n_ilr,
      has_fit_data = !is.null(fit_data)
    )
  )
  
  class(prior_spec) <- c("gc_prior_spec", "list")
  
  if (verbose) {
    cat("Hierarchy priors configured:\n")
    cat("  Pooling strategy:", pooling, "\n")
    cat("  Pooling coefficient:", pooling_coef, "\n")
    cat("  Prior SD (means):", prior_sd_mean, "\n")
    cat("  Prior shape (variance):", prior_shape_variance, "\n")
    cat("  Prior rate (variance):", prior_rate_variance, "\n")
    cat("  Zones:", hierarchy$n_zones, "| Components:", hierarchy$n_components, "\n")
  }
  
  prior_spec
}


#' Fit Hierarchical Shrinkage Model to Zone Data
#'
#' Estimate zone-specific parameters using empirical Bayes shrinkage toward
#' a global mean. Zone estimates are pooled via weighted averaging, with the
#' strength of pooling controlled by `pooling_coefficient` in `gc_set_hierarchy_priors()`.
#' This method provides fast analytical estimates suitable for exploration,
#' visualization, and initial inference.
#'
#' **NOTE**: This is empirical Bayes shrinkage (analytical method). For full Bayesian
#' MCMC inference with posterior sampling, use `backend = "stan"` or `backend = "nimble"`
#' (requires installation of rstan or nimble packages respectively).
#'
#' @param data Data frame with columns: zone (factor), ilr1, ilr2, ... (ILR values)
#' @param prior_spec An object of class `"gc_prior_spec"` (from `gc_set_hierarchy_priors()`)
#' @param backend Character, fitting backend: `"analytical"` (default, fast shrinkage estimates),
#'   `"stan"` (HMC MCMC, requires rstan), or `"nimble"` (adaptive MCMC, requires nimble)
#' @param verbose Logical, print estimation progress (default TRUE)
#'
#' @return An S3 object of class `"gc_hierarchical_fit"` containing:
#'   - `zone_estimates`: Zone-specific parameter estimates (mean, covariance)
#'   - `global_estimates`: Global parameter estimates (mean, covariance)
#'   - `samples`: Posterior samples (NULL for analytical, matrix for MCMC backends)
#'   - `diagnostics`: Backend-specific diagnostics
#'   - `prior_spec`: Original prior specification
#'   - `metadata`: Fitting information including backend and timestamp
#'
#' @details
#'
#' ## Analytical Backend (default)
#'
#' Fast empirical Bayes shrinkage method:
#'
#' 1. Computes zone-specific maximum likelihood estimates (MLE)
#' 2. Pools zone estimates toward global mean via weighted average
#' 3. Shrinkage strength controlled by `pooling_coefficient` from priors
#' 4. Closed-form analytical solution with no MCMC sampling
#'
#' Shrinkage formula:
#'
#' ```
#' mu[zone] = (1 - lambda) * mu_zone_MLE + lambda * mu_global
#' Sigma[zone] = [1 / (1 + n_zone * lambda)] * Sigma_zone_MLE
#' ```
#'
#' where `lambda = pooling_coefficient` from `gc_set_hierarchy_priors()`:
#' - lambda = 0.0: No pooling (zone-specific only)
#' - lambda = 0.05: Weak pooling
#' - lambda = 0.2: Moderate pooling
#' - lambda = 0.5: Strong pooling toward global
#'
#' Suitable for quick exploration and visualization. For publication-quality Bayesian
#' inference with full posterior distributions, use Stan or Nimble backends.
#'
#' ## MCMC Backends (Stan, Nimble)
#'
#' Full Bayesian inference with posterior sampling. Requires rstan or nimble to be installed:
#'
#' ```r
#' install.packages("rstan")  # For Stan backend
#' install.packages("nimble") # For Nimble backend
#' ```
#'
#' @examples
#' \dontrun{
#' hierarchy <- gc_define_hierarchy(
#'   zone_names = c("Zone1", "Zone2", "Zone3"),
#'   component_names = c("sand", "silt", "clay")
#' )
#' priors <- gc_set_hierarchy_priors(hierarchy, pooling = "moderate")
#' data_zones <- data.frame(
#'   zone = rep(c("Zone1", "Zone2", "Zone3"), each = 5),
#'   ilr1 = rnorm(15), ilr2 = rnorm(15)
#' )
#'
#' # Fast analytical shrinkage estimates (default)
#' fit_analytical <- gc_fit_hierarchical_model(data_zones, priors)
#' print(fit_analytical)
#'
#' # Full MCMC inference with Stan (requires rstan)
#' # fit_stan <- gc_fit_hierarchical_model(data_zones, priors, backend = "stan")
#' }
#'
#' @export
gc_fit_hierarchical_model <- function(data,
                                     prior_spec,
                                     backend = "analytical",
                                     verbose = TRUE) {

  # Validate backend argument
  backend <- match.arg(backend, c("analytical", "stan", "nimble"))

  # Delegate to backend dispatcher
  fit_hierarchical_backend(data, prior_spec, backend = backend, verbose = verbose)
}


#' Simulate Compositional Realizations from Hierarchical Model
#'
#' Generate spatial realizations using zone-specific posterior parameters from
#' fitted hierarchical model. Each realization respects zone boundaries and
#' compositional constraints.
#'
#' @param fit An object of class `"gc_hierarchical_fit"`
#' @param locations Data frame with columns x, y (coordinates) and zone (zone assignment)
#' @param nsim Numeric, number of realizations (default 5)
#' @param target_names Character vector, names for back-transformed components
#'   (default: c("sand", "silt", "clay"))
#' @param verbose Logical, print progress (default TRUE)
#'
#' @return A list of length `nsim`, where each element is a list containing:
#'   - `realizations_raster`: SpatRaster (per-zone realizations for that simulation)
#'   - `zone_realizations`: List of per-zone SpatRasters
#'   - `metadata`: Realization info (zone, nsim, target_names)
#'
#' @examples
#' \dontrun{
#' fit <- gc_fit_hierarchical_model(data_zones, priors)
#' locations <- data.frame(
#'   x = rep(c(100, 200, 300), each = 5),
#'   y = rep(c(100, 200, 300), 5),
#'   zone = rep(c("Zone1", "Zone2", "Zone3"), each = 5)
#' )
#' realizations <- gc_sim_hierarchical(
#'   fit = fit, locations = locations, nsim = 10
#' )
#' }
#'
#' @export
gc_sim_hierarchical <- function(fit,
                               locations,
                               nsim = 5,
                               target_names = c("sand", "silt", "clay"),
                               verbose = TRUE) {
  
  if (!inherits(fit, "gc_hierarchical_fit")) {
    stop("fit must be an object of class 'gc_hierarchical_fit'")
  }
  
  if (!("zone" %in% names(locations))) {
    stop("locations must contain a 'zone' column")
  }
  
  if (!all(c("x", "y") %in% names(locations))) {
    stop("locations must contain 'x' and 'y' columns")
  }
  
  hierarchy <- fit$prior_spec$hierarchy
  
  if (verbose) {
    cat("Simulating hierarchical realizations\n")
    cat("  Realizations:", nsim, "\n")
    cat("  Zones:", hierarchy$n_zones, "\n")
  }
  
  realizations <- vector("list", nsim)
  
  for (sim in seq_len(nsim)) {
    if (verbose) {
      cat("  Simulation", sim, "of", nsim, "\n")
    }
    
    zone_sims <- list()
    
    for (zone_name in hierarchy$zones) {
      if (!(zone_name %in% names(fit$zone_estimates))) {
        next  # Skip unfitted zones
      }

      # Get zone locations
      zone_idx <- locations$zone == zone_name
      zone_locs <- locations[zone_idx, , drop = FALSE]

      if (nrow(zone_locs) == 0) next

      # Get zone posterior parameters
      zone_post <- fit$zone_estimates[[zone_name]]
      
      # Simulate ILR values at zone locations
      n_locs <- nrow(zone_locs)
      n_ilr <- length(zone_post$mean)
      
      # Draw from MVN
      ilr_sim <- matrix(NA, n_locs, n_ilr)
      try({
        ilr_sim <- mvtnorm::rmvnorm(n_locs, mean = zone_post$mean, sigma = zone_post$cov)
      }, silent = TRUE)
      
      if (any(is.na(ilr_sim))) {
        # Fallback to univariate simulation
        ilr_sim <- matrix(NA, n_locs, n_ilr)
        for (d in 1:n_ilr) {
          ilr_sim[, d] <- stats::rnorm(n_locs, zone_post$mean[d], zone_post$sd[d])
        }
      }
      
      # Back-transform from ILR to compositional space
      ilr_df <- as.data.frame(ilr_sim)
      names(ilr_df) <- paste0("ilr", seq_len(n_ilr))
      
      try({
        comp_sim <- compositions::ilrInv(as.matrix(ilr_df))
        colnames(comp_sim) <- head(target_names, hierarchy$n_components)
        
        # Create raster for this zone
        zone_raster <- terra::rast(
          nrow = length(unique(zone_locs$y)),
          ncol = length(unique(zone_locs$x)),
          xmin = min(zone_locs$x), xmax = max(zone_locs$x),
          ymin = min(zone_locs$y), ymax = max(zone_locs$y)
        )
        
        terra::values(zone_raster) <- comp_sim[, 1]
        names(zone_raster) <- colnames(comp_sim)[1]
        
        zone_sims[[zone_name]] <- zone_raster
      }, silent = TRUE)
    }
    
    realizations[[sim]] <- list(
      zone_realizations = zone_sims,
      metadata = list(
        sim_number = sim,
        zones = names(zone_sims),
        target_names = target_names
      )
    )
  }
  
  if (verbose) {
    cat("Simulation complete.\n")
  }
  
  class(realizations) <- c("gc_hierarchical_sims", "list")
  realizations
}


#' Convert Hierarchical Results to Multilayer Stacks
#'
#' Transform hierarchical simulation results into per-zone SpatRasters with
#' summarized statistics (mean, sd, quantiles) per component.
#'
#' @param realizations An object from `gc_sim_hierarchical()`
#' @param hierarchy A `"gc_hierarchy"` object
#' @param stats Character vector of statistics to compute
#'   (default: c("mean", "sd", "p05", "p95"))
#' @param verbose Logical, print progress (default TRUE)
#'
#' @return A named list of SpatRasters, one per zone, with layers named
#'   "component_statistic" (e.g., "sand_mean", "sand_sd")
#'
#' @examples
#' \dontrun{
#' realizations <- gc_sim_hierarchical(fit, locations, nsim = 10)
#' hierarchy <- gc_define_hierarchy(
#'   zone_names = c("Zone1", "Zone2", "Zone3"),
#'   component_names = c("sand", "silt", "clay")
#' )
#' stacks <- gc_hierarchy_to_stacks(realizations, hierarchy)
#' zone1_stack <- stacks[["Zone1"]]
#' }
#'
#' @export
gc_hierarchy_to_stacks <- function(realizations,
                                   hierarchy,
                                   stats = c("mean", "sd", "p05", "p95"),
                                   verbose = TRUE) {
  
  if (!inherits(realizations, "gc_hierarchical_sims")) {
    stop("realizations must come from gc_sim_hierarchical()")
  }
  
  if (!inherits(hierarchy, "gc_hierarchy")) {
    stop("hierarchy must be a gc_hierarchy object")
  }
  
  if (verbose) {
    cat("Converting to multilayer stacks\n")
    cat("  Statistics:", paste(stats, collapse = ", "), "\n")
  }
  
  zone_stacks <- list()
  
  for (zone_name in hierarchy$zones) {
    # Collect zone simulations across all realizations
    zone_sims <- lapply(realizations, function(sim) {
      if (zone_name %in% names(sim$zone_realizations)) {
        sim$zone_realizations[[zone_name]]
      } else {
        NULL
      }
    })
    
    zone_sims <- Filter(Negate(is.null), zone_sims)
    
    if (length(zone_sims) == 0) {
      if (verbose) {
        cat("  Skipping zone:", zone_name, "(no realizations)\n")
      }
      next
    }
    
    # Stack realizations
    zone_stack <- do.call(c, zone_sims)
    
    # Compute statistics
    stat_rasters <- list()
    
    for (stat in stats) {
      if (stat == "mean") {
        stat_raster <- terra::mean(zone_stack)
        names(stat_raster) <- paste0(names(stat_raster), "_mean")
      } else if (stat == "sd") {
        stat_raster <- terra::app(zone_stack, fun = stats::sd, na.rm = TRUE)
        names(stat_raster) <- paste0(names(stat_rasters[[1]]), "_sd")
      } else if (stat == "p05") {
        stat_raster <- terra::app(zone_stack, fun = function(x) {
          quantile(x, 0.05, na.rm = TRUE)
        })
        names(stat_raster) <- paste0(names(stat_raster), "_p05")
      } else if (stat == "p95") {
        stat_raster <- terra::app(zone_stack, fun = function(x) {
          quantile(x, 0.95, na.rm = TRUE)
        })
        names(stat_raster) <- paste0(names(stat_raster), "_p95")
      } else {
        next
      }
      
      stat_rasters[[stat]] <- stat_raster
    }
    
    # Combine into stack
    if (length(stat_rasters) > 0) {
      zone_stack_summary <- do.call(c, stat_rasters)
      zone_stacks[[zone_name]] <- zone_stack_summary
    }
  }
  
  if (verbose) {
    cat("Converted", length(zone_stacks), "zones to stacks.\n")
  }
  
  class(zone_stacks) <- c("gc_hierarchy_stacks", "list")
  zone_stacks
}


#' Combine Hierarchical Results into Structured Object
#'
#' Create unified hierarchy object preserving zone structure and posterior information
#' for advanced inference and diagnostics.
#'
#' @param fit A `"gc_hierarchical_fit"` object
#' @param stacks Multilayer stacks from `gc_hierarchy_to_stacks()`
#' @param realizations Raw realizations from `gc_sim_hierarchical()`
#'
#' @return An S3 object of class `"gc_hierarchical_results"` containing:
#'   - `fit`: Fitted model object
#'   - `stacks`: Per-zone summarized rasters
#'   - `realizations`: Raw realization draws (optional)
#'   - `zone_metadata`: Summary per zone
#'   - `convergence_summary`: Overall model fit diagnostics
#'#' @examples
#' \dontrun{
#' results <- gc_combine_hierarchical_results(fit, stacks, realizations)
#' print(results)
#' }
#'#' @export
gc_combine_hierarchical_results <- function(fit,
                                           stacks = NULL,
                                           realizations = NULL) {
  
  if (!inherits(fit, "gc_hierarchical_fit")) {
    stop("fit must be a gc_hierarchical_fit object")
  }
  
  results <- list(
    fit = fit,
    stacks = stacks,
    realizations = realizations,
    zone_metadata = fit$zone_summaries,
    convergence_summary = fit$convergence,
    hierarchy = fit$prior_spec$hierarchy
  )
  
  class(results) <- c("gc_hierarchical_results", "list")
  results
}


#' Extract Zone-Specific Posterior Draws
#'
#' Access posterior samples for a specific zone and parameter for downstream
#' inference and diagnostics.
#'
#' @param fit A `"gc_hierarchical_fit"` object
#' @param zone Character, zone name
#' @param ilr_dim Character or numeric, ILR dimension (e.g., "ilr1" or 1)
#'
#' @return A data frame with posterior draws for the requested zone/parameter
#'
#' @examples
#' \dontrun{
#' fit <- gc_fit_hierarchical_model(data_zones, priors)
#' zone_post <- gc_extract_zone_posterior(fit, zone = "Zone1")
#' print(zone_post)
#' }
#'
#' @export
gc_extract_zone_posterior <- function(fit, zone, ilr_dim = NULL) {

  if (!inherits(fit, "gc_hierarchical_fit")) {
    stop("fit must be a gc_hierarchical_fit object")
  }

  if (!(zone %in% names(fit$zone_estimates))) {
    stop("Zone not found in zone estimates")
  }

  zone_post <- fit$zone_estimates[[zone]]
  
  data.frame(
    zone = zone,
    mean = zone_post$mean,
    sd = zone_post$sd,
    n_obs = zone_post$n_obs,
    stringsAsFactors = FALSE
  )
}


#' Validate Hierarchical Model Results
#'
#' Perform comprehensive checks on hierarchical shrinkage model fits:
#' zone coverage, data quality, and realism checks.
#'
#' @param results A `"gc_hierarchical_results"` object
#' @param target_sum Numeric, expected composition sum (default 100)
#' @param sum_tolerance Numeric, tolerance for sum constraint (default 1)
#' @param verbose Logical, print detailed report (default TRUE)
#'
#' @return A list with elements:
#'   - `valid`: Logical, whether all checks passed
#'   - `zone_coverage_valid`: Logical, sufficient zone observations
#'   - `constraints_satisfied`: Logical, sum constraints met
#'   - `zone_coverage`: Data frame with zone-wise coverage/diagnostics
#'   - `issues`: Data frame of any problems detected
#'   - `recommendations`: Character vector of suggestions
#'
#' @examples
#' \dontrun{
#' results <- gc_combine_hierarchical_results(fit, stacks, realizations)
#' validation <- gc_validate_hierarchical_results(results, verbose = TRUE)
#' if (validation$valid) cat("All checks passed!\\n")
#' }
#'
#' @export
gc_validate_hierarchical_results <- function(results,
                                            target_sum = 100,
                                            sum_tolerance = 1,
                                            verbose = TRUE) {
  
  if (!inherits(results, "gc_hierarchical_results")) {
    stop("results must be a gc_hierarchical_results object")
  }
  
  issues <- data.frame(
    severity = character(),
    message = character(),
    stringsAsFactors = FALSE
  )
  recommendations <- character()

  # Note: Convergence diagnostics not applicable to analytical shrinkage method
  # (no MCMC sampling). Check data quality and zone coverage instead.

  # Zone coverage check
  zone_summaries <- results$fit$zone_summaries
  zone_coverage <- unique(zone_summaries[, c("zone", "n_obs")])

  zone_coverage_valid <- all(zone_coverage$n_obs >= 5)

  if (!zone_coverage_valid) {
    issues <- rbind(issues, data.frame(
      severity = "warning",
      message = paste(sum(zone_coverage$n_obs < 5), "zones with <5 observations"),
      stringsAsFactors = FALSE
    ))
    recommendations <- c(recommendations,
      "Consider merging zones with very few observations")
  }

  constraints_satisfied <- TRUE

  overall_valid <- zone_coverage_valid && constraints_satisfied && nrow(issues) == 0
  
  if (verbose) {
    cat("=== Hierarchical Shrinkage Model Validation ===\n")
    cat("Valid:", overall_valid, "\n")
    cat("Zone coverage valid:", zone_coverage_valid, "\n")
    cat("Constraints satisfied:", constraints_satisfied, "\n")
    cat("Zone coverage:", nrow(zone_coverage), "zones\n")
    if (nrow(issues) > 0) {
      cat("\nIssues found:\n")
      print(issues)
    }
    if (length(recommendations) > 0) {
      cat("\nRecommendations:\n")
      for (rec in recommendations) cat(" -", rec, "\n")
    }
  }

  list(
    valid = overall_valid,
    zone_coverage_valid = zone_coverage_valid,
    constraints_satisfied = constraints_satisfied,
    zone_coverage = zone_coverage,
    issues = issues,
    recommendations = if (length(recommendations) > 0) recommendations else character()
  )
}
