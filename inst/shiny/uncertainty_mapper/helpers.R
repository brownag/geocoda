suppressPackageStartupMessages({
  library(geocoda)
  library(terra)
  library(sf)
  library(ggplot2)
  library(plotly)
})

load_example_data <- function() {
  set.seed(42)
  locations <- expand.grid(
    x = seq(0, 500, by = 125),
    y = seq(0, 500, by = 125),
    z = c(10, 30, 50, 70, 90)
  )
  nrows <- nrow(locations)
  locations$sand <- 40 + sin(locations$x / 200) * 8 + rnorm(nrows, sd = 4)
  locations$silt <- 35 - cos(locations$y / 200) * 6 + rnorm(nrows, sd = 4)
  locations$clay <- 100 - locations$sand - locations$silt
  locations$sand <- pmax(0, pmin(100, locations$sand))
  locations$silt <- pmax(0, pmin(100, locations$silt))
  locations$clay <- 100 - locations$sand - locations$silt
  locations$zone <- ifelse(
    locations$y > locations$x,
    ifelse(locations$x < 250, "Upslope_North", "Midslope"),
    ifelse(locations$y < 250, "Downslope", "Midslope")
  )
  return(locations)
}

validate_composition <- function(data, comp_cols) {
  issues <- c()
  if (!all(comp_cols %in% names(data))) {
    issues <- c(issues, "Some composition columns not found in data")
    return(issues)
  }
  comp_data <- data[, comp_cols, drop = FALSE]
  if (any(comp_data < 0, na.rm = TRUE)) {
    neg_count <- sum(comp_data < 0, na.rm = TRUE)
    issues <- c(issues, paste("Found", neg_count, "negative values"))
  }
  row_sums <- rowSums(comp_data, na.rm = TRUE)
  if (any(abs(row_sums - 100) > 1, na.rm = TRUE)) {
    bad_sums <- sum(abs(row_sums - 100) > 1, na.rm = TRUE)
    issues <- c(issues, paste(bad_sums, "rows don't sum to ~100%"))
  }
  if (any(is.na(comp_data))) {
    na_count <- sum(is.na(comp_data))
    issues <- c(issues, paste("Found", na_count, "missing values"))
  }
  return(issues)
}

validate_spatial <- function(data, x_col, y_col, z_col = NULL, zone_col = NULL) {
  issues <- c()
  required_cols <- c(x_col, y_col)
  if (!is.null(z_col)) required_cols <- c(required_cols, z_col)
  if (!is.null(zone_col)) required_cols <- c(required_cols, zone_col)
  missing <- setdiff(required_cols, names(data))
  if (length(missing) > 0) {
    issues <- c(issues, paste("Missing columns:", paste(missing, collapse = ", ")))
    return(issues)
  }
  if (!is.numeric(data[[x_col]]) || !is.numeric(data[[y_col]])) {
    issues <- c(issues, "X and Y coordinates must be numeric")
  }
  if (!is.null(z_col) && !is.numeric(data[[z_col]])) {
    issues <- c(issues, "Z coordinate must be numeric")
  }
  if (any(is.na(data[[x_col]])) || any(is.na(data[[y_col]]))) {
    issues <- c(issues, "Missing values in spatial coordinates")
  }
  if (!is.null(zone_col)) {
    zone_counts <- table(data[[zone_col]])
    if (any(zone_counts < 3)) {
      small_zones <- names(zone_counts)[zone_counts < 3]
      issues <- c(issues, paste("Zones with < 3 obs:", paste(small_zones, collapse = ", ")))
    }
  }
  return(issues)
}

prepare_hierarchical_data <- function(data, zone_col, comp_cols) {
  comp_data <- data[, comp_cols, drop = FALSE]
  ilr_data <- data.frame(
    zone = data[[zone_col]],
    ilr1 = log(comp_data[, 1] / sqrt(comp_data[, 2] * comp_data[, 3] + 1e-6) + 1e-6),
    ilr2 = log(comp_data[, 2] / (comp_data[, 3] + 1e-6) + 1e-6),
    x = data$x,
    y = data$y,
    z = if ("z" %in% names(data)) data$z else NA
  )
  return(ilr_data)
}

fit_zone_model <- function(data, backend = "analytical_3d", params = list()) {
  ilr_cols <- grep("^ilr", names(data), value = TRUE)
  zones <- unique(data$zone)
  hierarchy <- list(
    zones = as.character(zones),
    n_zones = length(zones),
    n_components = length(ilr_cols) + 1
  )
  class(hierarchy) <- "gc_hierarchy"
  prior_spec <- list(
    hierarchy = hierarchy,
    global_mean = colMeans(data[, ilr_cols, drop = FALSE], na.rm = TRUE),
    global_covariance = cov(data[, ilr_cols, drop = FALSE], use = "complete.obs"),
    pooling_coefficient = params$pooling_coef %||% 0.4,
    prior_sd_mean = rep(1, length(ilr_cols)),
    prior_shape_variance = 2,
    prior_rate_variance = 2
  )
  class(prior_spec) <- "gc_prior_spec"
  zone_centers <- aggregate(cbind(x, y, z) ~ zone, data, FUN = mean)
  fit_result <- geocoda::fit_hierarchical_backend(
    data = data,
    prior_spec = prior_spec,
    backend = backend,
    zone_centers = zone_centers,
    verbose = TRUE
  )
  return(list(
    fit = fit_result,
    prior_spec = prior_spec,
    ilr_cols = ilr_cols,
    zones = zones
  ))
}

simulate_realizations <- function(fit_result, n_sims = 100) {
  zone_estimates <- fit_result$fit$zone_estimates
  zones <- names(zone_estimates)
  sims <- list()
  for (z in zones) {
    est <- zone_estimates[[z]]
    sims[[z]] <- mvtnorm::rmvnorm(
      n = n_sims,
      mean = est$mean,
      sigma = est$cov
    )
  }
  return(sims)
}

ilr_to_composition <- function(ilr_sims, n_components = 3) {
  n_sims <- nrow(ilr_sims)
  ilr1 <- ilr_sims[, 1]
  ilr2 <- ilr_sims[, 2]
  e1 <- exp(ilr1)
  e2 <- exp(ilr2)
  denom <- 1 + e1 + e2
  sand <- 100 * e1 / denom
  silt <- 100 * e2 / denom
  clay <- 100 / denom
  sand <- pmax(0, pmin(100, sand))
  silt <- pmax(0, pmin(100, silt))
  clay <- 100 - sand - silt
  return(data.frame(sand = sand, silt = silt, clay = clay))
}

probability_map <- function(sims_by_zone, variable, threshold, operator = "gt") {
  prob_results <- list()
  for (z in names(sims_by_zone)) {
    sims <- sims_by_zone[[z]]
    if (operator == "gt") {
      exceed <- rowMeans(sims[, variable] > threshold)
    } else if (operator == "lt") {
      exceed <- rowMeans(sims[, variable] < threshold)
    } else {
      exceed <- rowMeans(sims[, variable] > threshold[1] & sims[, variable] < threshold[2])
    }
    prob_results[[z]] <- exceed
  }
  return(prob_results)
}

percentile_map <- function(sims_by_zone, variable, percentiles = c(10, 50, 90)) {
  pct_results <- list()
  for (z in names(sims_by_zone)) {
    sims <- sims_by_zone[[z]]
    pcts <- quantile(sims[, variable], probs = percentiles / 100)
    pct_results[[z]] <- pcts
  }
  return(pct_results)
}

plot_prob_map <- function(prob_data, zones_sf) {
  prob_df <- data.frame(
    zone = names(prob_data),
    probability = unlist(prob_data)
  )
  zones_sf$probability <- prob_df$probability[match(zones_sf$zone, prob_df$zone)]
  p <- ggplot() +
    geom_sf(data = zones_sf, aes(fill = probability), color = "black", size = 0.5) +
    scale_fill_viridis_c(name = "Probability", limits = c(0, 1)) +
    theme_minimal() +
    labs(title = "Probability Map")
  return(p)
}

plot_depth_profile <- function(data_by_depth, variable) {
  depth_stats <- data.frame(
    depth = names(data_by_depth),
    mean = sapply(data_by_depth, function(x) mean(x[[variable]])),
    sd = sapply(data_by_depth, function(x) sd(x[[variable]])),
    q10 = sapply(data_by_depth, function(x) quantile(x[[variable]], 0.1)),
    q90 = sapply(data_by_depth, function(x) quantile(x[[variable]], 0.9))
  )
  depth_stats$depth <- as.numeric(depth_stats$depth)
  p <- ggplot(depth_stats, aes(x = mean, y = -depth)) +
    geom_ribbon(aes(xmin = q10, xmax = q90), alpha = 0.3) +
    geom_line(size = 1) +
    geom_point(size = 3) +
    theme_minimal() +
    labs(
      title = paste("Depth Profile:", variable),
      x = paste(variable, "percentage"),
      y = "Depth (cm)"
    )
  return(p)
}

`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}

format_model_params <- function(backend, pooling, spatial_decay) {
  paste("Backend:", backend, "|", "Pooling:", round(pooling, 2), "|", "Decay:", spatial_decay)
}

zone_summary_table <- function(fit_result) {
  summaries <- fit_result$zone_summaries
  aggregate(
    cbind(posterior_mean, posterior_sd) ~ zone,
    data = summaries,
    FUN = mean
  )
}
