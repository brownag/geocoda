#' Compute Observation-Level Spatial Weights
#'
#' Weight observations based on proximity to zone boundaries and neighbor adjacency.
#' Observations closer to boundaries (where neighbor zones interact) receive higher weights,
#' enabling smoother spatial transitions while preserving zone interior characteristics.
#'
#' @param observations Simple feature (sf) with geometry, zone_id, mukey columns
#' @param mapunit_zones sf object with zone polygons
#' @param contrasts List with adjacency_strength and adjacency_matrix
#' @param distance_power Numeric, exponent for distance weighting (default 0.5)
#'   Higher values = stronger boundary emphasis (1.0 = linear, 0.5 = square-root)
#' @param direction_boost Numeric, bonus multiplier for neighbor-facing sides (default 1.3)
#' @param adjacency_scale Numeric, how much adjacency strength affects weights (default 0.5)
#' @param verbose Logical, print details
#'
#' @return sf object with new `spatial_weight` column
#'
#' @details
#' Weights combine three factors:
#' 1. **Distance to boundary**: Observations near zone edges get 1.5-2.0× boost
#' 2. **Direction alignment**: Observations facing neighbor zones get additional boost
#' 3. **Adjacency strength**: Stronger adjacency with neighbors increases weighting
#'
#' Combined weight formula:
#' weight = 1 + (boundary_distance_factor) × (1 + direction_factor + adjacency_factor)
#'
#' @examples
#' \dontrun{
#' obs_weighted <- compute_observation_spatial_weights(
#'   observations = obs_sf,
#'   mapunit_zones = zone_polygons,
#'   contrasts = adjacency_contrasts,
#'   distance_power = 0.5,
#'   direction_boost = 1.3
#' )
#' }
#'
#' @export
compute_observation_spatial_weights <- function(observations,
                                                mapunit_zones,
                                                contrasts,
                                                distance_power = 0.5,
                                                direction_boost = 1.3,
                                                adjacency_scale = 0.5,
                                                verbose = TRUE) {

  if (!inherits(observations, "sf")) {
    stop("observations must be sf object")
  }

  # Initialize weights to 1.0
  observations$spatial_weight <- 1.0

  if (verbose) {
    cat("Computing observation-level spatial weights...\n")
  }

  # For each zone, compute distance-to-boundary for observations
  unique_zones <- unique(observations$zone_id)

  for (zone_id in unique_zones) {
    zone_obs_idx <- observations$zone_id == zone_id
    zone_obs <- observations[zone_obs_idx, ]

    # Get zone polygon(s)
    zone_mukey <- unique(zone_obs$mukey)[1]
    zone_geom <- mapunit_zones[mapunit_zones$mukey == zone_mukey, ]

    if (nrow(zone_geom) == 0) next

    # Union geometries if multiple polygons
    zone_union <- sf::st_union(zone_geom)

    # Calculate distance to boundary for each observation
    distances_to_boundary <- sf::st_distance(
      zone_obs,
      sf::st_boundary(zone_union)
    )
    distances_to_boundary <- as.numeric(distances_to_boundary)

    # Normalize by zone size (max distance within zone)
    max_dist <- max(distances_to_boundary, na.rm = TRUE)
    if (max_dist > 0) {
      norm_distances <- distances_to_boundary / max_dist  # 0 at boundary, 1 at center
      boundary_boost <- 1 + (1 - norm_distances ^ distance_power)  # 1.0 at center, 2.0 at boundary
    } else {
      boundary_boost <- rep(1.0, length(distances_to_boundary))
    }

    # Direction weighting: calculate bearing from zone centroid to each neighbor
    zone_centroid <- sf::st_centroid(zone_union)
    zone_obs_bearings <- atan2(
      sf::st_coordinates(zone_obs)[, 2] - sf::st_coordinates(zone_centroid)[2],
      sf::st_coordinates(zone_obs)[, 1] - sf::st_coordinates(zone_centroid)[1]
    )

    # Find neighbors and their average direction
    neighbor_directions <- numeric(nrow(zone_obs))
    adjacency_strength_vals <- numeric(nrow(zone_obs))

    if (!is.null(contrasts$adjacency_matrix)) {
      # Find zone index in adjacency matrix
      zone_idx <- which(rownames(contrasts$adjacency_matrix) %in% zone_mukey)

      if (length(zone_idx) > 0 && zone_idx <= nrow(contrasts$adjacency_matrix)) {
        # Get neighbors
        neighbor_mask <- contrasts$adjacency_matrix[zone_idx, ] > 0
        neighbor_indices <- which(neighbor_mask)

        if (length(neighbor_indices) > 0) {
          # Average strength of adjacency to all neighbors
          mean_adjacency_strength <- mean(
            contrasts$adjacency_strength[zone_idx, neighbor_indices],
            na.rm = TRUE
          )
          adjacency_strength_vals[] <- mean_adjacency_strength

          # Calculate average neighbor direction
          neighbor_geoms <- mapunit_zones[neighbor_indices, ]
          neighbor_centroids <- sf::st_centroid(neighbor_geoms)
          avg_neighbor_bearing <- atan2(
            mean(sf::st_coordinates(neighbor_centroids)[, 2], na.rm = TRUE) -
              sf::st_coordinates(zone_centroid)[2],
            mean(sf::st_coordinates(neighbor_centroids)[, 1], na.rm = TRUE) -
              sf::st_coordinates(zone_centroid)[1]
          )

          # Angle difference (0 = facing neighbor, π = away from neighbor)
          angle_diffs <- abs(zone_obs_bearings - avg_neighbor_bearing)
          angle_diffs <- pmin(angle_diffs, 2*pi - angle_diffs)  # Shortest angle

          # Weight observations facing neighbors higher
          direction_boost_factors <- 1 - (angle_diffs / pi) * (direction_boost - 1)
          direction_boost_factors <- pmax(1.0, direction_boost_factors)
        } else {
          direction_boost_factors <- rep(1.0, nrow(zone_obs))
        }
      } else {
        direction_boost_factors <- rep(1.0, nrow(zone_obs))
      }
    } else {
      direction_boost_factors <- rep(1.0, nrow(zone_obs))
    }

    # Adjacency strength factor
    adjacency_factors <- 1 + adjacency_strength_vals * adjacency_scale

    # Combine all factors
    combined_weights <- boundary_boost * direction_boost_factors * adjacency_factors

    observations$spatial_weight[zone_obs_idx] <- combined_weights
  }

  if (verbose) {
    cat("  Weight range:", round(min(observations$spatial_weight), 3), "to",
        round(max(observations$spatial_weight), 3), "\n")
    cat("  Mean weight:", round(mean(observations$spatial_weight), 3), "\n")
  }

  observations
}


#' Compute Per-Zone Kriging Range Adjustments
#'
#' Calculate spatially-varying variogram range multipliers based on zone characteristics.
#' Zones with stronger neighbor connections get longer ranges (smoother spatial transitions),
#' while isolated zones maintain shorter ranges (preserve unique characteristics).
#'
#' @param zone_properties Data frame with zone statistics (n_neighbors, adjacency strength)
#' @param adjacency List with adjacency_matrix and adjacency_strength
#' @param base_range Numeric, baseline kriging range in degrees (default 0.01)
#' @param neighbor_count_effect Numeric, how much neighbor count influences range (default 0.3)
#'   Higher values = more range variation based on connectivity
#' @param adjacency_strength_effect Numeric, how much composition similarity affects range (default 0.2)
#'   Higher values = more range variation based on adjacency similarity
#' @param verbose Logical, print details
#'
#' @return Data frame with zone_id, range_multiplier, range_actual columns
#'
#' @details
#' Range multiplier formula:
#' multiplier = 1 + (neighbor_count_effect × n_neighbors / max_neighbors)
#'            + (adjacency_strength_effect × mean_adjacency_strength)
#'
#' Isolated zones: multiplier ≈ 1.0 (narrow range, sharp transitions)
#' Well-connected zones: multiplier ≈ 1.5-2.0 (wider range, smooth transitions)
#'
#' @examples
#' \dontrun{
#' adjustments <- compute_kriging_range_adjustments(
#'   zone_properties = zones,
#'   adjacency = adjacency_obj,
#'   base_range = 0.01,
#'   neighbor_count_effect = 0.3,
#'   adjacency_strength_effect = 0.2
#' )
#' }
#'
#' @export
compute_kriging_range_adjustments <- function(zone_properties,
                                              adjacency,
                                              base_range = 0.01,
                                              neighbor_count_effect = 0.3,
                                              adjacency_strength_effect = 0.2,
                                              verbose = TRUE) {

  if (is.null(adjacency$adjacency_matrix) || is.null(adjacency$adjacency_strength)) {
    warning("adjacency matrix or strength missing; returning base ranges")
    return(data.frame(
      zone_id = seq_len(nrow(zone_properties)),
      range_multiplier = 1.0,
      range_actual = base_range
    ))
  }

  n_zones <- nrow(adjacency$adjacency_matrix)
  multipliers <- numeric(n_zones)
  neighbor_counts <- numeric(n_zones)
  adjacency_strengths <- numeric(n_zones)

  # For each zone, compute multiplier based on neighbors
  for (i in seq_len(n_zones)) {
    # Count neighbors
    n_neighbors <- sum(adjacency$adjacency_matrix[i, ] > 0)
    neighbor_counts[i] <- n_neighbors

    # Mean adjacency strength to neighbors
    neighbor_idx <- which(adjacency$adjacency_matrix[i, ] > 0)
    if (length(neighbor_idx) > 0) {
      mean_strength <- mean(adjacency$adjacency_strength[i, neighbor_idx], na.rm = TRUE)
    } else {
      mean_strength <- 0
    }
    adjacency_strengths[i] <- mean_strength

    # Calculate multiplier
    max_neighbors <- max(neighbor_counts, na.rm = TRUE)
    if (max_neighbors > 0) {
      neighbor_contribution <- (n_neighbors / max_neighbors) * neighbor_count_effect
    } else {
      neighbor_contribution <- 0
    }

    strength_contribution <- mean_strength * adjacency_strength_effect
    multipliers[i] <- 1.0 + neighbor_contribution + strength_contribution
  }

  result <- data.frame(
    zone_id = seq_len(n_zones),
    n_neighbors = neighbor_counts,
    mean_adjacency_strength = adjacency_strengths,
    range_multiplier = multipliers,
    range_actual = base_range * multipliers
  )

  if (verbose) {
    cat("Kriging range adjustments:\n")
    cat("  Neighbor count range:", min(neighbor_counts), "to", max(neighbor_counts), "\n")
    cat("  Range multiplier:", round(min(multipliers), 3), "to",
        round(max(multipliers), 3), "\n")
    cat("  Actual range:", round(min(result$range_actual), 6), "to",
        round(max(result$range_actual), 6), "degrees\n")
  }

  result
}
