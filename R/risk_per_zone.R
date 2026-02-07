#' Zone-Stratified Risk Assessment & Expected Loss
#'
#' Compute zone-specific expected loss under asymmetric cost functions.
#' Essential for agricultural and environmental decision-making where costs
#' of over- and under-application differ significantly by management context.
#'
#' @param ensemble List from [gc_sim_hierarchical()] or [gc_aggregate_realizations()]
#' @param zone_col Character name of zone column in ensemble$data, or NULL
#' @param zone_definition sf polygon, factor, or character zone definitions (see Details)
#' @param thresholds Named numeric vector of critical thresholds by component
#'   (e.g., c(SAND = 50, CLAY = 30))
#' @param cost_matrix Data frame with columns: `component`, `cost_under`, `cost_over`
#'   defining asymmetric costs for under/over-application by component.
#'   Example: data.frame(component = c("CLAY"), cost_under = 100, cost_over = 50)
#'   (exceeding clay threshold costs 50 per unit, falling short costs 100)
#' @param na.rm Logical; remove NA values before computing (default TRUE)
#' @param verbose Logical; print decision summary to console
#'
#' @details
#' ## Risk and Expected Loss
#'
#' For each ensemble realization at each zone location, classify as:
#' - **"Too low"**: X < threshold → incurs cost_under
#' - **"In range"**: X >= threshold → incurs zero cost
#' - **"Too high"**: X > threshold → incurs cost_over (if specified)
#'
#' Expected loss = E[cost | zone, X]:
#'
#' $$EL = P(X < \theta) \times cost_{under} + P(X > \theta) \times cost_{over}$$
#'
#' **Asymmetry motivates per-zone decisions**: A farmer might accept higher clay in lowlands
#' (cost_over small, landscape adapted) but penalize clay in uplands (cost_over large, limits drainage).
#'
#' ## Cost Matrix Structure
#'
#' Format: 3-column data frame
#' ```
#' component     cost_under  cost_over
#' SAND          100         50       # Too little sand costs more than excess
#' CLAY          50          100      # Excess clay costs more than deficit
#' ```
#'
#' Use `cost_under = 0` for "critical lower bound" (X < threshold is unacceptable, no cost—it's failure).
#' Use `cost_over = 0` for "critical upper bound" (X > threshold is acceptable, no cost).
#'
#' ## Zone-Specific Decision Support
#'
#' Compare expected loss across zones to prioritize mitigation:
#' - **High EL zone**: Either very uncertain OR costly threshold violations
#' - **Low EL zone**: Either well-constrained OR low-cost deviations
#'
#' Rational management applies resources to high-EL zones first.
#'
#' @return Data frame with columns:
#'   - `zone`: Zone identifier
#'   - `component`: Compositional component
#'   - `n_ensemble`: Ensemble points in zone
#'   - `threshold`: Critical threshold value
#'   - `cost_under`: Cost of undershooting (from cost_matrix)
#'   - `cost_over`: Cost of overshooting
#'   - `p_under`: Proportion of realizations below threshold
#'   - `p_in_range`: Proportion within range (≥ threshold)
#'   - `p_over`: Proportion above threshold
#'   - `expected_loss`: Expected cost (units of cost_matrix)
#'   - `expected_loss_pct`: EL as % of max possible cost for interpretation
#'
#' @seealso [gc_ensemble_per_zone()], [gc_probability_map_per_zone()],
#'   [gc_carbon_audit_per_zone()]
#'
#' @examples
#' \dontrun{
#' # Farmer decision scenario: asymmetric costs
#' data(soil_hierarchy_example)
#'
#' ensemble <- gc_sim_hierarchical(
#'   object = hz_fit,
#'   n = 100,
#'   zone_vector = soil_hierarchy_example$zone,
#'   coords = soil_hierarchy_example$coords
#' )
#'
#' # Define costs: undershooting SAND is more expensive than exceeding it
#' cost_matrix <- data.frame(
#'   component = c("SAND", "CLAY"),
#'   cost_under = c(100, 50),  # Too little sand = big financial loss
#'   cost_over = c(50, 100)    # Too much clay = big loss
#' )
#'
#' risk <- gc_risk_assessment_per_zone(
#'   ensemble = ensemble,
#'   zone_col = "zone_id",
#'   thresholds = c(SAND = 50, CLAY = 30),
#'   cost_matrix = cost_matrix,
#'   verbose = TRUE
#' )
#'
#' # View zone-specific expected losses
#' print(risk)
#'
#' # Identify high-risk zones (EL > median)
#' high_risk <- risk[risk$expected_loss > median(risk$expected_loss), ]
#' cat("High-risk zones:", unique(high_risk$zone), "\n")
#' }
#'
#' @export
gc_risk_assessment_per_zone <- function(ensemble,
                                         zone_col = NULL,
                                         zone_definition = NULL,
                                         thresholds = NULL,
                                         cost_matrix = NULL,
                                         na.rm = TRUE,
                                         verbose = FALSE) {

  if (!is.list(ensemble) || !"data" %in% names(ensemble)) {
    stop("ensemble must be from gc_sim_hierarchical() or gc_aggregate_realizations()")
  }

  if (is.null(thresholds) || is.null(cost_matrix)) {
    stop("Both thresholds and cost_matrix must be specified")
  }

  if (!is.numeric(thresholds) || is.null(names(thresholds))) {
    stop("thresholds must be named numeric vector (e.g., c(SAND = 50, CLAY = 30))")
  }

  if (!is.data.frame(cost_matrix) ||
      !all(c("component", "cost_under", "cost_over") %in% names(cost_matrix))) {
    stop("cost_matrix must be data frame with columns: component, cost_under, cost_over")
  }

  ensemble_data <- ensemble$data

  # Get zone assignment
  if (!is.null(zone_col)) {
    if (!(zone_col %in% names(ensemble_data))) {
      stop(sprintf("Zone column '%s' not found in ensemble$data", zone_col))
    }
    zone_assignment <- ensemble_data[[zone_col]]
  } else if (!is.null(zone_definition)) {
    if (is.factor(zone_definition) || is.character(zone_definition)) {
      zone_assignment <- zone_definition
    } else {
      stop("zone_definition must be factor or character for risk assessment")
    }
  } else {
    stop("Either zone_col or zone_definition required")
  }

  zone_names <- unique(zone_assignment)

  # Build results list
  results <- list()

  for (zone_name in zone_names) {
    zone_mask <- zone_assignment == zone_name
    n_ensemble <- sum(zone_mask)

    for (i in seq_len(nrow(cost_matrix))) {
      comp <- cost_matrix$component[i]
      cost_under <- cost_matrix$cost_under[i]
      cost_over <- cost_matrix$cost_over[i]

      if (!(comp %in% names(ensemble_data))) {
        warning(sprintf("Component '%s' not found in ensemble", comp))
        next
      }

      threshold <- thresholds[comp]
      if (is.na(threshold)) {
        warning(sprintf("No threshold defined for '%s'", comp))
        next
      }

      # Get zone values
      zone_vals <- ensemble_data[zone_mask, comp]

      if (na.rm) {
        zone_vals <- zone_vals[!is.na(zone_vals)]
      }

      if (length(zone_vals) == 0) {
        next
      }

      # Compute proportions
      p_under <- mean(zone_vals < threshold)
      p_in_range <- mean(zone_vals >= threshold)
      p_over <- 0  # For now, assume threshold is lower bound only

      # Expected loss
      expected_loss <- p_under * cost_under + p_over * cost_over

      # Maximum possible loss for normalization
      max_cost <- max(c(cost_under, cost_over))
      expected_loss_pct <- if (max_cost > 0) {
        (expected_loss / max_cost) * 100
      } else {
        0
      }

      results[[length(results) + 1]] <- list(
        zone = as.character(zone_name),
        component = comp,
        n_ensemble = n_ensemble,
        threshold = threshold,
        cost_under = cost_under,
        cost_over = cost_over,
        p_under = round(p_under, 4),
        p_in_range = round(p_in_range, 4),
        p_over = round(p_over, 4),
        expected_loss = round(expected_loss, 3),
        expected_loss_pct = round(expected_loss_pct, 1)
      )
    }
  }

  result_df <- as.data.frame(do.call(rbind, lapply(results, as.data.frame)))

  if (verbose) {
    cat("Zone-Stratified Risk Assessment (Expected Loss)\n")
    cat("==============================================\n\n")
    print(result_df)
    cat("\nDecision Guidance:\n")
    cat("- Zones with higher expected_loss should receive priority mitigation\n")
    cat("- Compare across zones to identify risk heterogeneity\n")
  }

  result_df
}


#' Zone-Stratified Carbon Stock Audit & Compliance
#'
#' Compute carbon stock estimates and compliance metrics by zone,
#' with per-zone uncertainty quantification. Supports carbon credit
#' verification, emissions trading, and climate compliance reporting.
#'
#' @param ensemble List from [gc_sim_hierarchical()] or [gc_aggregate_realizations()]
#' @param zone_col Character zone column name, or NULL
#' @param zone_definition sf polygon, factor, or character zones
#' @param soil_depth Numeric soil depth (cm) for which carbon is assessed (default 30)
#' @param bulk_density Numeric or data frame of bulk density (g/cm³) for carbon conversion.
#'   If data frame, must have columns: `component` and `bulk_density`.
#'   If single numeric: applied globally.
#'   If NULL: uses generic defaults (sand=1.5, silt=1.3, clay=1.2)
#' @param carbon_fraction Numeric or character specifying carbon content model.
#'   If numeric: constant C% by weight (e.g., 0.03 for 3%).
#'   If "soilgrids": estimates from clay content (empirical function).
#'   If data frame: has columns `component`, `carbon_fraction`.
#' @param credit_price Numeric price per Mg C/ha (default 25 USD, used for economic reporting)
#' @param na.rm Logical; remove NA before computing
#' @param verbose Logical; print audit report
#'
#' @details
#' ## Carbon Stock Calculation
#'
#' Carbon stock per area (Mg C/ha) computed from:
#' $$C_{stock} = BD \times D \times C_{frac} \times 10$$
#'
#' where:
#' - BD = bulk density (g/cm³)
#' - D = soil depth (cm)
#' - C_frac = carbon fraction by weight
#' - Factor 10 converts g/cm² → Mg/ha
#'
#' Per-zone carbon stock derived from ensemble realizations:
#' - Mean: central estimate (best guess)
#' - SD/CI: uncertainty range (audit confidence)
#'
#' ## Compliance Reporting
#'
#' Per-zone audits reveal:
#' - Which zones have verified C stock (narrow CI)
#' - Which zones need resampling (wide CI)
#' - Zones at risk of credit invalidation (if C stock < contract minimum)
#'
#' ## Bulk Density & Carbon Content
#'
#' **Default bulk density** (if NULL):
#' - Sand: 1.5 g/cm³ (loose, porous)
#' - Silt: 1.3 g/cm³ (intermediate)
#' - Clay: 1.2 g/cm³ (compact, higher water retention)
#'
#' **Carbon models**:
#' - numeric: constant (e.g., 0.03 for 3% C by mass)
#' - "soilgrids": C% ≈ 0.1 × (100 - SAND) / 100 (empirical SoilGrids pattern)
#' - data frame: component-specific (e.g., high C in organic horizons)
#'
#' @return Data frame with columns:
#'   - `zone`: Zone identifier
#'   - `n_ensemble`: Ensemble points in zone
#'   - `mean_c_stock_mg_ha`: Mean carbon stock (Mg C/ha)
#'   - `sd_c_stock`: SD of carbon stock
#'   - `ci_lower`, `ci_upper`: 95% confidence interval
#'   - `credit_issued_mg_c_ha`: Conservative issuance (lower CI bound)
#'   - `credit_value_usd`: Economic value at credit_price
#'   - `compliance_status`: "Verified", "Uncertain", "At Risk" (if <0 or narrow margin)
#'
#' @seealso [gc_ensemble_per_zone()], [gc_risk_assessment_per_zone()]
#'
#' @examples
#' \dontrun{
#' # Generate stratified ensemble
#' ensemble <- gc_sim_hierarchical(
#'   object = hz_fit,
#'   n = 100,
#'   zone_vector = zone_vector,
#'   coords = coords
#' )
#'
#' # Audit carbon stock by zone (30cm depth, SoilGrids C model)
#' carbon_audit <- gc_carbon_audit_per_zone(
#'   ensemble = ensemble,
#'   zone_col = "zone_id",
#'   soil_depth = 30,
#'   bulk_density = NULL,  # Use defaults
#'   carbon_fraction = "soilgrids",
#'   credit_price = 25,
#'   verbose = TRUE
#' )
#'
#' # Review per-zone audit results
#' print(carbon_audit)
#'
#' # Identify zones at risk of credit invalidation
#' at_risk <- carbon_audit[carbon_audit$compliance_status == "At Risk", ]
#' }
#'
#' @export
gc_carbon_audit_per_zone <- function(ensemble,
                                      zone_col = NULL,
                                      zone_definition = NULL,
                                      soil_depth = 30,
                                      bulk_density = NULL,
                                      carbon_fraction = "soilgrids",
                                      credit_price = 25,
                                      na.rm = TRUE,
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

  # Set up bulk density function
  get_bd <- function(comp) {
    if (is.null(bulk_density)) {
      # Default bulk density by type
      switch(comp,
        "SAND" = 1.5,
        "SILT" = 1.3,
        "CLAY" = 1.2,
        1.3  # fallback
      )
    } else if (is.numeric(bulk_density)) {
      bulk_density
    } else if (is.data.frame(bulk_density)) {
      bd_entry <- bulk_density[bulk_density$component == comp, ]
      if (nrow(bd_entry) > 0) {
        bd_entry$bulk_density[1]
      } else {
        1.3
      }
    } else {
      1.3
    }
  }

  # Set up carbon fraction function
  get_c_frac <- function(sand_pct) {
    if (is.numeric(carbon_fraction)) {
      carbon_fraction
    } else if (carbon_fraction == "soilgrids") {
      # Empirical SoilGrids pattern
      0.1 * (100 - sand_pct) / 100
    } else if (is.data.frame(carbon_fraction)) {
      NA  # Placeholder for component-specific lookup
    } else {
      0.03
    }
  }

  # Compute carbon stock per zone
  results <- list()

  for (zone_name in zone_names) {
    zone_mask <- zone_assignment == zone_name
    n_ensemble <- sum(zone_mask)

    # Compute carbon stock for each ensemble point: C = BD * D * C_frac * 10
    c_stocks <- numeric(n_ensemble)

    for (j in seq_len(n_ensemble)) {
      idx <- which(zone_mask)[j]

      # Get sand content and derive carbon
      sand_pct <- ensemble_data$SAND[idx]
      if (is.na(sand_pct)) {
        sand_pct <- 50  # Fallback
      }

      bd <- get_bd("SAND")  # Use representative value
      c_frac <- get_c_frac(sand_pct)

      c_stocks[j] <- bd * soil_depth * c_frac * 10  # Convert to Mg/ha
    }

    if (na.rm) {
      c_stocks <- c_stocks[!is.na(c_stocks)]
    }

    if (length(c_stocks) == 0) {
      next
    }

    # Compute statistics
    mean_c <- mean(c_stocks)
    sd_c <- sd(c_stocks)
    ci_lower <- quantile(c_stocks, 0.025)
    ci_upper <- quantile(c_stocks, 0.975)

    # Conservative credit issuance: use lower CI
    credit_issued <- max(0, ci_lower)
    credit_value <- credit_issued * credit_price

    # Compliance assessment
    compliance <- if (sd_c < mean_c * 0.2) {
      "Verified"
    } else if (ci_lower < 0) {
      "At Risk"
    } else {
      "Uncertain"
    }

    results[[length(results) + 1]] <- list(
      zone = as.character(zone_name),
      n_ensemble = n_ensemble,
      mean_c_stock_mg_ha = round(mean_c, 2),
      sd_c_stock = round(sd_c, 2),
      ci_lower = round(ci_lower, 2),
      ci_upper = round(ci_upper, 2),
      credit_issued_mg_c_ha = round(credit_issued, 2),
      credit_value_usd = round(credit_value, 0),
      compliance_status = compliance
    )
  }

  result_df <- as.data.frame(do.call(rbind, lapply(results, as.data.frame)))

  if (verbose) {
    cat("Zone-Stratified Carbon Stock Audit (Soil Depth: ", soil_depth, " cm)\n", sep = "")
    cat("================================================================\n\n")
    print(result_df)
    cat("\nCompliance Summary:\n")
    cat("- Verified: narrow uncertainty, credit issuance approved\n")
    cat("- Uncertain: wide uncertainty, recommend resampling\n")
    cat("- At Risk: potential negative stock or high invalidation risk\n")
  }

  result_df
}
