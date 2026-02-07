suppressPackageStartupMessages({
  library(shiny)
  library(shinydashboard)
  library(plotly)
  library(ggplot2)
  library(terra)
  library(sf)
  library(mvtnorm)
  library(geocoda)
})

server <- function(input, output, session) {
  app_data <- reactiveValues(
    raw_data = NULL,
    validated_data = NULL,
    hierarchical_fit = NULL,
    simulations = NULL,
    zone_sf = NULL,
    sim_status = "Ready"
  )

  observeEvent(input$load_example, {
    app_data$raw_data <- load_example_data()
    updateSelectInput(session, "zone_col", choices = names(app_data$raw_data))
    updateSelectInput(session, "x_col", choices = names(app_data$raw_data))
    updateSelectInput(session, "y_col", choices = names(app_data$raw_data))
    updateSelectInput(session, "z_col", choices = names(app_data$raw_data))
    updateSelectInput(session, "comp_sand_col", choices = names(app_data$raw_data))
    updateSelectInput(session, "comp_silt_col", choices = names(app_data$raw_data))
    updateSelectInput(session, "comp_clay_col", choices = names(app_data$raw_data))
  })

  observeEvent(input$data_file, {
    tryCatch(
      {
        app_data$raw_data <- read.csv(input$data_file$datapath)
        col_names <- names(app_data$raw_data)
        updateSelectInput(session, "zone_col", choices = col_names)
        updateSelectInput(session, "x_col", choices = col_names)
        updateSelectInput(session, "y_col", choices = col_names)
        if ("z" %in% col_names || "depth" %in% col_names) {
          z_col <- if ("z" %in% col_names) "z" else "depth"
          updateSelectInput(session, "z_col", choices = col_names, selected = z_col)
        } else {
          updateSelectInput(session, "z_col", choices = col_names)
        }
        comp_candidates <- c("sand", "silt", "clay", "sand_pct", "silt_pct", "clay_pct")
        comp_cols <- intersect(comp_candidates, tolower(col_names))
        
        if (length(comp_cols) >= 3) {
          updateSelectInput(session, "comp_sand_col", choices = col_names, selected = comp_cols[1])
          updateSelectInput(session, "comp_silt_col", choices = col_names, selected = comp_cols[2])
          updateSelectInput(session, "comp_clay_col", choices = col_names, selected = comp_cols[3])
        } else {
          updateSelectInput(session, "comp_sand_col", choices = col_names)
          updateSelectInput(session, "comp_silt_col", choices = col_names)
          updateSelectInput(session, "comp_clay_col", choices = col_names)
        }
      },
      error = function(e) {
        showNotification(paste("Error loading file:", e$message), type = "error")
      }
    )
  })

  output$data_preview <- renderTable({
    if (is.null(app_data$raw_data)) return(NULL)
    head(app_data$raw_data, 10)
  })

  output$data_summary <- renderText({
    if (is.null(app_data$raw_data)) return("No data loaded")
    
    sprintf(
      "Rows: %d\nColumns: %d\nZones: %d\nColumns: %s",
      nrow(app_data$raw_data),
      ncol(app_data$raw_data),
      length(unique(app_data$raw_data[[input$zone_col]])),
      paste(names(app_data$raw_data), collapse = ", ")
    )
  })

  output$data_validation <- renderText({
    if (is.null(app_data$raw_data)) return("Load data first")
    comp_cols <- c(input$comp_sand_col, input$comp_silt_col, input$comp_clay_col)
    comp_issues <- validate_composition(app_data$raw_data, comp_cols)
    spatial_issues <- validate_spatial(
      app_data$raw_data, input$x_col, input$y_col, input$z_col, input$zone_col
    )
    all_issues <- c(comp_issues, spatial_issues)
    if (length(all_issues) == 0) {
      "[OK] Data validation passed - Ready for modeling"
    } else {
      paste("[Issues]\n", paste("- ", all_issues, collapse = "\n"))
    }
  })

  observeEvent(input$main_menu == "model_setup", {
    if (is.null(app_data$raw_data)) return()
    comp_cols <- c(input$comp_sand_col, input$comp_silt_col, input$comp_clay_col)
    tryCatch(
      {
        app_data$validated_data <- prepare_hierarchical_data(
          app_data$raw_data,
          zone_col = input$zone_col,
          comp_cols = comp_cols
        )
      },
      error = function(e) {
        showNotification(paste("Error preparing data:", e$message), type = "error")
      }
    )
  })

  observeEvent(input$run_simulation, {
    validate(
      need(!is.null(app_data$validated_data), "Prepare data first"),
      need(!is.null(input$backend), "Select backend")
    )
    app_data$sim_status <- "Initializing..."
    tryCatch(
      {
        app_data$sim_status <- "Fitting hierarchical model..."
        fit_result <- fit_zone_model(
          data = app_data$validated_data,
          backend = input$backend,
          params = list(pooling_coef = input$pooling_coef)
        )
        
        app_data$hierarchical_fit <- fit_result
        app_data$sim_status <- paste0("Generating ", input$nsim, " realizations...")
        zone_sims <- simulate_realizations(fit_result, n_sims = input$nsim)
        comp_sims <- lapply(zone_sims, function(z_ilr) {
          ilr_to_composition(z_ilr, n_components = 3)
        })
        app_data$simulations <- comp_sims
        app_data$sim_status <- "Simulation complete!"
        showNotification("Simulation finished successfully!", type = "message")
      },
      error = function(e) {
        app_data$sim_status <- paste("Error:", e$message)
        showNotification(paste("Simulation error:", e$message), type = "error")
      }
    )
  })

  output$sim_status <- renderText({
    paste("Status:", app_data$sim_status)
  })

  output$sim_log <- renderText({
    if (is.null(app_data$hierarchical_fit)) {
      "Run simulation to see logs"
    } else {
      paste(
        "Model:", app_data$hierarchical_fit$fit$metadata$backend,
        "\nZones:", length(app_data$hierarchical_fit$zones),
        "\nObservations:", app_data$hierarchical_fit$fit$metadata$n_obs_total,
        "\nPooling coefficient:", round(input$pooling_coef, 2)
      )
    }
  })

  observeEvent(input$compute_prob, {
    if (is.null(app_data$simulations)) {
      showNotification("Run simulation first!", type = "error")
      return()
    }
    prob_results <- probability_map(
      sims_by_zone = app_data$simulations,
      variable = input$prob_variable,
      threshold = input$prob_threshold,
      operator = input$prob_operator
    )
    prob_df <- data.frame(
      zone = names(prob_results),
      probability = unlist(prob_results)
    )

    p <- ggplot(prob_df, aes(x = reorder(zone, probability), y = probability)) +
      geom_col(fill = "steelblue", alpha = 0.7) +
      coord_flip() +
      theme_minimal() +
      labs(
        title = paste("Probability of", input$prob_variable, ">", input$prob_threshold, "%"),
        x = "Zone",
        y = "Probability"
      )

    output$prob_map <- renderPlotly({
      ggplotly(p)
    })

    output$prob_stats <- renderTable({
      prob_df
    })
  })

  observeEvent(input$compute_pct, {
    if (is.null(app_data$simulations)) {
      showNotification("Run simulation first!", type = "error")
      return()
    }
    pct_results <- percentile_map(
      sims_by_zone = app_data$simulations,
      variable = input$pct_variable,
      percentiles = as.numeric(input$pct_levels)
    )
    p <- plot_ly() %>%
      add_trace(
        type = "bar",
        x = names(pct_results),
        y = sapply(pct_results, `[`, 1),
        name = "P10",
        marker = list(color = "lightblue")
      ) %>%
      add_trace(
        x = names(pct_results),
        y = sapply(pct_results, function(x) {
          if (length(x) > 1) x[2] else NA
        }),
        name = "P50",
        marker = list(color = "steelblue")
      ) %>%
      add_trace(
        x = names(pct_results),
        y = sapply(pct_results, function(x) {
          if (length(x) > 2) x[3] else NA
        }),
        name = "P90",
        marker = list(color = "darkblue")
      ) %>%
      layout(
        title = paste("Percentile Map:", input$pct_variable),
        xaxis = list(title = "Zone"),
        yaxis = list(title = paste(input$pct_variable, "%")),
        barmode = "group"
      )

    output$pct_map <- renderPlotly({
      p
    })
  })

  observeEvent(input$compute_profile, {
    if (is.null(app_data$simulations) || is.null(input$depth_zone)) {
      showNotification("Run simulation first and select zone!", type = "error")
      return()
    }
    zone_sims <- app_data$simulations[[input$depth_zone]]
    depth_mean <- colMeans(zone_sims)
    depth_sd <- apply(zone_sims, 2, sd)
    depth_q10 <- apply(zone_sims, 2, quantile, 0.1)
    depth_q90 <- apply(zone_sims, 2, quantile, 0.9)

    profile_df <- data.frame(
      depth_index = seq_along(depth_mean),
      mean = depth_mean[[input$depth_variable]],
      sd = depth_sd[[input$depth_variable]],
      q10 = depth_q10[[input$depth_variable]],
      q90 = depth_q90[[input$depth_variable]]
    )

    p <- ggplot(profile_df, aes(x = mean, y = depth_index)) +
      geom_ribbon(aes(xmin = q10, xmax = q90), alpha = 0.3) +
      geom_line(color = "steelblue", size = 1) +
      geom_point(size = 3, color = "steelblue") +
      theme_minimal() +
      labs(
        title = paste("Depth Profile:", input$depth_zone),
        x = paste(input$depth_variable, "%"),
        y = "Depth Level"
      )

    output$depth_profile <- renderPlotly({
      ggplotly(p)
    })
  })

  output$scatter_3d <- renderPlotly({
    if (is.null(app_data$simulations)) {
      plot_ly() %>%
        add_text(
          text = "Run simulation first",
          textposition = "middle center"
        ) %>%
        layout(title = "3D Interactive Visualization")
    } else {
      all_sims <- do.call(rbind, app_data$simulations)

      p <- plot_ly(all_sims,
        x = ~get(input$var_x),
        y = ~get(input$var_y),
        size = ~get(input$var_z),
        color = ~get(input$var_z),
        type = "scatter",
        mode = "markers"
      ) %>%
        layout(
          title = "Compositional Space",
          xaxis = list(title = input$var_x),
          yaxis = list(title = input$var_y),
          showlegend = FALSE
        )

      p
    }
  })

  output$risk_curve <- renderPlotly({
    if (is.null(app_data$simulations)) {
      plot_ly() %>%
        add_text(text = "Run simulation first") %>%
        layout(title = "Risk Analysis")
    } else {
      thresholds <- seq(input$risk_type[1], input$risk_type[2], by = 5)
      prob_exceed <- sapply(thresholds, function(t) {
        all_sims <- do.call(rbind, app_data$simulations)
        mean(all_sims[[input$risk_variable]] > t)
      })

      p <- plot_ly() %>%
        add_trace(
          x = thresholds,
          y = prob_exceed,
          type = "scatter",
          mode = "lines+markers",
          fill = "tozeroy"
        ) %>%
        layout(
          title = paste("Risk Curve:", input$risk_variable),
          xaxis = list(title = "Threshold (%)"),
          yaxis = list(title = "Probability of Exceedance")
        )

      p
    }
  })

  output$download_prob_tif <- downloadHandler(
    filename = "probability_map.tif",
    content = function(file) {
      showNotification("Export not yet implemented", type = "message")
    }
  )

  output$download_config <- downloadHandler(
    filename = "model_config.json",
    content = function(file) {
      config <- list(
        backend = input$backend,
        pooling_coef = input$pooling_coef,
        nsim = input$nsim,
        timestamp = Sys.time()
      )
      writeLines(jsonlite::toJSON(config, pretty = TRUE), file)
    }
  )

  output$download_report <- downloadHandler(
    filename = "uncertainty_report.html",
    content = function(file) {
      showNotification("Report generation not yet implemented", type = "message")
    }
  )
}
