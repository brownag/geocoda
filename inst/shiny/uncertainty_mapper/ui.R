# Geocoda: Uncertainty Mapper Shiny UI

suppressPackageStartupMessages({
  library(shiny)
  library(shinydashboard)
  library(plotly)
})

header <- dashboardHeader(
  title = "Geocoda: Uncertainty Mapper",
  titleWidth = 300,
  tags$li(
    class = "dropdown",
    tags$a(
      href = "https://github.com/brownag/geocoda",
      target = "_blank",
      icon("github"),
      "View Source"
    )
  )
)

# Dashboard Sidebar
sidebar <- dashboardSidebar(
  width = 280,
  sidebarMenu(
    id = "main_menu",
    menuItem(
      "Home",
      tabName = "home",
      icon = icon("home")
    ),
    menuItem(
      "1. Data Upload",
      tabName = "data_upload",
      icon = icon("upload")
    ),
    menuItem(
      "2. Model Setup",
      tabName = "model_setup",
      icon = icon("cogs")
    ),
    menuItem(
      "3. Run Simulation",
      tabName = "simulation",
      icon = icon("play")
    ),
    menuItem(
      "4. Visualization",
      tabName = "visualization",
      icon = icon("chart-bar"),
      menuSubItem("Probability Maps", tabName = "prob_maps"),
      menuSubItem("Percentile Maps", tabName = "pct_maps"),
      menuSubItem("Depth Profiles", tabName = "depth_profiles"),
      menuSubItem("3D Interactive", tabName = "interactive_3d"),
      menuSubItem("Risk Curves", tabName = "risk_curves")
    ),
    menuItem(
      "5. Download Results",
      tabName = "download",
      icon = icon("download")
    ),
    hr(),
    menuItem(
      "About",
      tabName = "about",
      icon = icon("info-circle")
    )
  )
)

# Dashboard Body
body <- dashboardBody(
  # Custom CSS
  tags$head(
    tags$style(HTML("
      .main-header .logo {
        font-weight: bold;
        font-size: 16px;
      }
      .content-wrapper {
        background-color: #ecf0f5;
      }
      .box {
        box-shadow: 0 1px 1px rgba(0,0,0,0.1);
      }
      .progress-container {
        margin: 20px 0;
      }
    "))
  ),

  tabItems(
    # HOME TAB
    tabItem(
      tabName = "home",
      h1("Welcome to Geocoda: Uncertainty Mapper"),
      p(
        "Interactive tool for exploring 3D geostatistical uncertainties
        and hierarchical zone model predictions."
      ),
      fluidRow(
        box(
          title = "Quick Start",
          width = 12,
          p("Follow these 5 simple steps:"),
          ol(
            li(strong("Data Upload"), " - Load your compositional data (CSV)"),
            li(strong("Model Setup"), " - Configure hierarchical zones and parameters"),
            li(strong("Run Simulation"), " - Execute geostatistical simulation"),
            li(strong("Visualize"), " - Explore uncertainty in interactive maps"),
            li(strong("Download"), " - Export results as raster/shapefile")
          )
        )
      ),
      fluidRow(
        box(
          title = "Features",
          width = 6,
          tags$ul(
            tags$li("Support for compositional data (sand, silt, clay %)"),
            tags$li("3D spatial modeling with stratigraphic coordinates"),
            tags$li("Hierarchical zone models with multiple backends"),
            tags$li("100+ realizations for uncertainty quantification"),
            tags$li("Probability maps at custom thresholds")
          )
        ),
        box(
          title = "Data Requirements",
          width = 6,
          tags$ul(
            tags$li("CSV file with columns:"),
            tags$ul(
              tags$li("x, y, z (coordinates)"),
              tags$li("zone_id or zone_name"),
              tags$li("sand%, silt%, clay% (or other parts)")
            ),
            tags$li("Compositions must sum to ~100%"),
            tags$li("Minimum 5 observations per zone")
          )
        )
      )
    ),

    # DATA UPLOAD TAB
    tabItem(
      tabName = "data_upload",
      h2("Step 1: Upload Data"),
      fluidRow(
        box(
          title = "Load Data File",
          width = 6,
          fileInput(
            "data_file",
            "Choose CSV file:",
            accept = c("text/csv", "text/comma-separated-values", ".csv")
          ),
          p(tags$em("Or use example data:")),
          actionButton(
            "load_example",
            "Load Example Dataset",
            icon = icon("database"),
            class = "btn-primary"
          )
        ),
        box(
          title = "Data Preview",
          width = 6,
          tableOutput("data_preview")
        )
      ),
      fluidRow(
        box(
          title = "Data Summary",
          width = 12,
          h4("Dataset Overview"),
          verbatimTextOutput("data_summary"),
          h4("Data Quality Check"),
          verbatimTextOutput("data_validation")
        )
      )
    ),

    # MODEL SETUP TAB
    tabItem(
      tabName = "model_setup",
      h2("Step 2: Configure Model"),
      fluidRow(
        box(
          title = "Composition Definition",
          width = 6,
          selectInput(
            "comp_sand_col",
            "Sand column:",
            choices = NULL
          ),
          selectInput(
            "comp_silt_col",
            "Silt column:",
            choices = NULL
          ),
          selectInput(
            "comp_clay_col",
            "Clay column:",
            choices = NULL
          )
        ),
        box(
          title = "Spatial Setup",
          width = 6,
          selectInput(
            "zone_col",
            "Zone identifier column:",
            choices = NULL
          ),
          selectInput(
            "x_col",
            "Easting (X) column:",
            choices = NULL
          ),
          selectInput(
            "y_col",
            "Northing (Y) column:",
            choices = NULL
          ),
          selectInput(
            "z_col",
            "Depth/Elevation (Z) column:",
            choices = NULL
          )
        )
      ),
      fluidRow(
        box(
          title = "Model Parameters",
          width = 6,
          selectInput(
            "backend",
            "Fitting backend:",
            choices = c(
              "analytical (fast)" = "analytical",
              "analytical_3d (spatial)" = "analytical_3d",
              "stan_3d (Bayesian)" = "stan_3d",
              "nimble_3d (MCMC)" = "nimble_3d"
            ),
            selected = "analytical_3d"
          ),
          sliderInput(
            "pooling_coef",
            "Pooling strength (0=independent, 1=full pooling):",
            min = 0, max = 1, value = 0.4, step = 0.1
          ),
          sliderInput(
            "spatial_decay_power",
            "Spatial decay power (inverse distance):",
            min = 1, max = 3, value = 2, step = 0.5
          )
        ),
        box(
          title = "MCMC Settings (for Stan/Nimble)",
          width = 6,
          sliderInput(
            "mcmc_iter",
            "Iterations per chain:",
            min = 500, max = 5000, value = 2000, step = 500
          ),
          sliderInput(
            "mcmc_warmup",
            "Warmup iterations:",
            min = 100, max = 2000, value = 500, step = 100
          ),
          sliderInput(
            "mcmc_chains",
            "Number of chains:",
            min = 1, max = 8, value = 2, step = 1
          )
        )
      )
    ),

    # SIMULATION TAB
    tabItem(
      tabName = "simulation",
      h2("Step 3: Run Simulation"),
      fluidRow(
        box(
          title = "Simulation Parameters",
          width = 6,
          sliderInput(
            "nsim",
            "Number of realizations:",
            min = 10, max = 1000, value = 100, step = 50
          ),
          selectInput(
            "sim_grid_type",
            "Prediction grid type:",
            choices = c(
              "Regular grid" = "regular",
              "From observations" = "obs_only",
              "Custom coordinates" = "custom"
            )
          ),
          sliderInput(
            "grid_spacing",
            "Grid spacing (m):",
            min = 10, max = 200, value = 50, step = 10
          )
        ),
        box(
          title = "Computation",
          width = 6,
          actionButton(
            "run_simulation",
            "Run Simulation",
            icon = icon("play"),
            class = "btn-success btn-lg",
            width = "100%"
          ),
          br(), br(),
          checkboxInput(
            "show_progress",
            "Show detailed progress",
            value = TRUE
          ),
          checkboxInput(
            "save_intermediate",
            "Save intermediate results",
            value = FALSE
          )
        )
      ),
      fluidRow(
        box(
          title = "Simulation Status",
          width = 12,
          p(textOutput("sim_status")),
          progressBar(
            id = "sim_progress",
            value = 0,
            striped = TRUE,
            animated = TRUE
          ),
          verbatimTextOutput("sim_log")
        )
      )
    ),

    # PROBABILITY MAPS TAB
    tabItem(
      tabName = "prob_maps",
      h2("Probability Maps"),
      fluidRow(
        box(
          title = "Map Settings",
          width = 3,
          selectInput(
            "prob_variable",
            "Variable:",
            choices = c("sand", "silt", "clay")
          ),
          numericInput(
            "prob_threshold",
            "Threshold value (%):",
            value = 40, min = 0, max = 100, step = 5
          ),
          radioButtons(
            "prob_operator",
            "Show probability of:",
            choices = c(
              "Greater than (>)" = "gt",
              "Less than (<)" = "lt",
              "Within range" = "range"
            )
          ),
          actionButton(
            "compute_prob",
            "Update Map",
            class = "btn-primary",
            icon = icon("refresh")
          )
        ),
        box(
          title = "Probability Map",
          width = 9,
          plotlyOutput("prob_map", height = "600px")
        )
      ),
      fluidRow(
        box(
          title = "Statistics",
          width = 12,
          tableOutput("prob_stats")
        )
      )
    ),

    # PERCENTILE MAPS TAB
    tabItem(
      tabName = "pct_maps",
      h2("Percentile Maps"),
      fluidRow(
        box(
          title = "Map Settings",
          width = 3,
          selectInput(
            "pct_variable",
            "Variable:",
            choices = c("sand", "silt", "clay")
          ),
          checkboxGroupInput(
            "pct_levels",
            "Percentile levels:",
            choices = c("10" = 10, "50" = 50, "90" = 90),
            selected = c(10, 50, 90)
          ),
          actionButton(
            "compute_pct",
            "Update Map",
            class = "btn-primary",
            icon = icon("refresh")
          )
        ),
        box(
          title = "Percentile Map",
          width = 9,
          plotlyOutput("pct_map", height = "600px")
        )
      )
    ),

    # DEPTH PROFILES TAB
    tabItem(
      tabName = "depth_profiles",
      h2("Depth-Stratified Uncertainties"),
      fluidRow(
        box(
          title = "Profile Settings",
          width = 3,
          selectInput(
            "depth_variable",
            "Variable:",
            choices = c("sand", "silt", "clay")
          ),
          selectInput(
            "depth_zone",
            "Zone to profile:",
            choices = NULL
          ),
          actionButton(
            "compute_profile",
            "Update Profile",
            class = "btn-primary"
          )
        ),
        box(
          title = "Uncertainty by Depth",
          width = 9,
          plotlyOutput("depth_profile", height = "600px")
        )
      )
    ),

    # INTERACTIVE 3D TAB
    tabItem(
      tabName = "interactive_3d",
      h2("3D Interactive Visualization"),
      fluidRow(
        box(
          title = "3D Scatterplot Controls",
          width = 3,
          selectInput(
            "var_x",
            "X-axis:",
            choices = c("sand", "silt", "clay")
          ),
          selectInput(
            "var_y",
            "Y-axis:",
            choices = c("sand", "silt", "clay")
          ),
          selectInput(
            "var_z",
            "Z-axis (size):",
            choices = c("sand", "silt", "clay", "depth", "uncertainty")
          ),
          sliderInput(
            "sample_realizations",
            "Sample # of realizations to show:",
            min = 10, max = 500, value = 100, step = 50
          )
        ),
        box(
          title = "Composition Space",
          width = 9,
          plotlyOutput("scatter_3d", height = "600px")
        )
      )
    ),

    # RISK CURVES TAB
    tabItem(
      tabName = "risk_curves",
      h2("Risk Analysis Curves"),
      fluidRow(
        box(
          title = "Loss Function Settings",
          width = 3,
          selectInput(
            "risk_variable",
            "Variable:",
            choices = c("sand", "silt", "clay")
          ),
          selectInput(
            "risk_type",
            "Risk metric:",
            choices = c(
              "Probability of exceedance" = "exceed",
              "Expected loss" = "loss",
              "Regret" = "regret"
            )
          ),
          sliderInput(
            "threshold_range",
            "Threshold range (%):",
            min = 0, max = 100, value = c(20, 80), step = 5
          )
        ),
        box(
          title = "Risk Curve",
          width = 9,
          plotlyOutput("risk_curve", height = "600px")
        )
      )
    ),

    # DOWNLOAD TAB
    tabItem(
      tabName = "download",
      h2("Step 5: Download Results"),
      fluidRow(
        box(
          title = "Export Data",
          width = 6,
          h4("Raster Outputs (GeoTIFF)"),
          downloadButton("download_prob_tif", "Download Probability Map"),
          br(), br(),
          downloadButton("download_pct_tif", "Download Percentile Map"),
          br(), br(),
          h4("Vector Outputs (Shapefile)"),
          downloadButton("download_zones_shp", "Download Zone Results")
        ),
        box(
          title = "Metadata & Config",
          width = 6,
          h4("Model Configuration"),
          downloadButton("download_config", "Config (JSON)"),
          br(), br(),
          h4("Realizations Data"),
          downloadButton("download_sims_csv", "All Realizations (CSV)"),
          br(), br(),
          h4("Report"),
          downloadButton("download_report", "Full Report (HTML)")
        )
      )
    ),

    # ABOUT TAB
    tabItem(
      tabName = "about",
      h2("About Geocoda"),
      fluidRow(
        box(
          title = "Uncertainty Mapper",
          width = 12,
          p(
            "Interactive tool for 3D geostatistical uncertainty mapping
            with hierarchical zone models for compositional soil data."
          ),
          h4("Features:"),
          tags$ul(
            tags$li("3D spatial modeling with stratigraphic coordinates"),
            tags$li("Multiple fitting backends (analytical, Stan, Nimble)"),
            tags$li("Compositional data (ILR transform)"),
            tags$li("Probability and percentile maps"),
            tags$li("Risk-based decision support"),
            tags$li("Interactive visualizations")
          ),
          h4("Documentation:"),
          tags$ul(
            tags$li(
              tags$a(
                href = "https://brownag.github.io/geocoda/",
                target = "_blank",
                "Official Website"
              )
            ),
            tags$li(
              tags$a(
                href = "https://github.com/brownag/geocoda",
                target = "_blank",
                "Source Code (GitHub)"
              )
            )
          ),
          h4("Author:"),
          p("Andrew Brown, USDA-NRCS", br(), "Developed for soil survey and pedometric research")
        )
      )
    )
  )
)

# Create dashboard UI
ui <- dashboardPage(
  header = header,
  sidebar = sidebar,
  body = body,
  skin = "blue"
)

# Helper function for progress bar (not built-in to shinydashboard)
progressBar <- function(id, value = 0, striped = FALSE, animated = FALSE, ...) {
  shiny::tags$div(
    class = "progress progress-container",
    shiny::tags$div(
      class = paste(
        "progress-bar",
        if (striped) "progress-bar-striped",
        if (animated) "active"
      ),
      style = paste0("width: ", value, "%"),
      role = "progressbar"
    )
  )
}
