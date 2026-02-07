#' Launch Uncertainty Mapper Shiny Application
#'
#' Interactive web application for compositional geostatistical visualization.
#'
#' @param ... Arguments passed to shiny::shinyApp()
#' @return Invisible NULL. Launches Shiny application.
#' @examples
#' \dontrun{
#' run_uncertainty_mapper()
#' }
#' @export
run_uncertainty_mapper <- function(...) {
  source("ui.R", local = TRUE)
  source("server.R", local = TRUE)
  source("helpers.R", local = TRUE)
  shiny::shinyApp(ui = ui, server = server, ...)
}
