#' launches the shinyApp
#'
#' @export launchApp
#'
#' @return shiny application object

# wrapper for shiny::shinyApp()
launchApp <- function() {
  shinyApp(ui = shinyAppUi, server = shinyAppServer)
}
