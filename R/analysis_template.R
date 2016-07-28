#' User interface of template
#' @param id Character: namespace identifier
#' 
#' @importFrom shiny NS
#' @importFrom highcharter highchartOutput
#' 
#' @return HTML elements
templateUI <- function(id) {
    ns <- NS(id)
    highchartOutput(ns("plot"))
}

#' Server logic of template
#' 
#' @param input Shiny input
#' @param output Shiny output
#' @param session Shiny session
#' 
#' @importFrom highcharter renderHighchart
#' @importFrom stats runif
templateServer <- function(input, output, session) {
    num <- runif(20, 1, 100)
    output$plot <- renderHighchart( hchart(num) )
}

# attr(templateUI, "loader") <- "analysis"
# attr(templateUI, "name") <- "Template"
# attr(templateUI, "selectEvent") <- TRUE
# attr(templateServer, "loader") <- "analysis"