#' User interface of template
#' @importFrom shiny NS
#' @importFrom highcharter highchartOutput
templateUI <- function(id) {
    ns <- NS(id)
    highchartOutput(ns("plot"))
}

#' Server logic of template
#' @importFrom highcharter renderHighchart
templateServer <- function(input, output, session) {
    num <- runif(20, 1, 100)
    output$plot <- renderHighchart( hchart(num) )
}

# attr(templateUI, "loader") <- "analysis"
# attr(templateUI, "name") <- "Template"
# attr(templateServer, "loader") <- "analysis"