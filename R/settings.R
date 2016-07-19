#' User interface of the settings
#' 
#' @param id Character: identifier
#' @param tab Function to create tabs
#' 
#' @importFrom parallel detectCores
#' @importFrom shiny NS tagList sliderInput h4 helpText numericInput div icon
#' fluidRow column textOutput
#' 
#' @return HTML elements
settingsUI <- function(id, tab) {
    ns <- NS(id)
    
    if (requireNamespace("parallel", quietly = TRUE)) {
        cores <- parallel::detectCores()
        coresInput <- tagList(
            sliderInput(ns("cores"), h4("Number of cores"), value=1, min=1, 
                        step=1, max=cores, width="auto", post=" core(s)"),
            helpText("A total of", cores, "were detected.")
        )
    } else {
        coresInput <- numericInput(ns("cores"), h4("Number of cores"), value=1,
                                   min=1, step=1, width="auto")
    }
    
    tab(title=div(icon("wrench"), "Settings"),
        fluidRow(
            column(4, coresInput),
            column(4,
                   sliderInput(ns("precision"), h4("Numeric precision"),
                               value=3, min=0, max=10, step=1, width="auto",
                               post=" decimal(s)"),
                   textOutput(ns("precisionExample")),
                   helpText("Only applies to new calculations.")),
            column(4,
                   sliderInput(ns("significant"), h4("Significant digits"),
                               value=3, min=0, max=10, step=1, width="auto",
                               post=" digit(s)"),
                   textOutput(ns("significantExample")),
                   helpText("Only applies to new calculations."))
        )
    )
}

#' Server logic of the settings
#' 
#' @param session Shiny session
#' @param input Shiny input
#' @param output Shiny output
#' 
#' @importFrom shiny observe renderText
settingsServer <- function(input, output, session) {
    observe(setCores(input$cores))
    observe({
        setPrecision(input$precision)
        output$precisionExample <- renderText(
            paste("Example:",
                  formatC(283.5837243243332313139838387437323823823829,
                          digits=getPrecision(), format="f")))
    })
    
    observe({
        setSignificant(input$significant)
        output$significantExample <- renderText(
            paste("Example:", 
                  formatC(0.000005849839043444982905434543482092830943,
                          getSignificant(), format="g")))
    })
}

attr(settingsUI, "loader") <- "app"
attr(settingsServer, "loader") <- "app"