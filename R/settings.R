#' User interface
#' @importFrom parallel detectCores
settingsUI <- function(id, tab) {
    ns <- NS(id)
    
    if( require("parallel") ) {
        cores <- detectCores()
        coresInput <- tagList(
            sliderInput(ns("cores"), h4("Number of cores"), value=1, min=1, 
                        step=1, max=cores, width="auto"),
            helpText("You have a total of", cores, "cores to use")
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
                               value=3, min=0, max=20, step=1, width="auto"),
                   textOutput(ns("precisionExample")),
                   helpText("Only applies to new calculations.")),
            column(4,
                   sliderInput(ns("significant"), h4("Significant digits"),
                               value=3, min=0, max=20, step=1, width="auto"),
                   textOutput(ns("significantExample")),
                   helpText("Only applies to new calculations."))
        )
    )
}

#' Server logic
settingsServer <- function(input, output, session) {
    observe(setCores(input$cores))
    observe({
        setPrecision(input$precision)
        output$precisionExample <- renderText(
            paste(
                "Example:",
                formatC(283.5837243243332313139838387437323823823829283984784687284,
                        digits=getPrecision(), format="f"))
        )
    })
    observe({
        setSignificant(input$significant)
        output$significantExample <- renderText(
            paste(
                "Example:",
                formatC(0.000005849839043444982905434543482092830943294512474723566,
                        getSignificant(), format="g"))
        )
    })
}

attr(settingsUI, "loader") <- "app"
attr(settingsServer, "loader") <- "app"