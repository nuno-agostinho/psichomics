#' User interface
settingsUI <- function(id, tab) {
    ns <- NS(id)
    tab(title=div(icon("wrench"), "Settings"),
        numericInput(ns("cores"), h4("Number of cores"), value=1, min=1, 
                     step=1),
        
        numericInput(ns("precision"), h4("Numeric precision"), value=3, min=0,
                     step=1),
        textOutput(ns("precisionExample")),
        helpText("Only applies to new calculations."),
        
        numericInput(ns("significant"), h4("Significant digits"),
                     value=3, min=0, step=1),
        textOutput(ns("significantExample")),
        helpText("Only applies to new calculations.")
    )
}

#' Server logic
settingsServer <- function(input, output, session) {
    observe(setCores(input$cores))
    observe({
        setPrecision(input$precision)
        output$precisionExample <- renderText(
            formatC(283.5837243243332313139838387437323823823829283984784687284,
                    digits=getPrecision(), format="f"))
    })
    observe({
        setSignificant(input$significant)
        output$significantExample <- renderText(
            formatC(0.000005849839043444982905434543482092830943294512474723566,
                    getSignificant(), format="g"))
    })
}

attr(settingsUI, "loader") <- "app"
attr(settingsServer, "loader") <- "app"