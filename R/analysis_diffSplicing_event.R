#' @rdname appUI
#'
#' @importFrom highcharter highchartOutput
#' @importFrom shiny tagList uiOutput NS sidebarLayout numericInput h3 mainPanel
#' actionButton sidebarPanel
diffSplicingEventUI <- function(id) {
    ns <- NS(id)
    return(diffEventUI(id, ns, psi=TRUE))
}

#' @rdname appServer
#'
#' @importFrom highcharter renderHighchart
#' @importFrom shinyjs show hide
diffSplicingEventServer <- function(input, output, session) {
    ns <- session$ns
    diffEventServer(ns, input, output, session, psi=TRUE)

    observe({
        geneExpr <- getGeneExpression(input$geneExpr)
        gene     <- input$gene
        toggleRugBasedOnDataLength(session, geneExpr, gene)
    })

    # Toggle options only if required data is available
    observe({
        # Get splicing event's inclusion levels
        psi <- getInclusionLevels()
        if (is.null(psi)) {
            show("missingData")
            hide("singleEventOptions")
            hide("survivalButton")
            hide("singleEventInfo")
        } else {
            hide("missingData")
            show("singleEventOptions")
        }
    })
}

attr(diffSplicingEventUI, "loader") <- "diffSplicing"
attr(diffSplicingEventUI, "name") <- "Individual alternative splicing event"
attr(diffSplicingEventServer, "loader") <- "diffSplicing"
