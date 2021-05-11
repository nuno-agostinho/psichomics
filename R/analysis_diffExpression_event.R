#' @rdname appUI
#'
#' @importFrom highcharter highchartOutput
#' @importFrom shiny tagList uiOutput NS sidebarLayout numericInput h3 mainPanel
#' actionButton sidebarPanel
diffExpressionEventUI <- function(id) {
    ns <- NS(id)
    return(diffEventUI(id, ns, psi=FALSE))
}

#' @rdname appServer
#'
#' @importFrom highcharter renderHighchart
#' @importFrom shinyjs show hide
diffExpressionEventServer <- function(input, output, session) {
    ns <- session$ns
    diffEventServer(ns, input, output, session, psi=FALSE)

    observe({
        psi   <- getInclusionLevels()
        event <- getEvent()
        toggleRugBasedOnDataLength(session, psi, event)
    })

    # Update available gene choices depending on gene expression data loaded
    # Reactive avoids updating if the input remains the same
    updateGeneChoices <- reactive({
        geneExpr <- getGeneExpression(input$geneExpr)
        genes <- rownames(geneExpr)
        updateSelectizeInput(session, "gene", choices=genes, server=TRUE)
    })

    observe({
        geneExpr <- getGeneExpression()
        if (is.null(geneExpr)) {
            show("missingData")
        } else {
            updateSelectizeInput(session, "geneExpr",
                                 choices=rev(names(geneExpr)))
            hide("missingData")
        }
    })

    # Show options if gene expression data is available, update available gene
    # expression data choices and update available genes for selection
    observe({
        geneExpr <- getGeneExpression(input$geneExpr)
        if (is.null(geneExpr)) {
            hide("singleEventOptions")
            hide("survivalButton")
            hide("singleEventInfo")
            # Gene-related tasks
            hide("gene")
        } else {
            show("singleEventOptions")
            # Gene-related tasks
            hide("gene")
            updateGeneChoices()
            show("gene")
        }
    })
}

attr(diffExpressionEventUI, "loader") <- "diffExpression"
attr(diffExpressionEventUI, "name") <- "Individual gene"
attr(diffExpressionEventServer, "loader") <- "diffExpression"
