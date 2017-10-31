#' @rdname appUI
#' @importFrom shiny NS
diffExpressionUI <- function(id, tab) { 
    ns <- NS(id)
    uiList <- getUiFunctions(ns, "diffExpression", 
                             priority=c("diffExpressionTableUI",
                                        "diffExpressionEventUI"))
    return(uiList)
}

#' @rdname appServer
#' 
#' @importFrom shiny observe observeEvent renderPlot
#' @importFrom shinyjs hide show
diffExpressionServer <- function(input, output, session) {
    # Run server logic from the scripts
    server <- getServerFunctions("diffExpression",
                                 priority=c("diffExpressionTableServer",
                                            "diffExpressionEventServer"))
}

attr(diffExpressionUI, "loader") <- "analysis"
attr(diffExpressionUI, "name") <- "Differential expression analysis"
attr(diffExpressionServer, "loader") <- "analysis"