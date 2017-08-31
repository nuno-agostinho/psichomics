#' @rdname appUI
#' @importFrom shiny NS
diffSplicingUI <- function(id, tab) { 
    ns <- NS(id)
    uiList <- getUiFunctions(ns, "diffSplicing", 
                             priority=c("diffSplicingTableUI",
                                        "diffSplicingEventUI"))
    return(uiList)
}

#' @rdname appServer
#' 
#' @importFrom shiny observe observeEvent renderPlot
#' @importFrom shinyjs hide show
diffSplicingServer <- function(input, output, session) {
    # Run server logic from the scripts
    server <- getServerFunctions("diffSplicing",
                                 priority=c("diffSplicingTableServer",
                                            "diffSplicingEventServer"))
}

attr(diffSplicingUI, "loader") <- "analysis"
attr(diffSplicingUI, "name") <- "Differential splicing analysis"
attr(diffSplicingServer, "loader") <- "analysis"