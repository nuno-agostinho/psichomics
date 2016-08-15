#' User interface for the differential splicing analyses
#' 
#' @param id Character: identifier
#' @param tab Function to process HTML elements
#' 
#' @importFrom shinyBS bsTooltip
#' @importFrom shiny NS div icon fluidRow column selectizeInput conditionalPanel
#' 
#' @return HTML element as character
diffSplicingUI <- function(id, tab) { 
    ns <- NS(id)
    uiList <- getUiFunctions(ns, "diffSplicing", 
                             priority=c("diffSplicingTableUI",
                                        "diffSplicingEventUI"))
    
    ui <- lapply(uiList, function(ui) tabPanel(attr(ui, "name"), ui) )
    do.call(tabsetPanel, c(list(type="pills"), ui))
    
    # # Select analyses
    # do.call(navbarMenu, c(list(icon=icon("flask"), "Analyses"), ui))
}

#' Server logic for the differential splicing analyses
#' @param input Shiny input
#' @param output Shiny ouput
#' @param session Shiny session
#' 
#' @importFrom shiny observe observeEvent updateSelectizeInput
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