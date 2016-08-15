#' Missing information modal template
#'  
#' @param session Shiny session
#' @param dataType Character: type of data missing
#' @param buttonId Character: identifier of button to take user to load missing 
#' data
#' 
#' @examples
#' \dontrun{
#'  session <- session$ns
#'  buttonInput <- "takeMeThere"
#'  buttonId <- ns(buttonInput)
#'  dataType <- "Inclusion levels"
#'  missingDataModal(session, buttonId, dataType)
#'  observeEvent(input[[buttonInput]], missingDataGuide(dataType))
#' }
missingDataModal <- function(session, dataType, buttonId) {
    template <- function(buttonLabel) {
        errorModal(
            session, paste("Load", tolower(dataType)),
            "This analysis requires", tolower(dataType), "to proceed.",
            footer=actionButton(buttonId, buttonLabel, "data-dismiss"="modal",
                                class="btn-danger"))
    }
    
    switch(dataType,
           "Clinical data"=template("Load"),
           "Junction quantification"=template("Load"),
           "Inclusion levels"=template("Load or calculate"))
}

#' @rdname missingDataModal
loadRequiredData <- function(dataType) {
    panel <- switch(dataType,
                    "Clinical data"="TCGA",
                    "Junction quantification"="TCGA",
                    "Inclusion levels"="alternative"
    )
    
    return(sprintf("showDataPanel('%s');", panel))
}

#' @rdname missingDataModal
missingDataGuide <- function(dataType) {
    js <- loadRequiredData(dataType)
    runjs(js)
}

#' User interface for the data analyses
#' 
#' @param id Character: identifier
#' @param tab Function to process HTML elements
#' 
#' @importFrom shinyBS bsTooltip
#' @importFrom shiny NS div icon fluidRow column selectizeInput conditionalPanel
#' 
#' @return HTML element as character
analysesUI <- function(id, tab) { 
    ns <- NS(id)
    uiList <- getUiFunctions(ns, "analysis")
    
    # Load available analyses
    ui <- lapply(uiList, function(ui) tabPanel(attr(ui, "name"), ui) )
    do.call(navbarMenu, c(list(icon=icon("flask"), "Analyses"), ui))
}

#' Test the survival difference between two survival groups given a cutoff
#' 
#' @inheritParams processSurvTerms
#' @param cutoff Numeric: Cut-off of interest
#' @param data Numeric: elements of interest to test against the cut-off
#' @param group Pre-filled vector of missing values with the length of data
#' @param filter Boolean or numeric: interest of the data elements
#' @param ... Arguments to pass to \code{processSurvTerms}
#' @param session Shiny session
#' 
#' @importFrom survival survdiff
#' @return p-value of the survival difference
testSurvivalCutoff <- function(cutoff, data, filter, ..., group=NULL, 
                               session=NULL) {
    if (is.null(group)) groups <- rep(NA, nrow(clinical))
    group[filter] <- data >= cutoff
    
    # Assign a value based on the inclusion levels cut-off
    group[group == "TRUE"]  <- paste("Inclusion levels >=", cutoff)
    group[group == "FALSE"] <- paste("Inclusion levels <", cutoff)
    
    # Calculate survival curves
    if (!is.null(session)) {
        survTerms <- processSurvival(session, group, ...)
        if (is.null(survTerms)) return(NULL)
    } else {
        survTerms <- tryCatch(processSurvTerms(group, ...), error=return)
        if ("simpleError" %in% class(survTerms)) return(NA)
    }
    
    pvalue <- testSurvival(survTerms$form, data=survTerms$survTime)
    return(pvalue)
}

#' Server logic for the analyses
#' @param input Shiny input
#' @param output Shiny ouput
#' @param session Shiny session
#' 
#' @importFrom shiny observe observeEvent updateSelectizeInput
#' @importFrom shinyjs hide show
analysesServer <- function(input, output, session) {
    # Run server logic from the scripts
    server <- getServerFunctions("analysis")
}

attr(analysesUI, "loader") <- "app"
attr(analysesServer, "loader") <- "app"