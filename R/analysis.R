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
missingDataGuide <- function(dataType) {
    panel <- switch(dataType,
                    "Clinical data"="TCGA",
                    "Junction quantification"="TCGA",
                    "Inclusion levels"="alternative"
    )
    
    js <- sprintf("showDataPanel('%s');", panel)
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
    
    analysesSelectEvent <- sapply(uiList, attr, "selectEvent")
    names(analysesSelectEvent) <- sapply(uiList, attr, "name")
    sharedData$analysesSelectEvent <- analysesSelectEvent
    
    tab(div(icon("flask"), "Analyses"),
        # allows the user to choose which UI set is shown
        fluidRow(
            column(4, 
                   selectizeInput(ns("selectizeAnalysis"), 
                                  "Select analysis:",
                                  choices = NULL, options = list(
                                      placeholder = "Select an analysis type"),
                                  width="auto")),
            column(4, selectizeInput(ns("selectizeCategory"), 
                                     "Select data category:",
                                     choices = NULL, options = list(
                                         placeholder = "Select data category"),
                                     width="auto")),
            column(4, selectizeInput(
                ns("selectizeEvent"), "Select event:",
                choices = NULL, options = list(
                    placeholder = paste("Search splicing events by gene,",
                                        "chromosome and coordinates")),
                width="auto"))),
        lapply(uiList, function(ui) {
            conditionalPanel(
                condition=sprintf("input[id='%s'] == '%s'",
                                  ns("selectizeAnalysis"), attr(ui, "name")),
                ui)
        }),
        conditionalPanel(sprintf("input[id='%s']==''", ns("selectizeAnalysis")), 
                         h3(icon("hand-o-up"), "Select an analysis above"))
    )
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
#' @param modals Boolean: show error dialogs? TRUE by default
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
    
    # Update selectize input to show available analyses
    observe( updateSelectizeInput(
        session, "selectizeAnalysis", selected="",
        choices=names(sharedData$analysesSelectEvent)))
    
    # Update selectize input to show available categories
    observe({
        data <- getData()
        if (!is.null(data))
            updateSelectizeInput(session, "selectizeCategory",
                                 choices=names(data))
    })
    
    # Set the category of the data
    observeEvent(input$selectizeCategory, setCategory(input$selectizeCategory))
    
    # Updates selectize event to show available events
    observe({
        psi <- getInclusionLevels()
        if (!is.null(psi)) {
            choices <- rownames(psi)
            names(choices) <- gsub("_", " ", rownames(psi))
            choices <- sort(choices)
            updateSelectizeInput(session, "selectizeEvent", choices=choices,
                                 selected=list())
            
            # Set the selected alternative splicing event
            observeEvent(input$selectizeEvent, setEvent(input$selectizeEvent))
        } else {
            # Replace with empty list since NULLs are dropped
            updateSelectizeInput(session, "selectizeEvent", choices=list(),
                                 selected=list())
        }
    })
    
    # If showing exploratory analyses, hide selectizeEvent; otherwise, show it
    observe({
        vis <- function(func, ...)
            func("selectizeEvent", anim = TRUE, animType = "fade")
        
        specific <- sharedData$analysesSelectEvent
        if(input$selectizeAnalysis %in% names(specific)[specific])
            vis(show)
        else
            vis(hide)
    })
}

attr(analysesUI, "loader") <- "app"
attr(analysesServer, "loader") <- "app"