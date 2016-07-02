#' User interface
#' @importFrom shinyBS bsTooltip
analysesUI <- function(id, tab) {
    ns <- NS(id)
    uiList <- getUiFunctions(ns, "analysis")
    sharedData$names <- sapply(uiList, attr, "name")
    
    tab(div(icon("flask"), "Analyses"),
        # allows the user to choose which UI set is shown
        fluidRow(
            column(4, 
                   selectizeInput(ns("selectizeAnalysis"), 
                                  "Select analysis type:",
                                  choices = NULL, options = list(
                                      placeholder = "Select an analysis type"),
                                  width="auto")),
            column(4, selectizeInput(ns("selectizeCategory"), 
                                     "Select category:",
                                     choices = NULL, options = list(
                                         placeholder = "Select a category"),
                                     width="auto")),
            column(4, selectizeInput(ns("selectizeEvent"), "Select event:",
                                     choices = NULL, options = list(
                                         placeholder = "Select an event"),
                                     width="auto"))),
        # bsTooltip(ns("selectizeEvent"), placement="right",
        #           "Delete text and start typing to search events",
        #           options = list(container="body")),
        lapply(uiList, function(ui) {
            conditionalPanel(
                condition=sprintf("input[id='%s'] == '%s'",
                                  ns("selectizeAnalysis"), attr(ui, "name")),
                ui)
        })
    )
}

analysesServer <- function(input, output, session) {
    # Run server logic from the scripts
    getServerFunctions("analysis")
    
    # Update selectize input to show available analyses
    observe( updateSelectizeInput(session, "selectizeAnalysis",
                                  choices=sharedData$names) )
    
    # Update selectize input to show available categories
    observe({
        data <- getData()
        if (!is.null(data)) updateSelectizeInput(session, "selectizeCategory",
                                                 choices=names(data))
    })
    
    # Set the category of the data when possible
    observeEvent(input$selectizeCategory, setCategory(input$selectizeCategory))
    
    # Updates selectize event to show available events
    observe({
        psi <- getInclusionLevels()
        if (!is.null(psi)) {
            choices <- rownames(psi)
            names(choices) <- gsub("_", " ", rownames(psi))
            choices <- sort(choices)
            updateSelectizeInput(session, "selectizeEvent", choices=choices)
            
            # Set the selected alternative splicing event
            observeEvent(input$selectizeEvent, setEvent(input$selectizeEvent))
        }
        # } else {
        #     ## TODO(NunoA): Input doesn't seem to update when changing data...
        #     updateSelectizeInput(session, "selectizeEvent", choices=NULL,
        #                          options=list(
        #                              placeholder="No event available"))
        # }
    })
    
    # # If showing table, hide selectizeEvent; otherwise, show it
    # observe({
    #     vis <- function(func, ...)
    #         func(id("selectizeEvent"), anim = TRUE, animType = "fade")
    #     
    #     if(grepl("Inclusion levels", input[[id("selectizeAnalysis")]]))
    #         vis(shinyjs::hide)
    #     else
    #         vis(shinyjs::show)
    # })
}

attr(analysesUI, "loader") <- "app"
attr(analysesServer, "loader") <- "app"