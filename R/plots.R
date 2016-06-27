#' User interface
#' @importFrom shinyBS bsTooltip
plotsUI <- function(id, tab) {
    ns <- NS(id)
    uiList <- getUiFunctions(ns, "plots")
    sharedData$names <- sapply(uiList, attr, "name")
     
    tab("Plots",
        # allows the user to choose which UI set is shown
        fluidRow(
            column(4, selectizeInput(ns("selectizePlot"), "Select plot type:",
                                     choices = NULL, options = list(
                                         placeholder = "Select a plot type"),
                                     width="auto")),
            column(4, selectizeInput(ns("selectizeCategory"), "Select category:",
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
                                  ns("selectizePlot"), attr(ui, "name")), ui)
        })
    )
}

plotsServer <- function(input, output, session) {
    # Run server logic from the scripts
    getServerFunctions("plots")

    # Update selectize input to show available plots
    observe({
        updateSelectizeInput(session, "selectizePlot", choices=sharedData$names)
    })

    # Update selectize input to show available categories
    observe({
        data <- getData()
        if (!is.null(data))
            updateSelectizeInput(session, "selectizeCategory",
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
    
    # # If showing datatable, hide selectizeEvent; otherwise, show it
    # observe({
    #     vis <- function(func, ...)
    #         func(id("selectizeEvent"), anim = TRUE, animType = "fade")
    #     
    #     if(grepl("Inclusion levels", input[[id("selectizePlot")]]))
    #         vis(shinyjs::hide)
    #     else
    #         vis(shinyjs::show)
    # })
}

attr(plotsUI, "loader") <- "app"
attr(plotsServer, "loader") <- "app"