name <- "Plots"
id <- function(value) objectId(name, value)

# Loads valid scripts from the indicated folder
plotName <- "plot"
plotEnvs <- sourceScripts(folder = paste0(tabsFolder, "plots/"),
                          check = c(plotName, "ui"),
                          parentEnv = environment())
plotEnvs.server <- lapply(plotEnvs, "[[", "server")

# Get name of the loaded scripts
names <- sapply(plotEnvs, "[[", plotName)

ui <- function(tab)
    tab(name,
        # allows the user to choose which UI set is shown
        fluidRow(
            column(4, selectizeInput(id("selectizePlot"), "Select plot type:",
                                     choices = NULL,
                                     options = list(
                                         placeholder = "Select a plot type"))),
            column(4, selectizeInput(id("selectizeCategory"), "Select category:",
                                     choices = NULL,
                                     options = list(
                                         placeholder = "Select a category"))),
            column(4, selectizeInput(id("selectizeEvent"), "Select event:",
                                     choices = NULL,
                                     options = list(
                                         placeholder = "Select an event")))),
        bsTooltip(id("selectizeEvent"), placement="right",
                  "Delete text and start typing to search events",
                  options = list(container="body")),
        lapply(plotEnvs, function(env) {
            conditionalPanel(
                condition=sprintf("input[id='%s']=='%s'",
                                  id("selectizePlot"), env[[plotName]]), env$ui)
        }))

#' @importFrom shinyjs show hide
server <- function(input, output, session) {
    # Runs server logic from the scripts
    lapply(plotEnvs.server, do.call, list(input, output, session))
    
    # Updates selectize input to show available plots
    updateSelectizeInput(session, id("selectizePlot"), choices = names)
    
    # Updates selectize input to show available categories
    observe({
        data <- getData()
        if (!is.null(data))
            updateSelectizeInput(session, id("selectizeCategory"), 
                                 choices = names(data))
    })
    
    # Set the category of the data when possible
    observeEvent(input[[id("selectizeCategory")]], 
                 setCategory(input[[id("selectizeCategory")]]))
    
    # Updates selectize event to show available events
    observe({
        psi <- getInclusionLevels()
        if (!is.null(psi)) {
            choices <- rownames(psi)
            names(choices) <- gsub("_", " ", rownames(psi))
            updateSelectizeInput(session, id("selectizeEvent"),
                                 choices = choices)
        } else {
            ## TODO(NunoA): Input doesn't seem to update when changing data...
            updateSelectizeInput(session, id("selectizeEvent"),
                                 choices = NULL, options = list(
                                     placeholder = "No event available"))
        }
    })
    
    # Set the selected alternative splicing event
    observeEvent(input[[id("selectizeEvent")]],
                 setEvent(input[[id("selectizeEvent")]]))
    
    # If showing datatable, hide selectizeEvent; otherwise, show it
    observe({
        vis <- function(func, ...)
            func(id("selectizeEvent"), anim = TRUE, animType = "fade")
        
        if(grepl("Inclusion levels", input[[id("selectizePlot")]]))
            vis(hide)
        else
            vis(show)
    })
}