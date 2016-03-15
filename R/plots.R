name <- "Plots"

ui <- function(tab)
    tab(name,
        # allows the user to choose which UI set is shown
        fluidRow(
            column(2, selectizeInput("selectizePlot", "Select plot type:",
                                     choices = NULL,
                                     options = list(
                                         placeholder = "Select a plot type"))),
            column(2, selectizeInput("selectizeEvent", "Select event:",
                                     choices = NULL,
                                     options = list(
                                         placeholder = "Select an event")))),
        uiOutput("plots"))

server <- function(input, output, session) {
    # Loads valid scripts from the indicated folder
    plotEnvs <- sourceScripts(folder = paste0(tabsFolder, "plots/"),
                              check = c("name", "ui"))
    plotEnvs.server <- lapply(plotEnvs, "[[", "server")
    # Get name of the loaded scripts
    names <- sapply(plotEnvs, "[[", "name")
    
    # Runs server logic from the scripts
    lapply(plotEnvs.server, do.call, list(input, output, session))

    # Updates selectize input to show available plots
    updateSelectizeInput(session, "selectizePlot", choices = names)
    
    # Updates selectize event to show available events
    observe({
        psi <- getInclusionLevels()
        if (!is.null(psi)) {
            choices <- rownames(psi)
            names(choices) <- gsub("_", " ", rownames(psi))
            updateSelectizeInput(session, "selectizeEvent", choices = choices)
        }
    })
    
    # If showing datatable, hide selectizeEvent; otherwise, show it
    observe({
        vis <- function(func, ...)
            func("selectizeEvent", anim = TRUE, animType = "fade")
        
        if(grepl("Inclusion levels", input$selectizePlot)) vis(shinyjs::hide)
        else vis(shinyjs::show)
    })
    
    # Render the respective UI of the requested plot
    output$plots <- renderUI({
        lapply(plotEnvs, function(env) {
            conditionalPanel(
                condition = sprintf("input.selectizePlot=='%s'", env$name),
                env$ui)
        })
    })
}