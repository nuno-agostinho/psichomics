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
                                     choices = sort(rownames(mtcars)),
                                     options = list(
                                         placeholder = "Select an event")))),
        uiOutput("plots"))

server <- function(input, output, session) {
    # loads valid scripts from the indicated folder
    plotEnvs <- loadScripts(folder = paste0(tabsFolder, "plots/"),
                            vars = c("name", "ui"))
    plotEnvs.server <- lapply(plotEnvs, "[[", "server")
    # get name of the loaded scripts
    names <- sapply(plotEnvs, "[[", "name")
    
    # Runs server logic from the scripts
    lapply(plotEnvs.server, do.call, list(input, output, session))
    # Updates selectize input to show available plots
    updateSelectizeInput(session, "selectizePlot", choices = names)
    
    observe({
        # Hide selectizeEvent when showing datatable; in other cases, show it
        events <- function(f) f("selectizeEvent",
                                anim = TRUE, animType = "fade")
        if(input$selectizePlot == "datatable") events(hide)
        else events(show)
    })
    
    output$plots <- renderUI({
        lapply(plotEnvs, function(env) {
            conditionalPanel(
                condition = sprintf("input.selectizePlot=='%s'", env$name),
                env$ui)
        })
    })
}