name <- "Plots"

ui <- function()
    tabPanel(name,
             uiOutput("selectUI"),
             uiOutput("chosenUI"))

server <- function(input, output, session) {
    # loads valid scripts from the indicated folder
    envs <- loadScripts(folder = paste0(tabsFolder, "plots/"),
                        vars = c("name", "ui"))
    envs.server <- lapply(envs, "[[", "server")
    lapply(envs.server, do.call, list(input, output, session))
    
    output$selectUI <- renderUI({
        # get name of the loaded scripts
        names <- sapply(envs, "[[", "name")
        
        # allows the user to choose which UI set is shown
        selectizeInput("selectizePlot", "Select plot type:", choices = names,
                       options = list(placeholder = "Select a plot type"))
    })
    
    output$chosenUI <- renderUI({
        # if no option is avaliable, this section is not shown
        validate( need(input$selectizePlot, "No plots are available.") )
        # each UI set is loaded depending on the value of selectizePlot
        # WARNING: each script needs a unique name
        for (env in envs)
            if (input$selectizePlot == env$name) return(env$ui)
    })
}