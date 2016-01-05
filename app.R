library(shiny)
library(shinyBS)

tabsFolder <- "R/"

loadScripts <- function(folder, vars, exclude = "", ...){
    envs <- list()
    
    # Exclude unwanted files
    exclude.regex <- paste(exclude, collapse = "|")
    if (exclude != "")
        files <- grep(exclude.regex, files, invert = TRUE, value = TRUE)
    
    # Check if the script defines all the desired variables
    files <- list.files(folder, ...)
    envs <- lapply(files, function(file) {
        env <- new.env()
        sys.source(paste0(folder, file), env)
        if (all(sapply(vars, exists, envir = env))) {
            return(env)
        }
    })
    envs <- Filter(Negate(is.null), envs)
    return(envs)
}

#' Server function
#' 
#' This function has the instructions to build the Shiny app
#'
#' @param input 
#' @param output 
#' @param session 
server <- function(input, output, session) {
    tabs2 <- loadScripts(tabsFolder, c("name", "server"))
    tabs.server <- lapply(tabs2, "[[", "server")
    lapply(tabs.server, do.call, list(input, output, session))
    
    # Stop Shiny app when session ends (e.g. closing the window)
    # TODO(NunoA): maybe it'd be better if the app DIDN'T stop... either that
    # or make it easy to resume app
    session$onSessionEnded(function() {
        # Stop app and print message to console
        suppressMessages(stopped <- stopApp(returnValue=TRUE))
        if (stopped) cat("\nShiny app was exited ")
    })
}

#' The user-interface (ui) controls the layout and appearance of the app

# Loads the interface from each tab
tabs <- loadScripts(tabsFolder, c("name", "ui"))
tabs.ui <- lapply(tabs, "[[", "ui")

ui <- shinyUI(
    do.call(navbarPage, append(
        list(
            title = "spliced canceR", id = "nav",
            collapsible = TRUE, position = "fixed-top",
            header = list(
                # Avoids fixed-top navbar from obscuring content
                tags$style(type = "text/css", "body {padding-top: 70px;}"),
                # Alert appears fixed on the right at the top of other elements
                tags$style(type = "text/css", ".sbs-alert {position: fixed; right: 10px; z-index:9;}")
            )
        ),
        # Loads the interface from each tab
        tabs.ui
    ))
)

shinyApp(ui, server)