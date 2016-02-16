source("R/initial.R")

#' Server function
#' 
#' Instructions to build the Shiny app.
#'
#' @param input Input object
#' @param output Output object
#' @param session Session object
server <- function(input, output, session) {
    callScriptsFunction(tabsFolder, c("name", "server"), "server",
                        input, output, session)
    
    # Stop Shiny app when session ends (e.g. closing the window)
    session$onSessionEnded(function() {
        # Stop app and print message to console
        suppressMessages(stopped <- stopApp(returnValue=TRUE))
        if (stopped) cat("\nShiny app was closed ")
    })
}

header <- list(
    useShinyjs(),
    # Avoids fixed-top navbar from obscuring content
    tags$style(type = "text/css", "body {padding-top: 70px;}"),
    # Alert appears fixed on the top right above other elements
    tags$style(type = "text/css",
               ".sbs-alert{ position:fixed;
               right:10px;
               z-index:9;
               -webkit-filter: opacity(80%);
               filter: opacity(80%); }"))

# The user interface (ui) controls the layout and appearance of the app
ui <- shinyUI(
    do.call(navbarPage, c(
        list(title = "spliced canceR", id = "nav",
             collapsible = TRUE, position = "fixed-top",
             header = header),
        # Loads the interface of each tab
        callScriptsFunction(tabsFolder, c("name", "ui"), "ui")
    ))
)

shinyApp(ui, server)