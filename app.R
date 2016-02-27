source("R/initial.R")

#' Server function
#' 
#' Instructions to build the Shiny app.
#'
#' @param input Input object
#' @param output Output object
#' @param session Session object
server <- function(input, output, session) {
    callScriptsFunction(func = "server", input, output, session,
                        check = c("name", "server"))
    
    session$onSessionEnded(function() {
        # Stop app and print message to console
        suppressMessages(stopped <- stopApp(returnValue=TRUE))
        if (stopped) cat("\nShiny app was closed ")
    })
}

addCSS <- function(...) tags$style(type = "text/css", ...)

header <- list(
    useShinyjs(),
    ## TODO(NunoA): put all CSS info in a CSS file and import with "includeCSS"
    # Avoids fixed-top navbar from obscuring content
    addCSS("body { padding-top: 70px; }"),
    # Alert appears fixed on the top right above other elements
    addCSS(".sbs-alert { position:fixed;
                         right:10px;
                         z-index:9;
                         -webkit-filter: opacity(80%);
                         filter: opacity(80%); }"),
    addCSS(".shiny-progress-container { top: 48px; }"))

# The user interface (ui) controls the layout and appearance of the app
ui <- shinyUI(
    do.call(navbarPage, c(
        list(title = "spliced canceR", id = "nav", collapsible = TRUE,
             position = "fixed-top", header = header),
        # Loads the interface of each tab
        callScriptsFunction(func = "ui", check = c("name", "ui"), tabPanel)
    ))
)

shinyApp(ui, server)