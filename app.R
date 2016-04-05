#' @import shiny shinyBS shinyjs
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
    
    # session$onSessionEnded(function() {
    #     # Stop app and print message to console
    #     suppressMessages(stopped <- stopApp(returnValue="Shiny app was closed"))
    # })
}

# The user interface (ui) controls the layout and appearance of the app
# All the CSS modifications are in the file "www/styles.css"
ui <- shinyUI(
    do.call(navbarPage, c(
        list(title = "PSΨchomics", id = "nav", collapsible = TRUE,
             position = "fixed-top",
             header = list(includeCSS("www/styles.css"),
                           includeScript("www/functions.js"),
                           conditionalPanel(
                               condition="$('html').hasClass('shiny-busy')",
                               div(icon("cog", "fa-spin"), "Working...",
                                   class="text-right", id="loadmessage")
                           )),
             footer = useShinyjs()),
        # Loads the interface of each tab
        callScriptsFunction(func = "ui", check = c("name", "ui"), tabPanel)
    ))
)

shinyApp(ui, server)