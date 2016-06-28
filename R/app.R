#' @include begin.R
NULL

#' Get psichomics file inside a given directory
#' @param ... character vectors, specifying subdirectory and file(s) within some
#'  package. The default, none, returns the root of the package. Wildcards are
#'  not supported.
insideFile <- function(...) {
    return(system.file(..., package="psichomics"))
}

#' Check if a given function should be loaded by a 
#' @param loader Character: name of the file responsible to load such function 
#' @param child Function
#' @return Boolean vector
loadBy <- function(loader, FUN) {
    attribute <- attr(FUN, "loader")
    if (is.null(attribute))
        return(FALSE)
    else
        return(attribute == loader)
}

#' Matches server functions from a given loader
#' @param ... Extra arguments to pass to server functions
#' @inheritParams getUiFunctions
#' 
#' @return Invisible TRUE
getServerFunctions <- function(loader, ..., priority=NULL) {
    # Get all functions ending with "Server"
    server <- ls(getNamespace("psichomics"), all.names=TRUE, pattern="Server$")
    server <- c(priority, server[!server %in% priority])
    
    lapply(server, function(name) {
        # Parse function name to get the function itself
        FUN <- eval(parse(text=name))
        # Check if module should be loaded by app
        if (loadBy(loader, FUN)) {
            # Remove last "Server" from the name and use it as ID
            id <- gsub("Server$", "", name)
            callModule(FUN, id, ...)
        }
    })
    return(invisible(TRUE))
}

#' Matches user interface (UI) functions from a given loader
#' 
#' @param ns Shiny function to create namespaced IDs
#' @param loader Character: loader to run the functions
#' @param ... Extra arguments to pass to the user interface (UI) functions
#' @param priority Character: name of functions to prioritise by the given
#' order; for instance, c("data", "analyses") would load "data", then "analyses"
#' then remaining functions
#' 
#' @return List of functions related to the given loader
getUiFunctions <- function(ns, loader, ..., priority=NULL) {
    # Get all functions ending with "UI"
    ui <- ls(getNamespace("psichomics"), all.names=TRUE, pattern="UI$")
    ui <- c(priority, ui[!ui %in% priority])
    
    # Get the interface of each tab
    uiList <- lapply(ui, function(name) {
        # Parse function name to get the function itself
        FUN <- eval(parse(text=name))
        # Check if module should be loaded by app
        if (loadBy(loader, FUN)) {
            # Remove last "UI" from the name and use it as ID
            id  <- gsub("UI$", "", name)
            res <- FUN(ns(id), ...)
            # Give a name to the UI
            attr(res, "name") <- attr(FUN, "name")
            return(res)
        }
    })
    # Remove NULL elements from list
    uiList <- Filter(Negate(is.null), uiList)
    return(uiList)
}

#' The user interface (ui) controls the layout and appearance of the app
#' All the CSS modifications are in the file "shiny/www/styles.css"
appUI <- function() {
    uiList <- getUiFunctions(paste, "app", tabPanel, priority=c("dataUI",
                                                                "analysesUI"))
    
    header <- list(
        includeCSS(insideFile("shiny", "www", "styles.css")),
        includeScript(insideFile("shiny", "www", "functions.js")),
        includeScript(insideFile("shiny", "www", "fuzzy.min.js")),
        includeScript(insideFile("shiny", "www", "jquery.textcomplete.min.js")),
        conditionalPanel(
            condition="$('html').hasClass('shiny-busy')",
            div(icon("flask", "fa-spin"), "Working...",
                class="text-right", id="loadmessage")
        )
    )
    
    shinyUI(
        do.call(navbarPage, c(
            list(title = "PSΨchomics", id = "nav", collapsible = TRUE,
                 position = "fixed-top",
                 header = header,
                 footer = shinyjs::useShinyjs()),
            uiList)
        )
    )
}

#' Server function
#'
#' Instructions to build the Shiny app.
#'
#' @param input Input object
#' @param output Output object
#' @param session Session object
appServer <- function(input, output, session) {
    getServerFunctions("app", priority=c("dataServer", "analysesServer"))
    
    # session$onSessionEnded(function() {
    #     # Stop app and print message to console
    #     suppressMessages(stopped <- stopApp(returnValue="Shiny app was closed"))
    # })

    # Save selected groups
    observe({
        sharedData$selectedGroups <- input$selectedGroups
        sharedData$javascriptRead <- TRUE
    })
}

#' Start graphical interface of PSICHOMICS
#'
#' @param ... Parameters to pass to the function runApp
#' @param reload Boolean: reload package? FALSE by default
#'
#' @importFrom devtools load_all
#' @importFrom shiny shinyApp
#'
#' @export
psichomics <- function(..., reload = FALSE) {
    if (reload) load_all()
    app <- shinyApp(appUI(), appServer)
    runApp(app, launch.browser = TRUE, ...)
}