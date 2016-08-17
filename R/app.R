#' @include globalAccess.R
NULL

# TODO(NunoA): increase allowed size and warn the user to wait for large files
# Refuse files with size greater than the specified
MB = 5000 # File size in MB
options(shiny.maxRequestSize = MB * 1024^2)

#' Get psichomics file inside a given directory
#' @param ... character vectors, specifying subdirectory and file(s) within some
#'  package. The default, none, returns the root of the package. Wildcards are
#'  not supported.
insideFile <- function(...) {
    return(system.file(..., package="psichomics"))
}

#' Check if a given function should be loaded by the calling module
#' @param loader Character: name of the file responsible to load such function 
#' @param FUN Function
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
#' @importFrom shiny callModule
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
            # Pass all attributes and add identifier
            attributes(res) <- c(attributes(res), attributes(FUN)[-1])
            return(res)
        }
    })
    # Remove NULL elements from list
    uiList <- Filter(Negate(is.null), uiList)
    return(uiList)
}

#' Create a selectize input available from any page
#' @param id Character: input identifier
#' @param label Character: label
#' @return HTML element for a global selectize input
globalSelectize <- function(id, placeholder) {
    elem <- paste0(id, "Elem")
    select <- selectizeInput(
        elem, "", choices = NULL,
        options = list(
            onDropdownClose = I(
                paste0("function($dropdown) { $('#", id,
                       "')[0].style.display = 'none'; }")),
            onBlur = I(
                paste0("function() { $('#", id,
                       "')[0].style.display = 'none'; }")),
            placeholder = placeholder),
        width="auto")
    select[[3]][[1]] <- NULL
    select <- tagAppendAttributes(
        select, id=id,
        style=paste("width: 95%;", "position: absolute;", 
                    "margin-top: 5px !important;", "display: none;"))
    return(select)
}

#' Create a special selectize input in the navigatin bar
#' @inheritParams globalSelectize
#' @return HTML element to be included in a navigation bar
navSelectize <- function(id, label, placeholder=label) {
    value <- paste0(id, "Value")
    tags$li( tags$div(
        class="navbar-text",
        style="margin-top: 5px !important; margin-bottom: 0px !important;", 
        globalSelectize(id, placeholder),
        tags$small(
            tags$b(label),
            tags$a(
                href="#", "Change...",
                onclick=paste0(
                    '$("#', id, '")[0].style.display = "block";',
                    '$("#', id, ' > div > select")[0].selectize.clear();',
                    '$("#', id, ' > div > select")[0].selectize.focus();'))), 
        tags$br(), textOutput(value)))
}

#' The user interface (ui) controls the layout and appearance of the app
#' All the CSS modifications are in the file "shiny/www/styles.css"
#' @importFrom shinyjs useShinyjs
#' @importFrom shiny tabPanel includeCSS includeScript conditionalPanel div h4 
#' icon shinyUI navbarPage
#' @return HTML elements
appUI <- function() {
    uiList <- getUiFunctions(paste, "app", tabPanel,
                             priority=c("dataUI", "analysesUI"))
    
    header <- tagList(
        includeCSS(insideFile("shiny", "www", "styles.css")),
        includeCSS(insideFile("shiny", "www", "animate.min.css")),
        includeScript(insideFile("shiny", "www", "functions.js")),
        includeScript(insideFile("shiny", "www", "fuzzy.min.js")),
        includeScript(insideFile("shiny", "www", "jquery.textcomplete.min.js")),
        conditionalPanel(
            condition="$('html').hasClass('shiny-busy')",
            div(class="text-right", id="loadmessage",
                h4(tags$span(class="label", class="label-info",
                             icon("flask", "fa-spin"), "Working...")))))
    
    nav <- do.call(navbarPage, c(
        list(title="PS\u03A8chomics", id="nav", collapsible=TRUE, 
             header=header, position="fixed-top", footer=useShinyjs()),
        uiList))
    
    nav[[3]][[1]][[3]][[1]][[3]][[2]] <- shiny::tagAppendChild(
        nav[[3]][[1]][[3]][[1]][[3]][[2]], 
        tags$ul(class="nav navbar-nav navbar-right",
                navSelectize("selectizeCategory", "Selected data category",
                             "Select data category"),
                navSelectize("selectizeEvent", "Selected splicing event",
                             "Search by gene, chromosome and coordinates")))
    shinyUI(nav)
}

#' Server function
#'
#' Instructions to build the Shiny app.
#'
#' @param input Input object
#' @param output Output object
#' @param session Session object
#' 
#' @importFrom shiny observe
appServer <- function(input, output, session) {
    getServerFunctions("app", priority=c("dataServer", "analysesServer"))
    
    # Update selectize input to show available categories
    observe({
        data <- getData()
        if (!is.null(data)) {
            updateSelectizeInput(session, "selectizeCategoryElem",
                                 choices=names(data))
            
            # Set the category of the data
            observeEvent(input$selectizeCategoryElem, 
                         if (input$selectizeCategoryElem != "")
                             setCategory(input$selectizeCategoryElem))
        } else {
            updateSelectizeInput(session, "selectizeCategoryElem",
                                 choices=list(), selected=list())
        }
    })
    
    # Update selectize event to show available events
    observe({
        psi <- getInclusionLevels()
        if (!is.null(psi)) {
            choices <- rownames(psi)
            names(choices) <- gsub("_", " ", rownames(psi))
            choices <- sort(choices)
            updateSelectizeInput(session, "selectizeEventElem", choices=choices)
            
            # Set the selected alternative splicing event
            observeEvent(input$selectizeEventElem,
                         if (input$selectizeEventElem != "")
                             setEvent(input$selectizeEventElem))
        } else {
            # Replace with empty list since NULLs are dropped
            updateSelectizeInput(session, "selectizeEventElem", choices=list(),
                                 selected=list())
        }
    })
    
    # Show the selected category
    output$selectizeCategoryValue <- renderText({
        category <- getCategory()
        if (is.null(category))
            return("No data loaded")
        else if(category == "")
            return("No category selected")
        else
            return(category)
    })
    
    # Show the selected event
    output$selectizeEventValue <- renderText({
        event <- getEvent()
        if (is.null(event))
            return("No data loaded")
        else if (event == "")
            return("No event selected")
        else
            return(gsub("_", " ", event))
    })
    
    # session$onSessionEnded(function() {
    #     # Stop app and print message to console
    #     suppressMessages(stopped <- stopApp(returnValue="Shiny app was closed"))
    # })
}

#' Start graphical interface of PSICHOMICS
#'
#' @param ... Parameters to pass to the function runApp
#' @param reload Boolean: reload package? FALSE by default
#'
#' @importFrom devtools load_all
#' @importFrom shiny shinyApp runApp
#'
#' @export
psichomics <- function(..., reload = FALSE) {
    if (reload) load_all()
    app <- shinyApp(appUI(), appServer)
    runApp(app, launch.browser = TRUE, ...)
}