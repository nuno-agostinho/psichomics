#' @importFrom Rcpp sourceCpp
#' @useDynLib psichomics, .registration=TRUE
#' @include globalAccess.R
NULL

# TODO(NunoA): increase allowed size and warn the user to wait for large files
# Refuse files with size greater than the specified
MB = 20 # File size in GB
options(shiny.maxRequestSize = MB * 1024^5)

# Sanitize errors
options(shiny.sanitize.errors = TRUE)

#' psichomics article's link interface
#'
#' @importFrom shiny tags icon
#' 
#' @return HTML elements
#' @keywords internal
linkToArticle <- function() {
    authors <- c("Nuno Saraiva-Agostinho", "Nuno L Barbosa-Morais")
    title   <- paste("psichomics: graphical application for alternative",
                     "splicing quantification and analysis.")
    year    <- 2019
    journal <- "Nucleic Acids Research"
    volume  <- 47
    number  <- 2
    pages   <- "e7"
    
    tags$a(
        target="_blank", href="https://doi.org/10.1093/nar/gky888",
        tags$div(
            class="alert alert-info", role="alert",
            icon("paper-plane-o"), 
            sprintf("%s (%s).", paste(authors, collapse=" and "), year),
            tags$b(title), tags$i(paste0(journal, ".")),
            sprintf("%s(%s), %s", volume, number, pages)))
}

#' Check if a given function should be loaded by the calling module
#' @param loader Character: name of the file responsible to load such function 
#' @param FUN Function
#' @return Boolean vector
#' @keywords internal
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
#' @keywords internal
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
#' @param ns Shiny function to create IDs within a namespace
#' @param loader Character: loader to run the functions
#' @param ... Extra arguments to pass to the user interface (UI) functions
#' @param priority Character: name of functions to prioritise by the given
#' order; for instance, c("data", "analyses") would load "data", then "analyses"
#' then remaining functions
#' 
#' @return List of functions related to the given loader
#' @keywords internal
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
#' @param placeholder Character: input placeholder
#' 
#' @importFrom shiny selectizeInput tagAppendAttributes
#' 
#' @return HTML element for a global selectize input
#' @keywords internal
globalSelectize <- function(id, placeholder) {
    elem <- paste0(id, "Elem")
    hideElem <- sprintf("$('#%s')[0].style.display = 'none';", id)
    
    select <- selectizeInput(elem, "", choices=NULL, width="auto", options=list(
        onItemAdd=I(paste0("function(value, $item) {", hideElem, "}")),
        onBlur=I(paste0("function() {", hideElem, "}")),
        placeholder=placeholder))
    select[[3]][[1]] <- NULL
    select <- tagAppendAttributes(select, id=id, style=paste(
        "display: none;",
        "width: 95%;", "position: absolute;",  "margin-top: 5px !important;"))
    return(select)
}

#' Create a special selectize input in the navigation bar
#' 
#' @inheritParams globalSelectize
#' @param label Character: input label
#' 
#' @return HTML element to be included in a navigation bar
#' @keywords internal
navSelectize <- function(id, label, placeholder=label) {
    value <- paste0(id, "Value")
    tags$li( tags$div(
        class="navbar-text",
        style="margin-top: 5px !important; margin-bottom: 0px !important;", 
        globalSelectize(id, placeholder),
        tags$small(tags$b(label), tags$a(
            "Change...", onclick=paste0(
                '$("#', id, '")[0].style.display = "block";',
                '$("#', id, ' > div > select")[0].selectize.clear();',
                '$("#', id, ' > div > select")[0].selectize.focus();'))), 
        tags$br(), uiOutput(value)))
}

#' Modified \code{tabPanel} function to show icon and title
#' 
#' @note Icon is hidden at small viewports
#' 
#' @param title Character: title of the tab
#' @param icon Character: name of the icon
#' @param ... HTML elements to render
#' @param menu Boolean: create a dropdown menu-like tab? FALSE by default
#' 
#' @importFrom shiny navbarMenu tabPanel
#' 
#' @return HTML interface
#' @keywords internal
modTabPanel <- function(title, ..., icon=NULL, menu=FALSE) {
    if (is.null(icon))
        display <- title
    else
        display <- tagList(icon(class="hidden-sm", icon), title)
    
    if (menu)
        navbarMenu(display, ...)
    else
        tabPanel(display, ..., value=title)
}

#' User interface
#' 
#' The user interface (UI) controls the layout and appearance of the app. All
#' CSS modifications are in the file \code{shiny/www/styles.css}
#' 
#' @importFrom shinyjs useShinyjs
#' @importFrom shiny includeCSS includeScript conditionalPanel div h4 icon
#' shinyUI navbarPage tagAppendChild tagAppendAttributes
#' 
#' @return HTML elements
appUI <- function() {
    uiList <- getUiFunctions(paste, "app", modTabPanel,
                             priority=c("dataUI", "analysesUI"))
    
    header <- tagList(
        includeCSS(insideFile("shiny", "www", "styles.css")),
        includeCSS(insideFile("shiny", "www", "animate.min.css")),
        includeScript(insideFile("shiny", "www", "functions.js")),
        includeScript(insideFile("shiny", "www", "highcharts.ext.js")),
        includeScript(insideFile("shiny", "www", "fuzzy.min.js")),
        includeScript(insideFile("shiny", "www", "jquery.textcomplete.min.js")),
        conditionalPanel(
            condition="$('html').hasClass('shiny-busy')",
            div(class="text-right", id="loadmessage",
                h4(tags$span(class="label", class="label-info",
                             icon("flask", "fa-spin"), "Working...")))))
    
    nav <- do.call(navbarPage, c(
        list(title="psichomics", id="nav", collapsible=TRUE, 
             header=header, position="fixed-top", footer=useShinyjs()),
        uiList))

    # Hide the header from the navigation bar if the viewport is small
    nav[[3]][[1]][[3]][[1]][[3]][[1]] <- tagAppendAttributes(
        nav[[3]][[1]][[3]][[1]][[3]][[1]], class="hidden-sm")
    
    # Add global selectize input elements to navigation bar
    nav[[3]][[1]][[3]][[1]][[3]][[2]] <- tagAppendChild(
        nav[[3]][[1]][[3]][[1]][[3]][[2]], 
        tags$ul(class="nav navbar-nav navbar-right",
                navSelectize("selectizeCategory", "Selected dataset",
                             "Select dataset"),
                navSelectize("selectizeEvent", "Selected splicing event",
                             "Search by gene, chromosome and coordinates")))
    shinyUI(nav)
}

#' Enable history navigation
#' 
#' Navigate app according to the location given by the navigation bar. Code
#' and logic adapted from
#' \url{https://github.com/daattali/advanced-shiny/blob/master/navigate-history}
#' 
#' @param navId Character: identifier of the navigation bar
#' @param input Input object
#' @param session Session object
#' 
#' @importFrom shiny observe parseQueryString updateTabsetPanel
#' 
#' @return NULL (this function is used to modify the Shiny session's state)
#' @keywords internal
browserHistory <- function(navId, input, session) {
    # Update browser history when user changes the active tab
    observeEvent(input[[navId]], {
        autoNav <- getAutoNavigation()
        if (isTRUE(autoNav)) {
            setAutoNavigation(FALSE)
        } else {
            # Update browser history
            runjs(paste0("updateHistory({ page: '", input[[navId]], "'})"))
        }
    })
    
    # Navigate to a tab according to a given query string
    restorePage <- function(qs) {
        data <- parseQueryString(qs)
        if (!is.null(data$page)) {
            setAutoNavigation(TRUE)
            updateTabsetPanel(session, navId, data$page)
        }
    }
    
    # Navigate tabs while browsing history
    observeEvent(input$appLocation, { restorePage(input$appLocation) })
}

#' Server logic
#' 
#' Instructions to build the Shiny app
#'
#' @param input Input object
#' @param output Output object
#' @param session Session object
#' 
#' @importFrom shiny observe stopApp
#' 
#' @return NULL (this function is used to modify the Shiny session's state)
appServer <- function(input, output, session) {
    ns <- session$ns
    groupsServerOnce(input, output, session)
    getServerFunctions("app", priority=c("dataServer", "analysesServer"))
    browserHistory("nav", input, session)
    
    updateSelectizeChoices <- function(session, id, choices, server=FALSE) {
        if (!is.null(choices)) {
            selected <- choices[[1]]
        } else {
            choices  <- list()
            selected <- list()
        }
        updateSelectizeInput(session, id, choices=choices, selected=selected,
                             server=server)
    }
    
    # Update available categories
    observe(updateSelectizeChoices(session, "selectizeCategoryElem", 
                                   names(getData()), server=FALSE))
    
    # Set data category
    observeEvent(input$selectizeCategoryElem, {
        selected <- input$selectizeCategoryElem
        if (!is.null(selected) && selected != "") setCategory(selected)
    })
    
    # Update available events
    observe(updateSelectizeChoices(session, "selectizeEventElem", 
                                   getASevents(), server=TRUE))
    
    # Set alternative splicing event
    observeEvent(input[["selectizeEventElem"]], {
        selected <- input[["selectizeEventElem"]]
        if (!is.null(selected) && selected != "") setEvent(selected)
    })
    
    # Display selected category
    output$selectizeCategoryValue <- renderUI({
        category <- getCategory()
        if (is.null(category))
            return("No dataset loaded")
        else if(category == "")
            return("No dataset selected")
        else
            return(category)
    })
    
    # Display selected event
    output$selectizeEventValue <- renderUI({
        areEventsLoaded  <- !is.null(getASevents())
        
        selected <- getASevent()
        isSelectionValid <- !is.null(selected) && selected != ""
        
        if (!areEventsLoaded)
            return("No events quantified")
        else if (!isSelectionValid)
            return("No event is selected")
        else
            return(parseSplicingEvent(selected, char=TRUE))
    })
    
    session$onSessionEnded(function() {
        # Stop app and print message to console
        message("\n-- psichomics was closed --")
        suppressMessages(stopApp())
    })
}

#' Start graphical interface of psichomics
#'
#' @inheritDotParams shiny::runApp -appDir -launch.browser
#' @param reset Boolean: reset Shiny session? Requires package \code{devtools}
#' @param testData Boolean: auto-start with test data
#'
#' @importFrom shiny shinyApp runApp addResourcePath
#'
#' @return NULL (this function is used to modify the Shiny session's state)
#' @export
#'
#' @examples
#' \dontrun{
#' psichomics()
#' }
psichomics <- function(..., reset=FALSE, testData=FALSE) {
    # Add icons related to set operations
    addResourcePath("set-operations",
                    insideFile("shiny", "www", "set-operations"))
    
    if (reset) devtools::load_all()
    
    if (testData) {
        loadFile <- function(file) {
            if (!file.exists(file)) {
                # Fetch file online if not locally available
                link <- paste0("https://github.com/",
                               "nuno-agostinho/psichomics/raw/master/",
                               file)
                file <- url(link)
            }
            readRDS(file)
        }
        
        data <- NULL
        data[["Clinical data"]]    <- loadFile("vignettes/BRCA_clinical.RDS")
        data[["Gene expression"]]  <- loadFile("vignettes/BRCA_geneExpr.RDS")
        data[["Inclusion levels"]] <- loadFile("vignettes/BRCA_psi.RDS")
        data[["Sample metadata"]]  <- parseTcgaSampleInfo(colnames(
            data[["Inclusion levels"]]))
        setData(list("Test data"=data))
    }
    
    app <- shinyApp(appUI(), appServer)
    runApp(app, launch.browser = TRUE, ...)
}