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
linkToArticles <- function() {
    article <- list()
    article$description <- "Original article"
    article$authors <- c("N Saraiva-Agostinho", "NL Barbosa-Morais")
    article$title   <- paste(
        "psichomics: graphical application for alternative splicing",
        "quantification and analysis")
    article$year    <- 2019
    article$journal <- "Nucleic Acids Research"
    article$volume  <- 47
    article$number  <- 2
    article$pages   <- "e7"
    article$url     <- "https://doi.org/10.1093/nar/gky888"

    chapter <- list()
    chapter$description <- "Methods article"
    chapter$authors <- c("N Saraiva-Agostinho", "NL Barbosa-Morais")
    chapter$title   <- paste("Interactive Alternative Splicing Analysis of",
                             "Human Stem Cells Using psichomics")
    chapter$year    <- 2020
    chapter$journal <- "Methods in Molecular Biology"
    chapter$volume  <- 2117
    chapter$pages   <- "179-205"
    chapter$url     <- "https://doi.org/10.1007/978-1-0716-0301-7_10"

    prepareInfo <- function(data) {
        number <- data$number
        if (!is.null(number)) {
            number <- sprintf("(%s)", number)
        } else {
            number <- ""
        }
        link <- tags$a(target="_blank", href=data$url,
                       tags$b(data$title), sprintf("(%s)", data$year),
                       tags$i(paste0(data$journal, ".")),
                       sprintf("%s%s, %s", data$volume, number, data$pages))
        tagList(tags$dt(style="width: 110px", data$description),
                tags$dd(style="margin-left: 120px", link))
    }
    tags$div(class="alert alert-info", role="alert", style="padding: 10px",
             tags$dl(class="dl-horizontal",
                     style="margin-bottom: 0px !important;",
                     prepareInfo(article), prepareInfo(chapter)))
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
#' order; for instance, \code{c("data", "analyses")} would load \code{data},
#' then \code{analyses} and finally the remaining functions
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

#' Create a \code{selectize} input available from any page
#'
#' @param id Character: input identifier
#' @param placeholder Character: input placeholder
#' @param ASevent Boolean: select alternative splicing events?
#'
#' @importFrom shiny selectizeInput tagAppendAttributes
#'
#' @return HTML element for a global \code{selectize} input
#' @keywords internal
globalSelectize <- function(id, placeholder, ASevent=FALSE) {
    elem <- paste0(id, "Elem")
    hideElem <- sprintf("$('#%s')[0].style.display = 'none';", id)
    unmark   <- "this.$dropdown_content.unmark();"
    mark     <- "this.$dropdown_content.unmark()
                                       .mark(value, {exclude: ['text']});"

    onItemAdd <- I(paste("function(value, $item) {", hideElem, "}"))
    onBlur    <- I(paste("function() {", hideElem, "}"))
    onType    <- I(paste("function(value) {", mark, "}"))
    onLoad    <- I(paste(
        "function(data) { var value = this.currentResults.query;", mark, "}"))
    onOptionAdd <- I(paste(
        "function(value, data) {
            var tmp = data.label.split(\" __ \");
            data.svg = tmp[1];
            data.label = tmp[0];
            return(data);
        }"))
    onDropdownOpen <- I(paste("function($dropdown) {", unmark, "}"))

    render <- NULL
    if (ASevent) render <- I("{ option: renderEvent }")

    opts <- list(onItemAdd=onItemAdd, onBlur=onBlur, maxOptions=20,
                 placeholder=placeholder, render=render, highlight=FALSE,
                 onType=onType, onLoad=onLoad, onOptionAdd=onOptionAdd,
                 onDropdownOpen=onDropdownOpen)

    select <- selectizeInput(elem, "", choices=NULL, width="95%", options=opts)
    select[[3]][[1]] <- NULL
    select <- tagAppendAttributes(select, id=id, style=paste(
        "display: none;",
        "position: absolute;",  "margin-top: 5px !important;"))
    return(select)
}

#' Create a special \code{selectize} input in the navigation bar
#'
#' @inheritParams globalSelectize
#' @param label Character: input label
#'
#' @return HTML element to be included in a navigation bar
#' @keywords internal
navSelectize <- function(id, label, placeholder=label, ASevent=FALSE) {
    value <- paste0(id, "Value")
    tags$li( tags$div(
        class="navbar-text",
        style="margin-top: 5px !important; margin-bottom: 0px !important;",
        globalSelectize(id, placeholder, ASevent=ASevent),
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
#' @param menu Boolean: create a dropdown menu-like tab?
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

#' Replace a string with another in a list
#' @keywords internal
replaceStrInList <- function(tag, old, new) {
    FUN <- function(x) {
        res <- x
        if (grepl(old, x)) res <- gsub(old, new, x, fixed=TRUE)
        return(res)
    }
    rapply(tag, FUN, how="replace", classes="character")
}

#' Find an item in list of lists and return its coordinates
#' @keywords internal
traceInList <- function(ll, item) {
    if (is.list(ll)) {
        for (elem in seq(ll)) {
            res <- traceInList(ll[[elem]], item)
            if (!is.null(res)) return(c(elem, res))
        }
    } else if (is.character(ll)) {
        if (any(grepl(item, ll, fixed=TRUE))) return(numeric(0))
    }
}

#' User interface
#'
#' The user interface (UI) controls the layout and appearance of the app. All
#' CSS modifications are in the file \code{shiny/www/styles.css}
#'
#' @importFrom shinyjs useShinyjs
#' @importFrom shiny includeCSS includeScript conditionalPanel div h4 icon
#' shinyUI navbarPage tagAppendChild tagAppendAttributes
#' @importFrom purrr pluck pluck<-
#'
#' @return HTML elements
appUI <- function() {
    uiList <- getUiFunctions(paste, "app", modTabPanel,
                             priority=c("dataUI", "analysesUI"))

    header <- tagList(
        # Include CSS files
        includeCSS(insideFile("shiny", "www", "animate.compat.css")),
        includeCSS(insideFile("shiny", "www", "psichomics.css")),
        # Include JavaScript files
        includeScript(insideFile("shiny", "www", "jquery.mark.min.js")),
        includeScript(insideFile("shiny", "www", "highcharts.ext.js")),
        includeScript(insideFile("shiny", "www", "fuzzy.min.js")),
        includeScript(insideFile("shiny", "www", "jquery.textcomplete.min.js")),
        includeScript(insideFile("shiny", "www", "shinyBS.min.js")),
        includeScript(insideFile("shiny", "www", "psichomics.js")),
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
    nav <- replaceStrInList(nav, "navbar-header", "navbar-header hidden-sm")

    # Add global selectize input elements to navigation bar
    globalSelectizeElems <- tags$ul(
        class="nav navbar-nav navbar-right",
        navSelectize("selectizeCategory", "Selected dataset", "Select dataset"),
        navSelectize("selectizeEvent", "Selected splicing event",
                     "Search by gene and coordinates...", ASevent=TRUE))

    pos <- traceInList(nav, "navbar-nav")
    pos <- head(pos, -3)
    pluck(nav, !!!pos) <- tagList(pluck(nav, !!!pos, 1), globalSelectizeElems)
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
#' @inherit psichomics return
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

# Prepare representation of alternative splicing events
prepareASeventsRepresentation <- reactive({
    ASevent <- getASevents()
    if (!is.null(ASevent)) {
        diagram <- suppressWarnings(
            plotSplicingEvent(ASevent, class="pull-right"))
        parsed  <- parseSplicingEvent(ASevent, coords=TRUE, pretty=TRUE)
        coords  <- attr(diagram, "position")
        gene    <- prepareGenePresentation(parsed$gene)

        # Replace unsupported diagrams by text
        unsupported <- vapply(diagram, `==`, "", FUN.VALUE=logical(1))
        pos <- parsed$`full coordinates`
        if (is.null(pos)) {
            pos <- parsed$pos
            if (is.null(pos)) pos <- paste(parsed$start, parsed$end, sep=", ")
        }
        if (!is.null(pos)) {
            altText <- paste("altText:", prepareWordBreak(pos[unsupported]))
            diagram[unsupported] <- altText
            coords[unsupported]  <- paste("Full coordinates:", pos[unsupported])
        }
        id <- parsed$id
        if (is.null(parsed)) id <- ASevent
        info <- paste(sep=";", parsed$subtype,
                      sprintf("(chr%s, %s strand)", parsed$chr, parsed$strand),
                      id, gene, coords, ASevent)
        representation <- setNames(ASevent, paste(info, " __ ", diagram))
    } else {
        representation <- NULL
    }
    return(representation)
}, label="app_prepareASeventsRepresentation")

#' Server logic
#'
#' Instructions to build the Shiny app
#'
#' @param input Shiny input
#' @param output Shiny output
#' @param session Shiny session
#'
#' @importFrom shiny observe stopApp
#'
#' @inherit psichomics return
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

    # Update available alternative splicing events
    observe({
        representation <- prepareASeventsRepresentation()
        selected <- getASevent()
        if (!is.null(representation) && !is.null(selected)) {
            # Move the selected alternative splicing event to the top
            find <- match(selected, representation)
            if (!is.na(find)) {
                sort <- unique(c(find, seq(representation)))
                representation <- representation[sort]
            }
        }
        updateSelectizeChoices(session, "selectizeEventElem", representation,
                               server=TRUE)
    }, label="app_updateASevents")

    # Set alternative splicing event
    observeEvent(input[["selectizeEventElem"]], {
        selected <- input[["selectizeEventElem"]]
        if (!is.null(selected) && selected != "") {
            psi <- isolate(getInclusionLevels())
            attr(selected, "eventData") <- getSplicingEventData(psi)
            setEvent(selected)
        }
    })

    # Display selected category
    output$selectizeCategoryValue <- renderUI({
        category <- getCategory()
        if (is.null(category)) {
            return("No dataset loaded")
        } else if(category == "") {
            return("No dataset selected")
        } else {
            return(category)
        }
    })

    # Display selected event
    output$selectizeEventValue <- renderUI({
        areEventsLoaded  <- !is.null(getASevents())

        selected <- getASevent()
        isSelectionValid <- !is.null(selected) && selected != ""

        if (!areEventsLoaded) {
            return("No events quantified")
        } else if (!isSelectionValid) {
            return("No event is selected")
        } else {
            return(selected)
        }
    })

    if (!getOption("shinyproxy", FALSE)) {
        session$onSessionEnded(function() {
            # Stop app and print message to console
            message("\n-- psichomics was closed --")
            suppressMessages(stopApp())
        })
    }
}

#' Start graphical interface of psichomics
#'
#' @inheritParams shiny::runApp
#' @inheritDotParams shiny::runApp -appDir -launch.browser
#' @param reset Boolean: reset Shiny session? Requires package \code{devtools}
#' @param shinyproxy Boolean: prepare visual interface to run in Shinyproxy?
#' @param testData Boolean: load with test data
#' @param unparsableEvents Boolean: when testing data, load alternative splicing
#' quantification events that cannot be parsed?
#'
#' @importFrom shiny shinyApp runApp addResourcePath
#'
#' @return \code{NULL} (function is only used to modify the Shiny session's
#' state or internal variables)
#' @export
#'
#' @examples
#' \dontrun{
#' psichomics()
#' }
psichomics <- function(..., launch.browser=TRUE, reset=FALSE, shinyproxy=FALSE,
                       testData=FALSE, unparsableEvents=FALSE) {
    options(shinyproxy=shinyproxy)
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
        psi                        <- loadFile("vignettes/BRCA_psi.RDS")

        if (unparsableEvents) {
            rownames(psi) <- paste0("undefASevent", seq(nrow(psi)))
        }
        data[["Inclusion levels"]] <- psi
        data[["Sample metadata"]]  <- parseTCGAsampleInfo(colnames(psi))

        eventData <- suppressWarnings(
            parseSplicingEvent(rownames(psi), coords=TRUE))
        if (!is.null(eventData)) {
            class(eventData) <- c("eventData", class(eventData))
            attr(data[["Inclusion levels"]], "rowData") <- eventData
        }
        setData(list("Test data"=data))
    }
    app <- shinyApp(appUI(), appServer)
    runApp(app, launch.browser=launch.browser, ...)
}
