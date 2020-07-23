#' Filter alternative splicing quantification
#' 
#' @param psi Data frame or matrix: alternative splicing quantification
#' @param eventType Character: filter data based on event type; check all event
#' types available by using \code{getSplicingEventTypes(psi)}, where \code{psi}
#' is the alternative splicing quantification data; if \code{eventType = NULL},
#' events are not filtered by event type
#' @param eventSubtype Character: filter data based on event subtype; check all
#' event subtypes available in your data by using
#' \code{unique(getSplicingEventInformation(psi)$subtype)}, where \code{psi} is
#' the alternative splicing quantification data; if \code{eventSubtype = NULL},
#' events are not filtered by event subtype
#' @param minMedian Numeric: minimum of read count median per splicing event
#' @param maxMedian Numeric: maximum of read count median per splicing event
#' @param minLogVar Numeric: minimum log10(read count variance) per splicing
#' event
#' @param maxLogVar Numeric: maximum log10(read count variance) per splicing
#' event
#' @param minRange Numeric: minimum range of read counts across samples per 
#' splicing event
#' @param maxRange Numeric: maximum range of read counts across samples per 
#' splicing event
#' 
#' @family functions for PSI quantification
#' @return Boolean vector indicating which splicing events pass the thresholds
#' @export
#' 
#' @examples 
#' # Calculate PSI for skipped exon (SE) and mutually exclusive (MXE) events
#' annot <- readFile("ex_splicing_annotation.RDS")
#' junctionQuant <- readFile("ex_junctionQuant.RDS")
#' 
#' psi <- quantifySplicing(annot, junctionQuant, eventType=c("SE", "MXE"))
#' # Filter PSI
#' psi[filterPSI(psi, minMedian=0.05, maxMedian=0.95, minRange=0.15), ]
filterPSI <- function(psi, eventType=NULL, eventSubtype=NULL,
                      minMedian=-Inf, maxMedian=Inf,
                      minLogVar=-Inf, maxLogVar=Inf,
                      minRange=-Inf, maxRange=Inf,
                      minPSI=-Inf, maxPSI=Inf,
                      vasttoolsCoverage=NULL) {
    if (is.na(minMedian)) minMedian <- -Inf
    if (is.na(maxMedian)) maxMedian <- Inf
    if (is.na(minLogVar)) minLogVar <- -Inf
    if (is.na(maxLogVar)) maxLogVar <- Inf
    if (is.na(minRange)) minRange <- -Inf
    if (is.na(maxRange)) maxRange <- Inf
    if (is.na(minPSI)) minPSI <- -Inf
    if (is.na(maxPSI)) maxPSI <- Inf
    
    trueVector <- function(len) rep(TRUE, len)
    
    type <- getSplicingEventInformation(psi)$type
    if (!is.null(eventType) && !is.null(type)) {
        type <- type %in% eventType
        eventTypeText <- c("Type"=paste(eventType, collapse = ", "))
    } else {
        type <- trueVector(nrow(psi))
        eventTypeText <- NULL
    }
    
    subtype <- getSplicingEventInformation(psi)$subtype
    if (!is.null(eventSubtype) && !is.null(subtype)) {
        subtype <- subtype %in% eventSubtype
        eventSubtypeText <- c("Type"=paste(eventSubtype, collapse = ", "))
    } else {
        subtype <- trueVector(nrow(psi))
        eventSubtypeText <- NULL
    }
    
    if (!is.infinite(minMedian) & !is.infinite(maxMedian)) {
        medians <- customRowMedians(psi, na.rm=TRUE, fast=TRUE)
        medianThres <- medians >= minMedian & medians <= maxMedian
    } else {
        medianThres <- trueVector(nrow(psi))
    }
    
    if (!is.infinite(minLogVar) & !is.infinite(maxLogVar)) {
        vars <- log10(customRowVars(psi, na.rm=TRUE, fast=TRUE))
        varThres <- vars >= minLogVar & vars <= maxLogVar
    } else {
        varThres <- trueVector(nrow(psi))
    }
    
    if (!is.infinite(minPSI) & !is.infinite(maxPSI)) {
        mini <- customRowMins(psi, na.rm=TRUE, fast=TRUE)
        maxi <- customRowMaxs(psi, na.rm=TRUE, fast=TRUE)
        psiThres <- mini >= minPSI & maxi <= maxPSI
    } else {
        mini <- maxi <- NULL
        psiThres <- trueVector(nrow(psi))
    }
    
    if (!is.infinite(minRange) & !is.infinite(maxRange)) {
        if (!is.null(mini) && !is.null(maxi)) {
            ranges <- maxi - mini
        } else {
            ranges <- customRowRanges(psi, na.rm=TRUE, fast=TRUE)
        }
        rangeThres <- ranges >= minRange & ranges <= maxRange
    } else {
        rangeThres <- trueVector(nrow(psi))
    }
    
    thres <- which(type & subtype &
                       psiThres & medianThres & varThres & rangeThres)
    
    attr(thres, "filtered") <- c("Filter enabled"="Yes",
                                 eventTypeText,
                                 eventSubtypeText,
                                 "PSI >="=minPSI,
                                 "PSI <="=maxPSI,
                                 "Median PSI >="=minMedian,
                                 "Median PSI <="=maxMedian,
                                 "log10(PSI variance) >="=minLogVar,
                                 "log10(PSI variance) <="=maxLogVar,
                                 "PSI range >="=minRange,
                                 "PSI range <="=maxRange)
    return(thres)
}

#' Remove alternative splicing quantification values based on coverage
#' 
#' @param psi Data frame or matrix: alternative splicing quantification
#' @param minReads Currently this argument does nothing
#' @param vasttoolsCoverage Character: if you are using inclusion levels from
#' VAST-TOOLS, filter the data based on quality scores for read coverage, e.g.
#' use \code{vasttoolsCoverage = c("SOK", "OK", "LOW")} to only keep events with
#' good read coverage (by default, events are not filtered based on quality
#' scores); read \url{https://github.com/vastgroup/vast-tools} for more
#' information on VAST-TOOLS quality scores
#' 
#' @return Alternative splicing quantification data with missing values for any
#' values with insufficient coverage
#' @export
discardLowCoveragePSIvalues <- function(
    psi, minReads=10, vasttoolsScoresToDiscard=c("VLOW", "N")) {
    
    eventData <- getSplicingEventInformation(psi)
    if (is.null(eventData)) return(psi)
    
    # Remove events containing only missing values
    removeNAonlyEvents <- function(psi) psi[rowSums(!is.na(psi)) > 0, ]
    
    isVastTools <- unique(eventData$source) == "vast-tools"
    if (!is.null(vasttoolsScoresToDiscard) && isVastTools) {
        quality     <- eventData[ , endsWith(colnames(eventData), "-Q")]
        qualityVals <- unlist(quality)
        if ("OK" %in% vasttoolsScoresToDiscard) {
            # Ignore OK from score 4
            qualityVals <- gsub("OK@", "@", qualityVals)
        }
        scores <- sprintf(",%s,", vasttoolsScoresToDiscard)
        lowCvg <- Reduce(`|`, lapply(scores, grepl, qualityVals, fixed=TRUE))
        lowCvg <- matrix(lowCvg, ncol=ncol(quality))
        psi[lowCvg] <- NA
        psi <- removeNAonlyEvents(psi)
        attr(psi, "filtered") <- c(
            "VAST-TOOLS coverage"=vasttoolsScoresToDiscard)
    }
    return(psi)
}

getPSIsummaryStats <- function() {
    c("Mean PSI"="mean", "Median PSI"="median",
      "PSI variance"="var", "log10(PSI variance)"="log10(var)",
      "PSI range"="range")
}

#' Interface to filter alternative splicing
#'
#' @param ns Namespace function
#'
#' @importFrom shiny tagList uiOutput selectizeInput numericInput actionButton
#' @importFrom shinyBS bsTooltip
#' @importFrom shinyjs hidden disabled
#'
#' @return HTML elements
#' @keywords internal
inclusionLevelsFilterInterface <- function(ns) {
    filterGenesSelectize <- selectizeInput(
        ns("filterGenes"), label=NULL, selected=NULL, multiple=TRUE,
        width="100%", choices=c("Type to search for genes..."=""),
        options=list(plugins=list("remove_button")))
    
    sampleFiltering <- bsCollapsePanel(
        tagList(icon("filter"), "Sample filtering",
                contextUI(ns("sampleFilterText"))),
        value="Sample filtering",
        selectizeInput(ns("sampleFilter"), "Samples to discard", multiple=TRUE,
                       width="100%", choices=character(0)))
    
    numInputCol <- function(...) {
        column(6, numericInput(..., width="100%"))
    }
    
    psiFiltering <- bsCollapsePanel(
        tagList(icon("filter"), "PSI filtering",
                contextUI(ns("psiFilterText"))),
        value="PSI filtering",
        checkboxInput(ns("enablePSIfiltering"), value=FALSE, width="100%",
                      "Filter splicing events based on their PSI values"),
        fluidRow(
            numInputCol(ns("minMedian"), "Median >=",
                        min=0, max=1, value=0, step=0.1),
            numInputCol(ns("maxMedian"), "Median <=",
                        min=0, max=1, value=1, step=0.1)),
        fluidRow(
            numInputCol(ns("minLogVar"), "log10(variance) >=",
                        min=-10, max=0, value=-10, step=0.5),
            numInputCol(ns("maxLogVar"), "log10(variance) <=",
                        min=-10, max=0, value=0, step=0.5)),
        fluidRow(
            numInputCol(ns("minRange"), "Range >=",
                        min=0, max=1, value=0, step=0.1),
            numInputCol(ns("maxRange"), "Range <=",
                        min=0, max=1, value=1, step=0.1)))
    
    
    unparsableError <- function(id) {
        errorDialog(
            paste("Loaded splicing events could not be parsed. Data cannot be",
                  "filtered based on event type or cognate genes."),
            id=ns(id), style="margin: 10px;")
    }
    geneFiltering <- bsCollapsePanel(
        title=tagList(icon("filter"), "Cognate gene filtering",
                      contextUI(ns("geneFilterText"))),
        value="Filter by genes",
        hidden(unparsableError("unparsableGenes")),
        tags$span(id=ns("geneOptions"),
                  radioButtons(
                      ns("filter"), label=NULL,
                      c("Do not filter splicing events by genes"="noFilter",
                        "Filter by selected genes"="select",
                        "Filter by genes imported from a file"="file")),
                  conditionalPanel(
                      sprintf("input[id='%s'] == '%s'", ns("filter"), "select"),
                      filterGenesSelectize),
                  conditionalPanel(
                      sprintf("input[id='%s'] == '%s'", ns("filter"), "file"),
                      fileBrowserInput(
                          ns("filterGenesFile"), label=NULL,
                          clearable=TRUE, placeholder="No file selected"),
                      helpText(
                          "Provide a file containing gene symbols separated",
                          "by space, comma, tab or new line. For instance: ",
                          tags$code("BRCA1, BRAF, ABL")))))
    
    choicesX <- getPSIsummaryStats()
    choicesY <- c(choicesX, "Density"="none")
    options <- div(
        id=ns("options"),
        selectizeInput(ns("eventSubtype"), "Event types to keep", width="100%",
                       selected=NULL, choices=NULL, multiple=TRUE),
        hidden(unparsableError("unparsableEventSubtypes")),
        bsCollapse(sampleFiltering, psiFiltering, geneFiltering),
        checkboxInput(ns("preview"), value=FALSE,
                      "Preview plot (slow for large datasets)"),
        conditionalPanel(sprintf("input[id='%s']", ns("preview")),
                         fluidRow(
                             column(6, selectizeInput(
                                 ns("plotX"), "X axis", width="100%",
                                 choices=choicesX, selected="range")),
                             column(6, selectizeInput(
                                 ns("plotY"), "Y axis", width="100%",
                                 choices=choicesY, selected="none"))),
                         uiOutput(ns("alert")),
                         plotOutput(ns("plot"), height="200px"),
                         tags$small(helpText(class="pull-right", 
                                             textOutput(ns("summary"))))))
    tagList(
        uiOutput(ns("modal")),
        errorDialog(
            "No alternative splicing quantification loaded or quantified.",
            id=ns("missingData"), style="margin: 10px;"),
        hidden(options),
        disabled(processButton(ns("filterIncLevels"),
                               "Filter alternative splicing")))
}

#' @rdname appUI
inclusionLevelsFilterUI <- function(id, panel) {
    ns <- NS(id)
    title <- "Alternative splicing quantification filtering"
    panel(style="success", title=tagList(icon("filter"), title),
          value=title, inclusionLevelsFilterInterface(ns))
}

getCognateGenesToFilter <- function(inputFilter, inputFilterGenes,
                                    inputFilterGenesFile) {
    filter <- NULL
    isFileChosen <- !is.null(inputFilterGenesFile) && inputFilterGenesFile != ""
    if (inputFilter == "select") {
        # Get genes based on select input
        filter <- inputFilterGenes
        if (identical(filter, "")) filter <- NULL
    } else if (inputFilter == "file" && isFileChosen) {
        # Get genes provided in a file
        filter <- fread(inputFilterGenesFile, header=FALSE)
        if (nrow(filter) == 1) {
            filter <- as.character(filter)
        } else if (nrow(filter) > 1) {
            filter <- as.character(filter[[1]])
        } else {
            filter <- NULL
        }
    }
    return(filter)
}

#' @importFrom shiny updateCheckboxInput
#' @importFrom shinyjs show hide enable disable
psiFilteringOptionsSet <- function(session, input, output) {
    observeEvent(input$missing, missingDataGuide("Inclusion Levels"))
    prepareFileBrowser(session, input, "filterGenesFile")
    
    # Warn user if alternative splicing quantification is not loaded
    observe({
        if (is.null(getData()) || is.null(getInclusionLevels())) {
            hide("options")
            disable("filterIncLevels")
            show("missingData")
        } else {
            show("options")
            enable("filterIncLevels")
            hide("missingData")
        }
    })
    
    # Update event types
    observe({
        processed <- processPSI()
        if (!is.null(processed)) {
            types <- table(processed$subtype)
            types <- setNames(names(types), as.vector(types))
            updateSelectizeInput(
                session, "eventSubtype", choices=types, selected=types,
                options=list(
                    plugins=list("remove_button"),
                    render=I('{option: renderASeventTypeSelection,
                               item: renderASeventTypeSelection}')))
            hide("unparsableEventSubtypes")
            show("eventSubtype")
        } else {
            show("unparsableEventSubtypes")
            hide("eventSubtype")
        }
    })
    
    # Update sample filtering options
    observe({
        psi <- getInclusionLevels()
        if (!is.null(psi)) {
            updateSelectizeInput(
                session, "sampleFilter", server=TRUE, choices=colnames(psi),
                options=list(placeholder="Select samples to discard",
                             plugins=list("remove_button")))
        }
    })
    
    # Update gene for filtering based on alternative splicing quantification
    observe({
        filter    <- input$filter
        processed <- processPSI()
        
        if (!is.null(processed)) {
            # Avoid loading if already loaded
            if (filter == "select") {
                # Load genes based on alternative splicing quantification
                genes <- sort(unique(unlist(processed$gene)))
                updateSelectizeInput(session, "filterGenes", choices=genes,
                                     selected=character(0), server=TRUE)
            }
            
            hide("unparsableGenes")
            show("geneOptions")
        } else {
            show("unparsableGenes")
            hide("geneOptions")
        }
    })
    
    # Toggle PSI filtering options
    observe({
        if (input$enablePSIfiltering) {
            enable("minMedian")
            enable("maxMedian")
            enable("minLogVar")
            enable("maxLogVar")
            enable("minRange")
            enable("maxRange")
        } else {
            disable("minMedian")
            disable("maxMedian")
            disable("minLogVar")
            disable("maxLogVar")
            disable("minRange")
            disable("maxRange")
        }
    })
    
    # Update context
    output$sampleFilterText <- renderText({
        sampleFilter <- input$sampleFilter
        if (is.null(sampleFilter) || sampleFilter == "") {
            text <- "No samples to discard"
        } else {
            len  <- length(sampleFilter)
            text <- sprintf("%s sample%s to discard",
                            len, ifelse(len == 1, "", "s"))
        }
        return(text)
    })
    output$psiFilterText    <- renderText({
        if (input$enablePSIfiltering) {
            text <- "Enabled"
        } else {
            text <- "Disabled"
        }
        return(text)
    })
    output$geneFilterText <- renderText({
        filter <- getCognateGenesToFilter(input$filter, input$filterGenes,
                                          input$filterGenesFile)
        if (is.null(filter)) {
            text <- "Not filtering by genes"
        } else {
            len  <- length(filter)
            text <- sprintf("Filtering by %s gene%s",
                            len, ifelse(len == 1, "", "s"))
        }
        return(text)
    })
}

# Process PSI
processPSI <- reactive({
    psi <- getInclusionLevels()
    if (!is.null(psi)) {
        parsed <- parseSplicingEvent(rownames(psi), data=psi, pretty=TRUE)
        return(parsed)
    }
})

filterSplicingOperation <- function(session, psi, processed, input) {
    areThereInclusionLevels <- function(psi) {
        return(!is.null(psi) && nrow(psi) != 0 && ncol(psi) != 0)
    }
    filteredPSI  <- NULL
    eventSubtype <- input$eventSubtype
    if (eventSubtype == "") eventSubtype <- NULL
    if (!areThereInclusionLevels(psi)) return(NULL)

    suppressWarnings(updateProgress("Discarding samples"))
    sampleFilter <- input$sampleFilter
    if (!is.null(sampleFilter) && sampleFilter != "") {
        samplesToKeep    <- !colnames(filteredPSI) %in% sampleFilter
        filteredPSI      <- filteredPSI[ , samplesToKeep]
        sampleFilterText <- paste(sampleFilter, collapse=", ")
    } else {
        sampleFilterText <- "None"
    }
    sampleFilterSettings <- c("Discarded samples"=sampleFilterText)
    if (!areThereInclusionLevels(filteredPSI)) return(NULL)
    
    suppressWarnings(updateProgress("Filtering based on PSI values"))
    enablePSIfiltering <- input$enablePSIfiltering
    if (enablePSIfiltering) {
        filtered <- filterPSI(
            filteredPSI, eventSubtype=eventSubtype,
            minMedian=input$minMedian, maxMedian=input$maxMedian,
            minLogVar=input$minLogVar, maxLogVar=input$maxLogVar,
            minRange=input$minRange, maxRange=input$maxRange)
        filteredPSI    <- filteredPSI[filtered, ]
        filterSettings <- attr(filtered, "filtered")
    } else {
        filteredPSI    <- filteredPSI
        filterSettings <- c("Filter enabled"="No")
    }
    if (!areThereInclusionLevels(filteredPSI)) return(NULL)

    suppressWarnings(updateProgress("Filtering based on genes"))
    if (!is.null(processed$gene)) {
        filter <- getCognateGenesToFilter(input$filter, input$filterGenes,
                                          input$filterGenesFile)
        if (!is.null(filter)) {
            filteredPSI <- filteredPSI[processed$gene %in% filter, ]
        }
    } else {
        filter <- NULL
    }
    if (!areThereInclusionLevels(filteredPSI)) return(NULL)

    # Include settings used for alternative splicing quantification
    settings <- c(
        list(
            "Splicing event types"=input$eventSubtype,
            "Selected cognate genes"=if (is.null(
                filter)) "All available genes" else filter),
        sampleFilterSettings, filterSettings)
    attr(filteredPSI, "settings") <- settings
    attr(filteredPSI, "icon") <- list(symbol="calculator", colour="green")
    suppressWarnings(closeProgress())
    return(filteredPSI)
}

#' Set of functions to filter alternative splicing
#'
#' @importFrom shiny tags
#' @inherit inclusionLevelsFilterServer
#' @keywords internal
filterSplicingSet <- function(session, input, output) {
    ns <- session$ns
    
    # Operation to filter splicing based on user-defined options
    filterSplicingBasedOnInput <- reactive({
        psi       <- getInclusionLevels()
        processed <- processPSI()
        filterSplicingOperation(session, psi, processed, input)
    })
    
    # Filter inclusion levels
    filterSplicing <- function() {
        time <- startProcess("filterIncLevels")
        startProgress("Filtering alternative splicing", divisions=4)
        psi  <- filterSplicingBasedOnInput()
        if (!is.null(psi)) {
            setInclusionLevels(psi)
        } else {
            errorModal(session, "No splicing events returned",
                       "The filtering was too strict and returned no",
                       "alternative splicing events.",
                       caller="Alternative splicing quantification filtering")
            time <- NULL
        }
        endProcess("filterIncLevels", time)
    }
    
    # Show warnings if needed before filtering alternative splicing
    observeEvent(input$filterIncLevels, {
        if (is.null(getData()) || is.null(getInclusionLevels())) {
            missingDataModal(session, "Inclusion levels", ns("missing"))
        } else if (!is.null(getInclusionLevels()) &&
                   !is.null(getDifferentialSplicing())) {
            warningModal(
                session, "Differential splicing already performed",
                "Do you wish to replace the current alternative splicing",
                "quantification data and discard differential splicing",
                "results?", footer=actionButton(
                    ns("discard2"), "Replace and discard",
                    class="btn-warning", "data-dismiss"="modal"),
                caller="Alternative splicing quantification filtering")
        } else {
            filterSplicing()
        }
    })
    
    observeEvent(input$replace, filterSplicing())
    observeEvent(input$discard, {
        setDifferentialSplicing(NULL)
        setDifferentialSplicingSurvival(NULL)
        filterSplicing()
    })
    
    output$plot <- renderPlot({
        if (input$preview) {
            x <- input$plotX
            y <- input$plotY
            
            if (x == "none") return(NULL)
            stats  <- getPSIsummaryStats()
            xLabel <- names(stats[stats == x])
            
            if (y == "none") {
                y <- NULL
                yLabel <- "Density"
            } else {
                yLabel <- names(stats[stats == y])
            }
            
            data  <- getInclusionLevels()
            cache <- isolate(getInclusionLevelsSummaryStatsCache())
            data2 <- filterSplicingBasedOnInput()
            
            if (is.null(data2)) {
                errorAlert(
                    session, title="No splicing events returned",
                    "The filtering was too strict and returned no alternative",
                    "splicing events.", dismissible=FALSE,
                    caller="Alternative splicing quantification filtering")
                output$summary <- renderText(NULL)
            } else {
                output$alert   <- renderUI(NULL)
                
                keep  <- nrow(data2)
                total <- nrow(data)
                perc  <- round(keep / total * 100)
                msg   <- sprintf(
                    "Keeping %s alternative splicing events out of %s (%s%%)", 
                    keep, total, perc)
                output$summary <- renderText(msg)
            }
            
            legendLabels <- c("Original", "Filtered in")
            res <- plotRowStats(data, x, y, data2=data2, cache=cache,
                                legend=TRUE, legendLabels=legendLabels,
                                verbose=TRUE) +
                theme_light(14) +
                theme(legend.position="bottom") +
                labs(x=xLabel, y=yLabel)
            
            cache <- attr(res, "cache")
            if (!is.null(cache)) {
                setInclusionLevelsSummaryStatsCache(cache)
            }
            return(res)
        }
    })
}

#' @rdname appServer
inclusionLevelsFilterServer <- function(input, output, session) {
    observeEvent(getInclusionLevels(),
                 updateCheckboxInput(session, "preview", value=FALSE))
    
    psiFilteringOptionsSet(session, input, output)
    filterSplicingSet(session, input, output)
}

attr(inclusionLevelsFilterUI, "loader") <- "data"
attr(inclusionLevelsFilterServer, "loader") <- "data"