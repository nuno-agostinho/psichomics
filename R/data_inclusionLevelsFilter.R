#' Filter alternative splicing quantification
#'
#' @param psi Data frame or matrix: alternative splicing quantification
#' @param eventType Character: filter data based on event type; check all event
#' types available by using \code{getSplicingEventTypes(psi)}, where \code{psi}
#' is the alternative splicing quantification data; if \code{eventType = NULL},
#' events are not filtered by event type
#' @param eventSubtype Character: filter data based on event subtype; check all
#' event subtypes available in your data by using
#' \code{unique(getSplicingEventData(psi)$subtype)}, where \code{psi} is the
#' alternative splicing quantification data; if \code{eventSubtype = NULL},
#' events are not filtered by event subtype
#' @param minPSI Numeric: minimum PSI value
#' @param maxPSI Numeric: maximum PSI value
#' @param minMedian Numeric: minimum median PSI per splicing event
#' @param maxMedian Numeric: maximum median PSI per splicing event
#' @param minLogVar Numeric: minimum log10(PSI variance) per splicing event
#' @param maxLogVar Numeric: maximum log10(PSI variance) per splicing event
#' @param minRange Numeric: minimum PSI range across samples per splicing event
#' @param maxRange Numeric: maximum PSI range across samples per splicing event
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
                      minPSI=-Inf, maxPSI=Inf,
                      minMedian=-Inf, maxMedian=Inf,
                      minLogVar=-Inf, maxLogVar=Inf,
                      minRange=-Inf, maxRange=Inf) {
    if (is.na(minPSI)) minPSI <- -Inf
    if (is.na(maxPSI)) maxPSI <- Inf
    if (is.na(minMedian)) minMedian <- -Inf
    if (is.na(maxMedian)) maxMedian <- Inf
    if (is.na(minLogVar)) minLogVar <- -Inf
    if (is.na(maxLogVar)) maxLogVar <- Inf
    if (is.na(minRange)) minRange <- -Inf
    if (is.na(maxRange)) maxRange <- Inf

    trueVector <- function(len) rep(TRUE, len)

    type <- getSplicingEventData(psi)$type
    if (!is.null(eventType) && !is.null(type)) {
        type <- type %in% eventType
        eventTypeText <- c("Splicing event types"=paste(eventType,
                                                        collapse = ", "))
    } else {
        type <- trueVector(nrow(psi))
        eventTypeText <- c("Splicing event types"="All")
    }

    subtype <- getSplicingEventData(psi)$subtype
    if (!is.null(eventSubtype) && !is.null(subtype)) {
        parsed  <- parseSplicingEvent(rownames(psi), data=psi,
                                      pretty=TRUE)$subtype
        subtype1 <- subtype %in% eventSubtype
        subtype2 <- parsed %in% eventSubtype
        if (length(subtype1) != length(subtype2)) {
            warning("Unknown warning: subtype selection went wrong...")
            subtype          <- trueVector(nrow(psi))
            eventSubtypeText <- c("Splicing event subtypes"="All")
        } else {
            subtype  <- subtype1 | subtype2
            eventSubtypeText <- c("Splicing event subtypes"=paste(
                eventSubtype, collapse = ", "))
        }
    } else {
        subtype          <- trueVector(nrow(psi))
        eventSubtypeText <- c("Splicing event subtypes"="All")
    }

    if (!is.infinite(minMedian) || !is.infinite(maxMedian)) {
        medians         <- customRowMedians(psi, na.rm=TRUE, fast=TRUE)
        medianThres     <- medians >= minMedian & medians <= maxMedian
        medianThresText <- c("Median PSI >="=minMedian,
                             "Median PSI <="=maxMedian)
    } else {
        medianThres     <- trueVector(nrow(psi))
        medianThresText <- NULL
    }

    if (!is.infinite(minLogVar) || !is.infinite(maxLogVar)) {
        vars         <- log10(customRowVars(psi, na.rm=TRUE, fast=TRUE))
        varThres     <- vars >= minLogVar & vars <= maxLogVar
        varThresText <- c("log10(PSI variance) >="=minLogVar,
                          "log10(PSI variance) <="=maxLogVar)
    } else {
        varThres     <- trueVector(nrow(psi))
        varThresText <- NULL
    }

    if (!is.infinite(minPSI) || !is.infinite(maxPSI)) {
        mini         <- customRowMins(psi, na.rm=TRUE, fast=TRUE)
        maxi         <- customRowMaxs(psi, na.rm=TRUE, fast=TRUE)
        psiThres     <- mini >= minPSI & maxi <= maxPSI
        psiThresText <- c("PSI >="=minPSI, "PSI <="=maxPSI)
    } else {
        mini         <- maxi <- NULL
        psiThres     <- trueVector(nrow(psi))
        psiThresText <- NULL
    }

    if (!is.infinite(minRange) || !is.infinite(maxRange)) {
        if (!is.null(mini) && !is.null(maxi)) {
            ranges <- maxi - mini
        } else {
            ranges <- customRowRanges(psi, na.rm=TRUE, fast=TRUE)
        }
        rangeThres     <- ranges >= minRange & ranges <= maxRange
        rangeThresText <- c("PSI range >="=minRange, "PSI range <="=maxRange)
    } else {
        rangeThres     <- trueVector(nrow(psi))
        rangeThresText <- NULL
    }

    thres <- which(type & subtype &
                       psiThres & medianThres & varThres & rangeThres)

    attr(thres, "filtered") <- c(eventTypeText, eventSubtypeText, psiThresText,
                                 medianThresText, varThresText, rangeThresText)
    return(thres)
}

#' Remove alternative splicing quantification values based on coverage
#'
#' @param psi Data frame or matrix: alternative splicing quantification
#' @param minReads Currently this argument does nothing
#' @param vasttoolsScoresToDiscard Character: if you are using inclusion levels
#' from VAST-TOOLS, filter the data based on quality scores for read coverage,
#' e.g. use \code{vasttoolsScoresToDiscard = c("SOK", "OK", "LOW")} to only keep
#' events with good read coverage (by default, events are not filtered based on
#' quality scores); read \url{https://github.com/vastgroup/vast-tools} for more
#' information on VAST-TOOLS quality scores
#'
#' @return Alternative splicing quantification data with missing values for any
#' values with insufficient coverage
#' @export
discardLowCoveragePSIvalues <- function(
    psi, minReads=10, vasttoolsScoresToDiscard=c("VLOW", "N")) {

    eventData <- getSplicingEventData(psi)
    if (is.null(eventData)) return(psi)

    # Remove events containing only missing values
    removeNAonlyEvents <- function(psi) psi[rowSums(!is.na(psi)) > 0, ]

    isVastTools <- isTRUE(unique(eventData$source) == "vast-tools")
    if (!is.null(vasttoolsScoresToDiscard) && isVastTools) {
        qualityCol <- max(which(!endsWith(colnames(eventData), "-Q")))
        psi <- discardVastToolsByCvg(psi, eventData, qualityCol,
                                     vasttoolsScoresToDiscard)
        psi <- removeNAonlyEvents(psi)
        attr(psi, "filtered") <- c(
            "VAST-TOOLS coverage"=paste(vasttoolsScoresToDiscard,
                                        collapse=", "))
    }
    return(psi)
}

getPSIsummaryStats <- function() {
    c("Mean PSI"="mean", "Median PSI"="median",
      "PSI variance"="var", "log10(PSI variance)"="log10(var)",
      "PSI range"="range")
}

numericInputWithCheckbox <- function(ns, id, label, ..., check=FALSE) {
    checkbox <- checkboxInput(ns(paste0("enable", capitalize(id))), label,
                              value=check, width="100%")
    # Adjust margins
    checkbox[[2]][["style"]] <- paste(checkbox[[2]][["style"]],
                                      "margin-bottom: 5px;")
    checkbox[[3]][[1]][[2]]$style <- "margin-top: 0; margin-bottom: 0;"
    # Change label to bold
    checkbox[[3]][[1]][[3]][[1]][[3]][[2]][[2]]$style <- "font-weight: bold"

    checkbox <- span(id=ns(paste0("element", capitalize(id))), checkbox)
    column(6, checkbox, numericInput(ns(id), label=NULL, ..., width="100%"))
}

psiFilteringSetting <- function(ns, id, min=0, max=1, step=0.1, ..., label=id,
                                check=FALSE) {
    minValue <- paste0("min", id)
    maxValue <- paste0("max", id)
    minLabel <- paste(label, ">=")
    maxLabel <- paste(label, "<=")

    if (length(check) == 1) {
        minCheck <- maxCheck <- check
    } else if (length(check) == 2) {
        minCheck <- check[[1]]
        maxCheck <- check[[2]]
    }
    fluidRow(
        numericInputWithCheckbox(ns, minValue, minLabel, ..., check=minCheck,
                                 min=min, max=max, value=min, step=step),
        numericInputWithCheckbox(ns, maxValue, maxLabel, ..., check=maxCheck,
                                 min=min, max=max, value=max, step=step))
}

#' Interface to filter alternative splicing
#'
#' @param ns Namespace function
#'
#' @importFrom shiny tagList uiOutput selectizeInput numericInput actionButton
#' @importFrom shinyBS bsTooltip
#' @importFrom shinyjs hidden disabled
#' @importFrom R.utils capitalize
#'
#' @return HTML elements
#' @keywords internal
inclusionLevelsFilterInterface <- function(ns) {
    filterGenesSelectize <- selectizeInput(
        ns("filterGenes"), label=NULL, selected=NULL, multiple=TRUE,
        width="100%", choices=c("Type to search for genes..."=""),
        options=list(plugins=list("remove_button")))

    sampleFiltering <- bsCollapsePanel(
        tagList(icon("vial"), "Sample filtering",
                contextUI(ns("sampleFilterText"))),
        value="Sample filtering",
        selectizeInput(ns("sampleFilter"), "Samples to discard", multiple=TRUE,
                       width="100%", choices=character(0)))

    psiFiltering <- bsCollapsePanel(
        tagList(icon("sliders-h"), "PSI filtering",
                contextUI(ns("psiFilterText"))),
        value="PSI filtering",
        psiFilteringSetting(ns, "PSI"),
        psiFilteringSetting(ns, "Median"),
        psiFilteringSetting(ns, "LogVar", label="log10(var)",
                            min=-10, max=0, step=0.5),
        psiFilteringSetting(ns, "Range"))

    unparsableError <- function(id) {
        errorDialog(
            paste("Loaded splicing events could not be parsed. Data cannot be",
                  "filtered based on event type or cognate genes."),
            id=ns(id), style="margin: 10px;")
    }
    geneFiltering <- bsCollapsePanel(
        title=tagList(icon("dna"), "Gene filtering",
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
        hidden(selectizeInput(ns("vasttoolsScoresToDiscard"), width="100%",
                              "VAST-TOOLS: quality scores to discard",
                              choices=c("SOK", "OK", "LOW", "VLOW", "N"),
                              selected=c("VLOW", "N"), multiple=TRUE)),
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

#' @rdname appServer
#' @importFrom shiny updateCheckboxInput
inclusionLevelsFilterServer <- function(input, output, session) {
    observeEvent(getInclusionLevels(), priority=1,
                 updateCheckboxInput(session, "preview", value=FALSE))

    observeEvent(input$missing, missingDataGuide("Inclusion Levels"))
    prepareFileBrowser(session, input, "filterGenesFile")

    ns <- session$ns

    # Process PSI
    processPSI <- reactive({
        psi <- getInclusionLevels()
        if (!is.null(psi)) {
            parsed <- parseSplicingEvent(rownames(psi), data=psi, pretty=TRUE)
            return(parsed)
        }
    })

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

    # Toggle VAST-TOOLS-specific coverage options
    observe({
        psi       <- getInclusionLevels()
        eventData <- getSplicingEventData(psi)
        if (!is.null(eventData) &&
            isTRUE(unique(eventData$source) == "vast-tools")) {
            show("vasttoolsScoresToDiscard")
        } else {
            hide("vasttoolsScoresToDiscard")
        }
    })

    # Update event types
    observe({
        processed <- processPSI()
        if (!is.null(processed)) {
            types <- table(processed$subtype)
            label <- names(types)
            count <- paste(as.vector(types), "AS events")

            types <- setNames(label, paste(label, count, sep=" __ "))
            updateSelectizeInput(
                session, "eventSubtype", choices=types, selected=types,
                options=list(
                    highlight=FALSE, plugins=list("remove_button"),
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
    togglePSIsetting <- function(id) {
        checkbox <- paste0("enable", capitalize(id))
        observe(toggleState(id, input[[checkbox]]))
    }
    togglePSIsetting("minPSI")
    togglePSIsetting("maxPSI")
    togglePSIsetting("minMedian")
    togglePSIsetting("maxMedian")
    togglePSIsetting("minLogVar")
    togglePSIsetting("maxLogVar")
    togglePSIsetting("minRange")
    togglePSIsetting("maxRange")

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
        arePSIsettingsEnabled <- function(input) {
            elems    <- names(input)
            enablers <- elems[startsWith(elems, "enableM")]
            res      <- any(sapply(enablers, function(x) input[[x]]))
            return(res)
        }
        text <- ifelse(arePSIsettingsEnabled(input), "Enabled", "Disabled")
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

    # Discard low coverage
    discardLowCoverage <- function(psi, vasttoolsScoresToDiscard) {
        discardLowCvgReactive <- reactive({
            filtered <- discardLowCoveragePSIvalues(
                psi, vasttoolsScoresToDiscard)
            return(filtered)
        })
        discardLowCvgReactive()
    }

    filterSplicingOperation <- function(session, psi, processed, input) {
        eventSubtype <- input$eventSubtype
        sampleFilter <- input$sampleFilter
        vasttoolsScoresToDiscard <- input$vasttoolsScoresToDiscard

        checkPSIsetting <- function(id) {
            if (startsWith(id, "min")) {
                res <- -Inf
            } else if (startsWith(id, "max")) {
                res <- Inf
            }

            enabled <- input[[paste0("enable", capitalize(id))]]
            if (isTRUE(enabled)) res <- input[[id]]
            return(res)
        }
        minPSI    <- checkPSIsetting("minPSI")
        maxPSI    <- checkPSIsetting("maxPSI")
        minMedian <- checkPSIsetting("minMedian")
        maxMedian <- checkPSIsetting("maxMedian")
        minLogVar <- checkPSIsetting("minLogVar")
        maxLogVar <- checkPSIsetting("maxLogVar")
        minRange  <- checkPSIsetting("minRange")
        maxRange  <- checkPSIsetting("maxRange")

        filterEnabled   <- input$filter
        filterGenes     <- input$filterGenes
        filterGenesFile <- input$filterGenesFile

        areThereInclusionLevels <- function(psi) {
            return(!is.null(psi) && nrow(psi) != 0 && ncol(psi) != 0)
        }
        if (!areThereInclusionLevels(psi)) return(NULL)
        filteredPSI <- psi

        suppressWarnings(updateProgress("Filtering based on genes"))
        if (!is.null(processed$gene)) {
            filter <- getCognateGenesToFilter(filterEnabled, filterGenes,
                                              filterGenesFile)
            if (!is.null(filter)) {
                filteredPSI <- filteredPSI[processed$gene %in% filter, ]
            }
        } else {
            filter <- NULL
        }
        if (!areThereInclusionLevels(filteredPSI)) return(NULL)
        geneFilterSettings <- list("Selected cognate genes"=if (
            is.null(filter)) "All available genes" else filter)

        suppressWarnings(updateProgress("Discarding samples"))
        if (!is.null(sampleFilter) && sampleFilter != "") {
            samplesToKeep    <- !colnames(filteredPSI) %in% sampleFilter
            filteredPSI      <- filteredPSI[ , samplesToKeep]
            sampleFilterText <- paste(sampleFilter, collapse=", ")
        } else {
            sampleFilterText <- "None"
        }
        sampleFilterSettings <- c("Discarded samples"=sampleFilterText)
        if (!areThereInclusionLevels(filteredPSI)) return(NULL)

        suppressWarnings(updateProgress("Filtering based on event types"))
        if ((length(eventSubtype) == 1 && eventSubtype == "") ||
            is.null(eventSubtype)) {
            subtypeSettings <- NULL
        } else {
            filtered        <- filterPSI(filteredPSI, eventSubtype=eventSubtype)
            filteredPSI     <- filteredPSI[filtered, ]
            subtypeSettings <- attr(filtered, "filtered")[-c(1)]
        }
        if (!areThereInclusionLevels(filteredPSI)) return(NULL)

        eventData <- getSplicingEventData(filteredPSI)
        isVastTools <- isTRUE(unique(eventData$source) == "vast-tools")
        if ( !is.null(vasttoolsScoresToDiscard) && isVastTools ) {
            filteredPSI <- discardLowCoverage(filteredPSI,
                                              vasttoolsScoresToDiscard)
            vasttoolsFilterSettings <- attr(filteredPSI, "filtered")
            if (!areThereInclusionLevels(filteredPSI)) return(NULL)
        } else {
            vasttoolsFilterSettings <- NULL
        }

        suppressWarnings(updateProgress("Filtering based on PSI values"))
        filtered <- filterPSI(filteredPSI, minPSI=minPSI, maxPSI=maxPSI,
                              minMedian=minMedian, maxMedian=maxMedian,
                              minLogVar=minLogVar, maxLogVar=maxLogVar,
                              minRange=minRange, maxRange=maxRange)
        if (length(filtered) != nrow(filteredPSI)) {
            filteredPSI    <- filteredPSI[filtered, ]
            filterSettings <- attr(filtered, "filtered")[-c(1, 2)]
        } else {
            filterSettings <- NULL
        }
        if (!areThereInclusionLevels(filteredPSI)) return(NULL)

        # Include settings used for alternative splicing quantification
        settings <- c(geneFilterSettings, sampleFilterSettings, subtypeSettings,
                      vasttoolsFilterSettings, filterSettings)
        attr(filteredPSI, "settings") <- settings
        attr(filteredPSI, "icon") <- list(symbol="calculator", colour="green")
        suppressWarnings(closeProgress())
        return(filteredPSI)
    }

    # Operation to filter splicing based on user-defined options
    filterSplicingBasedOnInput <- reactive({
        psi       <- getInclusionLevels()
        processed <- processPSI()
        filterSplicingOperation(session, psi, processed, input)
    })

    # Filter inclusion levels
    filterSplicing <- function() {
        time <- startProcess("filterIncLevels")
        startProgress("Filtering alternative splicing", divisions=5)
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

attr(inclusionLevelsFilterUI, "loader") <- "data"
attr(inclusionLevelsFilterServer, "loader") <- "data"
