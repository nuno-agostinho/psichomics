#' List alternative splicing annotations
#'
#' @param species Character: filter results by species (regular expression)
#' @param assembly Character: filter results by assembly (regular expression)
#' @param date Character: filter results by date (regular expression)
#'
#' @family functions for PSI quantification
#' @return Named character vector with splicing annotation names
#'
#' @importFrom data.table data.table
#' @importFrom R.utils capitalize
#' @export
#'
#' @examples
#' listSplicingAnnotations() # Return all alternative splicing annotations
#' listSplicingAnnotations(assembly="hg19") # Search for hg19 annotation
#' listSplicingAnnotations(assembly="hg38") # Search for hg38 annotation
#' listSplicingAnnotations(date="201(7|8)") # Search for 2017 or 2018 annotation
listSplicingAnnotations <- function(species=NULL, assembly=NULL, date=NULL) {
    df <- data.table(
        species=c("human", "human", "human"),
        assembly=c("hg19/GRCh37", "hg19/GRCh37", "hg38"),
        date=c("2017-10-20", "2016-10-11", "2018-04-30"),
        filename=c("annotationHub_alternativeSplicingEvents.hg19_V2.rda",
                   "annotationHub_alternativeSplicingEvents.hg19.rda",
                   "annotationHub_alternativeSplicingEvents.hg38_V2.rda"))

    filter <- function(data, col, value) {
        if (!is.null(value)) {
            data <- data[grep(value, data[[col]], ignore.case=TRUE)]
        }
        return(data)
    }
    df <- filter(df, "species", species)
    df <- filter(df, "assembly", assembly)
    df <- filter(df, "date", date)

    ns <- sprintf("%s %s (%s)", capitalize(df$species), df$assembly, df$date)
    res <- setNames(df$filename, ns)
    return(res)
}


#' Load alternative splicing annotation from \code{AnnotationHub}
#'
#' @param annotation Character: annotation to load
#' @param cache Character: directory path of cache (if \code{NULL}, default
#' location is \code{AnnotationHub::getAnnotationHubOption("CACHE")})
#'
#' @importFrom BiocFileCache BiocFileCache
#' @importFrom AnnotationHub AnnotationHub query getAnnotationHubOption
#'
#' @family functions for PSI quantification
#' @return List of data frames containing the alternative splicing annotation
#' per event type
#' @export
#'
#' @examples
#' human <- listSplicingAnnotations()[[1]]
#' \dontrun{
#' annot <- loadAnnotation(human)
#' }
loadAnnotation <- function(annotation, cache=NULL) {
    if (is.null(cache)) cache <- getAnnotationHubOption("CACHE")

    if (!dir.exists(cache)) {
        BiocFileCache(cache=cache, ask=FALSE)
        message(sprintf(
            "The directory %s was created to store annotation data", cache))
    }
    ah <- tryCatch(AnnotationHub(cache=cache), error=function(e) {
        msg <- paste("Timeout reached while contacting AnnotationHub.",
                     "Resuming with local cache...")
        message(msg)
        AnnotationHub(cache=cache, localHub=TRUE)
    })
    annot <- gsub("^annotationHub_", "", annotation)
    annot <- ah[[names(query(ah, annot))]]
    return(annot)
}


#' List alternative splicing annotation files available, as well as custom
#' annotation
#'
#' @param ... Custom annotation loaded
#'
#' @return Named character vector with splicing annotation files available
#' @keywords internal
#'
#' @examples
#' psichomics:::listAllAnnotations()
listAllAnnotations <- function(...) {
    list("Available annotation files"=listSplicingAnnotations(),
         "Custom annotation"=c(
             ..., "Load annotation from file..."="loadAnnotation"))
}

#' Interface to quantify alternative splicing
#'
#' @param ns Namespace function
#'
#' @importFrom shiny tagList uiOutput selectizeInput numericInput actionButton
#' @importFrom shinyBS bsTooltip
#' @importFrom shinyjs hidden disabled
#'
#' @return HTML elements
#' @keywords internal
inclusionLevelsInterface <- function(ns) {
    eventTypes <- getSplicingEventTypes()
    names(eventTypes) <- sprintf("%s (%s)", names(eventTypes), eventTypes)

    filterGenesSelectize <- selectizeInput(
        ns("filterGenes"), label=NULL, selected=NULL, multiple=TRUE,
        width="100%", choices=c("Type to search for genes..."=""), options=list(
            create=TRUE, createOnBlur=TRUE, plugins=list("remove_button")))

    sampleFiltering <- bsCollapsePanel(
        tagList(icon("vial"), "Sample filtering",
                contextUI(ns("sampleFilterText"))),
        value="Sample filtering",
        selectizeInput(ns("sampleFilter"), "Samples to discard",
                       multiple=TRUE, width="100%", choices=character(0)))

    psiFiltering <- bsCollapsePanel(
        tagList(icon("sliders-h"), "PSI filtering",
                contextUI(ns("psiFilterText"))),
        value="PSI filtering",
        psiFilteringSetting(ns, "PSI"),
        psiFilteringSetting(ns, "Median"),
        psiFilteringSetting(ns, "LogVar", label="log10(var)",
                            min=-10, max=0, step=0.5),
        psiFilteringSetting(ns, "Range"))

    geneFiltering <- bsCollapsePanel(
        title=tagList(icon("dna"), "Gene filtering",
                      contextUI(ns("geneFilterText"))),
        value="Filter by genes",
        radioButtons(ns("filter"), label=NULL,
                     c("Do not filter splicing events by genes"="noFilter",
                       "Filter by selected genes"="select",
                       "Filter by genes imported from a file"="file")),
        conditionalPanel(
            sprintf("input[id='%s'] == '%s'", ns("filter"), "select"),
            div(id=ns("geneOptionsLoading"), class="progress",
                div(class="progress-bar progress-bar-striped active",
                    role="progressbar", style="width:100%",
                    "Loading genes from annotation")),
            hidden(div(
                id=ns("geneOptions"), filterGenesSelectize,
                div(id=ns("geneLoadingIcon"),
                    style="position: relative;",
                    helpText("Presented genes are based on the selected",
                             "alternative splicing annotation."))))),
        conditionalPanel(
            sprintf("input[id='%s'] == '%s'", ns("filter"), "file"),
            fileBrowserInput(
                ns("filterGenesFile"), label=NULL,
                clearable=TRUE, placeholder="No file selected"),
            helpText("Provide a file containing gene symbols separated",
                     "by space, comma, tab or new line. For instance: ",
                     tags$code("BRCA1, BRAF, ABL"))))

    options <- div(
        id=ns("options"),
        selectizeInput(ns("junctionQuant"), choices=NULL, width = "100%",
                       "Alternative splicing junction quantification"),
        selectizeInput(ns("annotation"), choices=listAllAnnotations(),
                       "Alternative splicing event annotation", width = "100%"),
        selectizeInput(ns("eventType"), "Event types to quantify",
                       selected = c("SE", "MXE", "A5SS", "A3SS", "AFE", "ALE"),
                       choices=eventTypes, multiple = TRUE, width = "100%",
                       options=list(plugins=list("remove_button"))),
        numericInput(ns("minReads"), width = "100%",
                     div("Minimum read counts' threshold",
                         icon("question-circle")), value = 10),
        bsTooltip(ns("minReads"), placement = "right",
                  options = list(container = "body"),
                  paste("Discard alternative splicing quantified using a",
                        "number of reads below this threshold.")),
        bsCollapse(sampleFiltering, psiFiltering, geneFiltering))

    tagList(
        uiOutput(ns("modal")),
        helpText("Exon inclusion levels are measured from exon-exon junction",
                 "quantification using the Percent Spliced-In (PSI) metric."),
        errorDialog("Junction quantification not loaded.",
                    id=ns("missingData"), style="margin: 10px;"),
        hidden(options),
        actionButton(ns("loadIncLevels"), "Load from file..."),
        disabled(processButton(ns("calcIncLevels"),
                               "Quantify alternative splicing")))
}

#' @rdname appUI
inclusionLevelsUI <- function(id, panel) {
    ns <- NS(id)
    title <- "Alternative splicing quantification"
    panel(style="success", title=tagList(icon("calculator"), title),
          value=title, inclusionLevelsInterface(ns))
}

#' Quantify alternative splicing events
#'
#' @param annotation List of data frames: annotation for each alternative
#' splicing event type
#' @param junctionQuant Data frame: junction quantification
#' @param eventType Character: splicing event types to quantify
#' @param minReads Integer: values whose number of total supporting read counts
#' is below \code{minReads} are returned as \code{NA}
#' @param genes Character: gene symbols for which to quantify splicing events
#' (if \code{NULL}, events from all genes are quantified)
#'
#' @importFrom fastmatch %fin%
#'
#' @family functions for PSI quantification
#' @return Data frame with the quantification of the alternative splicing events
#' @export
#'
#' @examples
#' # Calculate PSI for skipped exon (SE) and mutually exclusive (MXE) events
#' annot <- readFile("ex_splicing_annotation.RDS")
#' junctionQuant <- readFile("ex_junctionQuant.RDS")
#'
#' quantifySplicing(annot, junctionQuant, eventType=c("SE", "MXE"))
quantifySplicing <- function(annotation, junctionQuant,
                             eventType=c("SE", "MXE", "ALE", "AFE", "A3SS",
                                         "A5SS"),
                             minReads=10, genes=NULL) {
    if (!is.null(genes)) {
        # Filter for given gene symbols
        filterByGenes <- function(df, genes) {
            # Check which genes are desired by unlisting them all (register the
            # respective event's index for each gene)
            allGenes   <- df$Gene
            valid      <- as.vector(unlist(allGenes)) %fin% genes
            eventGenes <- vapply(allGenes, length, numeric(1), USE.NAMES=FALSE)
            eventIndex <- rep(seq(allGenes), eventGenes)
            return(df[unique(eventIndex[valid]), ])
        }

        genes <- unique(genes)
        annotation <- lapply(annotation, filterByGenes, genes)
    }

    # Convert data frame to matrix if needed (faster)
    mJunctionQuant <- junctionQuant
    if (!is(mJunctionQuant, "matrix"))
        mJunctionQuant <- as.matrix(junctionQuant)

    psi <- NULL
    eventData <- NULL
    for (acronym in eventType) {
        eventTypes <- getSplicingEventTypes()
        type <- names(eventTypes)[[match(acronym, eventTypes)]]
        thisAnnot <- annotation[[type]]
        
        if (!is.null(thisAnnot) && nrow(thisAnnot) > 0) {
            updateProgress("Calculating inclusion levels", type, value=acronym,
                           max=length(eventType))
            incLevels <- calculateInclusionLevels(acronym, mJunctionQuant,
                                                  thisAnnot, minReads)
            eventData <- rbind(eventData, attr(incLevels, "eventData"))
            psi <- rbind(psi, incLevels)
        }
    }

    # Convert matrix to data frame
    colns <- colnames(psi)
    psi   <- data.frame(psi)
    colnames(psi) <- colns

    if (is.null(psi)) psi <- data.frame(NULL)
    psi <- addObjectAttrs(
        psi, rowNames=TRUE,
        description="PSI values per alternative splicing events",
        dataType="Inclusion levels", tablename="Inclusion levels",
        rows="alternative splicing events", columns="samples")
    class(eventData) <- c("eventData", class(eventData))
    attr(psi, "rowData") <- eventData
    psi <- preserveAttributes(psi)
    return(psi)
}

#' Set of functions to load a custom alternative splicing annotation
#'
#' @importFrom shiny tags fileInput
#' @inherit inclusionLevelsServer
#'
#' @keywords internal
loadCustomSplicingAnnotationSet <- function(session, input, output) {
    # Show modal for loading custom splicing annotation
    observe({
        ns <- session$ns
        if (input$annotation == "loadAnnotation") {
            url <- paste0("https://nuno-agostinho.github.io/psichomics/",
                          "articles/AS_events_preparation.html")

            updateSelectizeInput(session, "annotation",
                                 selected=listSplicingAnnotations())
            infoModal(
                session, "Load alternative splicing annotation",
                helpText("To learn how to create and load custom alternative",
                         "splicing annotations,",
                         tags$a(href=url, target="_blank", "click here.")),
                fileInput(ns("customAnnot"), "Choose RDS file", accept=".rds"),
                uiOutput(ns("alert")),
                footer=actionButton(ns("loadCustom"), "Load annotation",
                                    class="btn-primary"))
        }
    })

    # Load custom splicing annotation
    observeEvent(input$loadCustom, {
        customAnnot <- input$customAnnot
        if (is.null(customAnnot)) {
            errorAlert(session, title="No file provided",
                       "Please select a RDS file.",
                       caller="Custom alternative splicing annotation")
        } else if (!grepl("\\.rds$", customAnnot$name, ignore.case=TRUE)) {
            errorAlert(session, title="File format not supported",
                       "Please select a RDS file.",
                       caller="Custom alternative splicing annotation")
        } else {
            custom <- customAnnot$datapath
            names(custom) <- customAnnot$name
            updateSelectizeInput(session, "annotation", selected=custom,
                                 choices=listAllAnnotations(custom))
            removeModal()
            removeAlert(output)
        }
    })
}

#' Set of functions to load splicing quantification
#'
#' @inherit inclusionLevelsServer
#'
#' @importFrom shiny tags
#' @importFrom shinyBS bsPopover
#'
#' @keywords internal
loadSplicingQuantificationSet <- function(session, input, output) {
    ns <- session$ns

    # Show modal for loading alternative splicing quantification
    observeEvent(input$loadIncLevels, {
        infoModal(
            session, "Load alternative splicing quantification",
            ASquantFileInput(ns("customASquant")),
            uiOutput(ns("alertIncLevels")),
            footer=processButton(ns("loadASquant"), "Load quantification"))
    })

    observeEvent(input$loadIncLevels, {
        prepareFileBrowser(session, input, "customASquant")
    }, once=TRUE)

    # Load alternative splicing quantification
    loadSplicing <- reactive({
        time <- startProcess("loadIncLevels")

        startProgress("Wait a moment", divisions=2)
        updateProgress("Loading alternative splicing quantification")

        allFormats <- loadFileFormats()
        formats <- allFormats[sapply(allFormats, "[[",
                                     "dataType") == "Inclusion levels"]

        psi <- tryCatch(loadFile(input$customASquant, formats),
                        warning=return, error=return)
        if (is(psi, "error")) {
            if (psi$message == paste("'file' must be a character string or",
                                     "connection"))
                errorAlert(session, title="No file provided",
                           "Please provide a file.", alertId="alertIncLevels",
                           caller="Alternative splicing quantification")
            else
                errorAlert(session, title="An error was raised",
                           psi$message, alertId="alertIncLevels",
                           caller="Alternative splicing quantification")
        } else if (is(psi, "warning")) {
            warningAlert(session, title="A warning was raised",
                         psi$message, alertId="alertIncLevels",
                         caller="Alternative splicing quantification")
        } else {
            removeAlert(output, "alertIncLevels")

            if ( is.null(getData()) ) {
                name <- file_path_sans_ext( basename(input$customASquant) )
                name <- gsub(" Inclusion levels.*$", "", name)
                if (name == "") name <- "Unnamed"

                data <- setNames(list(list("Inclusion levels"=psi)), name)
                data <- processDatasetNames(data)
                setData(data)
                setCategory(name)

                samples <- colnames(psi)
                parsed <- parseTCGAsampleInfo(samples)
                if ( !is.null(parsed) ) setSampleInfo(parsed)
            } else {
                setInclusionLevels(psi)
            }
            removeModal()
        }
        endProcess("loadIncLevels", time)
    })

    # Show warnings if needed before loading splicing quantification
    observeEvent(input$loadASquant, {
        if (!is.null(getInclusionLevels())) {
            if (!is.null(getDifferentialSplicing())) {
                warningModal(
                    session, "Differential splicing already performed",
                    "Do you wish to replace the loaded alternative splicing",
                    "quantification data and discard the differential splicing",
                    "results?", footer=actionButton(
                        ns("discard2"), "Replace and discard",
                        class="btn-warning", "data-dismiss"="modal"),
                    caller="Alternative splicing quantification")
            } else {
                warningModal(
                    session, "Alternative splicing already quantified",
                    "Do you wish to replace the loaded splicing quantification",
                    "data?", footer=actionButton(ns("replace2"), "Replace",
                                                 class="btn-warning",
                                                 "data-dismiss"="modal"),
                    caller="Alternative splicing quantification")
            }
        } else {
            loadSplicing()
        }
    })

    # Replace previous splicing quantification
    observeEvent(input$replace2, {
        setGroups("Samples", NULL)
        setGroups("AS events", NULL)
        loadSplicing()
    })

    # Discard differential analyses and replace previous splicing quantification
    observeEvent(input$discard2, {
        setDifferentialSplicing(NULL)
        setDifferentialSplicingSurvival(NULL)
        setGroups("Samples", NULL)
        setGroups("AS events", NULL)
        loadSplicing()
    })
}

psiFilteringSet <- function(session, input, output) {
    # Update sample filtering options
    observe({
        junctionQuant <- getJunctionQuantification()[[input$junctionQuant]]
        if (!is.null(junctionQuant)) {
            updateSelectizeInput(
                session, "sampleFilter", server=TRUE,
                choices=colnames(junctionQuant),
                options=list(placeholder="Select samples to discard",
                             plugins=list("remove_button")))
        }
    })

    # Update gene symbols for filtering based on selected annotation
    observe({
        annotation <- input$annotation
        filter <- input$filter

        # Avoid loading if already loaded
        if (filter == "select" && !is.null(annotation) &&
            !annotation %in% c("", "loadAnnotation") &&
            !identical(annotation, getAnnotationName())) {
            # Show loading bar
            show("geneOptionsLoading")
            hide("geneOptions")

            annotation <- input$annotation
            startProgress("Loading alternative splicing annotation",
                          divisions=2)
            annot <- readAnnot(session, annotation, showProgress=TRUE)
            updateProgress("Preparing gene list")
            genes <- sort(unique(unlist(lapply(annot, "[[", "Gene"))))
            updateSelectizeInput(session, "filterGenes", choices=genes,
                                 selected=character(0), server=TRUE)
            closeProgress("Gene list prepared")

            setAnnotationName(annotation)
            # Show gene options set
            hide("geneOptionsLoading")
            show("geneOptions")
        }
    })

    # Toggle visibility of loading icon
    observe({
        toggle("geneLoadingIcon",
               selector = paste0(
                   '$("#data-inclusionLevels-filterGenes").parent()',
                   '.children("div.selectize-control").hasClass("loading")'))
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
            return(isTRUE(res))
        }
        text <- ifelse(arePSIsettingsEnabled(input), "Enabled", "Disabled")
        return(text)
    })
    output$geneFilterText   <- renderText({
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

#' Read custom or remote annotation
#'
#' @inherit inclusionLevelsServer
#' @param annotation Character: chosen annotation
#' @param showProgress Boolean: show progress?
#'
#' @keywords internal
readAnnot <- function(session, annotation, showProgress=FALSE) {
    annot <- NULL
    if (grepl("^/var/folders/", annotation)) { # if custom annotation
        if (showProgress)
            updateProgress("Loading alternative splicing annotation")
        annot <- readRDS(annotation)
    } else if (grepl("^annotationHub_", annotation)) {
        if (showProgress)
            updateProgress("Downloading alternative splicing annotation")
        annot <- loadAnnotation(annotation)

        # Set species and assembly version
        allAnnot <- listSplicingAnnotations()
        annotID <- names(allAnnot)[match(annotation, allAnnot)]
        if (grepl("Human", annotID)) setSpecies("Human")
        if (grepl("hg19", annotID)) setAssemblyVersion("hg19")
    }
    return(annot)
}

#' Set of functions to quantify alternative splicing
#'
#' @importFrom shiny tags
#' @inherit inclusionLevelsServer
#' @keywords internal
quantifySplicingSet <- function(session, input) {
    ns <- session$ns

    # Calculate inclusion levels
    calcSplicing <- reactive({
        eventType  <- input$eventType
        minReads   <- input$minReads
        annotation <- input$annotation

        if (is.null(eventType) || is.null(minReads) || is.null(annotation)) {
            return(NULL)
        } else if (input$junctionQuant == "") {
            errorModal(session, "Select junction quantification",
                       "Select a junction quantification dataset",
                       caller="Alternative splicing quantification")
            endProcess("calcIncLevels")
            return(NULL)
        }

        time <- startProcess("calcIncLevels")
        startProgress("Quantifying alternative splicing", divisions=4)
        # Read annotation
        annot <- readAnnot(session, annotation, showProgress=TRUE)
        junctionQuant <- getJunctionQuantification()[[input$junctionQuant]]

        # Filter alternative splicing events based on cognate genes
        filter <- getCognateGenesToFilter(input$filter, input$filterGenes,
                                          input$filterGenesFile)
        # Discard samples
        sampleFilter <- input$sampleFilter
        if (!is.null(sampleFilter) && sampleFilter != "") {
            samplesToKeep <- !colnames(junctionQuant) %in% sampleFilter
            junctionQuant <- junctionQuant[ , samplesToKeep]
            sampleFilterText <- paste(sampleFilter, collapse=", ")
        } else {
            sampleFilterText <- "None"
        }
        sampleFilterSettings <- c("Discarded samples"=sampleFilterText)

        # Quantify splicing with splicing annotation and junction quantification
        updateProgress("Calculating inclusion levels")
        psi <- quantifySplicing(annot, junctionQuant, eventType, minReads,
                                genes=filter)

        if (nrow(psi) == 0) {
            errorModal(session, "No splicing events returned",
                       "The total reads of the alternative splicing events are",
                       "below the minium read counts' threshold.",
                       caller="Alternative splicing quantification")
            endProcess("calcIncLevels")
            return(NULL)
        }

        # Filter PSI values
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

        filtered <- filterPSI(psi, minPSI=minPSI, maxPSI=maxPSI,
                              minMedian=minMedian, maxMedian=maxMedian,
                              minLogVar=minLogVar, maxLogVar=maxLogVar,
                              minRange=minRange, maxRange=maxRange)
        if (length(filtered) != nrow(psi)) {
            filteredPSI    <- psi[filtered, ]
            filterSettings <- attr(filtered, "filtered")
        } else {
            filteredPSI    <- psi
            filterSettings <- NULL
        }

        # Include settings used for alternative splicing quantification
        allEventTypes <- getSplicingEventTypes()
        eventTypeName <- names(allEventTypes[allEventTypes %in% eventType])
        settings <- c(
            list(
                "Alternative splicing annotation"=annotation,
                "Exon-exon junction quantification (label)"=input$junctionQuant,
                "Exon-exon junction quantification (file)"=attr(junctionQuant,
                                                                "filename"),
                "Splicing event types"=eventTypeName,
                "Minimum read counts' threshold"=minReads,
                "Selected cognate genes"=if (is.null(
                    filter)) "All available genes" else filter),
            sampleFilterSettings, filterSettings)
        attr(filteredPSI, "settings") <- settings
        attr(filteredPSI, "icon") <- list(symbol="calculator", colour="green")

        setInclusionLevels(filteredPSI)
        endProcess("calcIncLevels", time)
    })

    # Show warnings if needed before quantifying alternative splicing
    observeEvent(input$calcIncLevels, {
        if (is.null(getData()) || is.null(getJunctionQuantification())) {
            missingDataModal(session, "Junction quantification", ns("missing"))
        } else if (!is.null(getInclusionLevels())) {
            if (!is.null(getDifferentialSplicing())) {
                warningModal(
                    session, "Differential splicing already performed",
                    "Do you wish to replace the current alternative splicing",
                    "quantification data and discard the differential splicing",
                    "results?", footer=actionButton(
                        ns("discard2"), "Replace and discard",
                        class="btn-warning", "data-dismiss"="modal"),
                    caller="Alternative splicing quantification")
            } else {
                warningModal(
                    session, "Alternative splicing already quantified",
                    "Do you wish to replace the loaded splicing quantification",
                    "data?", footer=actionButton(ns("replace"), "Replace",
                                                 class="btn-warning",
                                                 "data-dismiss"="modal"),
                    caller="Alternative splicing quantification")
            }
        } else {
            calcSplicing()
        }
    })

    observeEvent(input$replace, calcSplicing())
    observeEvent(input$discard, {
        setDifferentialSplicing(NULL)
        setDifferentialSplicingSurvival(NULL)
        calcSplicing()
    })
}

#' @rdname appServer
#'
#' @importFrom shiny reactive observeEvent helpText removeModal
#' @importFrom tools file_path_sans_ext
#' @importFrom shinyjs enable disable hide show
#' @importFrom data.table fread
inclusionLevelsServer <- function(input, output, session) {
    ns <- session$ns
    observeEvent(input$missing, missingDataGuide("Junction quantification"))

    prepareFileBrowser(session, input, "filterGenesFile")

    # Update available junction quantification according to loaded files
    observe({
        junctionQuant <- getJunctionQuantification()
        if (!is.null(junctionQuant)) {
            updateSelectizeInput(session, "junctionQuant",
                                 choices=c(names(junctionQuant),
                                           "Select junction quantification"=""))
        } else {
            updateSelectizeInput(
                session, "junctionQuant",
                choices=c("No junction quantification loaded"=""))
        }
    })

    # Warn user if junction quantification is not loaded
    observe({
        if (is.null(getData()) || is.null(getJunctionQuantification())) {
            hide("options")
            disable("calcIncLevels")
            show("missingData")
        } else {
            show("options")
            enable("calcIncLevels")
            hide("missingData")
        }
    })

    # Update default AS event annotation based on selected dataset
    observe({
        data <- getCategoryData()
        isRecountData <- !is.null(data) && !is.null(attr(data, "source")) &&
            attr(data, "source") == "recount"
        if (isRecountData) {
            selected <- grep("hg38", listSplicingAnnotations(), value=TRUE)[[1]]
        } else {
            selected <- grep("hg19", listSplicingAnnotations(), value=TRUE)[[1]]
        }
        updateSelectizeInput(session, "annotation", selected=selected)
    })

    quantifySplicingSet(session, input)
    loadCustomSplicingAnnotationSet(session, input, output)
    loadSplicingQuantificationSet(session, input, output)
    psiFilteringSet(session, input, output)
}

attr(inclusionLevelsUI, "loader") <- "data"
attr(inclusionLevelsServer, "loader") <- "data"
