#' List the alternative splicing annotation files available
#' 
#' @return Named character vector with splicing annotation files available
#' @export
#' 
#' @examples
#' listSplicingAnnotations()
listSplicingAnnotations <- function() {
    c("Human hg19/GRCh37 (2017-10-20)"=
          "annotationHub_alternativeSplicingEvents.hg19_V2.rda",
      "Human hg19/GRCh37 (2016-10-11)"=
          "annotationHub_alternativeSplicingEvents.hg19.rda",
      "Human hg38 (2017-10-20)"=
          "annotationHub_alternativeSplicingEvents.hg38.rda")
}

#' List alternative splicing annotation files available, as well as custom 
#' annotation
#' 
#' @param ... Custom annotation loaded
#' 
#' @return Named character vector with splicing annotation files available#' 
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
inclusionLevelsInterface <- function(ns) {
    eventTypes <- getSplicingEventTypes()
    names(eventTypes) <- sprintf("%s (%s)", names(eventTypes), eventTypes)
    
    filterGenesSelectize <- selectizeInput(
        ns("filterGenes"), label=NULL, selected=NULL, multiple=TRUE,
        width="100%", choices=c("Type to search for genes..."=""), options=list(
            # Allow to add new items
            create=TRUE, createOnBlur=TRUE, plugins=list("remove_button")))
    
    options <- div(
        id=ns("options"),
        selectizeInput(ns("junctionQuant"), choices=NULL, width = "100%",
                       "Alternative splicing junction quantification"),
        selectizeInput(ns("annotation"), choices=listAllAnnotations(),
                       "Alternative splicing event annotation", width = "100%"),
        selectizeInput(ns("eventType"), "Event type(s)", selected = "SE",
                       choices=eventTypes, multiple = TRUE, width = "100%",
                       options=list(plugins=list("remove_button"))),
        numericInput(ns("minReads"), width = "100%", 
                     div("Minimum read counts' threshold",
                         icon("question-circle")), value = 10),
        bsCollapse(
            bsCollapsePanel(
                title=tagList(icon("filter"),
                              "Filter splicing events by genes"), 
                value="Filter by genes",
                radioButtons(
                    ns("filter"), NULL,
                    c("Do not filter splicing events"="noFilter",
                      "Filter by selected genes"="select",
                      "Filter by genes imported from a file"="file")),
                conditionalPanel(
                    sprintf("input[id='%s'] == '%s'", ns("filter"), "select"),
                    div(id=ns("geneOptionsLoading"), class="progress",
                        div(class="progress-bar progress-bar-striped active",
                            role="progressbar", style="width:100%",
                            "Loading genes from annotation")),
                    hidden(div(id=ns("geneOptions"),
                        filterGenesSelectize,
                        div(id=ns("geneLoadingIcon"),
                            style="position: relative;",
                            # div(class="fa fa-spinner fa-spin",
                            #     style="position:absolute;", 
                            #     style="right:6px;",
                            #     style="bottom: 24px;", style="z-index: 2;")),
                        helpText("Presented genes are based on the selected",
                                 "alternative splicing annotation."))))),
                conditionalPanel(
                    sprintf("input[id='%s'] == '%s'", ns("filter"), "file"),
                    fileBrowserInput(ns("filterGenesFile"), NULL,
                                     placeholder="No file selected"),
                    helpText("Provide a file with gene symbols separated by a",
                             "space, comma, tab or new line. For instance: ",
                             tags$code("BRCA1, BRAF, ABL"))))),
        bsTooltip(ns("minReads"), placement = "right", 
                  options = list(container = "body"),
                  paste("Discard alternative splicing quantified using a",
                        "number of reads below this threshold.")))
    
    tagList(
        uiOutput(ns("modal")),
        helpText("Exon inclusion levels are measured from exon-exon junction",
                 "quantification using the Percent Spliced-In (PSI) metric."),
        errorDialog("No junction quantification is loaded.",
                    id=ns("missingData"), style="margin: 10px;"),
        hidden(options),
        actionButton(ns("loadIncLevels"), "Load from file"),
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
#' @param minReads Integer: discard alternative splicing quantified using a
#' number of reads below this threshold
#' @param genes Character: gene symbols for which the splicing quantification 
#' of associated splicing events is performed (by default, all splicing events 
#' undergo splicing quantification)
#' 
#' @return Data frame with the quantification of the alternative splicing events
#' @export
#' 
#' @importFrom fastmatch %fin%
#' 
#' @examples 
#' # Calculate PSI for skipped exon (SE) and mutually exclusive (MXE) events
#' annot <- readFile("ex_splicing_annotation.RDS")
#' junctionQuant <- readFile("ex_junctionQuant.RDS")
#' 
#' psi <- quantifySplicing(annot, junctionQuant, eventType=c("SE", "MXE"))
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
    for (acronym in eventType) {
        eventTypes <- getSplicingEventTypes()
        type <- names(eventTypes)[[match(acronym, eventTypes)]]
        thisAnnot <- annotation[[type]]
        updateProgress("Calculating inclusion levels", type, value=acronym, 
                       max=length(eventType))
        
        if (!is.null(thisAnnot) && nrow(thisAnnot) > 0) {
            psi <- rbind(psi, calculateInclusionLevels(
                acronym, mJunctionQuant, thisAnnot, minReads))
        }
    }
    
    # Convert matrix to data frame
    colns <- colnames(psi)
    psi   <- data.frame(psi)
    colnames(psi) <- colns
    
    if (is.null(psi)) psi <- data.frame(NULL)
    attr(psi, "rowNames")    <- TRUE
    attr(psi, "description") <- "PSI values per alternative splicing events"
    attr(psi, "dataType")    <- "Inclusion levels"
    attr(psi, "tablename")   <- "Inclusion levels"
    attr(psi, "rows")        <- "alternative splicing events"
    attr(psi, "columns")     <- "samples"
    return(psi)
}

#' Load alternative splicing annotation from \code{AnnotationHub}
#' 
#' @param annotation Character: annotation to load
#' 
#' @importFrom AnnotationHub AnnotationHub query
#' 
#' @return List of data frames containing the alternative splicing annotation
#' per event type
#' @export
#' @examples 
#' human <- listSplicingAnnotations()[[1]]
#' \dontrun{
#' annot <- loadAnnotation(human)
#' }
loadAnnotation <- function(annotation) {
    ah <- AnnotationHub()
    annot <- gsub("^annotationHub_", "", annotation)
    annot <- ah[[names(query(ah, annot))]]
    return(annot)
}

#' Set of functions to load a custom alternative splicing annotation
#' 
#' @importFrom shiny tags fileInput
#' @inherit inclusionLevelsServer
loadCustomSplicingAnnotationSet <- function(session, input, output) {
    # Show modal for loading custom splicing annotation
    observe({
        ns <- session$ns
        if (input$annotation == "loadAnnotation") {
            url <- "http://rpubs.com/nuno-agostinho/alt-splicing-annotation"
            
            updateSelectizeInput(session, "annotation", 
                                 selected=listSplicingAnnotations())
            infoModal(session, "Load alternative splicing annotation",
                      helpText("Load alternative splicing annotation from a",
                               "RDS file. To learn more on how to create a",
                               "custom splicing annotation,", 
                               tags$a(href=url, target="_blank",
                                      "click here.")),
                      fileInput(ns("customAnnot"), "Choose RDS file",
                                accept=".rds"),
                      selectizeInput(ns("customSpecies"), "Species", 
                                     choices="Human", 
                                     options=list(create=TRUE)),
                      selectizeInput(ns("customAssembly"), "Assembly",
                                     choices="hg19", options=list(create=TRUE)),
                      uiOutput(ns("alert")),
                      footer=actionButton(ns("loadCustom"), "Load annotation",
                                          class="btn-primary"))
        }
    })
    
    # Load custom splicing annotation
    observeEvent(input$loadCustom, {
        customAnnot <- input$customAnnot
        if (is.null(customAnnot)) {
            errorAlert(session, title="No file provided.",
                       "Please select a RDS file.")
        } else if (!grepl("\\.rds$", customAnnot$name, ignore.case=TRUE)) {
            errorAlert(session, title="File format not allowed.",
                       "Please select a RDS file.")
        } else {
            custom <- customAnnot$datapath
            names(custom) <- customAnnot$name
            updateSelectizeInput(session, "annotation", selected=custom,
                                 choices=listAllAnnotations(custom))
            setSpecies(input$customSpecies)
            setAssemblyVersion(input$customAssembly)
            removeModal()
            removeAlert(output)
        }
    })
}

#' Set of functions to load splicing quantification
#' 
#' @importFrom shiny tags
#' @importFrom shinyBS bsPopover
#' @inherit inclusionLevelsServer
loadSplicingQuantificationSet <- function(session, input, output) {
    ns <- session$ns
    
    # Show modal for loading alternative splicing quantification
    observeEvent(input$loadIncLevels, {
        infoModal(
            session, "Load alternative splicing quantification",
            ASquantFileInput(ns("customASquant"), ns("customSpecies2"),
                             ns("customAssembly2")),
            uiOutput(ns("alertIncLevels")),
            footer=processButton(ns("loadASquant"), "Load quantification"))
    })
    
    observeEvent(input$loadIncLevels, {
        prepareFileBrowser(session, input, "customASquant")
    }, once=TRUE)
    
    # Show warnings if needed before loading splicing quantification
    observeEvent(input$loadASquant, {
        if (!is.null(getInclusionLevels())) {
            if (!is.null(getDifferentialAnalyses())) {
                warningModal(session, "Warning",
                             "The calculated differential splicing analyses",
                             "will be discarded and the previously loaded",
                             "inclusion levels will be replaced.",
                             footer=actionButton(ns("discard2"),
                                                 "Discard and replace",
                                                 class="btn-warning",
                                                 "data-dismiss"="modal"))
            } else {
                warningModal(session, "Inclusion levels already quantified",
                             "Do you wish to replace the inclusion levels",
                             "loaded?",
                             footer=actionButton(ns("replace2"), "Replace",
                                                 class="btn-warning",
                                                 "data-dismiss"="modal"))
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
        setDifferentialAnalyses(NULL)
        setDifferentialAnalysesSurvival(NULL)
        setGroups("Samples", NULL)
        setGroups("AS events", NULL)
        loadSplicing()
    })
    
    # Load alternative splicing quantification
    loadSplicing <- reactive({
        time <- startProcess("loadIncLevels")
        
        startProgress("Wait a moment", divisions=2)
        updateProgress("Loading alternative splicing quantification")
        
        allFormats <- loadFileFormats()
        formats <- allFormats[sapply(allFormats, "[[", 
                                     "dataType") == "Inclusion levels"]
        
        psi <- tryCatch(parseValidFile(input$customASquant, formats),
                        warning=return, error=return)
        if (is(psi, "error")) {
            if (psi$message == paste("'file' must be a character string or",
                                     "connection"))
                errorAlert(session, title="Error", "No file was provided",
                           alertId="alertIncLevels")
            else
                errorAlert(session, title="Error", 
                           psi$message, alertId="alertIncLevels")
        } else if (is(psi, "warning")) {
            warningAlert(session, title="Warning", 
                         psi$message, alertId="alertIncLevels")
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
                parsed <- parseTcgaSampleInfo(samples) 
                if ( !is.null(parsed) ) setSampleInfo(parsed)
            } else {
                setInclusionLevels(psi)
            }
            setSpecies(input$customSpecies2)
            setAssemblyVersion(input$customAssembly2)
            removeModal()
        }
        endProcess("loadIncLevels", time)
    })
}

#' Read custom or remote annotation
#' @inherit inclusionLevelsServer
#' @param annotation Character: chosen annotation
#' @param showProgress Boolean: show progress? FALSE by default
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
quantifySplicingSet <- function(session, input) {
    ns <- session$ns
    
    # Show warnings if needed before calculating inclusion levels
    observeEvent(input$calcIncLevels, {
        if (is.null(getData()) || is.null(getJunctionQuantification())) {
            missingDataModal(session, "Junction quantification", ns("missing"))
        } else if (!is.null(getInclusionLevels())) {
            if (!is.null(getDifferentialAnalyses())) {
                warningModal(
                    session, "Warning",
                    "The calculated differential splicing analyses will be",
                    "discarded and the previously loaded inclusion levels will",
                    "be replaced.", footer=actionButton(
                        ns("discard"), "Discard and replace",
                        class="btn-warning", "data-dismiss"="modal"))
            } else {
                warningModal(
                    session, "Inclusion levels already quantified",
                    "Do you wish to replace the loaded inclusion levels?",
                    footer=actionButton(
                        ns("replace"), "Replace", class="btn-warning",
                        "data-dismiss"="modal"))
            }
        } else {
            calcSplicing()
        }
    })
    
    observeEvent(input$replace, calcSplicing())
    observeEvent(input$discard, {
        setDifferentialAnalyses(NULL)
        setDifferentialAnalysesSurvival(NULL)
        calcSplicing()
    })
    
    # Calculate inclusion levels
    calcSplicing <- reactive({
        eventType  <- input$eventType
        minReads   <- input$minReads
        annotation <- input$annotation
        
        if (is.null(eventType) || is.null(minReads) || is.null(annotation)) {
            return(NULL)
        } else {
            if (input$junctionQuant == "") {
                errorModal(session, "Select junction quantification",
                           "Select a junction quantification dataset")
                endProcess("calcIncLevels")
                return(NULL)
            }
        }
        time <- startProcess("calcIncLevels")
        startProgress("Quantifying alternative splicing", divisions=4)
        # Read annotation
        annot <- readAnnot(session, annotation, showProgress=TRUE)
        junctionQuant <- getJunctionQuantification()[[input$junctionQuant]]
        
        # Filter alternative splicing events based on their genes
        filter <- NULL
        if (input$filter == "select") {
            # Filter genes based on select input
            filter <- input$filterGenes
            if (identical(filter, "")) filter <- NULL
        } else if (input$filter == "file") {
            # Filter genes provided in a file
            filter <- fread(input$filterGenesFile, header=FALSE)
            if (nrow(filter) == 1) {
                filter <- as.character(filter)
            } else if (nrow(filter) > 1) {
                filter <- as.character(filter[[1]])
            } else {
                filter <- NULL
            }
        }
        
        # Quantify splicing with splicing annotation and junction quantification
        updateProgress("Calculating inclusion levels")
        psi <- quantifySplicing(annot, junctionQuant, eventType, minReads, 
                                genes=filter)
        if (nrow(psi) == 0) {
            errorModal(session, "No splicing events returned",
                       "The total reads of the alternative splicing events are",
                       "below the minium read counts' threshold.")
            endProcess("calcIncLevels")
            return(NULL)
        }
        setInclusionLevels(psi)
        
        endProcess("calcIncLevels", time)
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
    
    # # Toggle visibility of loading icon
    # observe({
    #     toggle("geneLoadingIcon",
    #            selector = paste0(
    #                '$("#data-inclusionLevels-filterGenes").parent()',
    #                '.children("div.selectize-control").hasClass("loading")'))
    # })
    
    quantifySplicingSet(session, input)
    loadCustomSplicingAnnotationSet(session, input, output)
    loadSplicingQuantificationSet(session, input, output)
}

attr(inclusionLevelsUI, "loader") <- "data"
attr(inclusionLevelsServer, "loader") <- "data"