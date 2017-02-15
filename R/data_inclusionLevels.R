#' Parse and prepare sample information from TCGA samples
#' 
#' @param samples Character: sample identifiers
#' @inheritParams getClinicalMatchFrom
#' 
#' @return Data frame containing metadata associated with each TCGA sample
parseTcgaSampleInfo <- function (samples, category=getCategory()) {
    parsed <- parseSampleGroups(samples)
    if ( all(is.na(parsed)) ) return(NULL)
    
    info <- data.frame(parsed)
    colnames(info) <- "Sample types"
    rownames(info) <- samples
    
    # Patient match
    patients <- getPatientId()
    match <- getClinicalMatchFrom("Inclusion levels", category)
    if ( !is.null(patients) ) {
        if (is.null(match))
            match <- getPatientFromSample(samples, patients)
        
        match <- match[samples]
        patients <- patients[match]
        info <- cbind(info, "Patient ID"=patients)
    }
    
    # Metadata
    attr(info, "rowNames") <- TRUE
    attr(info, "description") <- "Metadata for TCGA samples"
    attr(info, "dataType")  <- "Sample metadata"
    attr(info, "tablename") <- "Sample metadata"
    return(info)
}

#' Splicing event types available
#' @return Named character vector with splicing event types
#' @export
#' 
#' @examples 
#' getSplicingEventTypes()
getSplicingEventTypes <- function() {
    c("Skipped exon" = "SE",
      "Mutually exclusive exon" = "MXE",
      "Alternative 5' splice site" = "A5SS",
      "Alternative 3' splice site" = "A3SS",
      "Alternative first exon" = "AFE",
      "Alternative last exon" = "ALE")
}

#' List the alternative splicing annotation files available
#' 
#' @return Named character vector with splicing annotation files available
#' @export
#' 
#' @examples
#' listSplicingAnnotations()
listSplicingAnnotations <- function() {
    c("Human (hg19/GRCh37)"="annotationHub_alternativeSplicingEvents.hg19")
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
    list("Available annotation"=listSplicingAnnotations(),
         "Custom annotation"=c(
             ..., "Load annotation from file..."="loadAnnotation"))
}

#' Interface to quantify alternative splicing
#' 
#' @param ns Namespace function
#' 
#' @importFrom shiny tagList uiOutput selectizeInput numericInput actionButton
#' @importFrom shinyBS bsTooltip
#' 
#' @return HTML elements
inclusionLevelsInterface <- function(ns) {
    tagList(
        uiOutput(ns("modal")),
        helpText("Exon inclusion levels are measured from junction",
                 "quantification using the Percent Spliced-In (PSI) metric."),
        selectizeInput(ns("junctionQuant"), choices=NULL,
                       "Alternative splicing junction quantification"),
        selectizeInput(ns("annotation"), choices=listAllAnnotations(),
                       "Alternative splicing event annotation"),
        selectizeInput(ns("eventType"), "Event type(s)", selected = "SE",
                       choices=getSplicingEventTypes(), multiple = TRUE),
        numericInput(ns("minReads"), div("Minimum read counts threshold",
                                         icon("question-circle")), value = 10),
        bsTooltip(ns("minReads"), placement = "right", 
                  options = list(container = "body"),
                  paste("Inclusion levels calculated with a number of read",
                        "counts below this threshold are discarded.")),
        actionButton(ns("loadIncLevels"), "Load from file"),
        processButton(ns("calcIncLevels"), "Quantify events"))
}

#' Interface of the alternative splicing event quantification module
#' 
#' @param id Character: identifier
#' @param panel Function to process HTML elements
#' 
#' @return HTML elements
inclusionLevelsUI <- function(id, panel) {
    ns <- NS(id)
    title <- "Alternative splicing quantification"
    panel(style="info", title=tagList(icon("calculator"), title), value=title,
          inclusionLevelsInterface(ns))
}

#' Quantify alternative splicing events
#' 
#' @param annotation List of data frames: annotation for each alternative
#' splicing event type
#' @param junctionQuant Data frame: junction quantification
#' @param eventType Character: splicing event types to quantify
#' @param minReads Integer: minimum of read counts to consider a junction read 
#' in calculations
#' @param progress Function to track the progress
#' 
#' @return Data frame with the quantification of the alternative splicing events
#' @export
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
                             minReads=10, progress=echoProgress) {
    psi <- NULL
    for (acronym in eventType) {
        eventTypes <- getSplicingEventTypes()
        type <- names(eventTypes)[[match(acronym, eventTypes)]]
        
        if (!is.null(annotation[[type]])) {
            progress("Calculating inclusion levels", type, value=acronym, 
                     max=length(eventType))
            psi <- rbind(psi, calculateInclusionLevels(
                acronym, junctionQuant, annotation[[type]], minReads))
        }
    }
    if (!is.null(psi)) {
        attr(psi, "rowNames") <- TRUE
        attr(psi, "description") <- paste("Exon and intron inclusion levels",
                                          "for any given alternative splicing",
                                          "event.")
        attr(psi, "dataType")  <- "Inclusion levels"
        attr(psi, "tablename") <- "Inclusion levels"
    }
    return(psi)
}

#' Load alternative splicing annotation from AnnotationHub
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

#' Server logic of the alternative splicing event quantification module
#' 
#' @param input Shiny input
#' @param output Shiny ouput
#' @param session Shiny session
#' 
#' @importFrom shiny reactive observeEvent fileInput helpText removeModal
#' @importFrom tools file_path_sans_ext
#' 
#' @return NULL (this function is used to modify the Shiny session's state)
inclusionLevelsServer <- function(input, output, session) {
    ns <- session$ns
    
    observeEvent(input$takeMeThere, missingDataGuide("Junction quantification"))
    observeEvent(input$takeMeToClinical, missingDataGuide("Clinical data"))
    
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
    
    observeEvent(input$calcIncLevels, {
        if (is.null(getData()) || is.null(getJunctionQuantification())) {
            missingDataModal(session, "Junction quantification",
                             ns("takeMeThere"))
        } else if (!is.null(getInclusionLevels())) {
            if (!is.null(getDifferentialAnalyses())) {
                warningModal(session, "Warning",
                             "The calculated differential splicing analyses",
                             "will be discarded and the previously loaded",
                             "inclusion levels will be replaced.",
                             footer=actionButton(ns("discard"),
                                                 "Discard and replace",
                                                 class="btn-warning",
                                                 "data-dismiss"="modal"))
            } else {
                warningModal(session, "Inclusion levels already quantified",
                             "Do you wish to replace the inclusion levels",
                             "loaded?",
                             footer=actionButton(ns("replace"), "Replace",
                                                 class="btn-warning",
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
    
    calcSplicing <- reactive({
        eventType <- input$eventType
        minReads  <- input$minReads
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
        if (grepl("^/var/folders/", annotation)) { # if custom annotation
            updateProgress("Loading alternative splicing annotation")
            annot <- readRDS(annotation)
        } else if (grepl("^annotationHub_", annotation)) {
            updateProgress("Downloading alternative splicing annotation")
            annot <- loadAnnotation(annotation)
            
            # Set species and assembly version
            allAnnot <- listSplicingAnnotations()
            annotID <- names(allAnnot)[match(annotation, allAnnot)]
            if (grepl("Human", annotID)) setSpecies("Human")
            if (grepl("hg19", annotID)) setAssemblyVersion("hg19")
        }
        junctionQuant <- getJunctionQuantification()[[input$junctionQuant]]
        
        # Calculate inclusion levels with annotation and junction quantification
        updateProgress("Calculating inclusion levels")
        psi <- quantifySplicing(annot, junctionQuant, eventType, minReads, 
                                progress=updateProgress)
        setInclusionLevels(psi)
        
        samples <- colnames(psi)
        parsed <- parseTcgaSampleInfo(samples) 
        if ( !is.null(parsed) )
            setSampleInfo(parsed)
        
        endProcess("calcIncLevels", time)
    })
    
    # Show modal to load custom alternative splicing quantification
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
                                     choices="Human", options=list(create=TRUE)),
                      selectizeInput(ns("customAssembly"), "Assembly",
                                     choices="hg19", options=list(create=TRUE)),
                      uiOutput(ns("alert")),
                      footer=actionButton(ns("loadCustom"), "Load annotation",
                                          class="btn-primary"))
        }
    })
    
    # Load custom alternative splicing annotation
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
    
    # Show modal for loading alternative splicing quantification
    observeEvent(input$loadIncLevels, {
        ns <- session$ns
        
        infoModal(
            session, "Load alternative splicing quantification",
            helpText("A table containing the sample identifiers as columns and",
                     "the alternative splicing event identifiers as rows is",
                     "recommended. The event identifier should be similar to:"),
            tags$kbd(style="word-wrap: break-word;",
                     paste0("EventType_Chromosome_Strand_Coordinate1_",
                            "_Coordinate2_..._Gene")),
            tags$hr(),
            fileInput(ns("customASquant"), "Choose a file"),
            selectizeInput(ns("customSpecies2"), "Species", choices="Human",
                           options=list(create=TRUE)),
            selectizeInput(ns("customAssembly2"), "Assembly", choices="hg19",
                           options=list(create=TRUE)),
            uiOutput(ns("alertIncLevels")),
            footer=processButton(ns("loadASquant"), "Load quantification"))
    })
    
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
    
    observeEvent(input$replace2, {
        setSampleId(NULL)
        setGroupsFrom("Clinical data", NULL)
        loadSplicing()
    })
    observeEvent(input$discard2, {
        setDifferentialAnalyses(NULL)
        setDifferentialAnalysesSurvival(NULL)
        setSampleId(NULL)
        setGroupsFrom("Clinical data", NULL)
        loadSplicing()
    })
    
    # Load alternative splicing quantification
    loadSplicing <- reactive({
        time <- startProcess("loadIncLevels")
        
        startProgress("Wait a moment", divisions=2)
        updateProgress("Loading alternative splicing quantification")
        psi <- tryCatch(read.delim(input$customASquant$datapath, 
                                   row.names=1, check.names=FALSE),
                        error=return, warning=return)
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
            attr(psi, "rowNames") <- TRUE
            attr(psi, "description") <- paste("Exon and intron inclusion",
                                              "levels for any given",
                                              "alternative splicing event.")
            attr(psi, "dataType")  <- "Inclusion levels"
            attr(psi, "tablename") <- "Inclusion levels"
            
            if ( is.null(getData()) ) {
                name <- file_path_sans_ext( input$customASquant$name )
                name <- gsub(" Inclusion levels$", "", name)
                if (name == "") name <- "Unnamed"
                
                data <- setNames(list(list("Inclusion levels"=psi)), name)
                data <- processDatasetNames(data)
                setData(data)
                setCategory(name)
            }
            
            setInclusionLevels(psi)
            
            samples <- colnames(psi)
            parsed <- parseTcgaSampleInfo(samples) 
            if ( !is.null(parsed) )
                setSampleInfo(parsed)
            
            setSpecies(input$customSpecies2)
            setAssemblyVersion(input$customAssembly2)
            
            removeModal()
        }
        endProcess("loadIncLevels", time)
    })
}

attr(inclusionLevelsUI, "loader") <- "data"
attr(inclusionLevelsServer, "loader") <- "data"