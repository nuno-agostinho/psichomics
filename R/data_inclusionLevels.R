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
#' @param ... Custom annotation loaded
#' 
#' @return Named character vector with splicing annotation files available
#' @export
#' 
#' @examples
#' listSplicingAnnotation()
listSplicingAnnotation <- function(...) {
    list(
        "Available annotations"=c(
            "Human (hg19/GRCh37)"="hg19_splicingAnnotation.RDS"),
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
        helpText("Measure exon inclusion levels from junction quantification.",
                 "The Percent Spliced-In (PSI) metric is used."),
        selectizeInput(ns("junctionQuant"), choices=NULL,
                       "Alternative splicing junction quantification"),
        selectizeInput(ns("annotation"), choices=listSplicingAnnotation(),
                       "Alternative splicing event annotation"),
        selectizeInput(ns("eventType"), "Event type(s)", selected = "SE",
                       choices=getSplicingEventTypes(), multiple = TRUE),
        numericInput(ns("minReads"), div("Minimum read counts threshold",
                                         icon("question-circle")), value = 10),
        bsTooltip(ns("minReads"), placement = "right", 
                  options = list(container = "body"),
                  paste("Inclusion levels calculated with a number of read",
                        "counts below this threshold are discarded.")),
        processButton(ns("calcIncLevels"), "Calculate inclusion levels"))
}

#' Interface of the alternative splicing event quantification module
#' 
#' @param id Character: identifier
#' @param panel Function to process HTML elements
#' 
#' @return HTML elements
inclusionLevelsUI <- function(id, panel) {
    ns <- NS(id)
    title <- "Quantify alternative splicing events"
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
quantifySplicing <- function(annotation, junctionQuant, eventType="SE", 
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
        attr(psi, "dataType") <- "Inclusion levels"
    }
    return(psi)
}

#' Server logic of the alternative splicing event quantification module
#' 
#' @param input Shiny input
#' @param output Shiny ouput
#' @param session Shiny session
#' 
#' @importFrom shiny reactive observeEvent fileInput helpText removeModal
#' @return NULL (this function is used to modify the Shiny session's state)
inclusionLevelsServer <- function(input, output, session) {
    ns <- session$ns
    
    observeEvent(input$takeMeThere, missingDataGuide("Junction quantification"))
    
    observe({
        junctionQuant <- getJunctionQuantification()
        if (!is.null(junctionQuant)) {
            updateSelectizeInput(session, "junctionQuant",
                                 choices=c(names(junctionQuant),
                                           "Select junction quantification"=""))
        } else {
            updateSelectizeInput(session, "junctionQuant",
                                 choices=c("No junction quantification loaded"=""))
        }
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
        
        # Read annotation
        startProgress("Reading alternative splicing annotation", divisions=3)
        if (grepl("^/var/folders/", annotation)) { # if custom annotation
            annot <- readRDS(annotation)
        } else {
            annot <- readFile(annotation)
            
            # Set species and assembly version
            allAnnot <- listSplicingAnnotation()
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
        
        updateProgress("Matching clinical data")
        match <- getPatientFromSample(colnames(psi), getClinicalData())
        setClinicalMatchFrom("Inclusion levels", match)
        
        endProcess("calcIncLevels", time)
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
    
    observe({
        ns <- session$ns
        if (input$annotation == "loadAnnotation") {
            url <- "http://rpubs.com/nuno-agostinho/alt-splicing-annotation"
            
            updateSelectizeInput(session, "annotation", 
                                 selected=listSplicingAnnotation()[[1]])
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
                                 choices=listSplicingAnnotation(custom))
            setSpecies(input$customSpecies)
            setAssemblyVersion(input$customAssembly)
            removeModal()
            removeAlert(output)
        }
    })
}

attr(inclusionLevelsUI, "loader") <- "data"
attr(inclusionLevelsServer, "loader") <- "data"