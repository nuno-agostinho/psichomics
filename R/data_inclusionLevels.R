#' Splicing event types available
#' @return Named character vector with splicing event types
#' @export
#' 
#' @examples 
#' getSplicingEventTypes()
getSplicingEventTypes <- function() {
    c("Skipped exon (SE)" = "SE",
      "Mutually exclusive exons (MXE)" = "MXE",
      "Alternative 5' Splice Site (A5SS)" = "A5SS",
      "Alternative 3' Splice Site (A3SS)" = "A3SS",
      "Alternative first exon (AFE)" = "AFE",
      "Alternative last exon (ALE)" = "ALE")
}

#' List the alternative splicing annotation files available
#' @return Named character vector with splicing annotation files available
#' @export
#' 
#' @examples
#' listSplicingAnnotation()
listSplicingAnnotation <- function() {
    c("Human (hg19/GRCh37)"="hg19_splicingAnnotation.RDS")
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
        selectizeInput(ns("junctionQuant"), 
                       "Alternative splicing junction quantification",
                       choices=NULL),
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
#' eventType <- c("SE", "MXE")
#' annot <- readFile("ex_splicing_annotation.RDS")
#' junctionQuant <- readFile("ex_junctionQuant.RDS")
#' 
#' psi <- quantifySplicing(annot, junctionQuant, eventType=c("SE", "MXE"))
quantifySplicing <- function(annotation, junctionQuant, eventType="SE", 
                             minReads=10, progress=echoProgress) {
    psi <- NULL
    for (i in seq_along(eventType)) {
        type <- eventType[[i]]
        progress("Calculating inclusion levels",
                 names(getSplicingEventTypes())[[i]], value = i, 
                 max = length(eventType))
        
        if (i == "AFE") 
            annotation$AFE <- annotation$AFE[!is.na(annotation$AFE$C2.start), ]
        if (i == "ALE") 
            annotation$ALE <- annotation$ALE[!is.na(annotation$ALE$C1.end), ]
        psi <- rbind(psi, calculateInclusionLevels(
            type, junctionQuant, annotation[[type]], minReads))
    }
    attr(psi, "rowNames") <- TRUE
    attr(psi, "description") <- paste("Exon and intron inclusion levels",
                                      "for any given alternative splicing",
                                      "event.")
    attr(psi, "dataType") <- "Inclusion levels"
    return(psi)
}

#' Server logic of the alternative splicing event quantification module
#' 
#' @param input Shiny input
#' @param output Shiny ouput
#' @param session Shiny session
#' 
#' @importFrom shiny reactive observeEvent
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
        
        if (is.null(eventType) || is.null(minReads) || is.null(annotation))
            return(NULL)
        time <- startProcess("calcIncLevels")
        
        # Read annotation
        startProgress("Reading alternative splicing annotation", divisions=3)
        annot <- readFile(annotation)
        
        # Set species and assembly version
        allAnnot <- listSplicingAnnotation()
        annotID <- names(allAnnot)[match(annotation, allAnnot)]
        if (grepl("Human", annotID)) setSpecies("Human")
        if (grepl("hg19", annotID)) setAssemblyVersion("hg19")
        
        if (input$junctionQuant == "") {
            errorModal(session, "Select junction quantification",
                       "Select a junction quantification dataset")
            endProcess("calcIncLevels")
            return(NULL)
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
}

attr(inclusionLevelsUI, "loader") <- "data"
attr(inclusionLevelsServer, "loader") <- "data"