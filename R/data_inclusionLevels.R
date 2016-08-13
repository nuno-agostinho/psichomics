#' Splicing event types available
#' @return Named character vector with splicing event types
#' @export
#' 
#' @examples 
#' getSplicingEventTypes()
getSplicingEventTypes <- function() {
    c("Exon skipping (SE)" = "SE",
      "Mutually exclusive exons (MXE)" = "MXE",
      "Alternative 5' Splice Site (A5SS)" = "A5SS",
      "Alternative 3' Splice Site (A3SS)" = "A3SS",
      "Alternative first exon (AFE)" = "AFE",
      "Alternative last exon (ALE)" = "ALE")
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
                 "This is also known as Percentage Spliced In (PSI)."),
        selectizeInput(ns("junctionQuant"), "Junction quantification",
                       choices=NULL),
        selectizeInput(ns("annotation"),
                       "Alternative splicing event annotation",
                       choices=c("Human (hg19/GRCh37)"=
                                     "hg19_splicingAnnotation.RDS")),
        selectizeInput(ns("eventType"), "Event type(s)", selected = "SE",
                       choices=getSplicingEventTypes(), multiple = TRUE),
        numericInput(ns("minReads"), div("Minimum reads threshold",
                                         icon("question-circle")), value = 10),
        bsTooltip(ns("minReads"), placement = "right", 
                  options = list(container = "body"),
                  "Reads below this threshold will be discarded"),
        actionButton(ns("calcIncLevels"), class = "btn-primary",
                     "Calculate inclusion levels"))
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
#' @param progress Function to track the progess
#' 
#' @return Data frame with the quantification of the alternative splicing events
#' @export
#' 
#' @examples 
#' \dontrun{
#' annotation <- readRDS(system.file("extdata", "hg19_splicingAnnotation.RDS",
#'                       package="psichomics"))
#'                       
#' data <- loadFirehoseData(cohort = "ACC", data=c("Clinical", 
#'                                                 "junction_quantification"))
#' junctionQuant <- data[[1]]$`Junction quantification (Illumina HiSeq)`
#'                       
#' # Check splicing event types that can be calculated
#' getSplicingEventTypes()
#' 
#' eventType <- c("SE", "A5SS")
#' quantifySplicing(annotation, junctionQuant, eventType, minReads=10)
#' }
quantifySplicing <- function(annotation, junctionQuant, eventType="SE", 
                             minReads=10, progress=printPaste) {
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
    
    observeEvent(input$calcIncLevels, {
        if(is.null(getData()) || is.null(getJunctionQuantification())) {
            missingDataModal(session, "Junction quantification",
                             ns("takeMeThere"))
        } else {
            eventType <- input$eventType
            minReads  <- input$minReads
            
            if (is.null(eventType) || is.null(minReads)) return(NULL)
            
            # Read annotation
            startProgress("Reading alternative splicing annotation",
                          divisions = 3)
            annot <- readRDS(system.file("extdata", input$annotation, 
                                         package="psichomics"))
            
            # Set species and assembly version
            if (grepl("Human", "Human (hg19/GRCh37)")) setSpecies("Human")
            if (grepl("hg19", "Human (hg19/GRCh37)")) setAssemblyVersion("hg19")
            
            if (input$junctionQuant == "") return(
                errorModal(session, "Select junction quantification",
                           "Select a junction quantification dataset"))
            junctionQuant <- getJunctionQuantification()[[input$junctionQuant]]
            
            # Calculate inclusion levels with annotation and junction
            # quantification
            updateProgress("Calculating inclusion levels")
            psi <- quantifySplicing(annot, junctionQuant, eventType, 
                                    minReads, progress=updateProgress)
            setInclusionLevels(psi)
            
            updateProgress("Matching clinical data")
            match <- matchIdWithClinical(colnames(psi), getClinicalData())
            setClinicalMatchFrom("Inclusion levels", match)
            
            closeProgress()
        }
    })
}

attr(inclusionLevelsUI, "loader") <- "data"
attr(inclusionLevelsServer, "loader") <- "data"