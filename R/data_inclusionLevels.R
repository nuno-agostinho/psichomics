# createLink <- function(val) {
#     id <- gsub(" ", "_", val)
#     js <- 'document.getElementById(\'%s\').selectize.setValue(\'%s\')'
#
#     onclick <- paste(sprintf(js, objectId(name, "selectizePlot"), "plot4"),
#                      sprintf(js, objectId(name, "selectizeEvent"), id),
#                      sep = "; ")
#
#     html <- paste('<a id="%s" title="Get more information about this event"',
#                   'href="#" onclick="%s">%s</a>')
#     link <- sprintf(html, id, onclick, val)
#     return(link)
# }

#' Splicing event types available
#' @return Named character vector with splicing event types
getSplicingEventChoices <- function() {
    c("Skipping exon (SE)" = "SE",
      "Mutually exclusive exons (MXE)" = "MXE",
      "Alternative 5' Splice Site (A5SS)" = "A5SS",
      "Alternative 3' Splice Site (A3SS)" = "A3SS",
      "Alternative first exon (AFE)" = "AFE",
      "Alternative last exon (ALE)" = "ALE")
}

#' Interface
#' @param ns Namespace function
#' @return HTML elements
inclusionLevelsInterface <- function(ns) {
    tagList(
        uiOutput(ns("modal")),
        helpText("Measure exon inclusion levels from junction quantification.",
                 "This is also known as Percentage Spliced In (PSI)."),
        selectizeInput(ns("annotation"),
                       "Alternative splicing event annotation",
                       choices=c("Human (hg19/GRCh37)"=
                                     "hg19_splicingAnnotation.RDS")),
        selectizeInput(ns("eventType"), "Event type(s)", selected = "SE",
                       choices=getSplicingEventChoices(), multiple = TRUE),
        numericInput(ns("minReads"), "Minimum reads threshold", value = 10),
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
    panel(style="info", title=list(icon("calculator"), title), value=title,
          inclusionLevelsInterface(ns))
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
    levels <- reactive({
        eventType <- input$eventType
        minReads  <- input$minReads
        
        if (is.null(eventType) || is.null(minReads)) return(NULL)
        
        # Read annotation
        startProgress("Reading alternative splicing annotation", divisions = 3)
        annot <- readRDS(system.file("extdata", input$annotation,
                                     package = "psichomics"))
        
        # Set species and assembly version
        if (grepl("Human", "Human (hg19/GRCh37)")) setSpecies("Human")
        if (grepl("hg19", "Human (hg19/GRCh37)")) setAssemblyVersion("hg19")
        
        # Calculate inclusion levels with annotation and junction quantification
        updateProgress("Calculating inclusion levels")
        junctionQuant <- getJunctionQuantification()
        
        psi <- NULL
        for (i in seq_along(eventType)) {
            type <- eventType[[i]]
            updateProgress("Calculating inclusion levels", 
                           names(getSplicingEventChoices())[[i]], value = i, 
                           max = length(eventType))
            
            if (i == "AFE") annot$AFE <- annot$AFE[!is.na(annot$AFE$C2.start), ]
            if (i == "ALE") annot$ALE <- annot$ALE[!is.na(annot$ALE$C1.end), ]
            psi <- rbind(psi, calculateInclusionLevels(type, junctionQuant,
                                                       annot[[type]], minReads))
        }
        attr(psi, "rowNames") <- TRUE
        attr(psi, "description") <- paste("Exon and intron inclusion levels",
                                          "for any given alternative splicing",
                                          "event.")
        setInclusionLevels(psi)
        
        updateProgress("Matching clinical data")
        match <- matchIdWithClinical(colnames(psi), getClinicalData())
        match <- match[!is.na(match)] # remove non-matching IDs
        setClinicalMatchFrom("Inclusion levels", match)
        closeProgress()
    })
    
    observeEvent(input$takeMeThere, missingDataGuide("Junction quantification"))
    
    observeEvent(input$calcIncLevels, {
        if(is.null(getData()) || is.null(getJunctionQuantification()))
            missingDataModal(session, "Junction quantification",
                             ns("takeMeThere"))
        else
            levels()
    })
}

attr(inclusionLevelsUI, "loader") <- "data"
attr(inclusionLevelsServer, "loader") <- "data"