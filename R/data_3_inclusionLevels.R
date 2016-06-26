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

inclusionLevelsUI <- function(id, tab) {
    ns <- NS(id)
    choices <- c("Skipping exon (SE)" = "SE",
                 "Mutually exclusive exons (MXE)" = "MXE",
                 "Alternative 5' Splice Site (A5SS)" = "A5SS",
                 "Alternative 3' Splice Site (A3SS)" = "A3SS",
                 "Alternative first exon (AFE)" = "AFE",
                 "Alternative last exon (ALE)" = "ALE")
    
    tab("Inclusion levels",
        uiOutput(ns("modal")),
        helpText("Calculate exon inclusion levels. This is also",
                 "known as percentage spliced in (PSI or \u03A8)."),
        selectizeInput(ns("annotation"),
                       "Alternative splicing event annotation",
                       choices = c("Human (hg19/GRCh37)"="hg19_splicingAnnotation.RDS")),
        selectizeInput(ns("eventType"), "Event type(s)", selected = "SE",
                       choices = choices, multiple = TRUE),
        numericInput(ns("minReads"), "Minimum reads threshold", value = 10),
        actionButton(ns("calcIncLevels"), class = "btn-primary",
                     "Calculate inclusion levels"))
}

inclusionLevelsServer <- function(input, output, session) {
    choices <- c("Skipping exon (SE)" = "SE",
                 "Mutually exclusive exons (MXE)" = "MXE",
                 "Alternative 5' Splice Site (A5SS)" = "A5SS",
                 "Alternative 3' Splice Site (A3SS)" = "A3SS",
                 "Alternative first exon (AFE)" = "AFE",
                 "Alternative last exon (ALE)" = "ALE")
    
    levels <- reactive({
        eventType <- input$eventType
        minReads  <- input$minReads
        
        if (is.null(eventType) || is.null(minReads)) return(NULL)

        # Read annotation
        startProgress("Reading alternative splicing annotation", divisions = 3)
        annot <- readRDS(system.file("extdata", input$annotation,
                                     package = "psichomics"))

        # Calculate inclusion levels with annotation and junction quantification
        updateProgress("Calculating inclusion levels")
        junctionQuant <- getJunctionQuantification()

        psi <- NULL
        for (i in seq_along(eventType)) {
            type <- eventType[[i]]
            updateProgress("Calculating inclusion levels", names(choices)[[i]],
                           value = i, max = length(eventType))

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
    
    observeEvent(input$calcIncLevels, {
        if(is.null(getData()) || is.null(getJunctionQuantification()))
            errorModal(session, "Data missing",
                       "No junction quantification data loaded!")
        else
            levels()
    })
}

attr(inclusionLevelsUI, "loader") <- "data"
attr(inclusionLevelsServer, "loader") <- "data"