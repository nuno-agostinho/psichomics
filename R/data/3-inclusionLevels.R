## TODO(NunoA): Don't calculate inclusion levels with less than X reads (just
## put as NA since matrix calculation ignores NAs)

name <- "Exon/intron inclusion levels"

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

choices <- c("Skipping exon (SE)" = "SE",
             "Mutually exclusive exons (MXE)" = "MXE",
             "Alternative 5' Splice Site (A5SS)" = "A5SS",
             "Alternative 3' Splice Site (A3SS)" = "A3SS",
             "Alternative first exon (AFE)" = "AFE",
             "Alternative last exon (ALE)" = "ALE")

ui <- function() {
    tagList(
        helpText("Calculate exon and intron inclusion levels. This is also",
                 "known as percentage spliced in or PSI or even Ψ."),
        selectizeInput(id("eventType"), "Event type(s)", selected = "SE",
                       choices = choices, multiple = TRUE),
        numericInput(id("minReads"), "Minimum reads to consider", value = 10),
        actionButton(id("calcIncLevels"), class = "btn-primary",
                     "Calculate inclusion levels"))
}

server <- function(input, output, session) {
    levels <- reactive({
        eventType <- input[[id("eventType")]]
        minReads  <- input[[id("minReads")]]
        
        if (is.null(eventType) || is.null(minReads)) return(NULL)

        # Read annotation
        startProgress("Reading alternative splicing annotation", divisions = 3)
        annot <- readRDS("data/splicingAnnotation.RDS")
        
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
        setClinicalMatchFrom("Inclusion levels", match)
        closeProgress()
    })

    observeEvent(input[[id("calcIncLevels")]], {
        if(is.null(getData()) || is.null(getJunctionQuantification()))
            errorModal(session, "Data missing",
                       "No junction quantification data loaded!")
        else
            levels()
    })
}