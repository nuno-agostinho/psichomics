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

ui <- function() {
    tagList(
        helpText("Calculate exon and intron inclusion levels. This is also",
                 "known as percentage spliced in or PSI or even Ψ."),
        selectizeInput(id("eventType"), "Event type",
                       choices = c("Skipping exon" = "SE"), selected = "SE",
                       multiple = TRUE),
        numericInput(id("minReads"), "Minimum reads to consider", value = 10),
        actionButton(id("calcIncLevels"), class = "btn-primary",
                     "Calculate inclusion levels"))
}

server <- function(input, output, session) {
    levels <- reactive({
        eventType <- input[[id("eventType")]]
        minReads  <- input[[id("minReads")]]

        # Read annotation
        startProgress("Reading alternative splicing annotation", divisions = 3)
        annot <- readAnnotation(eventType)

        # Calculate inclusion levels with annotation and junction quantification
        updateProgress("Calculating inclusion levels")
        junctionQuant <- getJunctionQuantification()
        psi <- calculateInclusionLevels(eventType, junctionQuant, annot,
                                        minReads = minReads)
        attr(psi, "rowNames") <- TRUE
        attr(psi, "description") <- "Exon and intron inclusion levels for any given alternative splicing event."
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