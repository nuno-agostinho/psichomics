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
    list(
        selectizeInput(id("eventType"), "Event type",
                       choices = c("Skipping exon" = "SE"), selected = "SE",
                       multiple = TRUE),
        actionButton(id("calcIncLevels"), class = "btn-primary",
                     "Calculate inclusion levels"))
}

server <- function(input, output, session) {
    levels <- reactive({
        eventType <- input[[id("eventType")]]

        # Read annotation
        startProgress("Reading alternative splicing annotation",
                      divisions = 3)
        annot <- readAnnotation(eventType)

        # Calculate inclusion levels with annotation and junction quantification
        updateProgress("Calculating inclusion levels")
        junctionQuant <- getJunctionQuantification()
        psi <- calculateInclusionLevels(eventType, junctionQuant, annot)
        setInclusionLevels(psi)

        updateProgress("Matching clinical data")
        match <- matchIdWithClinical(colnames(psi), getClinicalData())
        # TODO: make sure this attribute allows to show row names in Data section
        # attr(match, "rowNames") <- TRUE
        setClinicalMatchFrom("Inclusion levels", match)
        closeProgress()
    })

    observeEvent(input[[id("calcIncLevels")]], {
        if(is.null(getData()) || is.null(getJunctionQuantification()))
            print("No junction quantification data!")
        else
            levels()
    })
}