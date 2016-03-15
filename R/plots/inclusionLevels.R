# The name used for the plot must be unique
name <- "Inclusion levels"

ui <- list(
    sidebarPanel(
        selectizeInput("eventType", "Event type",
                       choices = c("Skipping exon" = "SE"), selected = "SE",
                       multiple = TRUE),
        actionButton("calcIncLevels", "Calculate inclusion levels")),
    mainPanel( dataTableOutput("incLevels") )
)

createLink <- function(val) {
    id <- gsub(" ", "_", val)
    js <- 'document.getElementById(\'%s\').selectize.setValue(\'%s\')'
    onclick <- paste(sprintf(js, "selectizePlot", "plot4"),
                     sprintf(js, "selectizeEvent", id), sep = "; ")
    
    html <- paste('<a id="%s" title="Get more information about this event"',
                  'href="#" onclick="%s">%s</a>')
    link <- sprintf(html, id, onclick, val)
    return(link)
}

server <- function(input, output, session) {
    levels <- reactive({
        eventType <- input$eventType
        
        # Read annotation
        startProgress("Reading alternative splicing annotation",
                      divisions = 2)
        annot <- readAnnotation(eventType)
        
        # Calculate inclusion levels with annotation and junction quantification
        updateProgress("Calculating inclusion levels")
        junctionQuant <- getJunctionQuantification()
        setInclusionLevels(
            calculateInclusionLevels(eventType, junctionQuant, annot))
        
        updateProgress("Done!")
        closeProgress()
    })
    
    observeEvent(input$calcIncLevels, {
        if(is.null(getData()) || is.null(getJunctionQuantification()))
            print("No junction quantification data!")
        else
            levels()
    })
    
    output$incLevels <- renderDataTable({
        psi <- getInclusionLevels()
        if (is.null(psi)) return(NULL)
        cbind(Events = createLink(rownames(psi)), psi)
    }, escape = FALSE, options = list(pageLength = 10, scrollX = TRUE))
}