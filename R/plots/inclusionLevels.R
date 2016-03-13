# The name used for the plot must be unique
name <- "Inclusion levels"

ui <- list( uiOutput(name) )

createLink <- function(val) {
    id <- gsub(" ", "_", val)
    onclick <- sprintf(
        'document.getElementById(\'%s\').selectize.setValue(\'%s\')',
        c("selectizePlot", "selectizeEvent"),
        c(rep("plot4", length(id)), id))
    onclick <- paste(onclick, collapse = "; ")
    html <- paste('<a id="%s" title="Get more information about this event"',
                  'href="#" onclick="%s">%s</a>')
    link <- sprintf(html, id, onclick, val)
    return(link)
}

server <- function(input, output, session) {
    levels <- reactive({
        # Read annotation
        startProgress("Reading alternative splicing annotation",
                      divisions = 2)
        annot <- readAnnotation("SE")
        
        # Calculate inclusion levels with annotation and junction quantification
        updateProgress("Calculating inclusion levels")
        junctionQuant <- getJunctionQuantification()
        setInclusionLevels(
            calculateInclusionLevels("SE", junctionQuant, annot))
        
        updateProgress("Done!")
        closeProgress()
    })
    
    observe({
        if(is.null(getData()) || is.null(getJunctionQuantification()))
            disableTab("Plots")
        else
            enableTab("Plots")
        
        if(input$nav == 'Plots') levels()
    })
    
    output[[name]] <- renderUI( dataTableOutput("incLevels") )
    
    output$incLevels <- renderDataTable({
        psi <- getInclusionLevels()
        cbind(Events = createLink(rownames(psi)), psi)
    }, escape = FALSE, options = list(pageLength = 10, scrollX = TRUE))
}