# The name used for the plot must be unique
name <- "Inclusion levels 2"

ui <- list(
    uiOutput(name)
)

createLink <- function(val) {
    onclick <- sprintf('document.getElementById(\'%s\').selectize.setValue(\'%s\')',
                       c("selectizePlot", "selectizeEvent"),
                       c("plot4", "Valiant"))
    onclick <- paste(onclick, collapse = "; ")
    html <- paste('<a id="%s" title="Get more information about this event"',
                  'href="#" onclick="%s">%s</a>')
    link <- sprintf(html, gsub("_", " ", val), onclick, val)
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
        shared.data$psi <- calculateInclusionLevels(
            "SE", junctionQuant, annot)
        
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
        cbind(Events = createLink(rownames(shared.data$psi)), 
              shared.data$psi)
    }, escape = FALSE, options = list(pageLength = 10, scrollX = TRUE))
}