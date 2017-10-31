#' @rdname appUI
#' 
#' @importFrom shiny NS tagList sidebarPanel mainPanel sliderInput actionButton
#' uiOutput
#' @importFrom highcharter highchartOutput
templateUI <- function(id) {
    ns <- NS(id) # Identifier
    tagList(
        sidebarPanel( # Sidebar interface
            sliderInput(inputId=ns("num"), label="Number of samples",
                        min=1, max=20, value=4), 
            actionButton(inputId=ns("submit"), label="Plot")),
        mainPanel( # Main interface
            uiOutput(outputId=ns("text")),
            highchartOutput(outputId=ns("plot"))))
}

#' @rdname appServer
#' 
#' @importFrom shiny renderUI observeEvent isolate tagList tags
#' @importFrom highcharter renderHighchart
templateServer <- function(input, output, session) {
    # Wait for user to press the button
    observeEvent(input$submit, {
        # Break reactive dependence
        num <- isolate(input$num)
        
        # Retrieve random samples
        sample <- sort( sample(100, size=num, replace=TRUE) )
        
        # Print samples
        output$text <- renderUI(tagList(tags$b("Samples:"), 
                                        paste(sample, collapse=" ")))
        # Plot frequency
        output$plot <- renderHighchart( hchart(sample) )
    })
}

# attr(templateUI, "loader") <- "analysis"
# attr(templateUI, "name") <- "Template"
# attr(templateServer, "loader") <- "analysis"