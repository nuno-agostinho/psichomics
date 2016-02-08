# The name used for the plot must be unique
name <- "plot3"

ui <- list(sidebarLayout(
    sidebarPanel("Hey",
                 uiOutput("v")), 
    mainPanel(plotOutput(name))))

server <- function(input, output, session) {
    output[["v"]] <- renderUI({
        selectInput("var", "Choose the dataset variable to plot",
                    choices = names(inclusion.levels(shared.data$a)))
    })
        
    output[[name]] <- renderPlot({
        plot(inclusion.levels(shared.data$a)[[input$var]])
    })
}