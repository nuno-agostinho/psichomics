# The name used for the plot must be unique
name <- "plot1"

ui <- list(
    sidebarLayout(
        sidebarPanel(
            selectizeInput("x", "Pick x axis", choices = names(mtcars)),
            selectizeInput("y", "Pick y axis", choices = names(mtcars)),
            shiny::actionButton("change", "Change to plot2")
        ), 
        mainPanel( plotOutput(name) )
    )
)

server <- function(input, output, session) {
    output[[name]] <- renderPlot({
        ggplot(data=mtcars, aes_string(input$x, input$y)) + geom_bin2d()
    })
    
    observeEvent(input$change, {
        updateSelectizeInput(session, "selectizePlot", selected = "plot2")
    })
}