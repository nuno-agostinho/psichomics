# The name used for the plot must be unique
name <- "plot1"

ui <- list(
    sidebarLayout(
        sidebarPanel(
            selectizeInput("x",
                           "Pick x axis",
                           choices = names(h)),
            selectizeInput("y",
                           "Pick y axis",
                           choices = names(h))
        ), 
        mainPanel( plotOutput(name) )
    )
)

server <- function(input, output, session) {
    output[[name]] <- renderPlot({
        ggplot(data=h, aes_string(input$x, input$y)) + geom_bin2d()
    })
}