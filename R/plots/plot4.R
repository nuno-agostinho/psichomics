# The name used for the plot must be unique
name <- "plot4"

ui <- list(
    sidebarLayout(
        sidebarPanel(
            selectizeInput(
                "colours",
                "Pick a colour",
                choices = list("red", "orange", "blue", "pink", "green")
            )
        ), 
        mainPanel( plotOutput(name) )
    )
)

server <- function(input, output, session) {
    output[[name]] <- renderPlot(
        ggplot(data=mtcars, aes(gears, carb)) +
            geom_bar(stat = "identity", colour = input$colours)
    )
}