# The name used for the plot must be unique
name <- "plot5"

ui <- plotOutput(name)

server <- function(input, output, session) {
    output[[name]] <- renderPlot(
        ggplot(data=mtcars, aes(gear, carb)) + geom_dotplot()
    )
}