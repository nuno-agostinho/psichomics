# The name used for the plot must be unique
name <- "plot2"

ui <- plotOutput(name)

server <- function(input, output, session) {
    output[[name]] <- renderPlot({
        # validate(if(is.null(thisData$data)) return(NULL))
        ggplot(data=mtcars, aes(cyl, mpg)) + geom_bin2d()
    })
}