# The name used for the plot must be unique
plot <- "plot2"
id <- function(value) objectId(name, plot, value)

ui <- function() plotOutput(id(plot))

server <- function(input, output, session) {
    output[[id(plot)]] <- renderPlot({
        # validate(if(is.null(thisData$data)) return(NULL))
        ggplot(data=mtcars, aes(cyl, mpg)) + geom_bin2d()
    })
}