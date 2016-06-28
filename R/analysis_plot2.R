plots2UI <- function(id) {
    ns <- NS(id)
    plotOutput(ns("plot"))
}

plots2Server <- function(input, output, session) {
    output$plot <- renderPlot({
        # validate(if(is.null(thisData$data)) return(NULL))
        plot(mtcars$cyl, mtcars$mpg)
        # ggplot(data=mtcars, aes(cyl, mpg)) + geom_bin2d()
    })
}

# attr(plots2UI, "loader") <- "analysis"
# attr(plots2UI, "name") <- "Plots 2"
# attr(plots2Server, "loader") <- "analysis"