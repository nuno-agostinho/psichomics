templateUI <- function(id) {
    ns <- NS(id)
    plotOutput(ns("plot"))
}

templateServer <- function(input, output, session) {
    output$plot <- renderPlot({
        plot(mtcars$cyl, mtcars$mpg)
    })
}

# attr(plots2UI, "loader") <- "analysis"
# attr(plots2UI, "name") <- "Template"
# attr(plots2Server, "loader") <- "analysis"