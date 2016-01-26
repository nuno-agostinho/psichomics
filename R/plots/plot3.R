# The name used for the plot must be unique
name <- "plot3"

ui <- list(sidebarLayout(
    sidebarPanel("Hey"), 
    mainPanel(plotOutput(name))))

server <- function(input, output, session) {
    output[[name]] <- renderPlot(
        plot(inclusion.levels(data$a)$age)
    )
}