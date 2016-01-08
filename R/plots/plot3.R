# The name used for the plot must be unique
name <- "plot3"

ui <- list(sidebarLayout(
    sidebarPanel(selectizeInput("colours", "Pick a colour",
                                choices = list
                                ("red", "orange", "blue", "pink", "green"))), 
    mainPanel(plotOutput(name))))

server <- function(input, output, session) {
    output[[name]] <- renderPlot(
        ggplot(data=h, aes(Region, Life.Expectancy)) +
            geom_bar(stat = "identity", colour = "orange")
    )
}