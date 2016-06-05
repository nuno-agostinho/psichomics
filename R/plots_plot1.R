# The name used for the plot must be unique
plot <- "plot1"
id <- function(value) objectId(name, plot, value)

ui <- function() {
    list(
        sidebarLayout(
            sidebarPanel(
                selectizeInput(id("x"), "Pick x axis", choices = names(mtcars)),
                selectizeInput(id("y"), "Pick y axis", choices = names(mtcars)),
                actionButton(id("change"), "Change to plot2")
            ), 
            mainPanel( plotOutput(id(plot)) )
        )
    )
}

server <- function(input, output, session) {
    output[[id(plot)]] <- renderPlot({
        ggplot(data=mtcars, aes_string(input[[id("x")]], input[[id("y")]])) +
            geom_bin2d()
    })
    
    observeEvent(input[[id("change")]], {
        updateSelectizeInput(session, "Plots_selectizePlot", selected = "plot2")
    })
}