plots1UI <- function(id) {
    ns <- NS(id)
    tagList(
        sidebarLayout(
            sidebarPanel(
                selectizeInput(ns("x"), "Pick x axis", choices = names(mtcars)),
                selectizeInput(ns("y"), "Pick y axis", choices = names(mtcars))
                # actionButton(ns("change"), "Change to plot2")
            ), 
            mainPanel( plotOutput(ns("plot")) )
        )
    )
}

plots1Server <- function(input, output, session) {
    output$plot <- renderPlot({
        plot(mtcars[[input$x]], mtcars[[input$y]])
        # ggplot(data=mtcars, aes_string(input$x, input$y)) + geom_bin2d()
    })
    
    # observeEvent(input$change, {
    #     updateSelectizeInput(session, "Plots_selectizePlot", selected="plot2")
    # })
}

attr(plots1UI, "loader") <- "plots"
attr(plots1UI, "name") <- "Plots 1"
attr(plots1Server, "loader") <- "plots"