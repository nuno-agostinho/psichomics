#' User interface
#' @importFrom shiny column
#' @importFrom shinyBS bsAlert
widgetsUI <- function(id, tab) {
    ns <- NS(id)
    tab("Widgets",
        fluidRow(
            uiOutput(ns("localModal")),
            column(3, h3("Buttons"),
                          actionButton(ns("action"), label = "Action")),
            #submitButton("Submit")), # DOESN'T ALLOW TO UPDATE VALUES...
            column(3, h3("Single checkbox"),
                          checkboxInput(ns("checkbox"), label = "Choice A",
                                        value = TRUE)),
            column(3, checkboxGroupInput(
                ns("checkGroup"), label = h3("Checkbox group"), 
                choices = list("Choice 1" = 1, "Choice 2" = 2,
                               "Choice 3" = 3), selected = 1)),
            column(3, dateInput(ns("date"),
                                       label = h3("Date input"),
                                       value = "2014-01-01"))   
        ),
        fluidRow(
            column(3, dateRangeInput(ns("dates"), label = h3("Date range"))),
            column(3, fileInput(ns("file"), label = h3("File input"))),
            column(3, h3("Help text"),
                          helpText("Note: help text isn't a true widget,", 
                                   "but it provides an easy way to add text to",
                                   "accompany other widgets.")),
            column(3, numericInput(ns("num"), label = h3("Numeric input"),
                                          value = 1))   
        ),
        fluidRow(
            column(3, radioButtons(
                ns("radio"), label = h3("Radio buttons"),
                choices = list("Choice 1" = 1,
                               "Choice 2" = 2,
                               "Choice 3" = 3), selected = 1)),
            column(3, selectizeInput(
                ns("select"), label = h3("Selectize box"),
                choices = list("Primates" = c("Homo sapiens (human)" = 1),
                               "Rodents" = c("Mus musculus (mouse)" = 2),
                               "Carnivora" = c("Canis lupus (wolf)" = 3)),
                # more options: https://github.com/brianreavis/selectize.js/blob/master/docs/usage.md
                options = list(create = TRUE, createOnBlur = TRUE,
                               addPrecedence = TRUE),
                selected = 1, multiple = TRUE, size = 2)),
            column(3,
                          sliderInput(ns("slider1"), label = h3("Sliders"), min = 0,
                                      max = 100, value = 50),
                          sliderInput(ns("slider2"), "", min = 0, max = 100,
                                      value = c(25, 75))),
            column(3, textInput(ns("text"), label = h3("Text input"), 
                                       value = "Enter text..."))   
        )
    )
}

#' @importFrom shinyBS bsModal
widgetsServer <- function(input, output, session) {
    ns <- session$ns
    
    obsB <- observe({
        # Dismiss the alert or it won't update afterwards
        if (input$slider1 > 90) {
            output$localModal <- renderUI(
                bsModal(ns("modalExample"), "Awesome!", NULL, "Over 90!")
            )
        } else if (input$slider1 > 75) {
            output$localModal <- renderUI(
                bsModal(ns("modalExample"), "Great", NULL, "It's something.")
            )
        } else if (input$slider1 < 25) {
            output$localModal <- renderUI(
                bsModal(ns("modalExample"), "Too low", NULL, "Turn it up a bit")
            )
        } else return(NULL)
        toggleModal(session, "modalExample", toggle="open")
    })
}

attr(widgetsUI, "loader") <- "app"
attr(widgetsServer, "loader") <- "app"