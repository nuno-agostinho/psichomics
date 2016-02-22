name <- "Widgets"

ui <- function(tab) {
    tab(name,
        fluidRow(
            bsAlert(anchorId = "alert"),
            shiny::column(3, h3("Buttons"),
                          actionButton("action", label = "Action")),
            #submitButton("Submit")), # DOESN'T ALLOW TO UPDATE VALUES...
            shiny::column(3, h3("Single checkbox"),
                          checkboxInput("checkbox", label = "Choice A",
                                        value = TRUE)),
            shiny::column(3, checkboxGroupInput(
                "checkGroup", label = h3("Checkbox group"), 
                choices = list("Choice 1" = 1, "Choice 2" = 2,
                               "Choice 3" = 3), selected = 1)),
            shiny::column(3, dateInput("date", label = h3("Date input"),
                                       value = "2014-01-01"))   
        ),
        fluidRow(
            shiny::column(3,dateRangeInput("dates", label = h3("Date range"))),
            shiny::column(3,fileInput("file", label = h3("File input"))),
            shiny::column(3, h3("Help text"),
                          helpText("Note: help text isn't a true widget,", 
                                   "but it provides an easy way to add text to",
                                   "accompany other widgets.")),
            shiny::column(3, numericInput("num", label = h3("Numeric input"),
                                          value = 1))   
        ),
        fluidRow(
            shiny::column(3, radioButtons(
                "radio", label = h3("Radio buttons"),
                choices = list("Choice 1" = 1,
                               "Choice 2" = 2,
                               "Choice 3" = 3), selected = 1)),
            shiny::column(3, selectizeInput(
                "select", label = h3("Selectize box"),
                choices = list("Homo sapiens (human)" = 1,
                               "Mus musculus (mouse)" = 2,
                               "Canis lupus (wolf)" = 3),
                selected = 1,
                multiple = T)),
            shiny::column(3,
                          sliderInput("slider1", label = h3("Sliders"), min = 0,
                                      max = 100, value = 50),
                          sliderInput("slider2", "", min = 0, max = 100,
                                      value = c(25, 75))),
            shiny::column(3, textInput("text", label = h3("Text input"), 
                                       value = "Enter text..."))   
        )
    )
}

server <- function(input, output, session) {
    obsB <- observe({
        # Dismiss the alert or it won't update afterwards
        closeAlert(session, alertId = "exampleAlert")
        if (input$slider1 > 90)
            createAlert(session, anchorId = "alert", alertId = "exampleAlert",
                        title = "Awesome", content = "Over 90! Fantastic.",
                        style = "success")
        else if (input$slider1 > 75)
            createAlert(session, "alert", "exampleAlert", title = "Great",
                        content = "Above 75!")
        else if (input$slider1 < 25)
            createAlert(session, "alert", "exampleAlert", title = "Too low",
                        content = "Turn it up a bit, please.", style = "danger")
    })
}