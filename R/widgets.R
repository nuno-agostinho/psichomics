name <- "Widgets"
primary <- TRUE
id <- function(value) objectId(name, value)

#' @importFrom shiny column
ui <- function(tab) {
    tab(name,
        fluidRow(
            bsAlert(anchorId = id("alert")),
            column(3, h3("Buttons"),
                          actionButton(id("action"), label = "Action")),
            #submitButton("Submit")), # DOESN'T ALLOW TO UPDATE VALUES...
            column(3, h3("Single checkbox"),
                          checkboxInput(id("checkbox"), label = "Choice A",
                                        value = TRUE)),
            column(3, checkboxGroupInput(
                id("checkGroup"), label = h3("Checkbox group"), 
                choices = list("Choice 1" = 1, "Choice 2" = 2,
                               "Choice 3" = 3), selected = 1)),
            column(3, dateInput(id("date"),
                                       label = h3("Date input"),
                                       value = "2014-01-01"))   
        ),
        fluidRow(
            column(3, dateRangeInput(id("dates"), label = h3("Date range"))),
            column(3, fileInput(id("file"), label = h3("File input"))),
            column(3, h3("Help text"),
                          helpText("Note: help text isn't a true widget,", 
                                   "but it provides an easy way to add text to",
                                   "accompany other widgets.")),
            column(3, numericInput(id("num"), label = h3("Numeric input"),
                                          value = 1))   
        ),
        fluidRow(
            column(3, radioButtons(
                id("radio"), label = h3("Radio buttons"),
                choices = list("Choice 1" = 1,
                               "Choice 2" = 2,
                               "Choice 3" = 3), selected = 1)),
            column(3, selectizeInput(
                id("select"), label = h3("Selectize box"),
                choices = list("Primates" = c("Homo sapiens (human)" = 1),
                               "Rodents" = c("Mus musculus (mouse)" = 2),
                               "Carnivora" = c("Canis lupus (wolf)" = 3)),
                # more options: https://github.com/brianreavis/selectize.js/blob/master/docs/usage.md
                options = list(create = TRUE, createOnBlur = TRUE,
                               addPrecedence = TRUE),
                selected = 1, multiple = TRUE, size = 2)),
            column(3,
                          sliderInput(id("slider1"), label = h3("Sliders"), min = 0,
                                      max = 100, value = 50),
                          sliderInput(id("slider2"), "", min = 0, max = 100,
                                      value = c(25, 75))),
            column(3, textInput(id("text"), label = h3("Text input"), 
                                       value = "Enter text..."))   
        )
    )
}

server <- function(input, output, session) {
    obsB <- observe({
        # Dismiss the alert or it won't update afterwards
        closeAlert(session, alertId = id("exampleAlert"))
        if (input[[id("slider1")]] > 90)
            createAlert(session, anchorId = id("alert"),
                        alertId = id("exampleAlert"),
                        title = "Awesome", content = "Over 90! Fantastic.",
                        style = "success")
        else if (input[[id("slider1")]] > 75)
            createAlert(session, id("alert"), id("exampleAlert"), title = "Great",
                        content = "Above 75!")
        else if (input[[id("slider1")]] < 25)
            createAlert(session, id("alert"), id("exampleAlert"), title = "Too low",
                        content = "Turn it up a bit, please.", style = "danger")
    })
}