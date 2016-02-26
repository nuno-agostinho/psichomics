# The name used for the plot must be unique
name <- "plot5"

ui <- dataTableOutput(name)

server <- function(input, output, session) {
    output[[name]] <- renderDataTable(
        t(shared.data$data[[1]][["Junction quantification"]][1:10, 1:10])
    )
}