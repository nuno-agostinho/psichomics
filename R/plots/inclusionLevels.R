# The name used for the plot must be unique
name <- "Inclusion levels"

ui <- list(
    dataTableOutput(name)
)

createLink <- function(val) {
    a <- '<a id="%s" href="#">%s</button>'
    sprintf(a, gsub(" ", "_", val), val)
    #actionButton(val[1], val[1])
}

server <- function(input, output, session) {
    output[[name]] <- renderDataTable({
        cbind(cars = createLink(rownames(mtcars)), mtcars)
    }, escape = FALSE, options = list(pageLength = 10))
    
    lapply(gsub(" ", "_", rownames(mtcars)), function(i) {
        onclick(i, {
            updateSelectizeInput(session, "selectizePlot", selected = "plot4")
            updateSelectizeInput(session, "selectizeEvent",
                                 selected = gsub("_", " ", i))
        })
    })
}