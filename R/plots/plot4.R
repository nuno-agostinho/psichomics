# The name used for the plot must be unique
plot <- "plot4"
id <- function(value) objectId(name, plot, value)

ui <- list(
    sidebarLayout(
        sidebarPanel(
            selectizeInput(
                id("colours"),
                "Pick a colour",
                choices = list("red", "orange", "blue", "pink", "green")
            )
        ), 
        mainPanel( plotOutput(id(plot)) )
    )
)

server <- function(input, output, session) {
    output[[id(plot)]] <- renderPlot({
        f <- rep(NA, nrow(mtcars))
        f[match(input[[id("selectizeEvent")]], rownames(mtcars))] <- "color"
        ggplot(mtcars, aes(x = rownames(mtcars), y = hp, fill = f)) + 
            geom_bar(stat = "identity",
                     colour = input[[id("colours")]],
                     show.legend = FALSE)
    })
}