# The name used for the plot must be unique
plot <- "Differential analysis"
id <- function(value) objectId(name, plot, value)

ui <- tagList(
    sidebarLayout(
        sidebarPanel("Hi there!"),
        mainPanel(
            highchartOutput(id(plot))
        )
    )
)

hc_add_series_density <- function (hc, x, name, area = FALSE, ...) {
    type <- ifelse(area, "areaspline", "spline")
    data <- list.parse3(data.frame(cbind(x = x$x, y = x$y)))
    return(hc %>% hc_add_series(name = name, data = data, type = type))
}

server <- function(input, output, session) {
    observe({
        # Get selected event (if there is any)
        event <- getEvent()
        if (is.null(event) || event == "")
            return(NULL)
        
        # Get event's inclusion levels for all samples
        psi <- getInclusionLevels()
        ids <- names(psi)
        psi <- as.numeric(psi[event, ])
        
        output[[id(plot)]] <- renderHighchart({
            # Include X-axis zoom and hide markers without hovering
            hc <- highchart() %>%
                hc_chart(zoomType = "x") %>%
                hc_xAxis(min = 0, max = 1, title = list(
                    text = "Density of exon/intron inclusion levels")) %>%
                hc_plotOptions(series = list(marker = list(enabled = FALSE)))
            
            # Separate samples by their type
            typeList <- readRDS("data/TCGAsampleType.RDS")
            type <- gsub(".*?-([0-9]{2}).-.*", "\\1", ids, perl = TRUE)
            type <- typeList[type]
            
            for (group in unique(type)) {
                row <- psi[type == group]
                # Ignore data with low number of data points
                if (sum(!is.na(row)) >= 2) {
                    # Calculate the density of inclusion levels for each sample type
                    den <- density(row, na.rm = TRUE)
                    hc <- hc %>%
                        hc_add_series_density(den, group, area = TRUE)
                }
            }
            return(hc)
        })
    })
}