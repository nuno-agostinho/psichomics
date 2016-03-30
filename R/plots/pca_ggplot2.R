# The name used for the plot must be unique
plot <- "PCA ggplot2"
id <- function(value) objectId(name, plot, value)

ui <- list(
    sidebarPanel(
        checkboxGroupInput(id("preprocess"), "Preprocessing",
                           c("Center values" = "center",
                             "Scale values" = "scale"),
                           selected = c("center")),
        actionButton(id("calculate"), "Calculate PCA"),
        uiOutput(id("selectPCA"))
    ), mainPanel(
        plotOutput(id("scatterplot"),
                   hover = hoverOpts("plotHover", delay = 100,
                                     delayType = "throttle")),
        uiOutput(id("hoverInfo"))
    )
)

#' Calculate mouse position according to mouse event
#' 
#' @references 
#' Many thanks to Pawel Papkala for this bit of code. Adapted from
#' http://www.77dev.com/2016/03/custom-interactive-csshtml-tooltips.html
getMousePosition <- function(event) {
    # Calculate point position INSIDE the image as percent of total dimensions
    # from left (horizontal) and from top (vertical)
    left_pct <- (event$x - event$domain$left) / (event$domain$right - event$domain$left)
    top_pct <- (event$domain$top - event$y) / (event$domain$top - event$domain$bottom)
    
    # Calculate distance from left and bottom side of the picture in pixels
    left_px <- event$range$left + left_pct * (event$range$right - event$range$left)
    top_px <- event$range$top + top_pct * (event$range$bottom - event$range$top)
    return(list("left"=left_px, "top"=top_px))
}

#' Create tooltip when interacting with a plot
createTooltip <- function(event, ..., xDeviation = 0, yDeviation = 0) {
    pos <- getMousePosition(event)
    # Style of the tooltip
    style <- paste0("position:absolute; z-index:100; color:white;",
                    "background-color: rgba(0, 0, 0, .7);",
                    "left:", pos$left + xDeviation,
                    "px; top:", pos$top + yDeviation, "px;")
    return(wellPanel(class = "well-sm", style = style, ...))
}

server <- function(input, output, session) {
    ir.pca <- reactiveValues()
    
    observeEvent(input[[id("calculate")]], {
        log.ir <- log(iris[,1:4])
        preprocess <- input[[id("preprocess")]]
        ir.pca$data <- prcomp(log.ir,
                              center = "center" %in% preprocess,
                              scale. = "scale" %in% preprocess)
    })
    
    output[[id("selectPCA")]] <- renderUI({
        if (is.null(ir.pca$data)) 
            return(NULL)
        
        pcs <- colnames(ir.pca$data$x)
        list(
            hr(),
            fluidRow(
                column(4, selectizeInput(id("pcX"), "Choose X axis", choices = pcs)),
                column(4, selectizeInput(id("pcY"), "Choose Y axis", choices = pcs,
                                         selected = pcs[[2]]))),
            highchartOutput(id("variancePlot"))
        )
    })
    
    output[[id("scatterplot")]] <- renderPlot({
        if (is.null(ir.pca$data)) 
            return(NULL)
        
        xAxis <- input[[id("pcX")]]
        yAxis <- input[[id("pcY")]]
        
        if (!is.null(xAxis) & !is.null(yAxis)) {
            perc <- (ir.pca$data$sdev)^2/sum((ir.pca$data$sdev)^2) * 100
            names(perc) <- colnames(ir.pca$data$x)
            
            label <- sprintf("%s (%s%% explained variance)", 
                             names(perc[c(xAxis, yAxis)]), 
                             round(perc[c(xAxis, yAxis)], 2))
            
            ggplot(data.frame(ir.pca$data$x), aes_string(x=xAxis, y=yAxis)) + 
                geom_point() + labs(x=label[1], y=label[2])
        }
    })
    
    output[[id("hoverInfo")]] <- renderUI({
        if (is.null(ir.pca$data)) return(NULL)
        
        xAxis <- input[[id("pcX")]]
        yAxis <- input[[id("pcY")]]
        
        # Get nearest point to hover event within 5 pixels
        hover <- input$plotHover
        point <- nearPoints(data.frame(ir.pca$data$x),
                            hover, threshold = 5, maxpoints = 1, addDist = TRUE)
        if (nrow(point) == 0) return(NULL)
        
        # Create tooltip from hover information
        createTooltip(hover, xDeviation = 19, yDeviation = 2,
                      tags$b(xAxis), point[[xAxis]], br(),
                      tags$b(yAxis), point[[yAxis]])
    })
    
    output[[id("variancePlot")]] <- renderHighchart({
        highchart() %>%
            hc_title(text = "Explained variance by each Principal Component (PC)") %>%
            hc_add_series(name = "PCs", data = (ir.pca$data$sdev)^2,
                          type = "waterfall") %>%
            hc_plotOptions(series = list(dataLabels = list(
                align = "center",
                verticalAlign = "top",
                enabled = TRUE,
                formatter = JS(
                    paste0("function(){ var total = ",
                           sum((ir.pca$data$sdev)^2),
                           ";var perc = (this.y/total) * 100;",
                           "return (Highcharts.numberFormat(this.y) +'<br/>'+",
                           "Highcharts.numberFormat(perc) + '%')}"))))) %>%
            hc_xAxis(categories = colnames(ir.pca$data$x)) %>%
            hc_yAxis(title = list(text = "Explained variance")) %>%
            hc_legend(enabled = FALSE) %>%
            hc_tooltip(pointFormat = '{point.name} {point.y:.5f}') %>%
            hc_exporting(enabled = TRUE,
                         buttons = list(
                             contextButton = list(text = "Export",
                                                  verticalAlign = "bottom", 
                                                  y = -25)))
    })
}