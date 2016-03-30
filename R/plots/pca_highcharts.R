# The name used for the plot must be unique
plot <- "PCA highcharts"
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
        highchartOutput(id("scatterplot"))
    )
)

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
    
    output[[id("scatterplot")]] <- renderHighchart({
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

            df <- data.frame(ir.pca$data$x)
            highchart() %>%
                hc_add_series_scatter(df[[xAxis]], df[[yAxis]]) %>%
                hc_xAxis(title = list(text = label[1])) %>%
                hc_yAxis(title = list(text = label[2]))
        }
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