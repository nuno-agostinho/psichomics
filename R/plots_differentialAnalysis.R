## TODO(NunoA): plot using boxplots

#' @importFrom highcharter highchartOutput
diffAnalysisUI <- function(id) {
    ns <- NS(id)
    
    tagList(
        sidebarLayout(
            sidebarPanel(
                numericInput(ns("bandwidth"), "Density bandwidth", 0.01, step=0.01),
                h3("Non-parametric tests"),
                uiOutput(ns("basicStats")), hr(),
                # uiOutput(ns("spearman")), hr(),
                # uiOutput(ns("fisher")), hr(),
                uiOutput(ns("wilcox")), hr(),
                uiOutput(ns("kruskal")), hr(),
                uiOutput(ns("levene"))
            ), mainPanel(
                highchartOutput(ns("density"))
            )
        )
    )
}

#' @importFrom stats kruskal.test median wilcox.test var
#' @importFrom lawstat levene.test
diffAnalysisServer <- function(input, output, session) {
    observe({
        # Get selected event (if there is any)
        event <- getEvent()
        if (is.null(event) || event == "")
            return(NULL)
        
        # Get event's inclusion levels for all samples
        psi <- getInclusionLevels()
        ids <- names(psi)
        psi <- as.numeric(psi[event, ])
        
        # Separate samples by their type
        type <- getSampleTypes(ids)
        
        group <- unique(type)
        len <- length(group)
        
        # output$fisher <- renderUI({
        #     stat <- try(R.utils::evalWithTimeout(
        #         fisher.test(psi, factor(type)), 
        #         timeout = 1, 
        #         onTimeout = "error"))
        #     
        #     if (class(a) != "try-error") {
        #         tagList(
        #             h4(stat$method),
        #             tags$b("p-value: "), stat$p.value, br(),
        #             tags$b("Alternative hypothesis: "), stat$alternative
        #         )
        #     } else {
        #         tagList(
        #             h4("Fisher's Exact Test for Count Data"),
        #             "This test took too much to complete!"
        #         )
        #     }
        # })
        
        # output$spearman <- renderUI({
        #     group <- unique(type)
        #     len <- length(group)
        #     
        #     if (len != 2) {
        #         tagList(
        #             h4("Spearman's correlation"),
        #             "Can only perform this test on 2 groups.")
        #     } else {
        #         var <- var(psi[type == group[1]], psi[type == group[2]])
        #         cor <- cor(psi[type == group[1]], psi[type == group[2]])
        #         
        #         tagList(
        #             h4("Spearman's correlation"),
        #             tags$b("Variance: "), var, br(),
        #             tags$b("Correlation: "), cor)
        #     }
        # })
        
        output$wilcox <- renderUI({
            if (len > 2) {
                return(tagList(h4("Wilcoxon rank sum test"),
                               "Can only perform this test on 2 or 1 group."))
            } else if (len == 2) {
                stat <- tryCatch(list(stat=wilcox.test(psi[type == group[1]],
                                                       psi[type == group[2]])), 
                                 warning=function(w)
                                     return(list(
                                         stat=wilcox.test(psi[type == group[1]],
                                                          psi[type == group[2]]),
                                         warning=w)))
            } else if (len == 1) {
                stat <- tryCatch(list(stat=wilcox.test(psi)), warning=function(w)
                    return(list(stat=wilcox.test(psi), warning=w)))
            }
            
            if ("warning" %in% names(stat))
                warn <- tagList(
                    tags$code(paste("Warning:", stat$warning$message)), br())
            else
                warn <- NULL
            
            tagList(
                h4(stat$stat$method), warn,
                tags$b("Test value: "), stat$stat$statistic, br(),
                tags$b("p-value: "), stat$stat$p.value, br(),
                tags$b("Test parameters: "), stat$stat$parameter, br(),
                tags$b("Location parameter: "), stat$stat$null.value, br(),
                tags$b("Alternative hypothesis: "), stat$stat$alternative
            )
        })
        
        output$kruskal <- renderUI({
            if (len >= 2) {
                stat <- kruskal.test(psi, factor(type))
                tagList(
                    h4(stat$method),
                    tags$b("Test value (Chi squared): "), stat$statistic, br(),
                    tags$b("p-value: "), stat$p.value, br(),
                    tags$b("Degrees of freedom: "), stat$parameter
                )
            } else {
                tagList(
                    h4("Kruskal test"),
                    "Can only perform this test on 2 or more groups."
                )
            }
        })
        
        output$levene <- renderUI({
            if (len >= 2) {
                nas <- is.na(psi)
                stat <- levene.test(psi[!nas], factor(type[!nas]))
                tagList(
                    h4("Levene's Test for Homogeneity of Variance"),
                    tags$b("Test value: "), stat$statistic, br(),
                    tags$b("p-value: "), stat$p.value, br(),
                    tags$b("p-value without bootstrap: "), 
                    stat$non.bootstrap.p.value
                )
            } else {
                tagList(
                    h4("Levene's Test for Homogeneity of Variance"),
                    "Can only perform this test on 2 or more groups."
                )
            }
        })
        
        bandwidth <- input$bandwidth
        if (bandwidth <= 0 || is.na(bandwidth)) {
            errorModal(session, "Bandwidth must have a positive value",
                       "Insert a number higher than 0.")
            return(NULL)
        }
        
        output$density <- renderHighchart({
            # Include X-axis zoom and hide markers without hovering
            hc <- highchart() %>%
                hc_chart(zoomType = "x") %>%
                hc_xAxis(min = 0, max = 1, title = list(
                    text = "Density of exon/intron inclusion levels")) %>%
                hc_plotOptions(series = list(marker = list(enabled = FALSE)))
            
            count <- 0
            allRows <- list()
            plotLines <- list()
            for (group in unique(type)) {
                row <- psi[type == group]
                med <- roundDigits(median(row, na.rm = TRUE))
                var <- roundDigits(var(row, na.rm = TRUE))
                max <- roundDigits(max(row, na.rm = TRUE))
                min <- roundDigits(min(row, na.rm = TRUE))
                samples <- sum(!is.na(row))
                # Ignore data with low number of data points
                if (sum(!is.na(row)) >= 2) {
                    # Calculate the density of inclusion levels for each sample type
                    den <- density(row, bw = bandwidth, na.rm = TRUE)
                    hc <- hc %>%
                        hc_add_series_density(den, name=group, area=TRUE,
                                              median=med, var=var,
                                              samples=samples,
                                              max=max, min=min)
                    # Save plot line with information
                    plotLines[[count + 1]] <- list(
                        label = list(text = paste("Median:", med,
                                                  "/ Variance:", var)),
                        # Colour the same as the series
                        color=JS("Highcharts.getOptions().colors[",
                                 count, "]"),
                        dashStyle="shortdash",
                        width=2,
                        value=med,
                        zIndex = 7)
                    allRows[[count + 1]] <- row
                    count <- count + 1
                }
            }
            
            # Add plotLines with information
            hc <- hc %>% 
                hc_xAxis(plotLines = plotLines) %>%
                hc_tooltip(
                    headerFormat = paste(
                        span(style="color:{point.color}", "\u25CF "),
                        tags$b("{series.name}"), br()),
                    pointFormat = paste(
                        "Inclusion level: {point.x}", br(),
                        "Number of samples: {series.options.samples}", br(),
                        "Median: {series.options.median}", br(),
                        "Variance: {series.options.var}", br(),
                        "Range: {series.options.min} - {series.options.max}"))
            
            output$basicStats <- renderUI ({
                var <- vapply(allRows, var, numeric(1), na.rm = TRUE)
                med <- vapply(allRows, median, numeric(1), na.rm = TRUE)
                
                if (len == 2) {
                    deltaMedian <- tagList(
                        tags$b("|\u0394 Median|: "), abs(med[2] - med[1]), br())
                } else {
                    deltaMedian <- NULL   
                }
                
                tagList(h4("Basic statistics"), deltaMedian,
                        tags$b("Average variance: "), sum(var)/length(var))
            })
            return(hc)
        })
    })
}

attr(diffAnalysisUI, "loader") <- "plots"
attr(diffAnalysisUI, "name") <- "Differential analysis"
attr(diffAnalysisServer, "loader") <- "plots"