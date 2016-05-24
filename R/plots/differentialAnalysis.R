## TODO(NunoA): plot using boxplots

# The name used for the plot must be unique
plot <- "Differential analysis"
id <- function(value) objectId(name, plot, value)

ui <- tagList(
    sidebarLayout(
        sidebarPanel(
            h3("Non-parametric tests"),
            # uiOutput(id("spearman")), hr(),
            # uiOutput(id("fisher")), hr(),
            uiOutput(id("wilcox")), hr(),
            uiOutput(id("kruskal")), hr(),
            uiOutput(id("levene"))
        ), mainPanel(
            highchartOutput(id("density"))
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
        
        # Separate samples by their type
        typeList <- readRDS("data/TCGAsampleType.RDS")
        type <- gsub(".*?-([0-9]{2}).-.*", "\\1", ids, perl = TRUE)
        type <- typeList[type]
        
        group <- unique(type)
        len <- length(group)
        
        # output[[id("fisher")]] <- renderUI({
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
        
        # output[[id("spearman")]] <- renderUI({
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
        
        output[[id("wilcox")]] <- renderUI({
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
                tags$b("Alternative hypothesis: "), stat$stat$alternative, br()
            )
        })
        
        output[[id("kruskal")]] <- renderUI({
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
        
        output[[id("levene")]] <- renderUI({
            if (len >= 2) {
                stat <- car::leveneTest(psi, factor(type))
                tagList(
                    h4("Levene's Test for Homogeneity of Variance"),
                    HTML(tooltip_table(names(stat), stat))
                )
            } else {
                tagList(
                    h4("Levene's Test for Homogeneity of Variance"),
                    "Can only perform this test on 2 or more groups."
                )
            }
        })
        
        output[[id("density")]] <- renderHighchart({
            # Include X-axis zoom and hide markers without hovering
            hc <- highchart() %>%
                hc_chart(zoomType = "x") %>%
                hc_xAxis(min = 0, max = 1, title = list(
                    text = "Density of exon/intron inclusion levels")) %>%
                hc_plotOptions(series = list(marker = list(enabled = FALSE)))
            
            count <- 0
            for (group in unique(type)) {
                row <- psi[type == group]
                med <- round(median(row, na.rm = TRUE), 2)
                # Ignore data with low number of data points
                if (sum(!is.na(row)) >= 2) {
                    # Calculate the density of inclusion levels for each sample type
                    den <- density(row, na.rm = TRUE)
                    hc <- hc %>%
                        hc_add_series_density(den, group, area = TRUE) %>%
                        hc_xAxis(plotLines = list(
                            list(label = list(text = paste(
                                "Median:", med)),
                                # Colour the same as the series
                                color=JS("Highcharts.getOptions().colors[",
                                         count, "]"),
                                dashStyle="shortdash",
                                width=2,
                                value=med,
                                zIndex = 7)))
                    count <- count + 1
                }
            }
            return(hc)
        })
    })
}