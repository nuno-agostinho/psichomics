## TODO(NunoA): plot using boxplots

#' @importFrom highcharter highchartOutput
diffAnalysisUI <- function(id) {
    ns <- NS(id)
    
    tagList(
        uiOutput(ns("modal")),
        sidebarLayout(
            sidebarPanel(
                numericInput(ns("bandwidth"), "Density bandwidth", 0.01, 
                             step=0.01),
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

#' @importFrom highcharter highchart hc_chart hc_xAxis hc_plotOptions hc_tooltip
#' JS hc_add_series_scatter
prepareDensityPlot <- function(psi, type, bandwidth) {
    # js <- "Highcharts.SVGRenderer.prototype.symbols.line =
    #       function(x, y, w, h) {
    # return ['M', x, y, 'L', x, y]; };"
    # HTML("<script>", js, "</script>")
    
    # Include X-axis zoom and hide markers
    hc <- highchart() %>%
        hc_chart(zoomType = "x") %>%
        hc_xAxis(min = 0, max = 1, title = list(
            text = "Density of exon/intron inclusion levels")) %>%
        hc_plotOptions(series = list(fillOpacity=0.3,
                                     marker = list(enabled = FALSE))) %>%
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
    
    count <- 0
    allRows <- list()
    plotLines <- list()
    for (group in unique(type)) {
        row  <- psi[type == group]
        med  <- roundDigits(median(row, na.rm = TRUE))
        vari <- roundDigits(var(row, na.rm = TRUE))
        max  <- roundDigits(max(row, na.rm = TRUE))
        min  <- roundDigits(min(row, na.rm = TRUE))
        samples <- sum(!is.na(row))
        # Ignore data with low number of data points
        if (sum(!is.na(row)) >= 2) {
            color <- JS("Highcharts.getOptions().colors[", count, "]")
            
            # Calculate the density of inclusion levels for each sample type
            den <- density(row, bw = bandwidth, na.rm = TRUE)
            hc <- hc %>%
                hc_add_series_density(den, name=group, area=TRUE, median=med, 
                                      var=vari, samples=samples, max=max,
                                      min=min, fillColor=color) %>%
                hc_add_series_scatter(row, rep(0, length(row)),
                                      marker=list(enabled=TRUE, symbol="circle",
                                                  radius=4, fillColor=color))
            # Save plot line with information
            plotLines[[count + 1]] <- list(
                label = list(text = paste("Median:", med, "/ Variance:", vari)),
                # Colour the same as the series
                color=color,
                dashStyle="shortdash",
                width=2,
                value=med,
                zIndex = 7)
            allRows[[count + 1]] <- row
            count <- count + 1
        }
    }
    
    # Add plotLines with information
    hc <- hc %>% hc_xAxis(plotLines = plotLines)
    return(list(plot=hc, data=allRows))
}

basicStats <- function(data, len) {
    vari <- vapply(data, var, numeric(1), na.rm = TRUE)
    medi <- vapply(data, median, numeric(1), na.rm = TRUE)
    
    if (len == 2) {
        deltaMedian <- tagList(
            tags$b("|\u0394 Median|: "), abs(medi[2] - medi[1]), br())
    } else {
        deltaMedian <- NULL
    }
    
    ui <- tagList(h4("Basic statistics"), deltaMedian,
                  tags$b("Average variance: "),
                  sum(vari)/length(vari))
    return(ui)
}

wilcox <- function(psi, type) {
    group <- unique(type)
    len <- length(group)
    
    if (len > 2) {
        return(tagList(h4("Wilcoxon rank sum test"),
                       "Can only perform this test on 2 or 1 group."))
    } else if (len == 2) {
        psiA <- psi[type == group[1]]
        psiB <- psi[type == group[2]]
        stat <- tryCatch(list(stat=wilcox.test(psiA, psiB)), 
                         warning=function(w)
                             return(list(stat=wilcox.test(psiA, psiB),
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
}

levene <- function(psi, type) {
    len <- length(unique(type))
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
}

kruskal <- function(psi, type) {
    len <- length(unique(type))
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
}

fisher <- function(psi, type) {
    stat <- try(R.utils::evalWithTimeout(
        fisher.test(psi, factor(type)),
        timeout = 1,
        onTimeout = "error"))
    
    if (class(stat) != "try-error") {
        tagList(
            h4(stat$method),
            tags$b("p-value: "), stat$p.value, br(),
            tags$b("Alternative hypothesis: "), stat$alternative
        )
    } else {
        tagList(
            h4("Fisher's Exact Test for Count Data"),
            "This test took too much to complete!"
        )
    }
}

spearman <- function(psi, type) {
    group <- unique(type)
    len <- length(group)
    
    if (len != 2) {
        tagList(
            h4("Spearman's correlation"),
            "Can only perform this test on 2 groups.")
    } else {
        var <- var(psi[type == group[1]], psi[type == group[2]])
        cor <- cor(psi[type == group[1]], psi[type == group[2]])
        
        tagList(
            h4("Spearman's correlation"),
            tags$b("Variance: "), var, br(),
            tags$b("Correlation: "), cor)
    }
}

#' @importFrom stats kruskal.test median wilcox.test var
#' @importFrom lawstat levene.test
#' @importFrom highcharter renderHighchart
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
        
        bandwidth <- input$bandwidth
        if (bandwidth <= 0 || is.na(bandwidth)) {
            errorModal(session, "Bandwidth must have a positive value",
                       "Insert a number higher than 0.")
            return(NULL)
        }
        
        # Filter groups with less data points than required
        threshold <- 1
        names(psi) <- type
        
        psi <- lapply(unique(type),
                      function(t) {
                          psi <- psi[type == t]
                          if ( sum(!is.na(psi)) >= threshold )
                              return(psi)
                      })
        psi <- unlist(psi)
        type <- names(psi)
        len <- length(unique(type))
        
        prep <- prepareDensityPlot(psi, type, bandwidth)
        output$density <- renderHighchart(prep$plot)
        
        output$basicStats <- renderUI(basicStats(prep$data, len))
        output$wilcox   <- renderUI(wilcox(psi, type))
        output$kruskal  <- renderUI(kruskal(psi, type))
        output$levene   <- renderUI(levene(psi, type))
        # output$fisher   <- renderUI(fisher(psi, type))
        # output$spearman <- renderUI(spearman(psi, type))
    })
}

attr(diffAnalysisUI, "loader") <- "analysis"
attr(diffAnalysisUI, "name") <- "Differential analysis (per splicing event)"
attr(diffAnalysisServer, "loader") <- "analysis"