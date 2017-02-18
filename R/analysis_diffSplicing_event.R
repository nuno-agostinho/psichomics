## TODO(NunoA): plot using boxplots

#' Interface for the analysis of an alternative splicing event
#' @param id Character: identifier
#' @importFrom highcharter highchartOutput
#' @importFrom shiny tagList uiOutput NS sidebarLayout numericInput h3 mainPanel
#' actionButton
#' @return Character with the HTML interface
diffSplicingEventUI <- function(id) {
    ns <- NS(id)
    
    card <- function(id) {
        div(class="col-sm-6 col-md-4",
            div(class="thumbnail", style="background:#eee;",
                div(class="caption", uiOutput(ns(id)))))
    }
    
    # Take user to the survival analysis by PSI cut-off
    survival <- tagList(
        actionButton(ns("optimalSurv1"), onclick="showSurvCutoff()",
                     "Survival analysis by PSI cut-off", 
                     class="btn-info btn-md btn-block",
                     class="visible-lg visible-md"),
        actionButton(ns("optimalSurv2"), onclick="showSurvCutoff()",
                     "Survival analysis by PSI cut-off", 
                     class="btn-info btn-xs btn-block",
                     class="visible-sm visible-xs"))
    
    tagList(
        uiOutput(ns("modal")),
        sidebarLayout(
            sidebarPanel(
                selectGroupsUI(ns("diffGroups"),
                               label="Groups of samples to analyse",
                               noGroupsLabel="All samples as one group",
                               groupsLabel="Samples by selected groups"),
                numericInput(ns("bandwidth"), "Density smoothing bandwidth",
                             0.01, step=0.01, min=0.01),
                actionButton(ns("analyse"), "Perform analyses",
                             class="btn-primary"),
                uiOutput(ns("basicStats")), hr(),
                survival
            ), mainPanel(
                highchartOutput(ns("density")),
                h4("Parametric tests"),
                div(class="row",
                    card("ttest"),
                    card("levene")),
                h4("Non-parametric tests"),
                div(class="row",
                    card("wilcox"),
                    card("kruskal"),
                    card("fligner")))
        )
    )
}

#' Plot distribution through a density plot
#' 
#' The tooltip shows the median, variance, max, min and number of non-NA samples
#' of each data series.
#' 
#' @inheritParams wilcox
#' @param rug Boolean: include rug plot to better visualise data distribution 
#' (TRUE by default)
#' @param vLine Boolean: include vertical plot lines to indicate the mean and
#' median of each group even when those groups are omitted
#' @param ... Extra parameters passed to \code{density} to create the kernel
#' density estimates
#' @param title Character: plot title
#' 
#' @importFrom highcharter highchart hc_chart hc_xAxis hc_plotOptions hc_tooltip
#' JS
#' @importFrom stats median var density
#' 
#' @return Highcharter object with density plot
#' @export
#' @examples
#' data <- sample(20, rep=TRUE)/20
#' groups <- c(rep("A", 10), rep("B", 10))
#' plotDistribution(data, groups)
plotDistribution <- function(psi, groups, rug=TRUE, vLine=TRUE, ..., 
                             title=NULL) {
    # Include X-axis zoom and hide markers
    hc <- highchart() %>%
        hc_chart(zoomType = "x") %>%
        hc_xAxis(min = 0, max = 1, title = list(
            text = "Distribution of PSI values")) %>%
        hc_plotOptions(series = list(fillOpacity=0.3,
                                     marker = list(enabled = FALSE))) %>%
        hc_tooltip(
            headerFormat = paste(
                span(style="color:{point.color}", "\u25CF "),
                tags$b("{series.name}"), br()),
            pointFormat = paste(
                "Inclusion level: {point.x:.2f}", br(),
                "Number of samples: {series.options.samples}", br(),
                "Median: {series.options.median}", br(),
                "Variance: {series.options.var}", br(),
                "Range: {series.options.min} - {series.options.max}")) %>%
        export_highcharts()
    
    if (!is.null(title)) hc <- hc %>% hc_title(text=title)
    
    count <- 0
    plotLines <- list()
    for (group in sort(unique(groups))) {
        row  <- psi[groups == group]
        med  <- roundDigits(median(row, na.rm=TRUE))
        vari <- roundDigits(var(row, na.rm=TRUE))
        max  <- roundDigits(max(row, na.rm=TRUE))
        min  <- roundDigits(min(row, na.rm=TRUE))
        samples <- sum(!is.na(row))
        
        color <- JS("Highcharts.getOptions().colors[", count, "]")
        
        # Calculate the density of inclusion levels for each sample group
        den <- density(row, na.rm=TRUE, ...)
        hc <- hc %>%
            hc_add_series(den, type="area", name=group, median=med, var=vari,
                          samples=samples, max=max, color=color, min=min)
        # Rug plot
        if (rug) {
            hc <- hc_scatter(
                hc, row, rep(0, length(row)), name=group, marker=list(
                    enabled=TRUE, symbol="circle", radius=4, fillColor=color),
                median=med, var=vari, samples=samples, max=max, min=min)
        }
        # Save plot line with information
        if (vLine) {
            plotLines[[count + 1]] <- list(
                label = list(text = paste("Median:", med, "/ Variance:", vari)),
                # Colour the same as the series
                color=color,
                dashStyle="shortdash",
                width=2,
                value=med,
                zIndex = 7)
        }
        count <- count + 1
    }
    
    # Add plotLines with information
    if (vLine) hc <- hc %>% hc_xAxis(plotLines = plotLines)
    return(hc)
}

#' Basic statistics performed on data
#' 
#' Variance and median of each group. If data has 2 groups, also calculates the
#' delta variance and delta median.
#' 
#' @inheritParams wilcox
#' @importFrom shiny tagList br h4
#' 
#' @return HTML elements
basicStats <- function(psi, groups) {
    data <- lapply(unique(groups), function(g) psi[groups == g])
    
    len <- length(unique(groups))
    vari <- vapply(data, var, numeric(1), na.rm = TRUE)
    medi <- vapply(data, median, numeric(1), na.rm = TRUE)
    
    if (len == 2) {
        deltaMedian <- tagList(tags$b("|\u0394 Median|: "), 
                               roundDigits(abs(medi[2] - medi[1])), br())
        deltaVar <- tagList(tags$b("|\u0394 Variance|: "), 
                            roundDigits(abs(vari[2] - vari[1])), br())
    } else {
        deltaMedian <- NULL
        deltaVar <- NULL
    }
    
    avgMedian <- roundDigits( mean(medi) )
    avgVar <- roundDigits( mean(vari) )
    ui <- tagList(hr(), h4("Basic statistics"),
                  tags$b("Average median: "), avgMedian, br(), deltaMedian,
                  tags$b("Average variance: "), avgVar, br(), deltaVar)
    return(ui)
}

#' Perform Wilcoxon analysis and return interface to show the results
#' @param psi Numeric: quantification of one alternative splicing event
#' @param groups Character: group of each PSI index
#' @param stat Data frame or matrix: values of the analyses to be performed (if
#' NULL, the analyses will be performed)
#' 
#' @importFrom shiny tagList tags h4 br
#' @importFrom stats wilcox.test
#' @importFrom R.utils capitalize
#' 
#' @return HTML elements
wilcox <- function(psi, groups, stat=NULL) {
    warn <- NULL
    group <- unique(groups)
    len <- length(group)
    
    p.value <- NULL
    if (!is.null(stat)) {
        method      <- stat$`Wilcoxon method`
        statistic   <- stat$`Wilcoxon statistic`
        p.value     <- stat$`Wilcoxon p-value`
        null.value  <- stat$`Wilcoxon null value`
        alternative <- stat$`Wilcoxon alternative`
    }
            
    if (len != 2) {
        return(tagList(h4("Wilcoxon test"),
                       "Can only perform this test on 2 groups."))
    } else if (!is.null(p.value)) {
        adjusted <- grep("Wilcoxon p-value \\(.* adjusted\\)", colnames(stat), 
                         value=TRUE)
        if (length(adjusted) != 0) {
            adjustMethod <- gsub(".*\\((.* adjusted)\\).*", "\\1", adjusted)
            adjusted <- stat[ , adjusted]
            label <- paste0("p-value (", adjustMethod, "): ")
            adjusted <- tagList(tags$b(label), signifDigits(adjusted), br())
        } else {
            adjusted <- NULL
        }
    } else {
        psiA <- psi[groups == group[1]]
        psiB <- psi[groups == group[2]]
        stat <- tryCatch(list(stat=wilcox.test(psiA, psiB)), 
                         warning=function(w)
                             return(list(stat=wilcox.test(psiA, psiB),
                                         warning=w)))
        
        if ("warning" %in% names(stat))
            warn <- tags$div(class="alert alert-warning", role="alert",
                             capitalize(stat$warning$message))
        
        method      <- stat$stat$method
        statistic   <- stat$stat$statistic
        p.value     <- stat$stat$p.value
        adjusted    <- NULL
        null.value  <- stat$stat$null.value
        alternative <- stat$stat$alternative
    }
    
    tagList(
        h4(method), warn,
        tags$b("Test value: "), roundDigits(statistic), br(),
        tags$b("Location parameter: "), null.value, br(),
        tags$b("Alternative hypothesis: "), alternative,
        div(style="text-align:right",
            tags$b("p-value: "), signifDigits(p.value), br(), adjusted))
}

#' Perform unpaired t-test analysis and return interface to show the results
#' @param psi Numeric: quantification of one alternative splicing event
#' @param groups Character: group of each PSI index
#' @param stat Data frame or matrix: values of the analyses to be performed (if
#' NULL, the analyses will be performed)
#' 
#' @importFrom shiny tagList tags h4 br
#' @importFrom stats t.test
#' @importFrom R.utils capitalize
#' 
#' @return HTML elements
ttest <- function(psi, groups, stat=NULL) {
    warn <- NULL
    group <- unique(groups)
    len <- length(group)
 
    p.value <- NULL
    if (!is.null(stat)) {
        method      <- stat$`T-test method`
        statistic   <- stat$`T-test statistic`
        p.value     <- stat$`T-test p-value`
        null.value  <- stat$`T-test null value`
        alternative <- stat$`T-test alternative`
        parameter   <- stat$`T-test parameter`
        int1        <- stat$`T-test conf int1`
        int2        <- stat$`T-test conf int2`
    }
    
    if (len != 2) {
        return(tagList(h4("Unpaired t-test"),
                       "Can only perform this test on 2 groups."))
    } else if (!is.null(p.value)) {
        adjusted <- grep("T-test p-value \\(.* adjusted\\)", colnames(stat), 
                         value=TRUE)
        if (length(adjusted) != 0) {
            adjustMethod <- gsub(".*\\((.* adjusted)\\).*", "\\1", adjusted)
            adjusted <- stat[ , adjusted]
            label <- paste0("p-value (", adjustMethod, "): ")
            adjusted <- tagList(tags$b(label), signifDigits(adjusted), br())
        } else {
            adjusted <- NULL
        }
    } else {
        psiA <- psi[groups == group[1]]
        psiB <- psi[groups == group[2]]
        stat <- tryCatch(list(stat=t.test(psiA, psiB)), 
                         warning=function(w)
                             return(list(stat=t.test(psiA, psiB),
                                         warning=w)),
                         error=return)
        if (is(stat, "error")) {
            message <- stat$message
            check <- "not enough '%s' observations"
            checkX <- sprintf(check, "x")
            checkY <- sprintf(check, "y")
            
            fewObservations <- function(name)
                tagList("Not enough observations in group", tags$b(name),
                        "to perform this statistical test.")
            
            if (message == checkX)
                message <- fewObservations(group[1])
            else if (message == checkY)
                message <- fewObservations(group[2])
            else
                message <- capitalize(message)
            
            error <- tagList(h4("t-test"), tags$div(class="alert alert-danger",
                                                    role="alert", message))
            return(error)
        }
        
        if ("warning" %in% names(stat))
            warn <- tags$div(class="alert alert-warning", role="alert",
                             capitalize(stat$warning$message))
        
        method      <- stat$stat$method
        statistic   <- stat$stat$statistic
        p.value     <- stat$stat$p.value
        adjusted    <- NULL
        null.value  <- stat$stat$null.value
        alternative <- stat$stat$alternative
        parameter   <- stat$stat$parameter
        int1        <- stat$stat$conf.int[[1]]
        int2        <- stat$stat$conf.int[[2]]
    }
    
    tagList(
        h4(method), warn,
        tags$b("Test value: "), roundDigits(statistic), br(),
        tags$b("Test parameter: "), parameter, br(),
        tags$b("Difference in means: "), null.value, br(),
        tags$b("Alternative hypothesis: "), alternative, br(),
        tags$b("95\u0025 confidence interval: "), roundDigits(int1),
        roundDigits(int2),
        div(style="text-align:right",
            tags$b("p-value: "), signifDigits(p.value), br(), adjusted))
}


#' Perform Levene's test and return interface to show the results
#' @inheritParams wilcox
#' @importFrom shiny tagList tags h4 br
#' @return HTML elements
levene <- function(psi, groups, stat=NULL) {
    p.value <- NULL
    if (!is.null(stat)) {
        statistic <- stat$`Levene statistic`
        p.value   <- stat$`Levene p-value`
        non.bootstrap.p.value <- stat$`Levene non bootstrap p-value`
    }
    
    len <- length(unique(groups))
    if (len < 2) {
        return(tagList(h4("Levene's Test for Homogeneity of Variances"),
                       "Can only perform this test on 2 or more groups."))
    } else if (!is.null(p.value)) {
        adjusted <- grep("Levene .*p-value \\(.* adjusted\\)", colnames(stat), 
                         value=TRUE)
        if (length(adjusted) == 2) {
            adjustMethod <- gsub(".*\\((.* adjusted)\\).*", "\\1", adjusted)
            adjusted <- stat[ , adjusted]
            label1 <- paste0("p-value (", adjustMethod[[1]], "): ")
            label2 <- paste0("p-value without bootstrap (", adjustMethod[[2]],
                             "): ")
            adjustedNonBootstrap <- tagList(br(), tags$b(label2), 
                                            signifDigits(adjusted[[1]]), br())
            adjusted <- tagList(tags$b(label1), signifDigits(adjusted[[2]]), 
                                br())
        } else if (length(adjusted) == 1) {
            adjustMethod <- gsub(".*\\((.* adjusted)\\).*", "\\1", adjusted)
            adjusted <- stat[ , adjusted]
            label <- paste0("p-value (", adjustMethod, "): ")
            adjusted <- tagList(tags$b(label), signifDigits(adjusted), br())
        } else {
            adjusted <- NULL
            adjustedNonBootstrap <- NULL
        }
    } else {
        nas <- is.na(psi)
        stat <- leveneTest(psi[!nas], factor(groups[!nas]))
        statistic <- stat$statistic
        p.value   <- stat$p.value
        adjusted  <- NULL
        non.bootstrap.p.value <- stat$non.bootstrap.p.value
        adjustedNonBootstrap  <- NULL
    }
    
    if (!is.null(non.bootstrap.p.value)) {
        nonBootstrap <- tagList(tags$b("p-value without bootstrap: "), 
                                signifDigits(non.bootstrap.p.value), 
                                adjustedNonBootstrap)
    } else {
        nonBootstrap <- NULL
    }
        
    tagList(
        h4("Levene's Test for Homogeneity of Variance"),
        tags$b("Test value: "), roundDigits(statistic), br(),
        div(style="text-align:right",
            tags$b("p-value: "), signifDigits(p.value), br(), adjusted,
            nonBootstrap))
}

#' Perform Fligner-Killeen test and return interface to show the results
#' @inheritParams wilcox
#' @importFrom shiny tagList tags h4 br
#' @importFrom stats fligner.test
#' @return HTML elements
fligner <- function(psi, groups, stat=NULL) {
    len <- length(unique(groups))
    
    p.value <- NULL
    if (!is.null(stat)) {
        statistic <- stat$`Fligner-Killeen statistic`
        p.value   <- stat$`Fligner-Killeen p-value`
        parameter <- stat$`Fligner-Killeen parameter`
    }
    
    if (len < 2) {
        return(tagList(h4("Fligner-Killeen Test for Homogeneity of Variances"),
                       "Can only perform this test on 2 or more groups."))
    } else if (!is.null(p.value)) {
        adjusted <- grep("Fligner-Killeen .*p-value \\(.* adjusted\\)", 
                         colnames(stat), value=TRUE)
        if (length(adjusted) != 0) {
            adjustMethod <- gsub(".*\\((.* adjusted)\\).*", "\\1", adjusted)
            adjusted <- stat[ , adjusted]
            label <- paste0("p-value (", adjustMethod, "): ")
            adjusted <- tagList(tags$b(label), signifDigits(adjusted), br())
        } else {
            adjusted <- NULL
        }
    } else {
        nas <- is.na(psi)
        stat <- fligner.test(psi[!nas], factor(groups[!nas]))
        statistic <- stat$statistic
        p.value   <- stat$p.value
        adjusted  <- NULL
        parameter <- stat$parameter
    }
    
    tagList(
        h4("Fligner-Killeen's Test for Homogeneity of Variance"),
        tags$b("Test value: "), roundDigits(statistic), br(),
        tags$b("Test parameter: "), parameter, br(),
        div(style="text-align:right",
            tags$b("p-value: "), signifDigits(p.value), br(), adjusted))
}

#' Perform Kruskal's test and return interface to show the results
#' @inheritParams wilcox
#' @importFrom shiny tagList tags h4 br
#' @importFrom stats kruskal.test
#' @return HTML elements
kruskal <- function(psi, groups, stat=NULL) {
    len <- length(unique(groups))
    
    p.value <- NULL
    if (!is.null(stat)) {
        method    <- stat$`Kruskal method`
        statistic <- stat$`Kruskal statistic`
        p.value   <- stat$`Kruskal p-value`
        parameter <- stat$`Kruskal parameter`
    }
    
    if (len < 2) {
        return(tagList(h4("Kruskal test"),
                       "Can only perform this test on 2 or more groups."))
    } else if (!is.null(p.value)) {
        adjusted <- grep("Kruskal p-value \\(.* adjusted\\)", colnames(stat), 
                         value=TRUE)
        if (length(adjusted) != 0) {
            adjustMethod <- gsub(".*\\((.* adjusted)\\).*", "\\1", adjusted)
            adjusted <- stat[ , adjusted]
            label <- paste0("p-value (", adjustMethod, "): ")
            adjusted <- tagList(tags$b(label), signifDigits(adjusted), br())
        } else {
            adjusted <- NULL
        }
    } else {
        stat      <- kruskal.test(psi, factor(groups))
        method    <- stat$method
        statistic <- stat$statistic
        p.value   <- stat$p.value
        parameter <- stat$parameter
        adjusted  <- NULL
    }
    
    tagList(h4(method),
            tags$b("Test value \u03C7\u00B2: "), roundDigits(statistic), br(),
            tags$b("Degrees of freedom: "), parameter,
            div(style="text-align:right",
                tags$b("p-value: "), signifDigits(p.value), br(), adjusted))
}

#' Perform Fisher's exact test and return interface to show the results
#' @inheritParams wilcox
#' 
#' @importFrom shiny tagList tags h4 br
#' @importFrom stats fisher.test
#' @importFrom R.utils evalWithTimeout
#' 
#' @return HTML elements
fisher <- function(psi, groups) {
    stat <- try(evalWithTimeout(
        fisher.test(psi, factor(groups)),
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

#' Perform Spearman's test and return interface to show the results
#' @inheritParams wilcox
#' @importFrom shiny tagList tags h4 br
#' @importFrom stats var cor
#' @return HTML elements
spearman <- function(psi, groups) {
    group <- unique(groups)
    len <- length(group)
    
    if (len != 2) {
        tagList(
            h4("Spearman's correlation"),
            "Can only perform this test on 2 groups.")
    } else {
        var <- var(psi[groups == group[1]], psi[groups == group[2]])
        cor <- cor(psi[groups == group[1]], psi[groups == group[2]])
        
        tagList(
            h4("Spearman's correlation"),
            tags$b("Variance: "), var, br(),
            tags$b("Correlation: "), cor)
    }
}

#' Filter groups with less data points than the threshold
#' 
#' Groups containing a number of non-missing values less than the threshold are
#' discarded.
#' 
#' @param vector Unnamed elements
#' @param group Character: group of the elements
#' @param threshold Integer: number of valid non-missing values by group
#' 
#' @return Named vector with filtered elementes from valid groups. The group of 
#' the respective element is given in the name.
#' @export
#' @examples 
#' # Removes groups with less than two elements
#' filterGroups(1:4, c("A", "B", "B", "D"), threshold=2)
filterGroups <- function(vector, group, threshold=1) {
    names(vector) <- group
    vector <- lapply(unique(group), function(t) {
        vector <- vector[group == t]
        if ( sum(!is.na(vector)) >= threshold )
            return(vector)
    })
    return(unlist(vector))
}

#' Server logic for the analyses of a single alternative splicing event
#' @param input Shiny input
#' @param output Shiny ouput
#' @param session Shiny session
#' 
#' @importFrom highcharter renderHighchart
#' @importFrom shinyjs runjs
#' @return NULL (this function is used to modify the Shiny session's state)
diffSplicingEventServer <- function(input, output, session) {
    ns <- session$ns
    
    selectGroupsServer(session, "diffGroups")
    
    observeEvent(input$analyse, {
        # Get splicing event's inclusion levels
        psi <- getInclusionLevels()
        if (is.null(psi)) {
            missingDataModal(session, "Inclusion levels",
                             ns("missingInclusionLevels"))
            return(NULL)
        }
        
        # Get selected event
        event <- getEvent()
        if (is.null(event) || event == "") {
            errorModal(session, "No event selected",
                       "Please, select an alternative splicing event.")
            return(NULL)
        }
        
        # Check if bandwidth is valid
        bandwidth <- input$bandwidth
        if (bandwidth <= 0 || is.na(bandwidth)) {
            errorModal(session, "Bandwidth must have a positive value",
                       "The density smoothing bandwidth requires a number",
                       "higher than 0.")
            return(NULL)
        }
        
        # Prepare groups of samples to analyse
        groups <- getSelectedGroups(input, "diffGroups", samples=TRUE,
                                    filter=colnames(psi))
        if ( !is.null(groups) ) {
            attrGroups <- groups
            psi <- psi[ , unlist(groups), drop=FALSE]
            groups <- rep(names(groups), sapply(groups, length))
        } else {
            attrGroups <- "All samples"
            groups <- rep(attrGroups, ncol(psi))
        }
        
        # Check if analyses were already performed
        stats <- getDifferentialAnalyses()
        if (!is.null(stats) && identical(attrGroups, attr(stats, "groups")))
            stat <- stats[event, ]
        else
            stat <- NULL
        
        # Separate samples by their groups
        eventPSI <- as.numeric(psi[event, ])
        eventPSI <- filterGroups(eventPSI, groups)
        groups <- names(eventPSI)
        
        plot <- plotDistribution(eventPSI, groups, bw=input$bandwidth,
                                 title=gsub("_", " ", event))
        output$density <- renderHighchart(plot)
        
        output$basicStats <- renderUI(basicStats(eventPSI, groups))
        output$ttest      <- renderUI(ttest(eventPSI, groups, stat))
        output$wilcox     <- renderUI(wilcox(eventPSI, groups, stat))
        output$kruskal    <- renderUI(kruskal(eventPSI, groups, stat))
        output$levene     <- renderUI(levene(eventPSI, groups, stat))
        output$fligner    <- renderUI(fligner(eventPSI, groups, stat))
        # output$fisher   <- renderUI(fisher(eventPSI, groups))
        # output$spearman <- renderUI(spearman(eventPSI, groups))
    })
    
    observeEvent(input$missingInclusionLevels, 
                 missingDataGuide("Inclusion levels"))
}

attr(diffSplicingEventUI, "loader") <- "diffSplicing"
attr(diffSplicingEventUI, "name") <- "Single event"
attr(diffSplicingEventServer, "loader") <- "diffSplicing"