#' Interface for differential analyses on all splicing events
#' 
#' @param id Character: identifier
#' 
#' @importFrom shinyjs disabled hidden
#' @importFrom shiny downloadLink selectizeInput uiOutput actionButton tags
#' checkboxGroupInput helpText tagList sidebarLayout mainPanel sidebarPanel
#' hoverOpts plotOutput
#' @importFrom shinyBS bsCollapse bsCollapsePanel
#' @importFrom DT dataTableOutput
#' @importFrom highcharter highchartOutput
#' 
#' @return HTML elements
diffSplicingTableUI <- function(id) {
    ns <- NS(id)
    pvalueAdjust <- list("No p-value adjustment"="none",
                         "False Discovery Rate"=c(
                             "Benjamini-Hochberg's method"="BH",
                             "Benjamini-Yekutieli's method"="BY"),
                         "Family-wise Error Rate"=c(
                             "Bonferroni correction"="bonferroni",
                             "Holm's method"="holm",
                             "Hochberg's method"="hochberg",
                             "Hommel's method"="hommel"))
    
    statAnalysesOptions <- tagList(
        selectGroupsUI(ns("diffGroups"), label="Groups of samples to analyse",
                       noGroupsLabel="All samples as one group",
                       groupsLabel="Samples by selected groups"),
        checkboxGroupInput(
            ns("statsChoices"),
            "Choose statistical analyses to perform:",
            # Basic stats is on and disabled by JavaScript
            c("Variance and median"="basicStats",
              "Unpaired t-test (2 groups)"="ttest",
              "Wilcoxon rank sum test (2 groups)"="wilcoxRankSum",
              "Kruskal-Wallis rank sum test (2 or more groups)"="kruskal", 
              "Levene's test (2 or more groups)"="levene",
              "Fligner-Killeen test (2 or more groups)"="fligner",
              "Distribution of alternative splicing quantification per group"=
                  "density"),
            selected=c("basicStats", "kruskal", "levene", "density", "ttest",
                       "fligner", "wilcoxRankSum")),
        # Disable checkbox of basic statistics
        tags$script('$("[value=basicStats]").attr("disabled", true);'),
        helpText("For each alternative splicing event, groups with one or less",
                 "non-missing values are discarded."), hr(),
        selectizeInput(ns("pvalueAdjust"), selected="BH",
                       "Adjust p-values of statistical tests", pvalueAdjust),
        processButton(ns("startAnalyses"), "Perform analyses"))
    
    survivalOptions <- tagList(
        helpText("For each splicing event, find the PSI cut-off that maximizes",
                 "differences in survival."),
        radioButtons(ns("censoring"), "Data censoring", selected="right",
                     inline=TRUE, choices=c(Left="left", Right="right",
                                            Interval="interval", 
                                            "Interval 2"="interval2")),
        selectizeInput(ns("timeStart"), choices=character(), "Follow up time"),
        # If the chosen censoring contains the word 'interval', show this input
        conditionalPanel(
            sprintf("input[id='%s'].indexOf('interval') > -1", ns("censoring")),
            selectizeInput(ns("timeStop"), choices=character(), "Ending time")),
        helpText("For patients for which there is no event reported, time",
                 "to last follow up is used instead."),
        selectizeInput(ns("event"), choices=NULL, 
                       "Event of interest"),
        radioButtons(
            ns("selected"), "Perform survival analysis based on",
            choices=c(
                "Splicing events shown in the screen"="shown",
                "Filtered splicing events (may be a slow process)"="filtered",
                "All splicing events (slow process)"="all")),
        processButton(ns("survival"), "Plot survival curves")
    )
    
    sidebar <- sidebarPanel(
        bsCollapse(
            id=ns("diffAnalysesCollapse"), open="statAnalyses",
            bsCollapsePanel(
                list(icon("tasks"), "Perform statistical analyses"),
                value="statAnalyses", statAnalysesOptions, style="info"),
            bsCollapsePanel(
                list(icon("tasks"), "Event plot options and table filtering"),
                style="info", value="plotEvents", uiOutput(ns("eventOptions"))),
            bsCollapsePanel(
                list(icon("tasks"),
                     "Survival analyses by splicing quantification cut-off"),
                style="info", value="survivalOptionsPanel",
                hidden(div(id=ns("survivalOptions"), survivalOptions)),
                div(id=ns("noSurvivalOptions"), "No survival data available"))),
        hr(),
        disabled(div(id=ns("downloadStats"), class="btn-group",
                     tags$button(class="btn btn-default dropdown-toggle",
                                 type="button", "data-toggle"="dropdown",
                                 "aria-haspopup"="true",
                                 "aria-expanded"="false", 
                                 icon("download"), 
                                 "Save table", tags$span(class="caret")),
                     tags$ul(class="dropdown-menu", 
                             tags$li(downloadLink(ns("downloadAll"), 
                                                  "All data")),
                             tags$li(downloadLink(ns("downloadSubset"), 
                                                  "Filtered data"))))))
    
    td <- lapply(seq(10), function(i) {
        chart <- highchartOutput(paste0(ns("curves"), i), height="100px")
        chart[[1]] <- tagAppendAttributes(
            chart[[1]], style="min-width: 150px;")
        
        tags$td(style="word-wrap: break-word;", style="white-space: normal;", 
                style="padding: 5px;",
                uiOutput(paste0(ns("eventText"), i)), chart,
                uiOutput(paste0(ns("eventSurvStats"), i)))
    })
    
    tagList(uiOutput(ns("modal")),
            sidebarLayout(
                sidebar, 
                mainPanel(
                    plotOutput(ns("plot"), brush=ns("brush"),
                               hover=hoverOpts(ns("hover"), delay=50, 
                                               delayType="throttle")),
                    uiOutput(ns("tooltip")),
                    dataTableOutput(ns("statsTable")),
                    highchartOutput(ns("densitySparklines"), 0, 0),
                    hidden(
                        tags$table(id=ns("survTable"), class="table", 
                                   class="table-bordered",
                                   class="text-center",
                                   style="margin-top: 20px;",
                                   tags$tbody(
                                       do.call(tags$tr, td[1:5]),
                                       do.call(tags$tr, td[6:10])))))))
}

#' Optimal survival difference given an inclusion level cut-off for a specific
#' alternative splicing event
#' 
#' @importFrom shinyjs runjs show hide
#' @importFrom shiny renderText
#' 
#' @param session Shiny session
#' @param input Shiny input
#' @param output Shiny output
#' @return NULL (this function is used to modify the Shiny session's state)
optimSurvDiff <- function(session, input, output) {
    ns <- session$ns
    
    # Interface of survival analyses
    observe({
        clinical <- getClinicalData()
        if (!is.null(getDifferentialAnalyses()) && !is.null(clinical)) {
            show("survivalOptions")
            hide("noSurvivalOptions")
            updateClinicalParams(session, clinical)
        } else {
            hide("survivalOptions")
            show("noSurvivalOptions")
        }
    })
    
    # Update selectize input label depending on the chosen censoring type
    observe({
        if (!is.null(getDifferentialAnalyses()) &&
            !is.null(getClinicalData()) && !is.null(input$censoring)) {
            label <- "Follow up time"
            if (grepl("interval", input$censoring, fixed=TRUE))
                label <- "Starting time"
            updateSelectizeInput(session, "timeStart", label=label)
        }
    })
    
    #' Calculate optimal survival cut-off for the inclusion levels of a given
    #' alternative splicing event
    observeEvent(input$survival, {
        time <- startProcess("survival")
        isolate({
            clinical      <- getClinicalData()
            psi           <- getInclusionLevels()
            match         <- getClinicalMatchFrom("Inclusion levels")
            statsTable    <- getDifferentialAnalyses()
            statsFiltered <- getDifferentialAnalysesFiltered()
            optimSurv     <- getDifferentialAnalysesSurvival()
            # User input
            censoring <- input$censoring
            timeStart <- input$timeStart
            timeStop  <- input$timeStop
            event     <- input$event
            display   <- input$statsTable_rows_current
            filtered  <- input$statsTable_rows_all
            selected  <- input$selected
        })
        
        if ("shown" %in% selected) {
            if (!is.null(display)) {
                events <- rownames(statsFiltered)[display]
                subset <- psi[events, ]
            } else {
                errorModal(session, "Error with selected events",
                           "Unfortunately, it's not possible to get events",
                           "shown in the table.")
                endProcess("survival")
                return(NULL)
            }
        } else if ("filtered" %in% selected) {
            if (!is.null(filtered)) {
                events <- rownames(statsFiltered)[filtered]
                subset <- psi[events, ]
            } else {
                errorModal(session, "Error with selected events",
                           "Unfortunately, it's not possible to get the events",
                           "from the table.")
                endProcess("survival")
                return(NULL)
            }
        } else if ("all" %in% selected) {
            subset <- psi
        }
        startProgress("Performing survival analysis", nrow(subset))
        
        # Assign alternative splicing quantification to patients based on their
        # samples
        clinicalPSI <- getPSIperPatient(subset, match, clinical)
        
        opt <- apply(clinicalPSI, 1, function(eventPSI) {
            opt <- optimalPSIcutoff(clinical, eventPSI, censoring, event, 
                                    timeStart, timeStop)
            
            updateProgress("Survival analysis", console=FALSE)
            return(c("Optimal survival PSI cut-off"=opt$par,
                     "Minimal survival p-value"=opt$value))
        })
        
        if (length(opt) == 0) {
            errorModal(session, "No survival analyses",
                       "Optimal PSI cut-off for the selected alternative",
                       "splicing events returned no survival analyses.")
        } else {
            df <- data.frame(t(opt))
            if (is.null(optimSurv)) {
                # Prepare survival table
                nas <- rep(NA, nrow(statsTable))
                optimSurv <- data.frame(nas, nas)
                rownames(optimSurv) <- rownames(statsTable)
                colnames(optimSurv) <- colnames(df)
            }
            for (col in names(df)) optimSurv[rownames(df), col] <- df[ , col]
            setDifferentialAnalysesSurvival(optimSurv)
        }
        endProcess("survival", time)
    })
    
    # Hide survival curves if there's no data
    observe({
        if (is.null(getDifferentialAnalysesSurvival()))
            hide("survTable", anim=TRUE, animType="fade")
        else
            show("survTable", anim=TRUE, animType="fade")
    })
    
    # Render container and information for survival curves
    observe({
        optimSurv <- getDifferentialAnalysesSurvival()
        if (is.null(optimSurv)) {
            lapply(seq(10), function(i) {
                output[[paste0("eventText", i)]] <- renderText(NULL)
                output[[paste0("eventSurvStats", i)]] <- renderUI(NULL)
            })
            return(NULL)
        }
        
        statsTable    <- getDifferentialAnalyses()
        statsFiltered <- getDifferentialAnalysesFiltered()
        
        lapply(seq(10), function(i) {
            # Get matching index between filtered table and whole table
            row <- input$statsTable_rows_current
            if (is.null(row) || length(row) < i)
                index <- NULL
            else
                index <- match(rownames(statsFiltered)[row[i]],
                               rownames(statsTable))
            
            if (!is.null(index) && !is.na(index)) {
                cutoff <- optimSurv[index, 1]
                pvalue <- optimSurv[index, 2]
            }
            
            if (is.null(index) || is.na(index)) {
                stat <- NULL
                splicingEvent <- NULL
            } else if (is.na(cutoff)) {
                splicingEvent <- rownames(optimSurv)[index]
                stat <- div(class="panel panel-default", div(
                    class="panel-heading",
                    "Survival analysis not yet requested for this event"))
            } else if (is.na(pvalue)) {
                splicingEvent <- rownames(optimSurv)[index]
                stat <- tagList(
                    "No optimal PSI cut-off", tags$br(), "found for this event")
            } else {
                splicingEvent  <- rownames(optimSurv)[index]
                stat <- tags$table(class="table", class="table-condensed", 
                                   style="margin-bottom: 0",
                                   tags$tbody(
                                       tags$tr(tags$td(tags$b("PSI cut-off")),
                                               tags$td(roundDigits(cutoff))),
                                       tags$tr(tags$td(tags$b("p-value")),
                                               tags$td(pvalue))))
            }
            
            if (is.null(splicingEvent) || is.na(splicingEvent)) stat <- NULL
            output[[paste0("eventText", i)]] <- renderUI(
                tags$a(
                    parseSplicingEvent(splicingEvent, char=TRUE),
                    class="label label-default", style="display: inline-block;",
                    style="white-space: normal;", 
                    onclick=paste0("showSurvCutoff('", splicingEvent,
                                   "', true)")))
            output[[paste0("eventSurvStats", i)]] <- renderUI(tags$small(stat))
        })
    })
    
    # Render survival curves
    observe({
        optimSurv <- getDifferentialAnalysesSurvival()
        if (is.null(optimSurv)) {
            lapply(seq(10), function(i) 
                output[[paste0("curves", i)]] <- renderHighchart(NULL))
            return(NULL)
        }
        
        isolate({
            clinical      <- getClinicalData()
            psi           <- getInclusionLevels()
            match         <- getClinicalMatchFrom("Inclusion levels")
            # User input
            censoring <- input$censoring
            timeStart <- input$timeStart
            timeStop  <- input$timeStop
            event     <- input$event
        })
        
        statsTable    <- getDifferentialAnalyses()
        statsFiltered <- getDifferentialAnalysesFiltered()
        
        # Interface for the survival curves of 10 splicing events
        lapply(seq(10), function(i) {
            # Get matching index between filtered table and whole table
            row <- input$statsTable_rows_current
            if (is.null(row) || length(row) < i)
                index <- NULL
            else
                index <- match(rownames(statsFiltered)[row[i]], 
                               rownames(statsTable))
            if (is.null(index) || is.na(index)) {
                output[[paste0("curves", i)]] <- renderHighchart(NULL)
            } else {
                hc <- plotMiniSurvivalCurves(i, input, index, optimSurv,
                                             clinical, match, psi, censoring, 
                                             event, timeStart, timeStop)
                output[[paste0("curves", i)]] <- renderHighchart(hc)
            }
        })
    })
}

#' Perform and plot survival curves
#' 
#' @inheritParams processSurvTerms
#' @param survParams List of parameters to plot survival curves
#' @param i Integer: index of the survival curves plot of interest
#' @param input Shiny input
#' @param index Integer: index of event to draw
#' @param match Integer: samples matched with clinical patients
#' @param psi Data frame or matrix: alternative splicing quantification
#' 
#' @importFrom highcharter hc_legend hc_xAxis hc_yAxis hc_chart hc_plotOptions
#' 
#' @return A \code{"highchart"} object to plot
plotMiniSurvivalCurves <- function(i, input, index, survParams, clinical,
                                   match, psi, censoring, event, timeStart, 
                                   timeStop) {
    cutoff <- survParams[index, 1]
    
    if (!is.null(index) && !is.na(cutoff)) {
        show(paste0("curves", i), anim=TRUE)
        splicingEvent  <- rownames(survParams)[index]
        
        # Assign alternative splicing quantification to patients based on their
        # samples
        clinicalPSI <- getPSIperPatient(psi, match, clinical)
        eventPSI <- as.numeric(clinicalPSI[splicingEvent, ])
        
        # Assign a value based on the inclusion levels cut-off
        group <- labelBasedOnCutoff(eventPSI, cutoff, label="")
        
        survTerms <- processSurvTerms(clinical, censoring, event, timeStart, 
                                      timeStop, group)
        surv <- survfit(survTerms)
        
        hc <- plotSurvivalCurves(surv, mark = FALSE, title=NULL) %>%
            hc_legend(enabled=FALSE) %>%
            hc_xAxis(title=list(text=""), showLastLabel=TRUE) %>%
            hc_yAxis(title=list(text="")) %>%
            hc_chart(zoomType=NULL) %>%
            hc_chart(events=list(click=JS(paste0(
                "function(e) { showSurvCutoff('", splicingEvent,
                ", true') }"))))
    } else {
        hc <- NULL
        hide(paste0("curves", i), anim=TRUE)
    }
    return(hc)
}

#' Create plot for events
#' 
#' @param df Data frame
#' @param x Character: name of the variable used for the X axis
#' @param y Character: name of the variable used for the Y axis
#' @param params List of parameters to pass to \code{\link[ggplot2]{geom_point}}
#' related to most points
#' @param highlightX Integer: region of points in X axis to highlight
#' @param highlightY Integer: region of points in Y axis to highlight
#' @param highlightParams List of parameters to pass to
#' \code{\link[ggplot2]{geom_point}} related to highlighted points
#' @param selected Integer: index of rows/points to be coloured
#' @param selectedParams List of parameters to pass to 
#' \code{\link[ggplot2]{geom_point}} related to selected points
#' 
#' @importFrom ggplot2 ggplot aes_string geom_point theme_light
#' 
#' @return HTML elements
createEventPlotting <- function(df, x, y, params, highlightX, highlightY,
                                highlightParams, selected, selectedParams) {
    aes <- aes_string(paste0("`", x, "`"), paste0("`", y, "`"))
    
    # Get points highlighted in X and Y that were not selected
    if (!is.null(highlightX)) {
        highlightedX <- findInterval(df[[x]], highlightX, left.open=FALSE,
                                     rightmost.closed=TRUE) == 1
        if (attr(highlightX, "inverted")) highlightedX <- !highlightedX
        highlightedX <- which(highlightedX)
    } else {
        highlightedX <- seq(nrow(df))
    }
    
    if (!is.null(highlightY)) {
        highlightedY <- findInterval(df[[y]], highlightY, left.open=FALSE,
                                     rightmost.closed=TRUE) == 1
        if (attr(highlightY, "inverted")) highlightedY <- !highlightedY
        highlightedY <- which(highlightedY)
    } else {
        highlightedY <- seq(nrow(df))
    }
    
    if ( is.null(highlightX) && is.null(highlightY) ) {
        highlighted <- NULL
    } else {
        highlighted <- intersect(highlightedX, highlightedY)
    }
    setDifferentialAnalysesHighlightedEvents(highlighted)
    
    # Render remaining points
    plotted <- union(selected, highlighted)
    if (!is.null(plotted))
        remaining <- df[-plotted, ]
    else
        remaining <- df
    plot <- ggplot() + do.call("geom_point", c(
        list(data=remaining, aes, na.rm=TRUE), params))
    
    # Render highlighted points
    plot <- plot + do.call("geom_point", c(
        list(data=df[setdiff(highlighted, selected), ], aes, na.rm=TRUE), 
        highlightParams))
    
    # Render selected points
    plot <- plot + do.call("geom_point", c(
        list(data=df[selected, ], aes, na.rm=TRUE), selectedParams))
    return(plot + theme_light(16))
}

#' Create the interface for the tooltip of a plot
#' 
#' @inheritParams createEventPlotting
#' @param hover Mouse hover information for a given plot as retrieved from
#' \code{\link[shiny]{hoverOpts}}
#' 
#' @importFrom shiny tags nearPoints wellPanel
#' 
#' @return HTML elements
createTooltip <- function(df, hover, x, y) {
    point <- nearPoints(df, hover, threshold=10, maxpoints=1, addDist=TRUE,
                        xvar=x, yvar=y)
    if (nrow(point) == 0) return(NULL)
    
    # Calculate point position inside the image as percent of total 
    # dimensions from left (horizontal) and from top (vertical)
    xDomain  <- hover$domain$right - hover$domain$left
    right_pct <- (hover$domain$right - hover$x) / xDomain
    yDomain  <- hover$domain$top - hover$domain$bottom
    top_pct  <- (hover$domain$top - hover$y) / yDomain
    
    # Calculate distance from left and bottom in pixels
    xRange   <- hover$range$right - hover$range$left
    right_px <- right_pct * xRange + 25
    yRange  <- hover$range$bottom - hover$range$top
    top_px  <- hover$range$top + top_pct * yRange + 2
    
    # Tooltip
    wellPanel(
        class = "well-sm",
        style = paste0("position:absolute; z-index:100;",
                       "background-color: rgba(245, 245, 245, 0.85); ",
                       "right:", right_px, "px; top:", top_px, "px;"),
        tags$table(class="table table-condensed", style="margin-bottom: 0;",
                   tags$thead(
                       tags$tr(
                           tags$td(tags$b("Event")),
                           tags$td(
                               parseSplicingEvent(rownames(point), char=TRUE))
                       )
                   ),
                   tags$tbody(
                       tags$tr(
                           tags$td(tags$b(x)),
                           tags$td(point[[x]])
                       ),
                       tags$tr(
                           tags$td(tags$b(y)),
                           tags$td(point[[y]]))))
    )
}

#' Show variable transformation(s)
#' 
#' @param label Character: label to display
#' @param type Character: show the variable transformation for the chosen type;
#' NULL (by default) to show all variable transformations
#' 
#' @return Character labelling variable transformation(s)
transformOptions <- function(label, type=NULL) {
    transform <- c("No transformation"="no",
                   "|%s|"="abs",
                   "-%s"="inv",
                   "log10(%s)"="log10",
                   "log10(|%s|)"="log10abs",
                   "log10(-%s)"="log10-",
                   "-log10(%s)"="-log10",
                   "-log10(|%s|)"="-log10abs",
                   "-log10(-%s)"="-log10-")
    names(transform) <- sprintf(names(transform), label)
    
    if (!is.null(type)) {
        show <- names(transform)
        show[[1]] <- label
        show <- show[match(type, transform)]
        return(show)
    } else {
        return(transform)
    }
}

#' Transform values as per a given type of transformation
#' 
#' @param val Integer: values to transform
#' @param type Character: type of transformation
#' @param avoidZero Boolean: add the smallest non-zero number available to zero
#' values; avoids returning infinity values during Log transformation (which are
#' not plotted); useful for preserving p-values of 0, for instance; TRUE by
#' default
#' 
#' @return Integer containing transformed values
transformValues <- function(val, type, avoidZero=TRUE) {
    # Remove NAs
    if (avoidZero) {
        zeroes <- val == 0 & !is.na(val)
        val[zeroes] <- val[zeroes] + .Machine$double.xmin
    }
    
    trans <- switch(type,
                    "no"=val,
                    "abs"=abs(val),
                    "inv"=-val,
                    "log10"=log10(val),
                    "log10abs"=log10(abs(val)),
                    "log10-"=log10(-val),
                    "-log10"=-log10(val),
                    "-log10abs"=-log10(abs(val)),
                    "-log10-"=-log10(-val))
    return(trans)
}

#' Options for event plotting
#' 
#' @param ns Namespace function
#' @param df Data frame
#' 
#' @return HTML elements
eventPlotOptions <- function(ns, df) {
    # Only allow to select numeric columns    
    cols <- colnames(df)
    type <- sapply(cols, function(i) class(df[[i]]))
    numericCols <- cols[type == "numeric"]
    
    # Default option for X axis
    deltaMedian <- "Delta median"
    if (deltaMedian %in% numericCols)
        xSelected <- deltaMedian
    else
        xSelected <- NULL
    
    # Default option for Y axis
    pValue <- grepl("p-value", numericCols)
    if (any(pValue))
        ySelected <- numericCols[pValue][1]
    else
        ySelected <- NULL
    
    tabsetPanel(
        tabPanel(
            "X axis",
            selectizeInput(ns("xAxis"), "Select X axis", numericCols,
                           selected=xSelected),
            selectizeInput(ns("xTransform"), 
                           "Data transformation of X values",
                           transformOptions("x")),
            checkboxInput(ns("xHighlight"),
                          paste("Highlight points based on X values")),
            uiOutput(ns("xHighlightValues"))),
        tabPanel(
            "Y axis",
            selectizeInput(ns("yAxis"), "Select Y axis", numericCols,
                           selected=ySelected),
            selectizeInput(ns("yTransform"),
                           "Data transformation of Y values",
                           transformOptions("y")),
            checkboxInput(ns("yHighlight"),
                          paste("Highlight points based on Y values")),
            uiOutput(ns("yHighlightValues"))),
        navbarMenu(
            "Plot style",
            tabPanel("Base points",
                     plotPointsStyle(ns, "base", "Base points", size=2,
                                     colour="grey", alpha=0.3)),
            tabPanel("Highlighted points",
                     plotPointsStyle(ns, "highlighted", "Highlighted points",
                                     size=3, colour="orange", alpha=0.5)),
            tabPanel("Selected in the table",
                     plotPointsStyle(ns, "selected",
                                     "Selected in the table",
                                     help=paste("Click in a row of the table",
                                                "to emphasise the respective",
                                                "point in the plot."),
                                     size=8, colour="blue", alpha=0.5))))
}

#' Interface to modify the style of the plot points
#' 
#' @param ns Namespace function
#' @param id Character: identifier
#' @param description Character: display text for user
#' @param help Character: extra text to help the user
#' @param colour Character: default colour ("black" by default)
#' @param size Integer: default size (2 by default)
#' @param alpha Numeric: default transparency value; (opaque by default)
#' 
#' @importFrom shiny tagList h4 helpText sliderInput
#' @importFrom colourpicker colourInput
#' 
#' @return HTML elements
plotPointsStyle <- function(ns, id, description, help=NULL, size=2,
                            colour="black", alpha=1.0) {
    id2 <- function(att) ns(paste0(id, att))
    tagList(
        h4(description),
        if (!is.null(help)) helpText(help),
        sliderInput(id2("Size"), "Size", min=1, max=10, step=1, value=size),
        colourInput(id2("Colour"), "Colour", value=colour),
        sliderInput(id2("Alpha"), "Alpha", min=0, max=1, step=0.01, value=alpha)
    )
}

#' Transform data in data frame
#' 
#' @param input Shiny input
#' @param df Data frame
#' @param x Character: column name
#' @param y Character: column name
#' 
#' @return Data frame with transformed data in new columns and respective name
#' of created columns
transformData <- function(input, df, x, y) {
    xTrans <- input$xTransform
    xLabel <- transformOptions(x, xTrans)
    if (!x %in% colnames(df)) return(NULL)
    df[[xLabel]] <- transformValues(df[[x]], xTrans)
    
    yTrans <- input$yTransform
    yLabel <- transformOptions(y, yTrans)
    if (!y %in% colnames(df)) return(NULL)
    df[[yLabel]] <- transformValues(df[[y]], yTrans)
    
    return(list(data=df, xLabel=xLabel, yLabel=yLabel))
}

#' Server logic of the exploratory differential analyses
#' 
#' @param input Shiny input
#' @param output Shiny ouput
#' @param session Shiny session
#' 
#' @importFrom shiny checkboxGroupInput
#' @importFrom shinyBS updateCollapse
#' @importFrom shinyjs toggleState disable
#' @importFrom DT replaceData dataTableProxy
#' @importFrom utils write.table
#' @importFrom colourpicker colourInput
#' @importFrom grDevices palette
#' 
#' @return NULL (this function is used to modify the Shiny session's state)
diffSplicingTableServer <- function(input, output, session) {
    ns <- session$ns
    
    selectGroupsServer(session, "diffGroups")
    
    observeEvent(input$missingInclusionLevels, 
                 missingDataGuide("Inclusion levels"))
    
    performDiffAnalyses <- reactive({
        # Get event's inclusion levels
        psi <- getInclusionLevels()
        statsChoices <- input$statsChoices
        pvalueAdjust <- input$pvalueAdjust
        
        totalTime <- startProcess("startAnalyses")
        
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
        
        stats <- diffAnalyses(psi, groups, statsChoices,
                              pvalueAdjust=pvalueAdjust,
                              progress=updateProgress)
        attr(stats, "groups") <- attrGroups
        setDifferentialAnalyses(stats)
        setDifferentialAnalysesSurvival(NULL)
        updateCollapse(session, "diffAnalysesCollapse", "plotEvents")
        endProcess("startAnalyses", totalTime)
    })
    
    # Perform statistical analyses
    observeEvent(input$startAnalyses, {
        psi <- isolate(getInclusionLevels())
        diffAnalyses <- isolate(getDifferentialAnalyses())
        if ( is.null(psi) ) {
            missingDataModal(session, "Inclusion levels",
                             ns("missingInclusionLevels"))
        } else if ( !is.null(diffAnalyses) ) {
            warningModal(session, "Differential analyses already performed",
                         "Do you wish to replace the loaded analyses?",
                         footer=actionButton(ns("replace"), "Replace",
                                             class="btn-warning",
                                             "data-dismiss"="modal"))
        } else {
            performDiffAnalyses()
        }
    })
    
    observeEvent(input$replace, {
        performDiffAnalyses()
        # Reset previous results from differential analyses
        setDifferentialAnalysesFiltered(NULL)
        setDifferentialAnalysesBrushedEvents(NULL)
        setDifferentialAnalysesHighlightedEvents(NULL)
        setDifferentialAnalysesSurvival(NULL)
    })
    
    colours <- list("Base colours"=capitalize(palette()),
                    "Extra colours"=capitalize(setdiff(colours(), palette())))
    
    # Options to plot events
    output$eventOptions <- renderUI({
        stats <- getDifferentialAnalyses()
        if (is.null(stats))
            return("Perform differential splicing analyses, first.")
        eventPlotOptions(ns, stats)
    })
    
    lapply(c("x", "y"), function(axis) {
        observe({
            highlightUI <- function(label, min, max) {
                highlightId <- ns(paste0(label, "Highlight"))
                sliderId    <- ns(paste0(label, "Slider"))
                sliderInvId <- ns(paste0(label, "SliderInv"))
                
                conditionalPanel(
                    sprintf("input[id='%s']", highlightId),
                    sliderInput(sliderId, "Values to highlight", min=min,
                                max=max, value=c(min, max), dragRange=TRUE,
                                step=0.001, round=getPrecision(), sep=""),
                    checkboxInput(sliderInvId, "Invert highlighted values"),
                    helpText("The data in the table is also filtered",
                             "according to highlighted events."))
            }
            
            stats <- getDifferentialAnalyses()
            value <- input[[paste0(axis, "Axis")]]
            if (is.null(stats) || is.null(value)) {
                output[[paste0(axis, "HighlightValues")]] <- renderUI(NULL)
                return(NULL)
            }
            
            trans <- input[[paste0(axis, "Transform")]]
            label <- transformOptions(value, trans)
            if (!value %in% colnames(stats)) {
                output[[paste0(axis, "HighlightValues")]] <- renderUI(NULL)
                return(NULL)
            }
            vals  <- transformValues(stats[[value]], trans)
            
            rangeNos <- range(vals, na.rm=TRUE)
            minNo    <- min(rangeNos)
            maxNo    <- max(rangeNos)
            
            output[[paste0(axis, "HighlightValues")]] <- renderUI( 
                highlightUI(axis, minNo, maxNo) )
        })
    })
    
    # Plot events and render the plot tooltip
    observe({
        stats <- getDifferentialAnalyses()
        x <- input$xAxis
        y <- input$yAxis
        if (is.null(stats) || is.null(x) || is.null(y)) {
            output$plot <- renderPlot(NULL)
            output$tooltip <- renderUI(NULL)
            return(NULL)
        }
        
        res <- transformData(input, stats, x, y)
        if (is.null(res)) {
            output$plot <- renderPlot(NULL)
            output$tooltip <- renderUI(NULL)
            return(NULL)
        }
        
        stats  <- res$data
        xLabel <- res$xLabel
        yLabel <- res$yLabel
        
        output$plot <- renderPlot({
            if (input$xHighlight) {
                highlightX <- input$xSlider
                attr(highlightX, "inverted") <- input$xSliderInv
            } else {
                highlightX <- NULL
            }
            
            if (input$yHighlight) {
                highlightY <- input$ySlider
                attr(highlightY, "inverted") <- input$ySliderInv
            } else {
                highlightY <- NULL
            }
            
            # Check selected events
            selected <- input$statsTable_rows_selected
            events <- getDifferentialAnalysesHighlightedEvents()
            if (!is.null(selected) && !is.null(events))
                selected <- events[selected]
            
            params <- list(size=input$baseSize, col=input$baseColour,
                           alpha=input$baseAlpha)
            highlightParams <- list(size=input$highlightedSize,
                                    col=input$highlightedColour,
                                    alpha=input$highlightedAlpha)
            selectedParams  <- list(size=input$selectedSize,
                                    col=input$selectedColour,
                                    alpha=input$selectedAlpha)
            
            createEventPlotting(stats, xLabel, yLabel, params, highlightX, 
                                highlightY, highlightParams, selected, 
                                selectedParams)
        })
        output$tooltip <- renderUI( createTooltip(stats, input$hover, xLabel, 
                                                  yLabel) )
    })
    
    # Save brushed events in the plot
    observe({
        stats <- getDifferentialAnalyses()
        x <- input$xAxis
        y <- input$yAxis
        if (is.null(stats) || is.null(x) || is.null(y)) return(NULL)
        
        res <- transformData(input, stats, x, y)
        if (is.null(res)) return(NULL)
        
        stats  <- res$data
        xLabel <- res$xLabel
        yLabel <- res$yLabel
        
        brush  <- input$brush
        if (!is.null(brush)) {
            xBrush <- findInterval(
                stats[[xLabel]], vec = c(brush$xmin, brush$xmax),
                left.open=FALSE, rightmost.closed=TRUE) == 1
            yBrush <- findInterval(
                stats[[yLabel]], vec = c(brush$ymin, brush$ymax),
                left.open=FALSE, rightmost.closed=TRUE) == 1
            brushed <- intersect(which(xBrush), which(yBrush))
        } else {
            brushed <- NULL
        }
        setDifferentialAnalysesBrushedEvents(brushed)
    })
    
    # Render table with sparklines
    output$statsTable <- renderDataTableSparklines({
        stats <- getDifferentialAnalyses()
        if (!is.null(stats)) {
            # Remove columns of no interest
            colFilter <- !grepl("method|data.name", colnames(stats))
            
            # Filter by highlighted events
            events  <- getDifferentialAnalysesHighlightedEvents()
            brushed <- getDifferentialAnalysesBrushedEvents()
            
            if (!is.null(events) && !is.null(brushed)) {
                rowFilter <- intersect(events, brushed)
            } else if (!is.null(events)) {
                rowFilter <- events
            } else if (!is.null(brushed)) {
                rowFilter <- brushed
            } else {
                rowFilter <- TRUE
            }
            setDifferentialAnalysesFiltered(stats[rowFilter, ])
            
            # Properly display event identifiers
            rownames(stats) <- parseSplicingEvent(rownames(stats), char=TRUE)
            return(stats[rowFilter, colFilter])
        }
    }, style="bootstrap", filter="top", server=TRUE, # selection="none",
    extensions="Buttons", options=list(
        pageLength=10, rowCallback=JS("createDiffSplicingLinks"), 
        dom="Bfrtip", buttons=I("colvis"),
        columnDefs=list(list(targets=5, searchable=FALSE))))
    
    # Disable download button if statistical table is NULL
    observe({
        if ( is.null(getDifferentialAnalyses()) )
            disable("downloadStats")
        else
            enable("downloadStats")
    })
    
    # Download whole table
    output$downloadAll <- downloadHandler(
        filename=paste(getCategory(), "diff. splicing analyses"),
        content=function(file) {
            stats <- getDifferentialAnalyses()
            densityCol <- TRUE
            if (!is.null(stats)) {
                densityCol <- -match("PSI distribution", colnames(stats))
                if (is.na(densityCol)) densityCol <- TRUE
            }
            
            write.table(stats[densityCol], file, quote=FALSE, sep="\t",
                        row.names=FALSE)
        }
    )
    
    # Download filtered table
    output$downloadSubset <- downloadHandler(
        filename=paste(getCategory(), "Differential splicing analyses"),
        content=function(file) {
            stats <- getDifferentialAnalysesFiltered()
            densityCol <- TRUE
            if (!is.null(stats)) {
                densityCol <- -match("PSI distribution", colnames(stats))
                if (is.na(densityCol)) densityCol <- TRUE
            }
            
            write.table(stats[input$statsTable_rows_all, densityCol], 
                        file, quote=FALSE, sep="\t", row.names=FALSE)
        }
    )
    
    # Optimal survival difference given an inclusion level cut-off for a 
    # specific alternative splicing event
    optimSurvDiff(session, input, output)
}

attr(diffSplicingTableUI, "loader") <- "diffSplicing"
attr(diffSplicingTableUI, "name") <- "All events (table)"
attr(diffSplicingTableUI, "selectEvent") <- FALSE
attr(diffSplicingTableServer, "loader") <- "diffSplicing"