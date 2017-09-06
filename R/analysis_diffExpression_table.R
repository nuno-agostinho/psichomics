#' @rdname appUI
#' 
#' @importFrom shinyjs disabled hidden
#' @importFrom shiny downloadLink selectizeInput uiOutput actionButton tags
#' checkboxGroupInput helpText tagList sidebarLayout mainPanel
#' @importFrom shinyBS bsCollapse bsCollapsePanel
#' @importFrom DT dataTableOutput
#' @importFrom highcharter highchartOutput
diffExpressionTableUI <- function(id) {
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
    
    statAnalysesOptions <- div(
        id=ns("statAnalysesOptions"),
        selectizeInput(ns("geneExpr"), "Gene expression", choices=NULL,
                       width="100%"),
        bsCollapse(
            open=c("lmFit"),
            bsCollapsePanel(
                tagList(icon("compress"), "Gene-wise linear model fit"),
                value="lmFit",
                helpText("The", tags$code("limma::lmFit"), "function is used",
                         "to fit a linear model per gene, based on a design",
                         "matrix prepared from the two selected groups."),
                selectGroupsUI(ns("diffGroups"), label="Groups of samples",
                               maxItems=2)),
            bsCollapsePanel(
                tagList(icon("adjust"), "Differential expression statistics"), 
                value="eBayes",
                helpText(
                    "The", tags$code("limma::eBayes"), "function is used to",
                    "compute moderated t-tests and log-odds of differential",
                    "expression by empirical Bayes moderation of the standard",
                    "errors towards a common value."),
                sliderInput(
                    ns("ebayesProportion"), min=0, max=1, value=0.01, step=0.01,
                    width="100%",
                    "Assumed proportion of differentially expressed genes"),
                hr(),
                helpText("Assumed limit for the standard deviation of log2",
                         "fold-changes for differentially expressed genes:"),
                fluidRow(
                    column(
                        6, numericInput(ns("ebaysStdevMin"), "Lower limit",
                                        min=0, value=0.1, step=0.1, 
                                        width="100%")),
                    column(
                        6, numericInput(ns("ebaysStdevMax"), "Upper limit",
                                        min=0, value=4, step=0.1, 
                                        width="100%"))))),
        tags$b("Extra analyses to be performed:"),
        tags$ul(tags$li("Variance and median expression"),
                tags$li("Distribution of gene expression per group")),
        selectizeInput(ns("pvalueAdjust"), selected="BH", width="100%",
                       "P-value adjustment", pvalueAdjust),
        processButton(ns("startAnalyses"), "Perform analyses"))
    
    eventOptions <- div(
        id=ns("eventOptions"),
        tabsetPanel(
            tabPanel(
                "X axis",
                selectizeInput(ns("xAxis"), "Select X axis", choices=NULL,
                               width="100%"),
                selectizeInput(ns("xTransform"),
                               "Data transformation of X values",
                               transformOptions("x"), width="100%"),
                checkboxInput(ns("xHighlight"), width="100%",
                              paste("Highlight points based on X values")),
                uiOutput(ns("xHighlightValues"))),
            tabPanel(
                "Y axis",
                selectizeInput(ns("yAxis"), "Select Y axis", choices=NULL,
                               width="100%"),
                selectizeInput(ns("yTransform"), width="100%",
                               "Data transformation of Y values",
                               transformOptions("y")),
                checkboxInput(ns("yHighlight"), width="100%",
                              paste("Highlight points based on Y values")),
                uiOutput(ns("yHighlightValues"))),
            navbarMenu(
                "Plot style",
                tabPanel("Base points",
                         plotPointsStyle(
                             ns, "base", "Base points",
                             help=paste("These are points not highlighted or",
                                        "selected."),
                             size=2, colour="grey", alpha=0.3)),
                tabPanel("Highlighted points",
                         plotPointsStyle(
                             ns, "highlighted", "Highlighted points",
                             help=paste("Highlight points in the X and Y axes",
                                        "options."),
                             size=3, colour="orange", alpha=0.5)),
                tabPanel("Selected in the table",
                         plotPointsStyle(
                             ns, "selected", "Selected in the table",
                             help=paste("Click in a row of the table to",
                                        "emphasise the respective point in",
                                        "the plot."),
                             size=8, colour="blue", alpha=0.5))))
    )
    
    sidebar <- sidebar(
        bsCollapse(
            id=ns("diffExpressionCollapse"), open="statAnalyses",
            bsCollapsePanel(
                list(icon("tasks"), "Perform differential expression analysis"),
                value="statAnalyses", style="info",
                errorDialog(
                    paste("Gene expression data is required for",
                          "differential expression."),
                    id=ns("missingGeneExpr"),
                    buttonLabel="Load data",
                    buttonIcon="plus-circle",
                    buttonId=ns("loadGeneExpr")),
                hidden(statAnalysesOptions)),
            bsCollapsePanel(
                list(icon("tasks"), "Event plot options and table filtering"),
                style="info", value="plotEvents",
                errorDialog(
                    "Differential expression analysis not yet performed.",
                    id=ns("missingDiffExpression")),
                hidden(eventOptions))),
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
    
    tagList(
        uiOutput(ns("modal")),
        sidebarLayout(
            sidebar, 
            mainPanel(
                ggplotUI(ns("ge-volcano")),
                dataTableOutput(ns("statsTable")),
                highchartOutput(ns("highchartsSparklines"), 0, 0))))
}


#' Set of functions to perform differential analyses
#' 
#' @importFrom shinyBS updateCollapse
#' @importFrom limma eBayes lmFit toptable
#' 
#' @inherit diffExpressionTableServer
diffExpressionSet <- function(session, input, output) {
    ns <- session$ns
    
    observe({
        geneExpr <- getGeneExpression()
        if (is.null(geneExpr)) {
            show("missingGeneExpr")
            hide("statAnalysesOptions")
        } else {
            updateSelectizeInput(session, "geneExpr",
                                 choices=rev(names(geneExpr)))
            hide("missingGeneExpr")
            show("statAnalysesOptions")
        }
    })
    
    performDiffExpression <- reactive({
        geneExpr <- getGeneExpression()[[input$geneExpr]]
        
        totalTime <- startProcess("startAnalyses")
        
        # Prepare groups of samples to analyse and filter samples not available 
        # in the selected groups from the gene expression data
        groups <- getSelectedGroups(input, "diffGroups", samples=TRUE,
                                    filter=colnames(geneExpr))
        geneExpr     <- geneExpr[ , unlist(groups), drop=FALSE]
        attrGroups   <- groups
        isFromGroup1 <- colnames(geneExpr) %in% groups[[1]]
        design       <- cbind(1, ifelse(isFromGroup1, 1, 0))
        
        # Fit a gene-wise linear model based on selected groups
        fit <- lmFit(geneExpr, design)
        
        # Calculate moderated t-statistics and DE log-odds using limma::eBayes
        ebayesProportion <- input$ebayesProportion
        ebaysStdevMin    <- input$ebaysStdevMin
        ebaysStdevMax    <- input$ebaysStdevMax
        stats <- eBayes(fit, proportion=ebayesProportion, 
                        stdev.coef.lim=c(ebaysStdevMin, ebaysStdevMax))
        
        # Prepare data summary
        pvalueAdjust <- input$pvalueAdjust
        summary <- toptable(stats, number=nrow(fit), coef=2,
                            adjust.method=pvalueAdjust, confint=TRUE)
        names(summary) <- c("log2 Fold-Change", "conf. int1", "conf. int2",
                            "moderated t-statistics", "p-value", 
                            paste0("p-value (", pvalueAdjust, " adjusted)"),
                            "B-statistics")
        attr(summary, "groups") <- attrGroups
        
        # Calculate basic statistics and density plots
        groups <- rep(names(groups), sapply(groups, length))
        stats  <- diffAnalyses(geneExpr, groups, c("basicStats", "density"),
                               pvalueAdjust=NULL, progress=updateProgress,
                               geneExpr=input$geneExpr)
        final  <- cbind(stats[ , c(1, 5:6)], summary, stats[ , 7:ncol(stats)])
        
        setDifferentialExpression(final)
        # setDifferentialExpressionSurvival(NULL)
        updateCollapse(session, "diffExpressionCollapse", "plotEvents")
        endProcess("startAnalyses", totalTime)
    })
    
    # Perform statistical analyses
    observeEvent(input$startAnalyses, {
        isolate({
            ge <- getGeneExpression()
            diffAnalyses <- getDifferentialExpression()
            groups <- input$diffGroups
        })
        if ( is.null(ge) ) {
            missingDataModal(session, "Gene expression", ns("missingGeneExpr"))
        } else if ( is.null(groups) || length(input$diffGroups) != 2 ) {
            errorModal(session, "Select two groups",
                       "Currently, two groups are required for differential",
                       "expression analysis. Please select two groups.")
        } else if ( !is.null(diffAnalyses) ) {
            warningModal(session,
                         "Differential expression analyses already performed",
                         "Do you wish to replace the loaded analyses?",
                         footer=actionButton(ns("replace"), "Replace",
                                             class="btn-warning",
                                             "data-dismiss"="modal"))
        } else {
            performDiffExpression()
        }
    })
    
    # Replace previously performed differential analyses
    observeEvent(input$replace, {
        performDiffExpression()
        # Reset previous results from differential analyses
        setDifferentialExpressionFiltered(NULL)
        setZoom("ge-volcano", NULL)
        setSelectedPoints("ge-volcano", NULL)
        setHighlightedPoints("ge-volcano", NULL)
        # setDifferentialExpressionSurvival(NULL)
    })
}

#' Set of functions to plot differential analyses
#' 
#' @inherit diffExpressionTableServer
diffExpressionPlotSet <- function(session, input, output) {
    ns <- session$ns
    # Toggle visibility of elements regarding event options
    observe({
        stats <- getDifferentialExpression()
        if (is.null(stats)) {
            show("missingDiffExpression")
            hide("eventOptions")
        } else {
            hide("missingDiffExpression")
            show("eventOptions")
        }
    })
    
    # Update columns available to plot
    observe({
        stats <- getDifferentialExpression()
        # optimSurv <- getDifferentialExpressionSurvival()
        # if (!is.null(optimSurv)) {
        #     optimSurvCols <- c("Optimal GE cutoff", "Log rank p-value")
        #     stats[[optimSurvCols[1]]] <- optimSurv[[1]]
        #     stats[[optimSurvCols[2]]] <- optimSurv[[2]]
        # 
        #     # Show these columns at the end
        #     names <- colnames(stats)
        #     colsMatch <- match(optimSurvCols, names)
        #     stats <- stats[c(names[-colsMatch], names[colsMatch])]
        # }
        eventPlotOptions(session, stats, isolate(input$xAxis),
                         isolate(input$yAxis))
    })
    
    # Interface elements to highlight values in the plot
    lapply(c("x", "y"), function(axis) {
        observe({
            highlightUI <- function(label, min, max) {
                highlightId <- ns(paste0(label, "Highlight"))
                sliderId    <- ns(paste0(label, "Slider"))
                sliderInvId <- ns(paste0(label, "SliderInv"))
                
                # Round max and min numbers with two decimal points
                max <- ceiling(max*100)/100
                min <- floor(min*100)/100
                
                conditionalPanel(
                    sprintf("input[id='%s']", highlightId),
                    sliderInput(sliderId, "Values to highlight", min=min,
                                max=max, value=c(min, max), dragRange=TRUE,
                                step=0.01, round=getPrecision(), sep="",
                                width="100%"),
                    checkboxInput(sliderInvId, "Invert highlighted values"),
                    helpText("The data in the table is also filtered",
                             "according to highlighted events."))
            }
            
            stats <- getDifferentialExpression()
            # optimSurv <- getDifferentialExpressionSurvival()
            # if (!is.null(optimSurv)) {
            #     stats[["Optimal PSI cutoff"]] <- optimSurv[[1]]
            #     stats[["Log rank p-value"]]   <- optimSurv[[2]]
            # }
            
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
        stats    <- getDifferentialExpression()
        filtered <- getDifferentialExpressionFiltered()
        x <- input$xAxis
        y <- input$yAxis
        if (is.null(stats) || is.null(x) || is.null(y)) {
            output$plot <- renderPlot(NULL)
            output$tooltip <- renderUI(NULL)
            return(NULL)
        }
        
        # # Include survival data
        # optimSurv <- getDifferentialExpressionSurvival()
        # if (!is.null(optimSurv)) {
        #     stats[["Optimal PSI cutoff"]] <- optimSurv[[1]]
        #     stats[["Log rank p-value"]]   <- optimSurv[[2]]
        # }
        
        res <- transformData(input, stats, x, y)
        if (is.null(res)) {
            output$plot <- renderPlot(NULL)
            output$tooltip <- renderUI(NULL)
            return(NULL)
        }
        
        stats  <- res$data
        xLabel <- res$xLabel
        yLabel <- res$yLabel
        
        ggplotServer(input=input, output=output, id="ge-volcano", 
                     df=stats, x=xLabel, y=yLabel, plot={
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
                         selected <- getSelectedPoints("ge-volcano")
                         selected <- rownames(filtered)[selected]
                         selected <- which(rownames(stats) %in% selected)
                         if (length(selected) < 1) selected <- NULL
                         
                         events <- getHighlightedPoints("ge-volcano")
                         
                         params <- list(size=input$baseSize, 
                                        col=input$baseColour,
                                        alpha=input$baseAlpha)
                         highlightParams <- list(size=input$highlightedSize,
                                                 col=input$highlightedColour,
                                                 alpha=input$highlightedAlpha)
                         selectedParams  <- list(size=input$selectedSize,
                                                 col=input$selectedColour,
                                                 alpha=input$selectedAlpha)
                         
                         zoom <- getZoom("ge-volcano")
                         if (!is.null(zoom)) {
                             xlim <- c(zoom$xmin, zoom$xmax)
                             ylim <- c(zoom$ymin, zoom$ymax)
                         } else {
                             xlim <- NULL
                             ylim <- NULL
                         }
                         
                         createEventPlotting(
                             stats, xLabel, yLabel, params, highlightX, 
                             highlightY, highlightParams, selected,
                             selectedParams, xlim=xlim, ylim=ylim)
                     })
    })
    
    ggplotAuxServer(input, output, "ge-volcano")
}

#' Set of functions to render data table for differential analyses
#' 
#' @importFrom DT reloadData dataTableProxy dataTableAjax selectRows
#' @importFrom shinyjs enable disable
#' @importFrom utils write.table
#' 
#' @inherit diffExpressionTableServer
diffExpressionTableSet <- function(session, input, output) {
    ns <- session$ns
    
    # Save selected points in the table
    observe({
        selected <- input$statsTable_rows_selected
        setSelectedPoints("ge-volcano", selected)
    })
    
    # Render table with sparklines
    output$statsTable <- renderDataTableSparklines({
        stats <- getDifferentialExpression()
        if (!is.null(stats)) {
            # Discard columns of no interest
            cols <- colnames(stats)
            cols <- cols[!grepl("method|data.name", cols)]
            setDifferentialExpressionColumns(cols)
            return(stats[ , cols])
        }
    }, style="bootstrap", filter="top", server=TRUE, extensions="Buttons",
    options=list(pageLength=10, dom="Bfrtip", buttons=I("colvis"),
                 columnDefs=list(list(targets=5:6, visible=FALSE)),
                 columnDefs=list(list(targets=1, searchable=FALSE))))
    
    # Update table with filtered information
    proxy <- dataTableProxy(ns("statsTable"))
    observe({
        stats <- getDifferentialExpression()
        
        if (!is.null(stats)) {
            # # Bind preview of survival curves based on PSI cutoff
            # optimSurv <- getDifferentialExpressionSurvival()
            # if (!is.null(optimSurv)) {
            #     stats[["Optimal PSI cutoff"]] <- optimSurv[[1]]
            #     stats[["Log rank p-value"]]   <- optimSurv[[2]]
            #     stats[["Survival by PSI cutoff"]] <- optimSurv[[3]]
            # }
            
            # Filter by highlighted events and events in the zoomed area
            events  <- getHighlightedPoints("ge-volcano")
            zoom    <- getZoom("ge-volcano")
            
            zoomed <- NULL
            if (!is.null(zoom)) {
                x <- input$xAxis
                y <- input$yAxis
                if (!is.null(x) && !is.null(y)) {
                    res <- transformData(input, stats, x, y)
                    if (!is.null(res)) {
                        stats  <- res$data
                        xLabel <- res$xLabel
                        yLabel <- res$yLabel
                        
                        xStats <- stats[[xLabel]]
                        xZoom  <- zoom$xmin <= xStats & xStats <= zoom$xmax
                        yStats <- stats[[yLabel]]
                        yZoom  <- zoom$ymin <= yStats & yStats <= zoom$ymax
                        zoomed <- intersect(which(xZoom), which(yZoom))
                    }
                }
            }
            
            # Filter rows based on highlighted and/or zoomed in events
            if (!is.null(events) && !is.null(zoomed)) {
                rowFilter <- intersect(events, zoomed)  
            } else if (!is.null(events)) {
                rowFilter <- events
            } else if (!is.null(zoomed)) {
                rowFilter <- zoomed
            } else {
                rowFilter <- TRUE
            }
            stats <- stats[rowFilter, ]
            
            # Keep previously selected rows if possible
            before <- isolate(getDifferentialExpressionFiltered())
            selected <- isolate(input$statsTable_rows_selected)
            selected <- rownames(before)[isolate(selected)]
            selected <- which(rownames(stats) %in% selected)
            if (length(selected) < 1) selected <- NULL
            
            # Set new data
            setDifferentialExpressionFiltered(stats)
            
            # Properly display event identifiers
            rownames(stats) <- parseSplicingEvent(rownames(stats), char=TRUE)
            
            # Keep columns from data table (else, no data will be rendered)
            cols  <- getDifferentialExpressionColumns()
            stats <- stats[ , cols]
            
            # Check if paging should be reset
            resetPaging <- isolate(getDifferentialExpressionResetPaging())
            if (is.null(resetPaging)) resetPaging <- TRUE
            setDifferentialExpressionResetPaging(TRUE)
            
            # Round numbers based on significant digits
            cols <- colnames(stats)
            type <- sapply(cols, function(i) class(stats[[i]]))
            numericCols <- cols[type == "numeric"]
            
            # Round numbers based on significant digits
            if (nrow(stats) > 0) {
                for (col in numericCols) {
                    stats[ , col] <- suppressWarnings(
                        as.numeric(signifDigits(stats[ , col])))
                }
            }
            
            dataTableAjax(session, stats, outputId="statsTable")
            reloadData(proxy, resetPaging=resetPaging)
            if (!is.null(selected)) selectRows(proxy, selected)
        }
    })
    
    # Disable download button if statistical table is NULL
    observe({
        if ( is.null(getDifferentialExpression()) )
            disable("downloadStats")
        else
            enable("downloadStats")
    })
    
    # Discard columns from data frame containing information to render plots
    discardPlotsFromTable <- function(df) {
        plotCols <- TRUE
        if (!is.null(df)) {
            plotCols <- -match(c("GE distribution", "Survival by GE cutoff"),
                               colnames(df))
            plotCols <- plotCols[!is.na(plotCols)]
            if (length(plotCols) == 0) plotCols <- TRUE
        }
        return(df[ , plotCols])
    }
    
    # Download whole table
    output$downloadAll <- downloadHandler(
        filename=paste(getCategory(), "Differential splicing analyses"),
        content=function(file) {
            stats <- getDifferentialExpression()
            stats <- discardPlotsFromTable(stats)
            
            # # Include updated survival analyses
            # optimSurv <- getDifferentialExpressionSurvival()
            # stats[["Optimal PSI cutoff"]] <- optimSurv[[1]]
            # stats[["Log rank p-value"]]   <- optimSurv[[2]]
            
            write.table(stats, file, quote=FALSE, sep="\t", row.names=FALSE)
        }
    )
    
    # Download filtered table
    output$downloadSubset <- downloadHandler(
        filename=paste(getCategory(), "Differential splicing analyses"),
        content=function(file) {
            stats <- getDifferentialExpressionFiltered()
            stats <- discardPlotsFromTable(stats)
            stats <- stats[input$statsTable_rows_all, ]
            write.table(stats, file, quote=FALSE, sep="\t", row.names=FALSE)
        }
    )
}

#' @rdname appServer
diffExpressionTableServer <- function(input, output, session) {
    ns <- session$ns
    
    selectGroupsServer(session, "diffGroups")
    
    observeEvent(input$loadClinical, 
                 missingDataGuide("Clinical data"))
    observeEvent(input$loadGeneExpr, 
                 missingDataGuide("Gene expression"))
    observeEvent(input$missingGeneExpr, 
                 missingDataGuide("Gene expression"))
    
    diffExpressionSet(session, input, output)
    diffExpressionPlotSet(session, input, output)
    diffExpressionTableSet(session, input, output)
    
    # # Optimal survival difference given a gene expression cutoff per gene
    # optimSurvDiffSet(session, input, output)
}

attr(diffExpressionTableUI, "loader") <- "diffExpression"
attr(diffExpressionTableUI, "name") <- "Exploratory (multiple genes)"
attr(diffExpressionTableUI, "selectEvent") <- FALSE
attr(diffExpressionTableServer, "loader") <- "diffExpression"