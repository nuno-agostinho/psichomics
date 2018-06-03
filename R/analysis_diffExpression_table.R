#' @rdname appUI
#' 
#' @importFrom shinyjs disabled hidden
#' @importFrom shiny actionLink downloadLink selectizeInput uiOutput tags
#' actionButton checkboxGroupInput helpText tagList sidebarLayout mainPanel
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
        selectizeInput(ns("geneExpr"), "Gene expression dataset", choices=NULL,
                       width="100%"),
        bsCollapse(
            open=c("lmFit"),
            bsCollapsePanel(
                tagList(icon("compress"), "Gene-wise linear model fit"),
                value="lmFit",
                helpText("The", tags$code("limma::lmFit"), "function is used",
                         "to fit a linear model per gene based on a design",
                         "matrix prepared from the two selected groups."),
                selectGroupsUI(
                    ns("diffGroups"), maxItems=2,
                    label="Select two groups for differential expression")),
            bsCollapsePanel(
                tagList(icon("adjust"), "Differential expression statistics"), 
                value="eBayes",
                helpText(
                    "The", tags$code("limma::eBayes"), "function is used to",
                    "compute moderated t-tests and log-odds of differential",
                    "expression by empirical Bayes moderation of the standard",
                    "errors towards a common value."),
                radioButtons(
                    ns("ebayesPriorVar"), "Prior gene-wise variance modelling",
                    list("Constant pooled variance"="constant",
                         "Mean-variance trend (limma-trend)"="trend"),
                    selected="trend"),
                sliderInput(
                    ns("ebayesProportion"), min=0, max=1, value=0.01, step=0.01,
                    width="100%",
                    "Assumed proportion of differentially expressed genes"),
                hr(),
                helpText("Assumed limit for the standard deviation of log2",
                         "fold-changes for differentially expressed genes:"),
                fluidRow(
                    column(6, numericInput(ns("ebayesStdevMin"), "Lower limit",
                                           min=0, value=0.1, step=0.1, 
                                           width="100%")),
                    column(6, numericInput(ns("ebayesStdevMax"), "Upper limit",
                                           min=0, value=4, step=0.1, 
                                           width="100%"))))),
        tags$b("Extra analyses that are performed:"),
        tags$ul(tags$li("Variance and median expression"),
                tags$li("Distribution of gene expression per group")),
        selectizeInput(ns("pvalueAdjust"), selected="BH", width="100%",
                       "P-value adjustment", pvalueAdjust),
        processButton(ns("startAnalyses"), "Perform analyses"))
    
    labelsPanel <- tabPanel(
        "Labels",
        bsCollapse(
            bsCollapsePanel(
                list(icon("tasks"), "Label top differentially expressed genes"),
                value="top", checkboxInput(
                    ns("labelTopEnable"), width="100%",
                    "Enable labelling of top differentially spliced events"),
                div(id=ns("labelTopOptions"),
                    selectizeInput(ns("labelSortBy"), choices=NULL, 
                                   width="100%",
                                   "Sort differentially expressed genes by"),
                    radioButtons(ns("labelOrder"), "Sorting order", 
                                 choices=c("Decreasing order"="decreasing",
                                           "Increasing order"="increasing")),
                    sliderInput(ns("labelTop"), value=10, min=1, max=1000, 
                                width="100%", "Number of top genes to label"))),
            bsCollapsePanel(
                list(icon("tasks"), "Label selected genes"), value="genes",
                checkboxInput(
                    ns("labelGenesEnable"), width="100%",
                    "Enable labelling of selected genes"),
                selectizeInput(ns("labelGenes"), "Genes to label", width="100%",
                               choices=c("Type to search for a gene..."=""),
                               multiple=TRUE))),
        actionButton(ns("unlabelPoints"), "Remove labels"),
        processButton(ns("labelPoints"), "Label points"))
    eventOptions <- prepareEventPlotOptions(ns("eventOptions"), ns, labelsPanel)
    
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
                list(icon("tasks"), "Plot options and table filtering"),
                style="info", value="plotEvents",
                errorDialog(
                    "Differential expression analysis not yet performed.",
                    id=ns("missingDiffExpression")),
                hidden(eventOptions))))
    
    downloadTable <- div(
        class="btn-group dropup",
        tags$button(class="btn btn-default dropdown-toggle", type="button",
                    "data-toggle"="dropdown", "aria-haspopup"="true", 
                    "aria-expanded"="false", icon("download"), 
                    "Save table", tags$span(class="caret")),
        tags$ul(class="dropdown-menu", 
                tags$li(downloadLink(ns("downloadAll"), "All data")),
                tags$li(downloadLink(ns("downloadSubset"), "Filtered data"))))
    
    geneGroupCreation <- div(
        class="btn-group dropup",
        tags$button(class="btn btn-default dropdown-toggle", type="button",
                    "data-toggle"="dropdown", "aria-haspopup"="true",
                    "aria-expanded"="false", icon("object-group"),
                    "Create groups based on...", tags$span(class="caret")),
        tags$ul(class="dropdown-menu",
                disabled(tags$li(id=ns("groupBySelectedGenesContainer"),
                                 actionLink(ns("groupBySelectedGenes"),
                                            "Selected genes"))),
                tags$li(actionLink(ns("groupByDisplayedGenes"),
                                   "Genes displayed in the table"))))
    
    tagList(
        uiOutput(ns("modal")),
        sidebarLayout(
            sidebar, mainPanel(
                ggplotUI(ns("ge-volcano")),
                dataTableOutput(ns("statsTable")),
                hidden(div(id=ns("tableToolbar"), class="btn-toolbar",
                           role="toolbar", downloadTable, geneGroupCreation)),
                highchartOutput(ns("highchartsSparklines"), 0, 0))))
}


#' Set of functions to perform differential analyses
#' 
#' @importFrom shinyBS updateCollapse
#' @importFrom limma eBayes lmFit topTable
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
        groups <- getSelectedGroups(input, "diffGroups", "Samples",
                                    filter=colnames(geneExpr))
        geneExpr     <- geneExpr[ , unlist(groups), drop=FALSE]
        isFromGroup1 <- colnames(geneExpr) %in% groups[[1]]
        design       <- cbind(1, ifelse(isFromGroup1, 1, 0))
        
        # Fit a gene-wise linear model based on selected groups
        fit <- lmFit(geneExpr, design)
        
        # Calculate moderated t-statistics and DE log-odds using limma::eBayes
        ebayesProportion <- input$ebayesProportion
        ebayesStdevMin   <- input$ebayesStdevMin
        ebayesStdevMax   <- input$ebayesStdevMax
        ebayesPriorVar   <- input$ebayesPriorVar
        limmaTrend       <- identical(ebayesPriorVar, "trend")
        
        stats <- eBayes(fit, proportion=ebayesProportion, trend=limmaTrend,
                        stdev.coef.lim=c(ebayesStdevMin, ebayesStdevMax))
        
        # Prepare data summary
        pvalueAdjust <- input$pvalueAdjust
        summary <- topTable(stats, number=nrow(fit), coef=2, sort.by="none",
                            adjust.method=pvalueAdjust, confint=TRUE)
        names(summary) <- c(
            "log2 Fold-Change", "CI (low)", "CI (high)", "Average expression",
            "moderated t-statistics", "p-value", 
            paste0("p-value (", pvalueAdjust, " adjusted)"), "B-statistics")
        attr(summary, "groups") <- groups
        
        # Calculate basic statistics and density plots
        stats  <- diffAnalyses(geneExpr, groups, c("basicStats", "density"),
                               pvalueAdjust=NULL, geneExpr=input$geneExpr)
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
        setLabelledPoints("ge-volcano", NULL)
    })
}

#' Set of functions to plot differential analyses
#' 
#' @inherit diffExpressionTableServer
#' @importFrom shinyjs toggleState
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
                         isolate(input$yAxis), isolate(input$labelSortBy))
    })
    
    # Update genes available to label
    observe({
        diffExpr <- getDifferentialExpression()
        if (!is.null(diffExpr)) {
            updateSelectizeInput(session, "labelGenes", server=TRUE,
                                 choices=rownames(diffExpr), 
                                 selected=character(0))
        } else {
            updateSelectizeInput(session, "labelGenes", server=TRUE,
                                 choices=character(0), selected=character(0))
        }
    })
    
    # Interface elements to highlight values in the plot
    lapply(c("x", "y"), function(axis) {
        observe({
            highlightUI <- function(label, min, max) {
                highlightId <- ns(paste0(label, "Highlight"))
                sliderMinId <- ns(paste0(label, "SliderMin"))
                sliderMaxId <- ns(paste0(label, "SliderMax"))
                sliderInvId <- ns(paste0(label, "SliderInv"))
                
                # Round max and min numbers with two decimal points
                max <- ceiling(max*100)/100
                min <- floor(min*100)/100
                
                conditionalPanel(
                    sprintf("input[id='%s']", highlightId),
                    fluidRow(
                        column(6, numericInput(sliderMinId, "Lower limit",
                                               min=min, value=min, step=0.1,
                                               width="100%")),
                        column(6, numericInput(sliderMaxId, "Upper limit",
                                               max=max, value=max, step=0.1,
                                               width="100%"))),
                    checkboxInput(sliderInvId, "Invert highlighted values"),
                    helpText("The data in the table is also filtered",
                             "according to highlighted events."))
            }
            
            stats <- getDifferentialExpression()
            
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
    
    # Disable labelling elements as appropriate
    observe(toggleState("labelTopOptions", input$labelTopEnable))
    observe(toggleState("labelGenes", input$labelGenesEnable))
    
    # Prepare labelled points
    observeEvent(input$labelPoints, {
        isolate({
            labelTopEnable   <- input$labelTopEnable
            labelGenesEnable <- input$labelGenesEnable
            labelSortBy      <- input$labelSortBy
            labelTop         <- input$labelTop
            labelOrder       <- input$labelOrder == "decreasing"
            genes            <- input$labelGenes
            diffExpr         <- getDifferentialExpression()
        })
        
        labelled <- NULL
        
        # Label top genes
        if (labelTopEnable && !identical(labelSortBy, "")) {
            sorted   <- order(diffExpr[[labelSortBy]], decreasing=labelOrder)
            labelled <- head(sorted, labelTop)
        }
        
        # Label selected genes
        if (labelGenesEnable && !identical(genes, ""))
            labelled <- c(labelled, match(genes, rownames(diffExpr)))
        
        setLabelledPoints("ge-volcano", labelled)
    })
    
    # Unlabel points
    observeEvent(input$unlabelPoints, setLabelledPoints("ge-volcano", NULL))
    
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
        
        res <- transformData(input, stats, x, y)
        if (is.null(res)) {
            output$plot <- renderPlot(NULL)
            output$tooltip <- renderUI(NULL)
            return(NULL)
        }
        
        stats  <- res$data
        xLabel <- res$xLabel
        yLabel <- res$yLabel
        
        ggplotServer(
            input=input, output=output, id="ge-volcano", df=stats, x=xLabel, 
            y=yLabel, plot={
                diffType <- "ge-volcano"
                parseHighlight <- function(input, arg) {
                    argStr <- function(...) paste0(arg, ...)
                    
                    if (!input[[argStr("Highlight")]]) return(NULL)
                    
                    highlightMin <- input[[argStr("SliderMin")]]
                    highlightMax <- input[[argStr("SliderMax")]]
                    
                    noMin <- is.null(highlightMin) || is.na(highlightMin)
                    noMax <- is.null(highlightMax) || is.na(highlightMax)
                    if (noMin || noMax || highlightMin >= highlightMax)
                        return(NULL)
                    
                    highlight <- c(highlightMin, highlightMax)
                    attr(highlight, "inverted") <- input[[argStr("SliderInv")]]
                    return(highlight)
                }
                
                highlightX <- parseHighlight(input, "x")
                highlightY <- parseHighlight(input, "y")
                
                # Check selected events
                selected <- getSelectedPoints(diffType)
                selected <- rownames(filtered)[selected]
                selected <- which(rownames(stats) %in% selected)
                if (length(selected) < 1) selected <- NULL
                
                params <- list(size=input$baseSize, col=input$baseColour,
                               alpha=input$baseAlpha)
                highlightParams <- list(size=input$highlightedSize,
                                        col=input$highlightedColour,
                                        alpha=input$highlightedAlpha)
                selectedParams  <- list(size=input$selectedSize,
                                        col=input$selectedColour,
                                        alpha=input$selectedAlpha)
                labelledParams  <- list(size=input$labelledSize,
                                        col=input$labelledColour,
                                        alpha=input$labelledAlpha)
                
                zoom <- getZoom(diffType)
                if (!is.null(zoom)) {
                    xlim <- c(zoom$xmin, zoom$xmax)
                    ylim <- c(zoom$ymin, zoom$ymax)
                } else {
                    xlim <- NULL
                    ylim <- NULL
                }
                
                labelled <- getLabelledPoints(diffType)
                eventPlot <- createEventPlotting(
                    stats, xLabel, yLabel, params, highlightX, highlightY, 
                    highlightParams, selected, selectedParams, labelled, 
                    labelledParams, xlim=xlim, ylim=ylim)
                setHighlightedPoints(diffType, eventPlot$highlighted)
                eventPlot$plot[[1]]
            })
    })
    
    ggplotAuxServer(input, output, "ge-volcano")
}

#' Set of functions to render data table for differential analyses
#' 
#' @importFrom DT dataTableProxy selectRows replaceData
#' @importFrom shinyjs toggleElement
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
    proxy <- dataTableProxy("statsTable")
    observe({
        stats <- getDifferentialExpression()
        
        if (!is.null(stats)) {
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
            replaceData(proxy, stats, resetPaging=resetPaging, 
                        clearSelection="none")
        }
    })
    
    # Hide table tooltbar if statistical table is not displayed
    observe(toggleElement("tableToolbar",
                          condition=!is.null(getDifferentialExpression())))
    
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
            stats <- cbind("Gene"=rownames(stats), stats)
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
            stats <- cbind("Gene"=rownames(stats), stats)
            write.table(stats, file, quote=FALSE, sep="\t", row.names=FALSE)
        }
    )
    
    # Create groups based on a given filter
    groupsBasedoOnDifferentialExpression <- function(filter, description="") {
        stats <- getDifferentialExpressionFiltered()
        stats <- discardPlotsFromTable(stats)
        stats <- stats[filter, ]
        
        genes <- rownames(stats)
        
        ASevents <- getASevents()
        if ( !is.null(ASevents) )
            ASevents <- getSplicingEventFromGenes(genes, ASevents)
        else
            ASevents <- character(0)
        
        origin <- "Selection from differential expression analysis"
        group <- cbind("Names"="DFS selection", "Subset"=origin, "Input"=origin, 
                       "ASevents"=list(ASevents), "Genes"=list(genes))
        appendNewGroups("ASevents", group)
        infoModal(
            session, "New group created",
            "New group created", description, "and containing:",
            div(style="font-size: 22px;", length(ASevents), "splicing events"),
            div(style="font-size: 22px;", length(genes), "genes"))
    }
    
    # Create groups based on genes displayed in the table
    observeEvent(input$groupByDisplayedGenes,
                 groupsBasedoOnDifferentialExpression(
                     input$statsTable_rows_all,
                     "based on the genes shown in the table"))
    
    # Create groups based on selected genes
    observeEvent(input$groupBySelectedGenes,
                 groupsBasedoOnDifferentialExpression(
                     input$statsTable_rows_selected,
                     "based on selected genes"))
    
    # Disable groups based on selected genes when no groups are selected
    observe(toggleState("groupBySelectedGenesContainer",
                        !is.null(input$statsTable_rows_selected)))
}

#' @rdname appServer
diffExpressionTableServer <- function(input, output, session) {
    ns <- session$ns
    
    selectGroupsServer(session, "diffGroups", "Samples")
    
    observeEvent(input$loadClinical, missingDataGuide("Clinical data"))
    observeEvent(input$loadGeneExpr, missingDataGuide("Gene expression"))
    observeEvent(input$missingGeneExpr, missingDataGuide("Gene expression"))
    
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