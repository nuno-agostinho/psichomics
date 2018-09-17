#' @rdname appUI
#' 
#' @importFrom shinyjs hidden disabled
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
                "Label top differentially expressed genes", value="top",
                checkboxInput(
                    ns("labelTopEnable"), width="100%",
                    "Enable labelling of top differentially expressed genes"),
                div(id=ns("labelTopOptions"),
                    selectizeInput(
                        ns("labelSortBy"), choices=NULL, width="100%",
                        "Sort top differentially expressed genes by"),
                    radioButtons(ns("labelOrder"), "Sorting order", 
                                 choices=c("Decreasing order"="decreasing",
                                           "Increasing order"="increasing")),
                    sliderInput(
                        ns("labelTop"), value=10, min=1, max=1000, 
                        width="100%", "Number of top genes to label"))),
            bsCollapsePanel(
                "Label selected genes", 
                value="genes", checkboxInput(
                    ns("labelGeneEnable"), width="100%",
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
                list(icon("cogs"), "Perform differential expression analysis"),
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
                list(icon("sliders"), "Plot options and table filtering"),
                style="info", value="plotEvents",
                errorDialog(
                    "Differential expression analysis not yet performed.",
                    id=ns("missingDiffAnalyses")),
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
    
    groupCreation <- div(
        class="btn-group dropup",
        tags$button(class="btn btn-default dropdown-toggle", type="button",
                    "data-toggle"="dropdown", "aria-haspopup"="true",
                    "aria-expanded"="false", icon("object-group"),
                    "Create group based on...", tags$span(class="caret")),
        tags$ul(class="dropdown-menu",
                disabled(tags$li(id=ns("groupBySelectedContainer"),
                                 actionLink(ns("groupBySelected"),
                                            "Selected genes"))),
                tags$li(actionLink(ns("groupByDisplayed"),
                                   "Genes displayed in the table"))))
    
    tagList(
        uiOutput(ns("modal")),
        sidebarLayout(
            sidebar, mainPanel(
                ggplotUI(ns("ge-volcano")),
                dataTableOutput(ns("statsTable")),
                hidden(div(id=ns("tableToolbar"), class="btn-toolbar",
                           role="toolbar", downloadTable, groupCreation)),
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
            diffExpr <- getDifferentialExpression()
            groups <- input$diffGroups
        })
        if ( is.null(ge) ) {
            missingDataModal(session, "Gene expression", ns("missingGeneExpr"))
        } else if ( is.null(groups) || length(input$diffGroups) != 2 ) {
            errorModal(session, "Select two groups",
                       "Currently, two groups are required for differential",
                       "expression analysis. Please select two groups.",
                       caller="Differential expression analysis")
        } else if ( !is.null(diffExpr) ) {
            warningModal(session, "Differential expression already performed",
                         "Do you wish to discard the current results?",
                         footer=actionButton(
                             ns("replace"), "Discard", class="btn-warning",
                             "data-dismiss"="modal"),
                         caller="Differential expression analysis")
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

#' @rdname appServer
diffExpressionTableServer <- function(input, output, session) {
    selectGroupsServer(session, "diffGroups", "Samples")
    
    observeEvent(input$loadClinical, missingDataGuide("Clinical data"))
    observeEvent(input$loadGeneExpr, missingDataGuide("Gene expression"))
    observeEvent(input$missingGeneExpr, missingDataGuide("Gene expression"))
    
    diffExpressionSet(session, input, output)
    analysesPlotSet(
        session, input, output, "GE", "ge-volcano", getDifferentialExpression,
        getDifferentialExpressionFiltered, getDifferentialExpressionSurvival)
    analysesTableSet(
        session, input, output, "GE", "ge-volcano", getDifferentialExpression,
        getDifferentialExpressionFiltered, setDifferentialExpressionFiltered, 
        getDifferentialExpressionSurvival, getDifferentialExpressionColumns, 
        setDifferentialExpressionColumns, getDifferentialExpressionResetPaging,
        setDifferentialExpressionResetPaging)
    
    # # Optimal survival difference given a gene expression cutoff per gene
    # optimSurvDiffSet(session, input, output)
}

attr(diffExpressionTableUI, "loader") <- "diffExpression"
attr(diffExpressionTableUI, "name") <- "Exploratory (multiple genes)"
attr(diffExpressionTableUI, "selectEvent") <- FALSE
attr(diffExpressionTableServer, "loader") <- "diffExpression"