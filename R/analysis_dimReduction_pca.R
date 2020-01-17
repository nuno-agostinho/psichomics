## TODO(NunoA): add histogram in/above percentage of NAs per row to remove
## TODO(NunoA): logarithmic values
## TODO(NunoA): BoxCox transformation

#' Perform principal component analysis after processing missing values
#' 
#' @param ... Arguments passed on to \code{stats::prcomp}
#' @inheritParams stats::prcomp
#' @inheritParams reduceDimensionality
#' 
#' @family functions to analyse principal components
#' @return PCA result in a \code{prcomp} object
#' @export
#' 
#' @examples 
#' performPCA(USArrests)
performPCA <- function(data, center=TRUE, scale.=FALSE, 
                       missingValues=round(0.05 * nrow(data)), ...) {
    reduceDimensionality(data, "pca", missingValues=missingValues, 
                         center=center, scale.=scale., ...)
}

#' @rdname appUI
#' 
#' @importFrom highcharter highchartOutput
#' @importFrom shiny checkboxGroupInput tagList uiOutput hr downloadButton
#' sliderInput actionButton selectizeInput helpText textOutput
#' @importFrom shinyBS bsTooltip
#' @importFrom shinyjs hidden
#' @importFrom DT dataTableOutput
pcaUI <- function(id) {
    ns <- NS(id)
    
    pcaOptions <- div(
        id=ns("pcaOptions"),
        selectizeInput(ns("dataForPCA"), "Dataset to perform PCA on", 
                       width="100%", choices=NULL, options=list(
                           placeholder="No data available")),
        checkboxGroupInput(ns("preprocess"), "Preprocessing",
                           c("Center values"="center", "Scale values"="scale"),
                           selected=c("center"), width="100%"),
        selectGroupsUI(ns("dataGroups"), "Perform PCA on...",
                       noGroupsLabel="All samples",
                       groupsLabel="Samples from selected groups"),
        numericInput(ns("missingValues"), div(
            "Number of missing values to tolerate per event",
            icon("question-circle")), min=0, max=100, value=10, width="100%"),
        helpText(textOutput(ns("maxSamples"))),
        bsTooltip(ns("missingValues"), placement="right", paste(
            "For events with a tolerable number of missing values, the median",
            "value of the event across samples is used to replace those",
            "missing values. The remaining events are discarded."),
            options=list(container="body")),
        selectGroupsUI(
            ns("dataGroups2"), "Perform PCA on...",
            noGroupsLabel="All genes and splicing events",
            groupsLabel="Genes and splicing events from selected groups"),
        processButton(ns("calculate"), "Calculate PCA")
    )
    
    performPcaCollapse <- bsCollapsePanel(
        list(icon("sliders"), "Perform PCA"), value="Perform PCA", style="info",
        errorDialog(paste("No alternative splicing quantification or gene",
                          "expression data are available."),
                    id=ns("pcaOptionsDialog"), buttonLabel="Load data",
                    buttonIcon="plus-circle", buttonId=ns("loadData")),
        hidden(pcaOptions))
    
    varsToPlot <- c("all", "top100")
    names(varsToPlot) <- c("All variables",
                           paste("Top 100 variables that most contribute to",
                                 "selected principal components"))
    
    plotPcaCollapse <- bsCollapsePanel(
        list(icon("binoculars"), "Plot PCA"),
        value="Plot PCA", style="info",
        errorDialog("PCA has not yet been performed.", id=ns("noPcaPlotUI")),
        hidden(div(
            id=ns("pcaPlotUI"),
            selectizeInput(ns("pcX"), choices=NULL, width="100%",
                           "Principal component for the X axis"),
            selectizeInput(ns("pcY"), choices=NULL, width="100%",
                           "Principal component for the Y axis"),
            selectGroupsUI(ns("colourGroups"), "Sample colouring",
                           noGroupsLabel="Do not colour samples",
                           groupsLabel="Colour using selected groups"),
            radioButtons(
                ns("plotVariables"), "Variables to plot in loading plot", 
                varsToPlot, selected="top100", width="100%"),
            actionButton(ns("showVariancePlot"), "Show variance plot"),
            actionButton(ns("plot"), "Plot PCA", class="btn-primary"))))
    
    kmeansPanel <- conditionalPanel(
        sprintf("input[id='%s'] == '%s'", ns("clusteringMethod"), "kmeans"),
        sliderInput(ns("kmeansIterations"), 
                    "Maximum number of iterations",
                    min=10, max=100, value=20, width="100%"),
        sliderInput(ns("kmeansNstart"), 
                    "Number of initial random sets",
                    min=50, max=1000, value=100, width="100%"),
        selectizeInput(ns("kmeansMethod"), "K-means method", 
                       width="100%", c("Hartigan-Wong",
                                       "Lloyd-Forgy", "MacQueen")))
    pamPanel <- conditionalPanel(
        sprintf("input[id='%s'] == '%s'", ns("clusteringMethod"), "pam"),
        selectizeInput(ns("pamMetric"), width="100%",
                       "Metric to be used when calculating dissimilarities",
                       c("Euclidean", "Manhattan")))
    
    claraPanel <- conditionalPanel(
        sprintf("input[id='%s'] == '%s'", ns("clusteringMethod"), "clara"),
        selectizeInput(ns("claraMetric"), width="100%",
                       "Metric to be used when calculating dissimilarities",
                       c("Euclidean", "Manhattan")),
        sliderInput(
            ns("claraSamples"), "Samples to be randomly drawn",
            min=10, max=1000, value=50, step=10, width="100%"))
    
    clusteringCollapse <- bsCollapsePanel(
        list(icon("th-large"), "Partitioning clustering"),
        value="Partitioning clustering", style="info",
        errorDialog("PCA has not yet been plotted.",
                 id=ns("noClusteringUI")),
        hidden(
            div(id=ns("clusteringUI"),
                selectizeInput(
                    ns("clusteringMethod"),
                    "Partitioning algorithm", width="100%", selected="clara",
                    c("k-means"="kmeans", 
                      "Partitioning around medoids (PAM)"="pam", 
                      "Clustering Large Applications (CLARA)"="clara")),
                sliderInput(ns("clusterNumber"), "Number of clusters",
                            min=1, max=20, value=2, width="100%"),
                # bsCollapse(
                #     bsCollapsePanel(
                #         tagList(icon("plus-circle"), 
                #                 "Optimal number of clusters"),
                #         value="Optimal number of clusters",
                #         selectizeInput(
                #             ns("estimatationOptimalClusters"), width="100%",
                #             "Method to estimate optimal number of clusters",
                #             c("Within cluster sums of squares"="wss",
                #               "Average silhouette"="silhouette",
                #               "Gap statistics"="gap_stat")),
                #         highchartOutput(ns("optimalClusters")))),
                kmeansPanel, pamPanel, claraPanel,
                actionButton(ns("saveClusters"), "Create groups from clusters"),
                processButton(ns("plotClusters"), "Plot clusters"))))
    
    tagList(
        uiOutput(ns("modal")),
        sidebar(
            bsCollapse(
                id=ns("pcaCollapse"), open="Perform PCA",
                performPcaCollapse,
                plotPcaCollapse,
                clusteringCollapse)
        ), mainPanel(
            highchartOutput(ns("scatterplot")),
            highchartOutput(ns("scatterplotLoadings")),
            hidden(dataTableOutput(ns("varContrTable"))),
            hidden(downloadButton(ns("saveVarContr"), "Save table", "btn-info"))
        )
    )
}

#' Create the explained variance plot from a PCA
#' 
#' @param pca \code{prcomp} object
#' 
#' @importFrom highcharter highchart hc_chart hc_title hc_add_series 
#' hc_plotOptions hc_xAxis hc_yAxis hc_legend hc_tooltip hc_exporting
#' @importFrom shiny tags
#' 
#' @family functions to analyse principal components
#' @return Plot variance as an \code{highchart} object
#' @export
#' @examples 
#' pca <- prcomp(USArrests)
#' plotVariance(pca)
plotVariance <- function(pca) {
    # Get a proportional value to eigenvalues based on standard deviation
    eigenvalue <- unname( pca$sdev ^ 2 )
    variance <- eigenvalue * 100 / sum(eigenvalue)
    cumvar <- cumsum(variance)
    
    # Prepare data
    data <- lapply(seq(eigenvalue), function(i) {
        return(list(y=variance[i], eigenvalue=eigenvalue[i], cumvar=cumvar[i]))
    })
    
    hc <- highchart() %>%
        hc_chart(zoomType="xy", backgroundColor=NULL) %>%
        hc_title(text=paste("Variance explained by each",
                            "Principal Component (PC)")) %>%
        hc_add_series(data=data, type="waterfall", cumvar=cumvar) %>%
        hc_plotOptions(series=list(dataLabels=list(
            format=paste0("{point.eigenvalue:.2f}", tags$br(),
                          "{point.y:.2f}%"),
            align="center", verticalAlign="top", enabled=TRUE))) %>%
        hc_xAxis(title=list(text="Principal Components"), 
                 categories=seq(length(data)), crosshair=TRUE) %>%
        hc_yAxis(title=list(text="Percentage of variance"), min=0, max=100) %>%
        hc_legend(enabled=FALSE) %>%
        hc_tooltip(
            headerFormat=paste(tags$b("Principal component {point.x}"),
                               tags$br()),
            pointFormat=paste0(
                "Eigenvalue: {point.eigenvalue:.2f}", tags$br(),
                "Variance: {point.y:.2f}%", tags$br(),
                "Cumulative variance: {point.cumvar:.2f}%")) %>%
        export_highcharts()
    return(hc)
}

#' Calculate the contribution of PCA loadings to the selected principal
#' components
#'
#' Total contribution of a variable is calculated as per
#' \code{((Cx * Ex) + (Cy * Ey))/(Ex + Ey)}, where:
#' \itemize{
#'   \item{\code{Cx} and \code{Cy} are the contributions of a variable to
#'   principal components \code{x} and \code{y}}
#'   \item{\code{Ex} and \code{Ey} are the eigenvalues of principal components
#'   \code{x} and \code{y}}
#' }
#'
#' @inheritParams plotPCA
#' 
#' @source
#' \url{http://www.sthda.com/english/articles/31-principal-component-methods-in-r-practical-guide/112-pca-principal-component-analysis-essentials/}
#'
#' @family functions to analyse principal components
#' @return Data frame containing the correlation between variables and selected 
#' principal components and the contribution of variables to the selected 
#' principal components (both individual and total contribution)
#' @export
#' 
#' @examples 
#' pca <- performPCA(USArrests)
#' calculateLoadingsContribution(pca)
calculateLoadingsContribution <- function(pca, pcX=1, pcY=2) {
    loadings <- data.frame(pca$rotation)[, c(pcX, pcY)]
    sdev <- pca$sdev[c(pcX, pcY)]
    # Get a proportional value to eigenvalues based on standard deviation
    eigenvalue <- sdev ^ 2
    # Correlation between variables and principal components
    varCorr <- t(loadings) * sdev
    quality <- varCorr ^ 2
    # Total contribution of the variables for the selected PCs
    contr <- quality * 100 / rowSums(quality)
    totalContr <- colSums(contr * eigenvalue) / sum(eigenvalue)
    
    table <- cbind(loadings, t(contr)/colSums(t(contr))*100, 
                   totalContr/sum(totalContr)*100)
    values <- sprintf("PC%s loading", c(pcX, pcY))
    colnames(table) <- c(
        values,
        sprintf("Contribution to PC%s (%%)", c(pcX, pcY)),
        sprintf("Contribution to PC%s and PC%s (%%)", pcX, pcY))
    
    # Parse alternative splicing events or genes
    if ( areSplicingEvents(rownames(table)) ) {
        extra <- parseSplicingEvent(rownames(table), pretty=TRUE)
        extra$gene <- sapply(extra$gene, paste0, collapse=", ")
        extra$pos  <- sapply(extra$pos,  paste0, collapse=", ")
        colnames(extra) <- c("Event type", "Chromosome", "Strand", "Gene",
                             "Event position")
        extra <- extra[ , c(4, 1, 2, 3, 5)]
        table <- cbind(extra, table)
    } else {
        table <- cbind("Genes"=rownames(table), table)
    }
    
    # Sort by total contribution to principal components
    table <- table[order(table[ , ncol(table)], decreasing=TRUE), ]
    table <- cbind("Rank"=seq(nrow(table)), table)
    
    attr(table, "xValues") <- values[1]
    attr(table, "yValues") <- values[2]
    return(table)
}

#' Create a scatterplot from a PCA object
#' 
#' @param pca \code{prcomp} object
#' @param pcX Character: name of the X axis of interest from the PCA
#' @param pcY Character: name of the Y axis of interest from the PCA
#' @param groups Matrix: groups to plot indicating the index of interest of the
#' samples (use clinical or sample groups)
#' @param individuals Boolean: plot PCA individuals
#' @param loadings Boolean: plot PCA loadings/rotations
#' @param nLoadings Integer: Number of variables to plot, ordered by those that 
#' most contribute to selected principal components (this allows for faster 
#' performance as only the most contributing variables are rendered); if 
#' \code{NULL}, all variables are plotted
#' 
#' @importFrom highcharter highchart hc_chart hc_xAxis hc_yAxis hc_tooltip %>%
#' tooltip_table
#' 
#' @family functions to analyse principal components
#' @return Scatterplot as an \code{highchart} object
#' @export
#' 
#' @examples
#' pca <- prcomp(USArrests, scale=TRUE)
#' plotPCA(pca)
#' plotPCA(pca, pcX=2, pcY=3)
#' 
#' # Plot both individuals and loadings
#' plotPCA(pca, pcX=2, pcY=3, loadings=TRUE)
plotPCA <- function(pca, pcX=1, pcY=2, groups=NULL, individuals=TRUE, 
                    loadings=FALSE, nLoadings=NULL) {
    if (is.character(pcX)) pcX <- as.numeric(gsub("[A-Z]", "", pcX))
    if (is.character(pcY)) pcY <- as.numeric(gsub("[A-Z]", "", pcY))
    
    imp <- summary(pca)$importance[2, ]
    perc <- as.numeric(imp)
    
    label <- sprintf("%s (%s%% explained variance)",
                     names(imp[c(pcX, pcY)]), 
                     roundDigits(perc[c(pcX, pcY)]*100))
    
    hc <- highchart() %>%
        hc_chart(zoomType="xy") %>%
        hc_xAxis(title=list(text=label[1]), crosshair=TRUE) %>%
        hc_yAxis(title=list(text=label[2]), gridLineWidth=0,
                 minorGridLineWidth=0, crosshair=TRUE) %>%
        hc_tooltip(pointFormat="{point.sample}") %>%
        export_highcharts()
    
    if (individuals) {
        df <- data.frame(pca$x)
        if (is.null(groups)) {
            hc <- hc_scatter(hc, df[[pcX]], df[[pcY]], sample=rownames(df))
        } else {
            # Colour data based on the selected groups
            for (group in names(groups)) {
                rows <- groups[[group]]
                colour <- attr(groups, "Colour")[[group]]
                values <- df[rows, ]
                if (!all(is.na(values))) {
                    hc <- hc_scatter(
                        hc, values[[pcX]], values[[pcY]], name=group, 
                        sample=rownames(values), showInLegend=TRUE,
                        color=colour)
                }
            }
        }
    }
    if (loadings) {
        contr      <- calculateLoadingsContribution(pca, pcX, pcY)
        if (!is.null(nLoadings)) contr <- head(contr, nLoadings)
        xValues    <- contr[ , attr(contr, "xValues")]
        yValues    <- contr[ , attr(contr, "yValues")]
        contrPCx   <- contr[ , ncol(contr) - 2]
        contrPCy   <- contr[ , ncol(contr) - 1]
        contrTotal <- contr[ , ncol(contr)]
        
        names <- parseSplicingEvent(rownames(contr), char=TRUE)
        dfX <- c(paste0("PC", pcX, " loading"),
                 paste0("PC", pcY, " loading"),
                 paste0("Contribution to PC", pcX),
                 paste0("Contribution to PC", pcY),
                 paste0("Contribution to PC", pcX, " and PC", pcY) )
        dfY <- c(sprintf(" {point.x:.%sf}", getPrecision()),
                 sprintf(" {point.y:.%sf}", getPrecision()),
                 sprintf(" {point.contrPCx:.%sf}%%", getPrecision()),
                 sprintf(" {point.contrPCy:.%sf}%%", getPrecision()),
                 sprintf(" {point.contr:.%sf}%%", getPrecision()))
        ## TODO(NunoA): color points with a gradient; see colorRampPalette()
        # For loadings, add series (but don't add to legend)
        hc <- hc_scatter(hc, xValues, yValues, unname(contrTotal), 
                         name="Loadings", sample=names, contr=contrTotal,
                         contrPCx=contrPCx, contrPCy=contrPCy) %>%
            hc_subtitle(text=sprintf(
                "Bubble size ~ relative contribution to PC%s and PC%s",
                pcX, pcY)) %>%
            hc_tooltip(useHTML=TRUE, headerFormat="", pointFormat=paste0(
                tags$b(style="text-align: center; white-space:pre-wrap;",
                       "{point.sample}"), tags$br(), "<small>",
                tooltip_table(dfX, dfY), "</small>"))
    }
    return(hc)
}

#' Server logic for clustering PCA data
#' 
#' @inheritParams appServer
#' 
#' @importFrom stats kmeans
#' @importFrom cluster pam clara silhouette
#' @importFrom shiny renderTable tableOutput
#' 
#' @inherit psichomics return
#' @keywords internal
clusterSet <- function(session, input, output) {
    clusterPCA <- reactive({
        algorithm <- input$clusteringMethod
        clusters  <- input$clusterNumber
        pca <- getPCA()
        pcX <- input$pcX
        pcY <- input$pcY
        
        if ( !is.null(pca$x) )
            groups <- getSelectedGroups(input, "colourGroups", "Samples",
                                        filter=rownames(pca$x))
        else
            groups <- NULL
        
        if (is.null(pca) || is.null(pcX) || is.null(pcY)) return(NULL)
        pcaScores <- pca$x[ , c(pcX, pcY)]
        
        clustering <- NULL
        if (algorithm == "kmeans") {
            isolate({
                iterations <- input$kmeansIterations
                nstart     <- input$kmeansNstart
                method     <- input$kmeansMethod
            })
            
            if (method == "Lloyd-Forgy") method <- "Lloyd"
            clustering <- kmeans(pcaScores, clusters, iter.max=iterations, 
                                 nstart=nstart, algorithm=method)
            clustering <- clustering$cluster
        } else if (algorithm == "pam") {
            metric     <- tolower(isolate(input$pamMetric))
            clustering <- pam(pcaScores, clusters, metric=metric, 
                              cluster.only=TRUE)
        } else if (algorithm == "clara") {
            isolate({
                metric  <- tolower(input$claraMetric)
                samples <- input$claraSamples
            })
            
            clustering <- clara(pcaScores, clusters, metric=metric,
                                samples=samples, medoids.x=FALSE, 
                                keep.data=FALSE, pamLike=TRUE)
            clustering <- clustering$clustering
        }
        return(clustering)
    })
    
    observeEvent(input$plotClusters, {
        isolate({
            pca <- getPCA()
            pcX <- input$pcX
            pcY <- input$pcY
            
            if ( !is.null(pca$x) )
                groups <- getSelectedGroups(input, "colourGroups", "Samples",
                                            filter=rownames(pca$x))
            else
                groups <- NULL
        })
        
        if (is.null(pca) || is.null(pcX) || is.null(pcY)) return(NULL)
        
        startProcess("plotClusters")
        clustering <- clusterPCA()
        
        hc <- plotPCA(pca, pcX, pcY, groups) %>% 
            plotClusters(pca$x[ , c(pcX, pcY)], clustering) %>% 
            hc_title(text="Clinical samples (PCA scores)") %>%
            hc_legend(symbolHeight=8, symbolWidth=8)
        output$scatterplot <- renderHighchart(hc)
        endProcess("plotClusters")
    })
    
    # # Render optimal clusters
    # output$optimalClusters <- renderHighchart({
    #     algorithm <- input$clusteringMethod
    #     pca <- getPCA()
    #     pcX <- input$pcX
    #     pcY <- input$pcY
    #     
    #     if ( !is.null(pca$x) )
    #         groups <- getSelectedGroups(input, "colourGroups", "Samples",
    #                                     filter=rownames(pca$x))
    #     else
    #         groups <- NULL
    #     
    #     if (is.null(pca) || is.null(pcX) || is.null(pcY)) return(NULL)
    #     pcaScores <- pca$x[ , c(pcX, pcY)]
    #     
    #     clusters <- 1:20
    #     estimation <- input$estimatationOptimalClusters
    #     if (algorithm == "kmeans") {
    #         iterations <- input$kmeansIterations
    #         nstart     <- input$kmeansNstart
    #         method     <- input$kmeansMethod
    #         
    #         if (method == "Lloyd-Forgy") method <- "Lloyd"
    # 
    #         res <- lapply(clusters, function(n) {
    #             kmeans(pcaScores, n, iter.max=iterations, nstart=nstart, 
    #                    algorithm=method)
    #         })
    #     } else if (algorithm == "pam") {
    #         metric <- tolower(input$pamMetric)
    #         res <- lapply(clusters, function(n) {
    #             pam(pcaScores, n, metric=metric, cluster.only=TRUE)
    #         })
    #     } else if (algorithm == "clara") {
    #         metric  <- tolower(input$claraMetric)
    #         samples <- input$claraSamples
    #         
    #         res <- lapply(clusters, function(n) {
    #             clara(pcaScores, n, metric=metric, samples=samples, 
    #                   medoids.x=FALSE, keep.data=FALSE, pamLike=TRUE)
    #         })
    #     }
    #     
    #     if (estimation == "wss") {
    #         withinss <- sapply(res, "[[", "tot.withinss")
    #         hc <- highchart() %>% hc_add_series(withinss) %>%
    #             hc_xAxis(categories=clusters) %>% hc_legend(enabled=FALSE)
    #         return(hc)
    #     } else if (estimation == "silhouette") {
    #         sil     <- silhouette(res)
    #         cluster <- sil[ , 2]
    #         width   <- sil[ , 3]
    #         names(width) <- cluster
    #         hc      <- highchart()
    #         for (i in sort(unique(cluster))) {
    #             hc <- hc %>% 
    #                 hc_add_series(unname(width[names(width) == i]), 
    #                               type="bar") %>%
    #                 hc_xAxis(categories=clusters) %>% 
    #                 hc_legend(enabled=FALSE)
    #         }
    #         return(hc)
    #     }
    # })
    
    # Create data groups from clusters
    observeEvent(input$saveClusters, {
        clustering <- clusterPCA()
        if (!is.null(clustering)) {
            new <- split(names(clustering), clustering)
            names <- paste("Cluster", names(new))
            groups <- cbind("Names"=names, 
                            "Subset"="PCA clustering", "Input"="PCA clustering",
                            "Samples"=new)
            rownames(groups) <- names
            
            # Match samples with subjects (if loaded)
            subjects <- isolate(getSubjectId())
            if (!is.null(subjects)) {
                indiv <- lapply(new, function(i)
                    unname(getSubjectFromSample(i, patientId=subjects)))
                groups <- cbind(groups[ , seq(3), drop=FALSE], "Patients"=indiv,
                                groups[ ,      4, drop=FALSE])
            }
            
            if (!is.null(groups)) appendNewGroups("Samples", groups)
            infoModal(
                session, "Groups successfully created",
                "The following groups were created based on the selected",
                "clustering options.", hr(),
                tableOutput(session$ns("clusteringTable")),
                footer=actionButton(session$ns("goToGroups"), "Show groups",
                                    class="btn-info", "data-dismiss"="modal"))
            
            # Render as table for user
            colnames(groups)[1] <- "Group"
            groups[ , "Samples"]  <- sapply(groups[ , "Samples"], length)
            cols <- c(1, 4)
            if (!is.null(subjects)) {
                groups[ , "Patients"] <- sapply(groups[ , "Patients"], length)
                cols <- c(cols, 5)
            }
            output$clusteringTable <- renderTable(groups[ , cols], digits=0,
                                                  align="c")
        }
    })
    
    observeEvent(input$goToGroups, runjs("showGroups('Samples');"))
}

#' @rdname appServer
#' 
#' @importFrom shiny downloadHandler
#' @importFrom shinyjs runjs hide show
#' @importFrom highcharter %>% hc_chart hc_xAxis hc_yAxis hc_tooltip
#' @importFrom stats setNames
#' @importFrom DT renderDataTable
pcaServer <- function(input, output, session) {
    ns <- session$ns
    
    selectGroupsServer(session, "dataGroups", "Samples")
    selectGroupsServer(session, "dataGroups2", "ASevents")
    selectGroupsServer(session, "colourGroups", "Samples")
    
    observe({
        dataForPCA <- NULL
        selectedDataForPCA <- input$dataForPCA
        if (selectedDataForPCA == "Inclusion levels")
            dataForPCA <- isolate(getInclusionLevels())
        else if (grepl("^Gene expression", selectedDataForPCA))
            dataForPCA <- isolate(getGeneExpression(selectedDataForPCA))
        if (is.null(dataForPCA)) NULL
        
        groups <- getSelectedGroups(input, "dataGroups", "Samples",
                                    filter=colnames(dataForPCA))
        if ( !is.null(groups) ) 
            dataForPCA <- dataForPCA[ , unlist(groups), drop=FALSE]
        
        samples    <- ncol(dataForPCA)
        defaultVal <- round(samples * 0.05) # default: 5% of samples
        updateNumericInput(session, "missingValues", max=samples, 
                           value=defaultVal)
        
        observe({
            missing <- input$missingValues
            text <- sprintf(
                "%s available samples (the selected %s represent %s%%)",
                samples, missing, round(missing / samples * 100))
            output$maxSamples <- renderText(text)
        })
    })
    
    observe({
        incLevels <- getInclusionLevels()
        geneExpr  <- getGeneExpression()
        if (is.null(incLevels) && is.null(geneExpr)) {
            hide("pcaOptions")
            show("pcaOptionsDialog")
        } else {
            show("pcaOptions")
            hide("pcaOptionsDialog")
        }
    })
    
    observe({
        if (!is.null(getPCA())) {
            hide("noPcaPlotUI", animType="fade")
            show("pcaPlotUI", animType="fade")
        } else {
            show("noPcaPlotUI", animType="fade")
            hide("pcaPlotUI", animType="fade")
            
            show("noClusteringUI", animType="fade")
            hide("clusteringUI", animType="fade")
        }
    })
    
    # Update available data input
    observe({
        geneExpr  <- getGeneExpression()
        incLevels <- getInclusionLevels()
        if (!is.null(incLevels) || !is.null(geneExpr)) {
            choices <- c(attr(incLevels, "dataType"), rev(names(geneExpr)))
            updateSelectizeInput(session, "dataForPCA", choices=choices)
        }
    })
    
    observeEvent(input$loadData, missingDataGuide("Inclusion levels"))
    observeEvent(input$takeMeThere, missingDataGuide("Inclusion levels"))
    
    # Perform principal component analysis (PCA)
    observeEvent(input$calculate, {
        selectedDataForPCA <- input$dataForPCA
        if (selectedDataForPCA == "Inclusion levels") {
            dataForPCA  <- isolate(getInclusionLevels())
            dataType    <- "Inclusion levels"
            groups2Type <- "ASevents"
        } else if (grepl("^Gene expression", selectedDataForPCA)) {
            dataForPCA <- isolate(getGeneExpression(selectedDataForPCA))
            dataType   <- "Gene expression"
            groups2Type <- "Genes"
        } else {
            missingDataModal(session, "Inclusion levels", ns("takeMeThere"))
            return(NULL)
        }
        
        if (is.null(dataForPCA)) {
            missingDataModal(session, "Inclusion levels", ns("takeMeThere"))
        } else {
            time <- startProcess("calculate")
            isolate({
                groups <- getSelectedGroups(input, "dataGroups", "Samples",
                                            filter=colnames(dataForPCA))
                groups2 <- getSelectedGroups(input, "dataGroups2", groups2Type, 
                                             filter=rownames(dataForPCA))
                preprocess <- input$preprocess
                missingValues <- input$missingValues
            })
            
            # Subset data based on the selected groups
            if ( !is.null(groups) ) 
                dataForPCA <- dataForPCA[ , unlist(groups), drop=FALSE]
            if ( !is.null(groups2) )
                dataForPCA <- dataForPCA[unlist(groups2), , drop=FALSE]
            
            # Raise error if data has no rows
            if (nrow(dataForPCA) == 0) {
                errorModal(session, "No data returned by PCA",
                           "PCA returned nothing. Check if everything is as",
                           "expected and try again.",
                           caller="Principal component analysis")
                endProcess("calculate", closeProgressBar=FALSE)
                return(NULL)
            }
            
            # Transpose the data to have individuals as rows
            dataForPCA <- t(dataForPCA)
            
            # Perform principal component analysis (PCA) on the subset data
            pca <- performPCA(dataForPCA, missingValues=missingValues,
                              center="center" %in% preprocess,
                              scale.="scale" %in% preprocess)
            if (is.null(pca)) {
                errorModal(session, "No individuals to plot PCA", 
                           "Try increasing the tolerance of missing values",
                           "per event.", caller="Principal component analysis")
            } else if (inherits(pca, "error")) {
                ## TODO(NunoA): what to do in this case?
                errorModal(
                    session, "PCA calculation error", 
                    "Constant/zero columns cannot be resized to unit variance",
                    caller="Principal component analysis")
            } else {
                attr(pca, "dataType") <- dataType
                attr(pca, "firstPCA") <- is.null(getPCA())
                setPCA(pca)
                
                # Clear previously plotted charts
                output$scatterplot <- renderHighchart(NULL)
                output$scatterplotLoadings <- renderHighchart(NULL)
                hide("varContrTable")
                hide("saveVarContr")
            }
            updateCollapse(session, "pcaCollapse", "Plot PCA")
            endProcess("calculate", closeProgressBar=FALSE)
        }
    })
    
    # Update select inputs of the principal components
    observe({
        pca <- getPCA()
        if (is.null(pca)) {
            choices <- c("PCA has not yet been performed"="")
            updateSelectizeInput(session, "pcX", choices=choices)
            updateSelectizeInput(session, "pcY", choices=choices)
            return(NULL)
        }
        
        imp <- summary(pca)$importance[2, ]
        perc <- as.numeric(imp)
        names(perc) <- names(imp)
        
        # Update inputs to select principal components
        label <- sprintf("%s (%s%% explained variance)", 
                         names(perc), roundDigits(perc * 100))
        choices <- setNames(names(perc), label)
        choices <- c(choices, "Select a principal component"="")
        
        updateSelectizeInput(session, "pcX", choices=choices)
        updateSelectizeInput(session, "pcY", choices=choices, 
                             selected=choices[[2]])
    })
    
    # Show variance plot
    observeEvent(input$showVariancePlot,
                 infoModal(session, size="large", "Variance plot",
                           highchartOutput(ns("variancePlot"))))
    
    # Plot the explained variance plot
    output$variancePlot <- renderHighchart({
        pca <- getPCA()
        if (is.null(pca)) {
            if (input$plot > 0) {
                errorModal(session, "PCA has not yet been performed",
                           "Perform a PCA and plot it afterwards.",
                           caller="Principal component analysis")
            }
            return(NULL)
        }
        plotVariance(pca)
    })
    
    # Plot the principal component analysis
    observeEvent(input$plot, {
        isolate({
            pca <- getPCA()
            pcX <- input$pcX
            pcY <- input$pcY
            plotVariables <- input$plotVariables
            
            if ( !is.null(pca$x) )
                groups <- getSelectedGroups(input, "colourGroups", "Samples",
                                            filter=rownames(pca$x))
            else
                groups <- NULL
        })
        
        output$scatterplot <- renderHighchart({
            if (!is.null(pcX) && !is.null(pcY)) {
                plotPCA(pca, pcX, pcY, groups) %>% 
                    hc_title(text="Clinical samples (PCA scores)")
            }
        })
        
        output$scatterplotLoadings <- renderHighchart({
            if (!is.null(pcX) && !is.null(pcY)) {
                dataType <- attr(pca, "dataType")
                if (dataType == "Inclusion levels") {
                    title <- "Alternative splicing events (PCA loadings)"
                    onClick <- sprintf(
                        "function() {
                            sample = this.options.sample;
                            sample = sample.replace(/ /g, '_');
                            showDiffSplicing(sample, %s); }",
                        toJSarray(isolate(names(groups))))
                } else if (dataType == "Gene expression") {
                    title <- "Genes (PCA loadings)"
                    
                    onClick <- sprintf(
                        "function() {
                            sample = this.options.sample;
                            showDiffExpression(sample, %s, '%s'); }",
                        toJSarray(isolate(names(groups))),
                        isolate(input$dataForPCA))
                }
                
                if (plotVariables == "all") nLoadings <- NULL
                else if (plotVariables == "top100") nLoadings <- 100
                
                plotPCA(pca, pcX, pcY, individuals=FALSE, loadings=TRUE,
                        nLoadings=nLoadings) %>%
                    hc_title(text=title) %>%
                    hc_plotOptions(series=list(cursor="pointer", 
                                               point=list(events=list(
                                                   click=JS(onClick)))))
            }
        })
        
        if (is.character(pcX)) pcX <- as.numeric(gsub("[A-Z]", "", pcX))
        if (is.character(pcY)) pcY <- as.numeric(gsub("[A-Z]", "", pcY))
        data <- calculateLoadingsContribution(pca, pcX, pcY)
        
        show("varContrTable")
        output$varContrTable <- renderDataTable(
            data, style="bootstrap", server=TRUE, rownames=FALSE, 
            selection="none", options=list(scrollX=TRUE))
        
        show("saveVarContr")
        output$saveVarContr <- downloadHandler(
            filename=function() {
                paste(getCategory(), "PCA variable contribution")
            }, content=function(con) {
                write.table(data, con, quote=FALSE, sep="\t", row.names=FALSE)
            }
        )
        
        hide("noClusteringUI", animType="fade")
        show("clusteringUI", animType="fade")
        
        updateSliderInput(session, "kmeansNstart", max=nrow(pca$x), value=100)
        updateSliderInput(session, "claraSamples", max=nrow(pca$x), value=50)
    })
    
    clusterSet(session, input, output)
}

attr(pcaUI, "loader") <- "dimReduction"
attr(pcaUI, "name") <- "Principal Component Analysis (PCA)"
attr(pcaServer, "loader") <- "dimReduction"