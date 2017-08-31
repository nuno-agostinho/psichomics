## TODO(NunoA): add histogram in/above percentage of NAs per row to remove
## TODO(NunoA): add brushing capabilities (brush function from Shiny? only
## rectangle selection available?)
## TODO(NunoA): logarithmic values
## TODO(NunoA): BoxCox transformation

## TODO(NunoA): create clusters and use those clusters as groups of data
##
##  km <- kmeans(scores, centers=7, nstart=5)
##  ggdata <- data.frame(scores, Cluster=km$cluster, Species=df$Species)
##  ggplot(ggdata) +
##      geom_point(aes(x=PC1, y=PC2, color=factor(Cluster)), size=5, shape=20) +
##      stat_ellipse(aes(x=PC1,y=PC2,fill=factor(Cluster)),
##          geom="polygon", level=0.95, alpha=0.2) +
##      guides(color=guide_legend("Cluster"), fill=guide_legend("Cluster"))
## 
## Source: http://stackoverflow.com/questions/20260434

#' Perform principal component analysis after processing missing values from 
#' data frame
#' 
#' @inheritParams stats::prcomp
#' @param naTolerance Integer: percentage of NA tolerance
#' @param data Data frame: data
#' 
#' @importFrom stats prcomp
#' @importFrom miscTools rowMedians colMedians
#' 
#' @return PCA result in a \code{prcomp} object
#' @export
#' 
#' @examples 
#' performPCA(USArrests)
performPCA <- function(data, center=TRUE, scale.=FALSE, naTolerance=0) {
    # # Get individuals (rows) with less than a given percentage of NAs
    # nas <- rowSums(is.na(data))
    # # hist(nas/ncol(data)*100)
    # data <- data[nas/ncol(data)*100 <= naTolerance, , drop=FALSE]
    # if (nrow(data) == 0) return(NULL)
    
    # # Replace NAs with the medians for each individual (row)
    # medians <- rowMedians(data, na.rm=TRUE)
    # data[is.na(data)] <- rep(medians, sum(is.na(data)))
    
    # Get loadings (columns) with less than a given percentage of NAs
    nas <- colSums(is.na(data))
    data <- data[, nas/nrow(data) * 100 <= naTolerance, drop=FALSE]
    if (ncol(data) == 0) return(NULL)
    
    # Replace NAs with the medians for each loading (column)
    medians <- colMedians(data, na.rm=TRUE)
    nas <- colSums(is.na(data))
    data[is.na(data)] <- rep(medians, nas)
    
    # Perform principal component analysis
    pca <- tryCatch(prcomp(data, center=center, scale.=scale.), error=return)
    
    # PCA is useless if done with only one point
    if ("x" %in% names(pca) && nrow(pca$x) == 1)
        return(NULL)
    else
        return(pca)
}

#' @rdname appUI
#' 
#' @importFrom highcharter highchartOutput
#' @importFrom shinyBS bsTooltip
#' @importFrom shiny checkboxGroupInput tagList uiOutput hr
#' sliderInput actionButton selectizeInput
#' @importFrom shinyjs hidden
pcaUI <- function(id) {
    ns <- NS(id)
    
    pcaOptions <- div(
        id=ns("pcaOptions"),
        selectizeInput(ns("dataForPCA"), "Data to perform PCA on", width="100%",
                       choices=NULL, options=list(
                           placeholder="No data available")),
        checkboxGroupInput(ns("preprocess"), "Preprocessing",
                           c("Center values"="center", "Scale values"="scale"),
                           selected=c("center"), width="100%"),
        sliderInput(ns("naTolerance"), div(
            "Percentage of missing values to tolerate per event",
            icon("question-circle")),
            min=0, max=100, value=0, post="%", width="100%"),
        bsTooltip(ns("naTolerance"), placement="right", paste(
            "For events with a tolerable percentage of missing",
            "values, the median value of the event across",
            "samples is used to replace those missing values.",
            "The remaining events are discarded."),
            options=list(container="body")),
        selectGroupsUI(ns("dataGroups"), "Samples to use for PCA",
                       noGroupsLabel="All samples",
                       groupsLabel="Samples from selected groups"),
        processButton(ns("calculate"), "Calculate PCA")
    )
    
    tagList(
        uiOutput(ns("modal")),
        sidebar(
            bsCollapse(
                id=ns("pcaCollapse"), open="Perform PCA",
                bsCollapsePanel(
                    list(icon("tasks"), "Perform PCA"),
                    value="Perform PCA", style="info",
                    errorDialog(
                        paste(
                            "No alternative splicing quantification or gene",
                            "expression data are available."),
                        id=ns("pcaOptionsDialog"),
                        buttonLabel="Load data",
                        buttonIcon="plus-circle",
                        buttonId=ns("loadIncLevels")),
                    hidden(pcaOptions)),
                bsCollapsePanel(
                    list(icon("binoculars"), "Plot PCA"),
                    value="Plot PCA", style="info",
                    errorDialog("PCA has not yet been performed.",
                                id=ns("noPcaPlotUI")),
                    hidden(
                        div(id=ns("pcaPlotUI"),
                            selectizeInput(
                                ns("pcX"), choices=NULL, width="100%",
                                "Choose principal component for the X axis"),
                            selectizeInput(
                                ns("pcY"), choices=NULL, width="100%",
                                "Choose principal component for the Y axis"),
                            selectGroupsUI(
                                ns("colourGroups"), "Sample colouring",
                                noGroupsLabel="Do not colour samples",
                                groupsLabel="Colour using selected groups"),
                            actionButton(ns("showVariancePlot"), 
                                         "Show variance plot"),
                            actionButton(ns("plot"), "Plot PCA",
                                         class="btn-primary")))))
        ), mainPanel(uiOutput(ns("pcaPlots")))
    )
}

#' Create the explained variance plot
#' 
#' @param pca PCA values
#' 
#' @importFrom highcharter highchart hc_chart hc_title hc_add_series 
#' hc_plotOptions hc_xAxis hc_yAxis hc_legend hc_tooltip hc_exporting
#' @importFrom shiny tags
#' 
#' @return Plot variance as an Highcharter object
#' @export
#' @examples 
#' pca <- prcomp(USArrests)
#' plotVariance(pca)
plotVariance <- function(pca) {
    eigenvalue <- unname( pca$sdev ^ 2 )
    variance <- eigenvalue * 100 / sum(eigenvalue)
    cumvar <- cumsum(variance)
    ns <- paste("PC", seq_along(eigenvalue))
    
    # Prepare data
    data <- lapply(seq(eigenvalue), function(i) {
        return(list(y=variance[i], eigenvalue=eigenvalue[i], cumvar=cumvar[i]))
    })
    
    hc <- highchart() %>%
        hc_chart(zoomType="xy", backgroundColor=NULL) %>%
        hc_title(text=paste("Explained variance by each",
                            "Principal Component (PC)")) %>%
        hc_add_series(data=data, type="waterfall", cumvar=cumvar) %>%
        hc_plotOptions(series=list(dataLabels=list(
            format=paste0("{point.eigenvalue:.2f}", tags$br(),
                          "{point.y:.2f}%"),
            align="center", verticalAlign="top", enabled=TRUE))) %>%
        hc_xAxis(title=list(text="Principal Components"), 
                 categories=seq(length(data)), crosshair=TRUE) %>%
        hc_yAxis(title=list(text="Percentage of variances"), min=0, max=100) %>%
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

#' Create a scatterplot from a PCA object
#' 
#' @param pca \code{prcomp} object
#' @param pcX Character: name of the X axis of interest from the PCA
#' @param pcY Character: name of the Y axis of interest from the PCA
#' @param groups Matrix: groups to plot indicating the index of interest of the
#' samples (use clinical or sample groups)
#' @param individuals Boolean: plot PCA individuals (TRUE by default)
#' @param loadings Boolean: plot PCA loadings/rotations (FALSE by default)
#' 
#' @importFrom highcharter highchart hc_chart hc_xAxis hc_yAxis hc_tooltip %>%
#' @return Scatterplot as an \code{highcharter} object
#' 
#' @export
#' @examples
#' pca <- prcomp(USArrests, scale=TRUE)
#' plotPCA(pca)
#' plotPCA(pca, pcX=2, pcY=3)
#' 
#' # Plot both individuals and loadings
#' plotPCA(pca, pcX=2, pcY=3, loadings=TRUE)
plotPCA <- function(pca, pcX=1, pcY=2, groups=NULL, individuals=TRUE, 
                    loadings=FALSE) {
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
            # Colour data by the selected clinical groups
            for (group in names(groups)) {
                rows <- groups[[group]]
                values <- df[rows, ]
                if (!all(is.na(values))) {
                    hc <- hc_scatter(
                        hc, values[[pcX]], values[[pcY]], name=group, 
                        sample=rownames(values), showInLegend=TRUE)
                }
            }
        }
    }
    if (loadings) {
        sdev <- pca$sdev[c(pcX, pcY)]
        eigenvalue <- sdev ^ 2
        loadings <- data.frame(pca$rotation)[, c(pcX, pcY)]
        # Correlation between variables and principal components
        varCoor <- t(loadings) * sdev
        quality <- varCoor ^ 2
        # Total contribution of the variables for the selected PCs
        contr <- quality * 100 / rowSums(quality)
        totalContr <- colSums(contr * eigenvalue)
        
        names <- parseSplicingEvent(rownames(loadings), char=TRUE)
        ## TODO(NunoA): color points with a gradient; see colorRampPalette()
        # For loadings, add series (but don't add to legend)
        hc <- hc_scatter(hc, varCoor[1, ], varCoor[2, ], unname(totalContr), 
                         name="Loadings", sample=names) %>%
            hc_subtitle(text=paste("Bubble size: contribution of a variable",
                                   "to the selected principal components"))
    }
    return(hc)
}

#' @rdname appServer
#' 
#' @importFrom shinyjs runjs hide show
#' @importFrom highcharter %>% hc_chart hc_xAxis hc_yAxis hc_tooltip
#' @importFrom stats setNames
pcaServer <- function(input, output, session) {
    ns <- session$ns
    
    selectGroupsServer(session, "dataGroups")
    selectGroupsServer(session, "colourGroups")
    
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
    
    observeEvent(input$loadIncLevels, missingDataGuide("Inclusion levels"))
    observeEvent(input$takeMeThere, missingDataGuide("Inclusion levels"))
    
    # Perform principal component analysis (PCA)
    observeEvent(input$calculate, {
        selectedDataForPCA <- input$dataForPCA
        if (selectedDataForPCA == "Inclusion levels") {
            dataForPCA <- isolate(getInclusionLevels())
            dataType   <- "Inclusion levels"
        } else if (grepl("^Gene expression", selectedDataForPCA)) {
            dataForPCA <- isolate(getGeneExpression()[[selectedDataForPCA]])
            dataType   <- "Gene expression"
        } else {
            missingDataModal(session, "Inclusion levels", ns("takeMeThere"))
            return(NULL)
        }
        
        if (is.null(dataForPCA)) {
            missingDataModal(session, "Inclusion levels", ns("takeMeThere"))
        } else {
            time <- startProcess("calculate")
            isolate({
                groups <- getSelectedGroups(input, "dataGroups", samples=TRUE,
                                            filter=colnames(dataForPCA))
                preprocess <- input$preprocess
                naTolerance <- input$naTolerance
            })
            
            # Subset data by the selected clinical groups
            if ( !is.null(groups) ) 
                dataForPCA <- dataForPCA[ , unlist(groups), drop=FALSE]
            
            # Raise error if data has no rows
            if (nrow(dataForPCA) == 0) {
                errorModal(session, "No data!", paste(
                    "PCA returned nothing. Check if everything is as",
                    "expected and try again."))
                endProcess("calculate", closeProgressBar=FALSE)
                return(NULL)
            }
            
            # Transpose the data to have individuals as rows
            dataForPCA <- t(dataForPCA)
            
            # Perform principal component analysis (PCA) on the subset data
            pca <- performPCA(dataForPCA, naTolerance=naTolerance,
                              center="center" %in% preprocess,
                              scale.="scale" %in% preprocess)
            if (is.null(pca)) {
                errorModal(session, "No individuals to plot PCA", 
                           "Try increasing the tolerance of NAs per event")
            } else if (inherits(pca, "error")) {
                ## TODO(NunoA): what to do in this case?
                errorModal(
                    session, "PCA calculation error", 
                    "Constant/zero columns cannot be resized to unit variance")
            } else {
                attr(pca, "dataType") <- dataType
                attr(pca, "firstPCA") <- is.null(getPCA())
                setPCA(pca)
            }
            updateCollapse(session, "pcaCollapse", "Plot PCA")
            endProcess("calculate", closeProgressBar=FALSE)
        }
    })
    
    # Show variance plot
    observeEvent(input$showVariancePlot,
                 infoModal(session, size="large", "Variance plot",
                           highchartOutput(ns("variancePlot"))))
    
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
    
    # Plot the explained variance plot
    output$variancePlot <- renderHighchart({
        pca <- getPCA()
        if (is.null(pca)) {
            if (input$plot > 0) {
                errorModal(session, "PCA has not yet been performed",
                           "Perform a PCA and plot it afterwards.")
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
            
            if ( !is.null(pca$x) )
                groups <- getSelectedGroups(input, "colourGroups", samples=TRUE,
                                            filter=rownames(pca$x))
            else
                groups <- NULL
        })
        
        output$pcaPlots <- renderUI({
            if (attr(pca, "firstPCA")) {
                tagList(
                    highchartOutput(ns("scatterplot")),
                    highchartOutput(ns("scatterplotLoadings")))
            } else {
                tagList(
                    fluidRow(
                        column(6, highchartOutput(ns("scatterplot"))),
                        column(6, highchartOutput(ns("scatterplot2")))),
                    fluidRow(
                        column(6, highchartOutput(ns("scatterplotLoadings"))),
                        column(6, highchartOutput(ns("scatterplotLoadings2")))))
            }
        })
        
        scatterplot <- reactive({
            if (!is.null(pcX) & !is.null(pcY)) {
                plotPCA(pca, pcX, pcY, groups) %>% 
                    hc_chart(plotBackgroundColor="#FCFCFC") %>%
                    hc_title(text="Clinical samples (PCA scores)")
            } else {
                return(NULL)
            }
        })
        
        scatterplotLoadings <- reactive({
            if (!is.null(pcX) & !is.null(pcY)) {
                dataType <- attr(pca, "dataType")
                if (dataType == "Inclusion levels") {
                    title <- "Alternative splicing events (PCA loadings)"
                    onClick <- sprintf(
                        "function() {
                            sample = this.options.sample;
                            sample = sample.replace(/ /g, '_');
                            showDiffSplicing(sample, %s); }",
                        toJSarray(isolate(groups)))
                } else if (dataType == "Gene expression") {
                    title <- "Genes (PCA loadings)"
                    
                    onClick <- sprintf(
                        "function() {
                            sample = this.options.sample;
                            showDiffExpression(sample, %s, '%s'); }",
                        toJSarray(isolate(groups)),
                        isolate(input$dataForPCA))
                }
                
                plotPCA(pca, pcX, pcY, individuals=FALSE, loadings=TRUE) %>% 
                    hc_chart(plotBackgroundColor="#FCFCFC") %>%
                    hc_title(text=title) %>%
                    hc_plotOptions(series=list(cursor="pointer", 
                                               point=list(events=list(
                                                   click=JS(onClick)))))
            } else {
                return(NULL)
            }
        })
        
        if (attr(pca, "firstPCA")) {
            output[["scatterplot"]] <- renderHighchart(scatterplot())
            output[["scatterplotLoadings"]] <- renderHighchart(
                scatterplotLoadings())
        } else {
            output[["scatterplot2"]] <- renderHighchart(scatterplot())
            output[["scatterplotLoadings2"]] <- renderHighchart(
                scatterplotLoadings())
        }
    })
}

attr(pcaUI, "loader") <- "analysis"
attr(pcaUI, "name") <- "Principal Component Analysis (PCA)"
attr(pcaServer, "loader") <- "analysis"