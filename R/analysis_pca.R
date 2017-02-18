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
## Source:
## http://stackoverflow.com/questions/20260434/test-significance-of-clusters-on-a-pca-plot

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

#' User interface of the principal component analysis
#' 
#' @param id Character: identifier
#' 
#' @importFrom highcharter highchartOutput
#' @importFrom shinyBS bsTooltip
#' @importFrom shiny checkboxGroupInput sidebarPanel tagList uiOutput hr
#' sliderInput actionButton selectizeInput
#' @importFrom shinyjs hidden
#' 
#' @return HTML element
pcaUI <- function(id) {
    ns <- NS(id)
    tagList(
        uiOutput(ns("modal")),
        sidebarPanel(
            selectizeInput(ns("dataForPCA"), "Data to perform PCA on",
                           choices=NULL,
                           options=list(placeholder="No data available")),
            checkboxGroupInput(ns("preprocess"), "Preprocessing",
                               c("Center values"="center",
                                 "Scale values"="scale"),
                               selected=c("center")),
            sliderInput(ns("naTolerance"), div(
                "Percentage of missing values to tolerate per event",
                icon("question-circle")),
                min=0, max=100, value=0, post="%"),
            bsTooltip(ns("naTolerance"), placement="right", 
                      paste("For events with a tolerable percentage of missing",
                            "values, the median value of the event across",
                            "samples is used to replace those missing values.",
                            "The remaining events are discarded."),
                      options=list(container="body")),
            selectGroupsUI(ns("dataGroups"), "Samples to use for PCA",
                           noGroupsLabel="All samples",
                           groupsLabel="Samples from selected groups"),
            processButton(ns("calculate"), "Calculate PCA"),
            hidden(
                div(id=ns("pcaPlotUI"),
                    hr(),
                    selectizeInput(ns("pcX"), choices=NULL,
                                   "Choose principal component for the X axis"),
                    selectizeInput(ns("pcY"), choices=NULL,
                                   "Choose principal component for the Y axis"),
                    selectGroupsUI(ns("colourGroups"), "Sample colouring",
                                   noGroupsLabel="Do not colour samples",
                                   groupsLabel="Colour using selected groups"),
                    actionButton(ns("showVariancePlot"), "Show variance plot"),
                    actionButton(ns("plot"), "Plot PCA", class="btn-primary")
                )
            )
        ), mainPanel(
            highchartOutput(ns("scatterplot")),
            highchartOutput(ns("scatterplotLoadings"))
        )
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
#' @param pcX Character: name of the xAxis of interest from the PCA
#' @param pcY Character: name of the yAxis of interest from the PCA
#' @param groups Matrix: groups to plot indicating the index of interest of the
#' samples (use clinical or sample groups)
#' @param individuals Boolean: plot PCA individuals (TRUE by default)
#' @param loadings Boolean: plot PCA loadings/rotations (FALSE by default)
#' 
#' @importFrom highcharter highchart hc_chart hc_xAxis hc_yAxis hc_tooltip %>%
#' @return Scatterplot as an Highcharter object
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
        
        names <- gsub("_", " ", rownames(loadings))
        ## TODO(NunoA): color points with a gradient; see colorRampPalette()
        # For loadings, add series (but don't add to legend)
        hc <- hc_scatter(hc, varCoor[1, ], varCoor[2, ], unname(totalContr), 
                         name="Loadings", sample=names) %>%
            hc_subtitle(text=paste("Bubble size: contribution of a variable",
                                   "to the selected principal components"))
    }
    return(hc)
}

#' Server logic for the principal component analysis
#' 
#' @param input Shiny input
#' @param output Shiny output
#' @param session Shiny session
#' 
#' @importFrom shinyjs runjs hide show
#' @importFrom highcharter %>% hc_chart hc_xAxis hc_yAxis hc_tooltip
#' @importFrom stats setNames
#' @return NULL (this function is used to modify the Shiny session's state)
pcaServer <- function(input, output, session) {
    ns <- session$ns
    
    selectGroupsServer(session, "dataGroups")
    selectGroupsServer(session, "colourGroups")
    
    # Update available data input
    observe({
        inclusionLevels <- getInclusionLevels()
        if (!is.null(inclusionLevels)) {
            updateSelectizeInput(session, "dataForPCA",
                                 choices=attr(inclusionLevels, "dataType"))
        }
    })
    
    observeEvent(input$takeMeThere, missingDataGuide("Inclusion levels"))
    
    # Perform principal component analysis (PCA)
    observeEvent(input$calculate, {
        if (input$dataForPCA == "Inclusion levels")
            psi <- isolate(getInclusionLevels())
        else {
            missingDataModal(session, "Inclusion levels", ns("takeMeThere"))
            return(NULL)
        }
            
        
        if (is.null(psi)) {
            missingDataModal(session, "Inclusion levels", ns("takeMeThere"))
        } else {
            time <- startProcess("calculate")
            isolate({
                groups <- getSelectedGroups(input, "dataGroups", samples=TRUE,
                                            filter=colnames(psi))
                preprocess <- input$preprocess
                naTolerance <- input$naTolerance
            })
            
            # Subset data by the selected clinical groups
            if ( !is.null(groups) ) psi <- psi[ , unlist(groups), drop=FALSE]
            
            # Raise error if data has no rows
            if (nrow(psi) == 0) {
                errorModal(session, "No data!", paste(
                    "PCA returned nothing. Check if everything is as",
                    "expected and try again."))
                endProcess("calculate", closeProgressBar=FALSE)
                return(NULL)
            }
            
            # Transpose the data to have individuals as rows
            psi <- t(psi)
            
            # Perform principal component analysis (PCA) on the subset data
            pca <- performPCA(psi, naTolerance=naTolerance,
                              center="center" %in% preprocess,
                              scale.="scale" %in% preprocess)
            if (is.null(pca)) {
                errorModal(session, "No individuals to plot PCA", 
                           "Try increasing the tolerance of NAs per event")
            } else if (inherits(pca, "error")) {
                ## TODO(NunoA): what to do in this case?
                errorModal(session, "PCA calculation error", 
                           "Constant/zero columns cannot be resized to unit",
                           "variance")
            } else {
                setInclusionLevelsPCA(pca)
            }
            endProcess("calculate", closeProgressBar=FALSE)
        }
    })
    
    # Show variance plot
    observeEvent(input$showVariancePlot,
                 infoModal(session, size="large", "Variance plot",
                           highchartOutput(ns("variancePlot"))))
    
    # Update select inputs of the principal components
    observe({
        pca <- getInclusionLevelsPCA()
        if (is.null(pca)) {
            choices <- c("No PCA performed yet"="")
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
        pca <- getInclusionLevelsPCA()
        if (is.null(pca)) {
            if (input$plot > 0) {
                errorModal(session, "No PCA performed",
                           "Perform a PCA and plot it afterwards.")
            }
            return(NULL)
        }
        plotVariance(pca)
    })
    
    # Plot the principal component analysis
    observeEvent(input$plot, {
        isolate({
            pca <- getInclusionLevelsPCA()
            pcX <- input$pcX
            pcY <- input$pcY
            
            if ( !is.null(pca$x) )
                groups <- getSelectedGroups(input, "colourGroups", samples=TRUE,
                                            filter=rownames(pca$x))
            else
                groups <- NULL
        })
        
        output$scatterplot <- renderHighchart(
            if (!is.null(pcX) & !is.null(pcY)) {
                plotPCA(pca, pcX, pcY, groups) %>% 
                    hc_chart(plotBackgroundColor="#FCFCFC") %>%
                    hc_title(text="Clinical samples (PCA individuals)")
            } else {
                return(NULL)
            })
        
        output$scatterplotLoadings <- renderHighchart(
            if (!is.null(pcX) & !is.null(pcY)) {
                plotPCA(pca, pcX, pcY, individuals=FALSE, loadings=TRUE) %>% 
                    hc_chart(plotBackgroundColor="#FCFCFC") %>%
                    hc_title(
                        text="Alternative splicing events (PCA loadings)") %>%
                    hc_plotOptions(
                        series=list(cursor="pointer", point=list(events=list(
                            click=JS("function() {
                                         sample = this.options.sample;
                                         sample = sample.replace(/ /g, '_');
                                         showDiffSplicing(sample);
                                      }")))))
            } else {
                return(NULL)
            })
    })
    
    observe(
        if (!is.null(getInclusionLevelsPCA()))
            show("pcaPlotUI", animType="fade")
    )
}

attr(pcaUI, "loader") <- "analysis"
attr(pcaUI, "name") <- "Principal Component Analysis (PCA)"
attr(pcaServer, "loader") <- "analysis"