## TODO(NunoA): project either individuals (as default) or events
## TODO(NunoA): add histogram in/above percentage of NAs per row to remove
## TODO(NunoA): add brushing capabilities (brush function from Shiny? only
## rectangle selection available?)
## TODO(NunoA): logarithmic values
## TODO(NunoA): BoxCox transformation

## TODO(NunoA): create clusters and use those clusters as groups of data
##
##  km <- kmeans(scores, centers = 7, nstart = 5)
##  ggdata <- data.frame(scores, Cluster=km$cluster, Species=df$Species)
##  ggplot(ggdata) +
##      geom_point(aes(x=PC1, y=PC2, color=factor(Cluster)), size=5, shape=20) +
##      stat_ellipse(aes(x=PC1,y=PC2,fill=factor(Cluster)),
##          geom="polygon", level=0.95, alpha=0.2) +
##      guides(color=guide_legend("Cluster"),fill=guide_legend("Cluster"))
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
    # data <- data[nas/ncol(data)*100 <= naTolerance, , drop = FALSE]
    # if (nrow(data) == 0) return(NULL)
    
    # # Replace NAs with the medians for each individual (row)
    # medians <- rowMedians(data, na.rm=TRUE)
    # data[is.na(data)] <- rep(medians, sum(is.na(data)))
    
    # Get loadings (columns) with less than a given percentage of NAs
    nas <- colSums(is.na(data))
    data <- data[, nas/nrow(data)*100 <= naTolerance, drop = FALSE]
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
                               c("Center values" = "center",
                                 "Scale values" = "scale"),
                               selected = c("center")),
            sliderInput(ns("naTolerance"), div(
                "Percentage of missing values to tolerate per event",
                icon("question-circle")),
                min = 0, max=100, value=0, post="%"),
            bsTooltip(ns("naTolerance"), placement="right", 
                      paste("For events with a tolerable percentage of missing",
                            "values, the median value of the event across",
                            "samples is used to replace those missing values.",
                            "The remaining events are discarded."),
                      options=list(container="body")),
            selectGroupsUI(ns("dataGroups"), "Filter data groups"),
            actionButton(ns("calculate"), class = "btn-primary", 
                         "Calculate PCA"),
            hidden(
                div(id=ns("pcaPlotUI"),
                    hr(),
                    selectizeInput(ns("pcX"), "Choose X axis", choices=NULL),
                    selectizeInput(ns("pcY"), "Choose Y axis", choices=NULL),
                    selectGroupsUI(ns("colourGroups"),
                                   "Clinical groups to colour the PCA"),
                    checkboxGroupInput(ns("plotShow"), "Show in plot",
                                       c("Samples (scores)"="individuals",
                                         "Splicing event (loadings)"="events"),
                                       selected="individuals"),
                    actionButton(ns("showVariancePlot"), "Show variance plot"),
                    actionButton(ns("plot"), class = "btn-primary", "Plot PCA")
                )
            )
        ), mainPanel(
            highchartOutput(ns("scatterplot"))
        )
    )
}

#' Create the explained variance plot
#' 
#' @param pca PCA values
#' 
#' @importFrom highcharter highchart hc_chart hc_title hc_add_series 
#' hc_plotOptions hc_xAxis hc_yAxis hc_legend hc_tooltip hc_exporting
#' 
#' @return Plot variance as an Highcharter object
#' @export
#' @examples 
#' pca <- princomp(USArrests)
#' plotVariance(pca)
plotVariance <- function(pca) {
    sdevSq <- pca$sdev ^ 2
    ns <- paste("PC", seq_along(sdevSq))
    
    hc <- highchart() %>%
        hc_chart(zoomType = "xy", backgroundColor = NULL) %>%
        hc_title(text = paste("Explained variance by each",
                              "Principal Component (PC)")) %>%
        hc_add_series(data = unname(sdevSq), type = "waterfall") %>%
        hc_plotOptions(series = list(dataLabels = list(
            align = "center", verticalAlign = "top", enabled = TRUE,
            formatter = JS(
                "function() {",
                "var total = ", sum(sdevSq), ";",
                "var perc = (this.y/total) * 100;",
                "return (Highcharts.numberFormat(this.y) +'<br/>'+",
                "Highcharts.numberFormat(perc) + '%')}")))) %>%
        hc_xAxis(categories = ns, crosshair = TRUE) %>%
        hc_yAxis(title = list(text = "Explained variance")) %>%
        hc_legend(enabled = FALSE) %>%
        hc_tooltip(pointFormat = '{point.name} {point.y:.2f} {point.perc}') %>%
        hc_exporting(enabled=TRUE, buttons=list(contextButton=list(
            text="Export", y=-50, verticalAlign="bottom", theme=list(fill=NULL)
        )))
    return(hc)
}

#' Create a scatterplot from a PCA object
#' 
#' @param pca \code{prcomp} object
#' @param pcX Character: name of the xAxis of interest from the PCA
#' @param pcY Character: name of the yAxis of interest from the PCA
#' @param clinicalGroups Matrix: groups to plot indicating the index of interest
#' @param individuals Boolean: plot PCA individuals (TRUE by default)
#' @param loadings Boolean: plot PCA loadings/rotations (FALSE by default)
#' 
#' @importFrom highcharter highchart hc_chart hc_xAxis hc_yAxis hc_tooltip %>%
#' hc_add_series_scatter
#' @return Scatterplot as an Highcharter object
#' 
#' @export
#' @examples 
#' pca <- prcomp(USArrests)
#' plotPCA(pca)
#' plotPCA(pca, pcX="PC2", pcY="PC3")
plotPCA <- function(pca, pcX="PC1", pcY="PC2", clinicalGroups=NULL, 
                    individuals=TRUE, loadings=FALSE) {
    imp <- summary(pca)$importance[2, ]
    perc <- as.numeric(imp)
    names(perc) <- names(imp)
    
    label <- sprintf("%s (%s%% explained variance)",
                     names(perc[c(pcX, pcY)]), 
                     roundDigits(perc[c(pcX, pcY)]*100))
    
    hc <- highchart() %>%
        hc_chart(zoomType = "xy") %>%
        hc_xAxis(title = list(text=label[1])) %>%
        hc_yAxis(title = list(text=label[2])) %>%
        hc_tooltip(pointFormat="{point.sample}")
    
    if (individuals) {
        df <- data.frame(pca$x)
        if (is.null(clinicalGroups)) {
            hc <- hc_add_series_scatter(hc, df[[pcX]], df[[pcY]], 
                                        sample=rownames(df))
        } else {
            # Colour data by the selected clinical groups
            for (group in names(clinicalGroups)) {
                rows <- clinicalGroups[[group]]
                values <- df[rows, ]
                if (!all(is.na(values))) {
                    hc <- hc_add_series_scatter(
                        hc, values[[pcX]], values[[pcY]], name=group, 
                        sample=rownames(values), showInLegend=TRUE)
                }
            }
        }
    }
    if (loadings) {
        m <- data.frame(pca$rotation)
        # For loadings, add series (but don't add to legend)
        hc <- hc_add_series_scatter(hc, m[[pcX]], m[[pcY]])
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
pcaServer <- function(input, output, session) {
    ns <- session$ns
    
    selectGroupsServer(session, "dataGroups", "Clinical data")
    selectGroupsServer(session, "colourGroups", "Clinical data")
    
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
        else
            return(NULL)
        
        if (is.null(psi)) {
            missingDataModal(session, "Inclusion levels", ns("takeMeThere"))
        } else {
            isolate({
                selected <- input$dataGroups
                clinicalGroups <- getGroupsFrom("Clinical data")
                clinical <- getClinicalData()
                
                preprocess <- input$preprocess
                naTolerance <- input$naTolerance
            })
            
            # Subset data by the selected clinical groups
            if (!is.null(selected)) {
                ns <- getMatchingSamples(clinicalGroups[selected],
                                         samples=colnames(psi), clinical)
                psi <- psi[ , unlist(ns)]
            }
            
            # Raise error if data has no rows
            if (nrow(psi) == 0) {
                errorModal(session, "No data!", paste(
                    "Calculation returned nothing. Check if everything is as",
                    "expected and try again."))
            }
            
            # Transpose the data to have individuals as rows
            psi <- t(psi)
            
            # Perform principal component analysis (PCA) on the subset data
            pca <- performPCA(psi, naTolerance = naTolerance,
                              center = "center" %in% preprocess,
                              scale. = "scale" %in% preprocess)
            if (is.null(pca)) {
                errorModal(session, "No individuals to plot PCA", 
                           "Try increasing the tolerance of NAs per event")
            } else if (inherits(pca, "error")) {
                ## TODO(NunoA): what to do in this case?
                errorModal(session, "PCA calculation error", 
                           "Constant/zero columns cannot be resized to unit",
                           "variance")
            }
            sharedData$inclusionLevelsPCA <- pca
        }
    })
    
    # Show variance plot
    observeEvent(input$showVariancePlot,
                 infoModal(session, "Variance plot", 
                           highchartOutput(ns("variancePlot")),
                           size = "large"))
    
    # Update select inputs of the principal components
    observe({
        pca <- sharedData$inclusionLevelsPCA
        if (is.null(pca)) {
            updateSelectizeInput(session, "pcX",
                                 choices=c("No PCA performed yet"=""))
            updateSelectizeInput(session, "pcY",
                                 choices=c("No PCA performed yet"=""))
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
        pca <- sharedData$inclusionLevelsPCA
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
            pca <- sharedData$inclusionLevelsPCA
            pcX <- input$pcX
            pcY <- input$pcY
            selected <- input$colourGroups
            show <- input$plotShow
            clinical <- getClinicalData()
            clinicalGroups <- getGroupsFrom("Clinical data")
            psi <- getInclusionLevels()
        })
        
        output$scatterplot <- renderHighchart(
            if (!is.null(pcX) & !is.null(pcY)) {
                if (is.null(selected))
                    groups <- NULL
                else {
                    groups <- clinicalGroups[selected]
                    groups <- getMatchingSamples(groups, colnames(psi),
                                                 clinical)
                }
                plotPCA(pca, pcX, pcY, groups, "individuals" %in% show, 
                        "events" %in% show)
            } else {
                return(NULL)
            })
    })
    
    observe({
        pca <- sharedData$inclusionLevelsPCA
        if (!is.null(pca)) shinyjs::show("pcaPlotUI", animType="fade")
    })
}

attr(pcaUI, "loader") <- "analysis"
attr(pcaUI, "name") <- "Principal Component Analysis (PCA)"
attr(pcaServer, "loader") <- "analysis"