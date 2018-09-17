#' Perform independent component analysis after processing missing values
#' 
#' @param ... Arguments passed on to \code{fastICA::fastICA}
#' @inheritParams stats::prcomp
#' @inheritParams reduceDimensionality
#' @inheritParams fastICA::fastICA
#' 
#' @return ICA result in a \code{prcomp} object
#' @export
#' 
#' @seealso \code{\link{plotICA}}, \code{\link{performPCA}} and
#' \code{\link{plotPCA}}
#' 
#' @examples 
#' performICA(USArrests)
performICA <- function(data, n.comp=min(5, ncol(data)), center=TRUE, 
                       scale.=FALSE, missingValues=round(0.05 * nrow(data)),
                       alg.typ=c("parallel", "defaltion"),
                       fun=c("logcosh", "exp"), alpha=1.0, ...) {
    alg.typ <- match.arg(alg.typ)
    fun <- match.arg(fun)
    reduceDimensionality(data, "ica", missingValues=missingValues, 
                         center=center, scale.=scale., n.comp=n.comp, 
                         alg.typ=alg.typ, fun=fun, alpha=alpha, ...)
}

#' Create multiple scatterplots from ICA
#' 
#' @param ica Object resulting from \code{\link{performICA}}
#' @param components Numeric: independent components to plot
#' @param groups Matrix: groups to plot indicating the index of interest of the
#' samples (use clinical or sample groups)
#' @inheritDotParams pairsD3::pairsD3 -x
#' 
#' @importFrom pairsD3 pairsD3
#' @return Multiple scatterplots as a \code{pairsD3} object
#' 
#' @export
#' @examples
#' data <- scale(USArrests)
#' ica  <- fastICA::fastICA(data, n.comp=4)
#' plotICA(ica)
#' 
#' # Colour by groups
#' groups <- NULL
#' groups$sunny <- c("California", "Hawaii", "Florida")
#' groups$ozEntrance <- c("Kansas")
#' groups$novel <- c("New Mexico", "New York", "New Hampshire", "New Jersey")
#' plotICA(ica, groups=groups)
plotICA <- function(ica, components=seq(10), groups=NULL, ...) {
    ica <- ica$S
    colnames(ica) <- paste0("IC", seq(ncol(ica)))
    if (!is.null(groups)) {
        ica <- ica[unlist(groups), ]
        groups <- rep(names(groups), sapply(groups, length))
    }
    components <- components[components <= ncol(ica)]
    return(pairsD3(ica[ , components], group=groups, ...))
}

#' Create a scatterplot for ICA
#' 
#' @param ica Object containing an ICA
#' @param icX Character: name of the X axis
#' @param icY Character: name of the Y axis
#' @param groups Matrix: groups to plot indicating the index of interest of the
#' samples (use clinical or sample groups)
#' 
#' @importFrom highcharter highchart hc_chart hc_xAxis hc_yAxis hc_tooltip %>%
#' @return Scatterplot as an \code{highcharter} object
#' 
#' @examples
#' ica <- performICA(USArrests, scale=TRUE)
#' psichomics:::plotSingleICA(ica)
#' psichomics:::plotSingleICA(ica, icX=2, icY=3)
#' 
#' # Colour by groups
#' groups <- NULL
#' groups$sunny <- c("California", "Hawaii", "Florida")
#' groups$ozEntrance <- c("Kansas")
#' groups$novel <- c("New Mexico", "New York", "New Hampshire", "New Jersey")
#' psichomics:::plotSingleICA(ica, groups=groups)
plotSingleICA <- function(ica, icX=1, icY=2, groups=NULL) {
    if (is.character(icX)) icX <- as.numeric(gsub("[A-Z]", "", icX))
    if (is.character(icY)) icY <- as.numeric(gsub("[A-Z]", "", icY))
    
    df <- data.frame(ica$S)
    label <- colnames(df[ , c(icX, icY)])
    
    hc <- highchart(height="400px") %>%
        hc_chart(zoomType="xy") %>%
        hc_xAxis(title=list(text=label[1]), crosshair=TRUE) %>%
        hc_yAxis(title=list(text=label[2]), gridLineWidth=0,
                 minorGridLineWidth=0, crosshair=TRUE) %>%
        hc_tooltip(pointFormat="{point.sample}") %>%
        export_highcharts()
    
    if (is.null(groups)) {
        hc <- hc_scatter(hc, df[[icX]], df[[icY]], sample=rownames(df))
    } else {
        # Colour data based on the selected groups
        for (group in names(groups)) {
            rows <- groups[[group]]
            colour <- attr(groups, "Colour")[[group]]
            values <- df[rows, ]
            if (!all(is.na(values))) {
                hc <- hc_scatter(
                    hc, values[[icX]], values[[icY]], name=group, 
                    sample=rownames(values), showInLegend=TRUE,
                    color=colour)
            }
        }
    }
    return(hc)
}

#' @rdname appUI
#' 
#' @importFrom shinyBS bsTooltip
#' @importFrom shiny checkboxGroupInput tagList uiOutput hr sliderInput 
#' actionButton selectizeInput
#' @importFrom shinyjs hidden
#' @importFrom pairsD3 pairsD3Output
icaUI <- function(id) {
    ns <- NS(id)
    
    icaOptions <- div(
        id=ns("icaOptions"),
        selectizeInput(ns("dataForICA"), "Data to perform ICA on", width="100%",
                       choices=NULL, options=list(
                           placeholder="No data available")),
        sliderInput(ns("componentNumber"), "Number of components",
                    width="100%", value=5, min=2, max=10),
        checkboxGroupInput(ns("preprocess"), "Preprocessing",
                           c("Center values"="center", "Scale values"="scale"),
                           selected=c("center"), width="100%"),
        numericInput(ns("missingValues"), div(
            "Number of missing values to tolerate per event",
            icon("question-circle")), min=0, max=100, value=10, width="100%"),
        bsTooltip(ns("missingValues"), placement="right", paste(
            "For events with a tolerable percentage of missing",
            "values, the median value of the event across",
            "samples is used to replace those missing values.",
            "The remaining events are discarded."),
            options=list(container="body")),
        selectGroupsUI(ns("dataGroups"), "Perform ICA on...",
                       noGroupsLabel="All samples",
                       groupsLabel="Samples from selected groups"),
        selectGroupsUI(
            ns("dataGroups2"), "Perform ICA on...",
            noGroupsLabel="All genes and splicing events",
            groupsLabel="Genes and splicing events from selected groups"),
        processButton(ns("calculate"), "Calculate ICA")
    )
    
    performIcaCollapse <- bsCollapsePanel(
        list(icon("sliders"), "Perform ICA"), value="Perform ICA", style="info",
        errorDialog(paste("No alternative splicing quantification or gene",
                          "expression data are available."),
                    id=ns("icaOptionsDialog"), buttonLabel="Load data",
                    buttonIcon="plus-circle", buttonId=ns("loadData")),
        hidden(icaOptions))
    
    plotIcaCollapse <- bsCollapsePanel(
        list(icon("binoculars"), "Plot ICA"),
        value="Plot ICA", style="info",
        errorDialog("ICA has not yet been performed.", id=ns("noIcaPlotUI")),
        hidden(div(
            id=ns("icaPlotUI"),
            selectizeInput(
                ns("plotComponents"), choices=NULL, width="100%", multiple=TRUE,
                "Independent components to plot (10 maximum)", options=list(
                    maxItems=10, plugins=list('remove_button', 'drag_drop'))),
            selectGroupsUI(ns("colourGroups"), "Sample colouring",
                           noGroupsLabel="Do not colour samples",
                           groupsLabel="Colour using selected groups"),
            bsCollapse(
                bsCollapsePanel(
                    list(icon("paint-brush"), "Plot style"), value="Plot style",
                    sliderInput(ns("plotCex"), "Point size", min=1, max=10, 
                                step=1, value=3, width="100%"),
                    sliderInput(ns("plotOpacity"), "Point opacity", min=0,
                                max=1, step=0.01, value=0.9, width="100%"))),
            actionButton(ns("showVariancePlot"), "Show variance plot"),
            actionButton(ns("plot"), "Plot ICA", class="btn-primary"))))
    
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
        errorDialog("ICA has not yet been plotted.",
                    id=ns("noClusteringUI")),
        hidden(
            div(id=ns("clusteringUI"),
                selectizeInput(
                    ns("clusteringComponents"), choices=NULL, width="100%", 
                    "Indepedent components to cluster (exactly 2)",
                    multiple=TRUE, options=list(
                        maxItems=2, 
                        plugins=list('remove_button', 'drag_drop'))),
                selectizeInput(
                    ns("clusteringMethod"),
                    "Partitioning algorithm", width="100%", selected="clara",
                    c("k-means"="kmeans", 
                      "Partitioning around medoids (PAM)"="pam", 
                      "Clustering Large Applications (CLARA)"="clara")),
                sliderInput(ns("clusterNumber"), "Number of clusters",
                            min=1, max=20, value=2, width="100%"),
                kmeansPanel, pamPanel, claraPanel,
                processButton(ns("saveClusters"), "Create groups from clusters")
            )))
    
    tagList(
        uiOutput(ns("modal")),
        sidebar(
            bsCollapse(
                id=ns("icaCollapse"), open="Perform ICA",
                performIcaCollapse,
                plotIcaCollapse,
                clusteringCollapse) ),
        mainPanel( pairsD3Output(ns("scatterplot"), height="600px") )
    )
}

#' Server logic for clustering ICA data
#' 
#' @inheritParams appServer
#' 
#' @importFrom stats kmeans
#' @importFrom cluster pam clara silhouette
#' @importFrom shiny renderTable tableOutput
#' @importFrom pairsD3 renderPairsD3
#' 
#' @return NULL (this function is used to modify the Shiny session's state)
clusterICAset <- function(session, input, output) {
    clusterICA <- reactive({
        algorithm <- input$clusteringMethod
        clusters  <- input$clusterNumber
        ica <- getICA()
        clusteringComponents <- input$clusteringComponents
        
        if ( !is.null(ica$S) )
            groups <- getSelectedGroups(input, "colourGroups", "Samples",
                                        filter=rownames(ica$S))
        else
            groups <- NULL
        
        if (is.null(ica) || is.null(clusteringComponents) || 
            length(clusteringComponents) < 2) return(NULL)
        icaScores <- ica$S[ , clusteringComponents]
        
        clustering <- NULL
        if (algorithm == "kmeans") {
            isolate({
                iterations <- input$kmeansIterations
                nstart     <- input$kmeansNstart
                method     <- input$kmeansMethod
            })
            
            if (method == "Lloyd-Forgy") method <- "Lloyd"
            clustering <- kmeans(icaScores, clusters, iter.max=iterations, 
                                 nstart=nstart, algorithm=method)
            clustering <- clustering$cluster
        } else if (algorithm == "pam") {
            metric     <- tolower(isolate(input$pamMetric))
            clustering <- pam(icaScores, clusters, metric=metric, 
                              cluster.only=TRUE)
        } else if (algorithm == "clara") {
            isolate({
                metric  <- tolower(input$claraMetric)
                samples <- input$claraSamples
            })
            
            clustering <- clara(icaScores, clusters, metric=metric,
                                samples=samples, medoids.x=FALSE, 
                                keep.data=FALSE, pamLike=TRUE)
            clustering <- clustering$clustering
        }
        return(clustering)
    })
    
    # Create data groups from clusters
    observeEvent(input$saveClusters, {
        clustering <- clusterICA()
        if (!is.null(clustering)) {
            ica <- getICA()
            if ( !is.null(ica$S) )
                selectedGroups <- getSelectedGroups(
                    input, "colourGroups", "Samples",
                    filter=rownames(ica$S))
            else
                selectedGroups <- NULL
            ics <- input$clusteringComponents
            icX <- ics[1]
            icY <- ics[2]

            hc <- plotSingleICA(ica, icX, icY, selectedGroups) %>%
                plotClusters(ica$S[ , c(icX, icY)], clustering) %>%
                hc_legend(symbolHeight=8, symbolWidth=8)
            
            infoModal(
                session, "Groups successfully created", size="medium",
                fluidRow(
                    column(
                        6, "The following groups were created based on the",
                        "selected clustering options. They are available for",
                        "selection and modification from any group selection", 
                        "input.", hr(),
                        tableOutput(session$ns("clusteringTable"))),
                    column(6, hc)))
            
            new <- split(names(clustering), clustering)
            names <- paste("Cluster", names(new))
            groups <- cbind("Names"=names, 
                            "Subset"="ICA clustering", "Input"="ICA clustering",
                            "Samples"=new)
            rownames(groups) <- names
            
            # Match samples with patients (if loaded)
            patients <- isolate(getPatientId())
            if (!is.null(patients)) {
                indiv  <- lapply(new, function(i)
                    unname(getPatientFromSample(i, patientId=patients)))
                groups <- cbind(groups[ , 1:3, drop=FALSE], "Patients"=indiv, 
                                groups[ ,   4, drop=FALSE])
            }
            
            if (!is.null(groups)) appendNewGroups("Samples", groups)
            
            # Render as table for user
            colnames(groups)[1] <- "Group"
            groups[ , "Samples"]  <- sapply(groups[ , "Samples"], length)
            cols <- c(1, 4)
            if (!is.null(patients)) {
                groups[ , "Patients"] <- sapply(groups[ , "Patients"], length)
                cols <- c(cols, 5)
            }
            output$clusteringTable <- renderTable(groups[ , cols], digits=0)
        }
    })
}

#' @rdname appServer
#' 
#' @importFrom shinyjs runjs hide show
#' @importFrom pairsD3 renderPairsD3
#' @importFrom stats setNames
icaServer <- function(input, output, session) {
    ns <- session$ns
    
    selectGroupsServer(session, "dataGroups", "Samples")
    selectGroupsServer(session, "dataGroups2", "ASevents")
    selectGroupsServer(session, "colourGroups", "Samples")
    
    observe({
        incLevels <- getInclusionLevels()
        geneExpr  <- getGeneExpression()
        if (is.null(incLevels) && is.null(geneExpr)) {
            hide("icaOptions")
            show("icaOptionsDialog")
        } else {
            show("icaOptions")
            hide("icaOptionsDialog")
        }
    })
    
    observe({
        if (!is.null(getICA())) {
            hide("noIcaPlotUI", animType="fade")
            show("icaPlotUI", animType="fade")
        } else {
            show("noIcaPlotUI", animType="fade")
            hide("icaPlotUI", animType="fade")
            
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
            updateSelectizeInput(session, "dataForICA", choices=choices)
        }
    })
    
    observeEvent(input$dataForICA, {
        selectedDataForICA <- input$dataForICA
        if (selectedDataForICA == "Inclusion levels")
            dataForICA <- isolate(getInclusionLevels())
        else if (grepl("^Gene expression", selectedDataForICA))
            dataForICA <- isolate(getGeneExpression()[[selectedDataForICA]])
        else return(NULL)
        
        val <- ncol(dataForICA)
        maximum <- min(5, val)
        updateSliderInput(session, "componentNumber", value=maximum, max=val)
    })
    
    observeEvent(input$loadData, missingDataGuide("Inclusion levels"))
    observeEvent(input$takeMeThere, missingDataGuide("Inclusion levels"))
    
    # Perform independent component analysis (ICA)
    observeEvent(input$calculate, {
        selectedDataForICA <- input$dataForICA
        if (selectedDataForICA == "Inclusion levels") {
            dataForICA  <- isolate(getInclusionLevels())
            dataType    <- "Inclusion levels"
            groups2Type <- "ASevents"
        } else if (grepl("^Gene expression", selectedDataForICA)) {
            dataForICA  <- isolate(getGeneExpression()[[selectedDataForICA]])
            dataType    <- "Gene expression"
            groups2Type <- "Genes"
        } else {
            missingDataModal(session, "Inclusion levels", ns("takeMeThere"))
            return(NULL)
        }
        
        if (is.null(dataForICA)) {
            missingDataModal(session, "Inclusion levels", ns("takeMeThere"))
        } else {
            time <- startProcess("calculate")
            isolate({
                groups <- getSelectedGroups(input, "dataGroups", "Samples",
                                            filter=colnames(dataForICA))
                groups2 <- getSelectedGroups(input, "dataGroups2", groups2Type, 
                                             filter=rownames(dataForICA))
                preprocess      <- input$preprocess
                componentNumber <- input$componentNumber
                missingValues   <- input$missingValues
            })
            
            # Subset data based on the selected groups
            if ( !is.null(groups) ) 
                dataForICA <- dataForICA[ , unlist(groups), drop=FALSE]
            if ( !is.null(groups2) )
                dataForICA <- dataForICA[unlist(groups2), , drop=FALSE]
            
            # Raise error if data has no rows
            if (nrow(dataForICA) == 0) {
                errorModal(session, "No data from ICA", 
                           "ICA returned nothing. Check if everything is as",
                           "expected and try again.",
                           caller="Independent component analysis")
                endProcess("calculate", closeProgressBar=FALSE)
                return(NULL)
            }
            
            # Transpose the data to have individuals as rows
            dataForICA <- t(dataForICA)
            
            # Perform independent component analysis (ICA) on the subset data
            ica <- performICA(dataForICA, n.comp=componentNumber, 
                              missingValues=missingValues,
                              center="center" %in% preprocess,
                              scale.="scale" %in% preprocess)
            if (is.null(ica)) {
                errorModal(session, "No individuals to plot ICA", 
                           "Try increasing the tolerance of missing values per",
                           "event", caller="Independent component analysis")
            } else if (inherits(ica, "error")) {
                ## TODO(NunoA): what to do in this case?
                errorModal(
                    session, "ICA calculation error", 
                    "Constant/zero columns cannot be resized to unit variance",
                    caller="Independent component analysis")
            } else {
                attr(ica, "dataType") <- dataType
                attr(ica, "firstICA") <- is.null(getICA())
                setICA(ica)
            }
            updateCollapse(session, "icaCollapse", "Plot ICA")
            endProcess("calculate", closeProgressBar=FALSE)
        }
    })
    
    # Update select inputs of the indepedent components
    observe({
        ica <- getICA()
        if (is.null(ica)) {
            choices <- c("ICA has not yet been performed"="")
            updateSelectizeInput(session, "plotComponents", choices=choices)
            updateSelectizeInput(session, "clusteringComponents", 
                                 choices=choices)
            return(NULL)
        }
        
        choices <- setNames(colnames(ica$S), colnames(ica$S))
        choices <- c(choices, "Select indepedent components"="")
        
        updateSelectizeInput(session, "plotComponents", choices=choices,
                             selected=head(choices, 10))
        updateSelectizeInput(session, "clusteringComponents", choices=choices,
                             selected=head(choices, 2))
    })
    
    # Plot the independent component analysis
    observeEvent(input$plot, {
        isolate({
            ica <- getICA()
            components <- input$plotComponents
            if ( !is.null(ica$S) ) {
                groups <- getSelectedGroups(input, "colourGroups", "Samples",
                                            filter=rownames(ica$S))
                colour <- attr(groups, "Colour")
            } else {
                groups <- NULL
                colour <- "black"
            }
            
            size  <- input$plotCex
            alpha <- input$plotOpacity
        })
        
        output$scatterplot <- renderPairsD3({
            plotICA(ica, components, groups, col=unname(colour), big=TRUE, 
                    leftmar=15, opacity=alpha, cex=size)
        })
        
        hide("noClusteringUI", animType="fade")
        show("clusteringUI", animType="fade")
        
        updateSliderInput(session, "kmeansNstart", max=nrow(ica$S), value=100)
        updateSliderInput(session, "claraSamples", max=nrow(ica$S), value=50)
    })
    
    clusterICAset(session, input, output)
}

attr(icaUI, "loader") <- "dimReduction"
attr(icaUI, "name") <- "Independent Component Analysis (ICA)"
attr(icaServer, "loader") <- "dimReduction"