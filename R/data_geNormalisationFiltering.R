#' Interface to normalise and filter gene expression
#'
#' @param ns Namespace function
#'
#' @importFrom shiny numericInput div column fluidRow tags helpText
#' selectizeInput checkboxInput
#' @importFrom shinyjs hidden
#'
#' @return HTML elements
#' @keywords internal
geNormalisationFilteringInterface <- function(ns) {
    filters <- div(
        id=ns("filteringInterface"),
        fluidRow(
            column(6, numericInput(ns("minMean"), "Mean >=",
                                   min=-1, max=100, value=0, width="100%")),
            column(6, numericInput(ns("minVar"), "Variance >=",
                                   min=-1, max=100, value=0, width="100%"))),
        # fluidRow(
        #     column(6, numericInput(ns("maxMean"), "Max mean",
        #                           min=-1, max=100, value=100, width="100%"))),
        #     column(6, numericInput(ns("maxVar"), "Max variance",
        #                           min=-1, max=100, value=100, width="100%"))),
        fluidRow(
            column(6, numericInput(ns("minCounts"), "Counts >=",
                                   min=0, max=100, value=10, width="100%")),
            column(6, numericInput(ns("minTotalCounts"), "Total counts >=",
                                   min=0, max=100, value=15, width="100%"))),
        helpText(textOutput(ns("filteredGenes"))))

    filteringAssistant <- NULL
    # filteringAssistant <- div(
    #     id=ns("assistantInterface"),
    #     hr(), h4("Filtering assistant"),
    #     selectizeInput(ns("assistantPlot"), "Plot type", width="100%",
    #                    c("Boxplot of the mean expression"="mean",
    #                      "Boxplot of the variance expression"="var")),
    #     highchartOutput(ns("filteringAssistant"), height="150px")
    # )

    options <- div(
        id=ns("options"),
        selectizeInput(ns("geneExpr"), "Gene expression dataset", width="100%",
                       choices=NULL),
        bsCollapse(
            bsCollapsePanel(
                tagList(icon("filter"), "Sample filtering"),
                value="Sample filtering",
                selectizeInput(ns("sampleFilter"), "Samples to discard",
                               multiple=TRUE, width="100%",
                               choices=character(0))),
            bsCollapsePanel(
                tagList(icon("filter"), "Gene filtering"), value="Filtering",
                checkboxInput(ns("enableFiltering"), value=TRUE, width="100%",
                              "Enable gene-wise filtering"),
                filters, filteringAssistant),
            bsCollapsePanel(
                tagList(icon("balance-scale"), "Normalisation"),
                value="Normalisation",
                helpText("Scale raw library sizes using the function",
                         tags$code("edgeR::calcNormFactors"), ", unless the",
                         tags$code("quantile"), "method is selected."),
                selectizeInput(
                    ns("normalisation"), "Normalisation method", width="100%",
                    c("Weighted trimmed mean of M-values (TMM)"="TMM",
                      "Relative log expression (RLE)"="RLE",
                      "Upper-quartile normalisation"="upperquartile",
                      "No normalisation"="none",
                      "Quantile"="quantile"),
                    options=list(render=I('{ option: renderGEnormOptions }'))),
                conditionalPanel(
                    sprintf("input[id='%s'] == '%s'", ns("normalisation"),
                            "upperquartile"),
                    sliderInput(ns("upperquartilePercentile"), width="100%",
                                paste("Percentile of the counts used to",
                                      "calculate scale factors"),
                                min=0, max=1, value=0.75, step=0.01)),
                checkboxInput(ns("voom"), width="100%",
                              "Perform mean-variance modelling using voom"),
                helpText("If library sizes are very different,",
                         tags$code("limma::voom"), "should be more powerful",
                         "and preferred.")),
            bsCollapsePanel(
                tagList(icon("retweet"), "Compute CPM and log-transform"),
                value="Log-transformation",
                helpText("Compute log2-transformed counts per million",
                         "(log2CPM) using", tags$code("edgeR::cpm"), "(or",
                         tags$code("limma::voom"), ", if selected)."),
                numericInput(
                    ns("priorCount"), value=0.25, width="100%",
                    paste("Average count to add to each observation to avoid",
                          "zeroes after log-transformation")),
                helpText())),
        checkboxInput(
            ns("convertToGeneSymbol"), width="100%",
            paste("Replace unambiguous ENSEMBL gene identifiers with their",
                  "gene symbols"), value=TRUE))

    tagList(
        uiOutput(ns("modal")),
        errorDialog("No gene expression data is loaded.",
                    id=ns("missingData"), style="margin: 10px;"),
        hidden(options),
        actionButton(ns("loadGeneExpr"), "Load from file..."),
        disabled(processButton(ns("processGeneExpr"),
                               "Filter and normalise gene expression")))
}

#' @rdname appUI
geNormalisationFilteringUI <- function(id, panel) {
    ns <- NS(id)
    title <- "Gene expression filtering and normalisation"
    panel(style="success", title=tagList(icon("cogs"), title),
          value=title, geNormalisationFilteringInterface(ns))
}

#' Filter and normalise gene expression
#'
#' @param geneExpr Matrix or data frame: gene expression
#' @param geneFilter Boolean: filtered genes
#' @param method Character: normalisation method, including \code{TMM},
#' \code{RLE}, \code{upperquartile}, \code{none} or \code{quantile} (see
#' Details)
#' @inheritParams edgeR::calcNormFactors
#' @param log2transform Boolean: perform log2-transformation?
#' @param priorCount Average count to add to each observation to avoid zeroes
#' after log-transformation
#' @param performVoom Boolean: perform mean-variance modelling (voom)?
#'
#' @details \code{edgeR::calcNormFactors} will be used to normalise gene
#' expression if one of the followin methods is set: \code{TMM}, \code{RLE},
#' \code{upperquartile} or \code{none}. However, \code{limma::voom} will be
#' used for normalisation if \code{performVoom = TRUE} and the selected method
#' is \code{quantile}.
#'
#' @importFrom edgeR DGEList [.DGEList calcNormFactors cpm
#' @importFrom limma voom
#'
#' @return Filtered and normalised gene expression
#' @export
#'
#' @examples
#' geneExpr <- readFile("ex_gene_expression.RDS")
#' normaliseGeneExpression(geneExpr)
normaliseGeneExpression <- function(geneExpr, geneFilter=NULL, method="TMM",
                                    p=0.75, log2transform=TRUE,
                                    priorCount=0.25, performVoom=FALSE) {
    updateProgress("Processing gene expression", divisions=3)

    updateProgress("Filtering gene expression")
    if (is.null(geneFilter)) geneFilter <- TRUE
    else if (!any(geneFilter)) return(NULL)

    originalGeneExpr <- geneExpr
    geneExpr <- DGEList(geneExpr)
    geneExprNorm <- geneExpr[geneFilter, , keep.lib.sizes=TRUE]

    updateProgress("Normalising gene expression")
    if (!performVoom && method == "quantile") method <- "none"
    if (method != "quantile")
        geneExprNorm <- calcNormFactors(geneExprNorm, method=method, p=p)

    if (!performVoom) {
        geneExprNorm <- cpm(geneExprNorm, log=log2transform,
                            prior.count=priorCount)
    } else {
        norm <- if (method == "quantile") "quantile" else "none"
        geneExprNorm <- voom(geneExprNorm, normalize.method=norm)
    }

    updateProgress("Preparing gene expression data")
    if (!is(geneExprNorm, "EList")) geneExprNorm <- data.frame(geneExprNorm)
    colnames(geneExprNorm) <- colnames(geneExpr)
    
    geneExprNorm <- inheritAttrs(geneExprNorm, originalGeneExpr)
    if (is(geneExprNorm, "EList"))
        geneExprNorm$E <- inheritAttrs(geneExprNorm$E, originalGeneExpr)
    return(geneExprNorm)
}

#' Set of functions to load splicing quantification
#'
#' @inherit geNormalisationFilteringServer
#'
#' @importFrom shiny tags
#' @importFrom shinyBS bsPopover
#'
#' @keywords internal
loadGeneExpressionSet <- function(session, input, output) {
    ns <- session$ns

    # Show modal for loading gene expression data
    observeEvent(input$loadGeneExpr, {
        infoModal(
            session, "Load gene expression data",
            geneExprFileInput(ns("customGeneExpr")),
            uiOutput(ns("alertGeneExpr")),
            footer=processButton(ns("loadCustomGE"), "Load quantification"))
    })

    observeEvent(input$loadGeneExpr, {
        prepareFileBrowser(session, input, "customGeneExpr")
    }, once=TRUE)

    observeEvent(input$loadCustomGE, loadGeneExpression())

    # Load alternative splicing quantification
    loadGeneExpression <- reactive({
        time <- startProcess("loadGeneExpr")

        updateProgress("Wait a moment", divisions=2)
        updateProgress("Loading gene expression")

        allFormats <- loadFileFormats()
        formats <- allFormats[sapply(allFormats, "[[",
                                     "dataType") == "Gene expression"]

        geneExpr <- tryCatch(parseValidFile(input$customGeneExpr, formats),
                             warning=return, error=return)
        if (is(geneExpr, "error")) {
            if (geneExpr$message == paste("'file' must be a character string",
                                          "or connection"))
                errorAlert(session, title="No file provided",
                           "Please provide a file", alertId="alertGeneExpr",
                           caller="Gene expression normalisation and filtering")
            else
                errorAlert(session, title="An error was raised",
                           geneExpr$message, alertId="alertGeneExpr",
                           caller="Gene expression normalisation and filtering")
        } else if (is(geneExpr, "warning")) {
            warningAlert(session, title="A warning was raised",
                         geneExpr$message, alertId="alertGeneExpr",
                         caller="Gene expression normalisation and filtering")
        } else {
            removeAlert(output, "alertGeneExpr")

            if ( is.null(getData()) ) {
                name <- file_path_sans_ext( basename(input$customGeneExpr) )
                name <- gsub(" Gene expression.*$", "", name)
                if (name == "") name <- "Unnamed"

                data <- setNames(list(list("Gene expression"=geneExpr)), name)
                data <- processDatasetNames(data)
                setData(data)
                setCategory(name)

                samples <- colnames(geneExpr)
                parsed <- parseTcgaSampleInfo(samples)
                if ( !is.null(parsed) ) setSampleInfo(parsed)
            } else {
                name <- renameDuplicated("Gene expression",
                                         names(getCategoryData()))
                setDataTable(name, geneExpr)
            }
            removeModal()
        }
        endProcess("loadGeneExpr", time)
    })
}

#' Convert gene identifiers
#'
#' @param annotation OrgDb: genome wide annotation for an organism, e.g.
#' \code{org.Hs.eg.db}
#' @param genes Character: genes to be converted
#' @param key Character: type of identifier used, e.g. \code{ENSEMBL}; read
#' \code{?AnnotationDbi::columns}
#' @param target Character: type of identifier to convert to; read
#' \code{?AnnotationDbi::columns}
#' @param ignoreDuplicatedTargets Boolean: if \code{TRUE}, identifiers that
#' share targets with other identifiers will not be converted
#'
#' @importFrom AnnotationDbi select
#' @importFrom data.table data.table
#' @importFrom org.Hs.eg.db org.Hs.eg.db
#'
#' @return Character vector of the respective targets of gene identifiers. The
#' previous identifiers remain other identifiers have the same target (in case
#' \code{ignoreDuplicatedTargets = TRUE}) or if no target was found.
#' @export
#'
#' @examples
#' if ( require("org.Hs.eg.db") ) {
#'     columns(org.Hs.eg.db)
#'
#'     genes <- c("ENSG00000012048", "ENSG00000083093", "ENSG00000141510",
#'                "ENSG00000051180")
#'     convertGeneIdentifiers(org.Hs.eg.db, genes,
#'                            key="ENSEMBL", target="SYMBOL")
#' }
convertGeneIdentifiers <- function(annotation, genes, key="ENSEMBL",
                                   target="SYMBOL",
                                   ignoreDuplicatedTargets=TRUE) {
    stopifnot(is(annotation, "OrgDb"))

    if (key == "ENSEMBL") {
        # Remove ENSEMBL identifiers
        genesClean <- gsub("\\..*", "", genes)
        # Keep version for gene identifier containing the string "PAR_Y"
        par_y <- grep("PAR", genes)
        genesClean[par_y] <- genes[par_y]
    } else {
        genesClean <- genes
    }

    match <- tryCatch(
        suppressMessages(select(annotation, genesClean, target, key)),
        error=return)

    if (is(match, "error")) return(setNames(genes, genes))
    match <- data.table(match, key=key)

    # Ignore missing values
    match <- match[!is.na(match[[target]]), ]

    # Collapse genes with more than one matching target
    colnames(match)[2] <- "target"
    collapsed <- match[
        , list(target=paste(unique(target), collapse="/")), by=key]

    if (ignoreDuplicatedTargets) {
        # Ignore genes sharing the same target
        geneTargets <- collapsed[["target"]]
        collapsed   <- collapsed[
            !geneTargets %in% unique(geneTargets[duplicated(geneTargets)]), ]
    }

    # Replace identifiers by their matching targets (if possible)
    converted <- collapsed[["target"]][match(genesClean, collapsed[[key]])]
    genes[!is.na(converted)] <- converted[!is.na(converted)]
    names(genes) <- genesClean
    return(genes)
}

#' Filter genes based on their expression
#'
#' @param geneExpr Data frame or matrix: gene expression
#' @param minMean Numeric: minimum of read count mean per gene
#' @param maxMean Numeric: maximum of read count mean per gene
#' @param minVar Numeric: minimum of read count variance per gene
#' @param maxVar Numeric: maximum of read count variance per gene
#' @param minCounts Numeric: minimum number of read counts per gene for at least
#' some samples
#' @param minTotalCounts Numeric: minimum total number of read counts per gene
#'
#' @importFrom edgeR filterByExpr
#'
#' @return Boolean vector indicating which genes have sufficiently large counts
#' @export
#' 
#' @examples 
#' geneExpr <- readFile("ex_gene_expression.RDS")
#' 
#' # Add some genes with low expression
#' geneExpr <- rbind(geneExpr, 
#'                   lowReadGene1=c(rep(4:5, 10)),
#'                   lowReadGene2=c(rep(5:1, 10)),
#'                   lowReadGene3=c(rep(10:1, 10)),
#'                   lowReadGene4=c(rep(7:8, 10)))
#' 
#' # Filter out genes with low reads across samples
#' geneExpr[filterGeneExpr(geneExpr), ]
filterGeneExpr <- function(geneExpr, minMean=0, maxMean=Inf, minVar=0,
                           maxVar=Inf, minCounts=10, minTotalCounts=15) {
    geneExprMean <- rowMeans(geneExpr)
    geneExprVar  <- rowVars(geneExpr)

    varMeanFilter <- geneExprMean >= minMean & geneExprMean <= maxMean &
        geneExprVar >= minVar & geneExprVar <= maxVar

    lowCountFilter <- filterByExpr(geneExpr[varMeanFilter, ],
                                   min.count=minCounts,
                                   min.total.count=minTotalCounts)
    filteredGenes <- varMeanFilter
    filteredGenes[names(lowCountFilter[!lowCountFilter])] <- FALSE
    return(filteredGenes)
}

#' Plot distribution of gene expression per sample
#'
#' @param geneExpr Data frame or matrix: gene expression
#' @inheritDotParams renderBoxplot
#'
#' @importFrom highcharter %>% hc_yAxis
#'
#' @return Gene expression distribution plots
#' @export
#'
#' @examples
#' df <- data.frame(geneA=c(2, 4, 5),
#'                  geneB=c(20, 3, 5),
#'                  geneC=c(5, 10, 21))
#' colnames(df) <- paste("Sample", 1:3)
#' plotGeneExprPerSample(df)
plotGeneExprPerSample <- function(geneExpr, ...) {
    if (is(geneExpr, "EList")) geneExpr <- geneExpr$E
    renderBoxplot(geneExpr, ...) %>%
        hc_yAxis(title=list(text="Gene expression"))
}

#' Sum columns using an \code{\link{EList-class}} object
#' @inheritParams base::colSums
#' 
#' @return Numeric vector with the sum of the columns
#' @export
setMethod("colSums", signature="EList", function(x, na.rm=FALSE, dims=1) {
    colSums(x$E, na.rm=na.rm, dims=dims)
})

#' @rdname appServer
#'
#' @importFrom shiny reactive observeEvent helpText removeModal
#' updateNumericInput
#' @importFrom tools file_path_sans_ext
#' @importFrom shinyjs enable disable hide show
#' @importFrom data.table fread
#' @importFrom highcharter hcboxplot hc_plotOptions hc_xAxis hc_chart
#'
#' @keywords internal
geNormalisationFilteringServer <- function(input, output, session) {
    ns <- session$ns
    observeEvent(input$missing, missingDataGuide("Gene expression"))

    # Warn user if gene expression is not loaded
    observe({
        if (is.null(getGeneExpression())) {
            hide("options")
            disable("processGeneExpr")
            show("missingData")
        } else {
            show("options")
            enable("processGeneExpr")
            hide("missingData")
        }
    })

    # Update available gene expression data according to loaded files
    observe({
        geneExpr <- getGeneExpression()
        if (!is.null(geneExpr)) {
            updateSelectizeInput(session, "geneExpr",
                                 choices=c(names(geneExpr),
                                           "Select gene expression data"=""))
        } else {
            updateSelectizeInput(
                session, "geneExpr",
                choices=c("No gene expression data loaded"=""))
        }
    })

    getFilter <- reactive({
        geneExpr <- isolate(input$geneExpr)
        if (is.null(geneExpr) || geneExpr == "") return(NULL)
        geneExpr <- isolate(getGeneExpression(geneExpr))

        minMean        <- input$minMean
        maxMean        <- Inf
        minVar         <- input$minVar
        maxVar         <- Inf
        minCounts      <- input$minCounts
        minTotalCounts <- input$minTotalCounts

        sampleFilter   <- input$sampleFilter

        if (is.na(minMean) || is.na(maxMean) ||
            is.na(minVar) || is.na(maxVar) ||
            is.na(minCounts) || is.na(minTotalCounts)) {
            return(NULL)
        } else {
            if (!is.null(sampleFilter) && sampleFilter != "") {
                samplesToKeep <- !colnames(geneExpr) %in% sampleFilter
                geneExpr <- geneExpr[ , samplesToKeep]
            }
            
            filtered <- filterGeneExpr(geneExpr, minMean, maxMean, minVar,
                                       maxVar, minCounts, minTotalCounts)
            return(filtered)
        }
    })

    output$filteredGenes <- renderText({
        geneExpr <- input$geneExpr
        if (is.null(geneExpr) || geneExpr == "") return(NULL)
        geneExpr <- isolate(getGeneExpression(geneExpr))

        filter <- sum(getFilter())
        total  <- nrow(geneExpr)
        ratio  <- filter/total * 100

        if (input$enableFiltering) {
            msg <- sprintf("Selecting %s genes (%s%%) out of %s.",
                           filter, round(ratio), total)
        } else {
            msg <- sprintf("Selecting all %s genes.", total)
        }
        return(msg)
    })

    # Update sample filtering options
    observeEvent(input$geneExpr, {
        geneExpr <- isolate(input$geneExpr)
        if (is.null(geneExpr) || geneExpr == "") return(NULL)
        geneExpr <- isolate(getGeneExpression(geneExpr))

        updateSelectizeInput(
            session, "sampleFilter", server=TRUE,
            choices=colnames(geneExpr),
            options=list(placeholder="Select samples to discard", 
                         plugins=list("remove_button")))
    })

    # Update filtering options based on selected gene expression data
    observeEvent(input$geneExpr, {
        geneExpr <- isolate(input$geneExpr)
        if (is.null(geneExpr) || geneExpr == "") return(NULL)
        geneExpr     <- isolate(getGeneExpression(geneExpr))

        sampleFilter <- isolate(input$sampleFilter)
        if (!is.null(sampleFilter) && sampleFilter != "") {
            samplesToKeep <- !colnames(geneExpr) %in% sampleFilter
            geneExpr <- geneExpr[ , samplesToKeep]
        }

        # Update mean range
        geneExprMean <- rowMeans(geneExpr)
        maxMean      <- max(geneExprMean, na.rm=TRUE)
        updateNumericInput(session, "minMean", max=maxMean)
        # updateNumericInput(session, "maxMean", max=maxMean, value=maxMean)

        # Update variance range
        geneExprVar <- rowVars(geneExpr)
        maxVar      <- max(geneExprVar, na.rm=TRUE)
        updateNumericInput(session, "minVar", max=maxVar)
        # updateNumericInput(session, "maxVar", max=maxVar, value=maxVar)

        # output$filteringAssistant <- renderHighchart({
        #     type <- input$assistantPlot
        #     if (type == "") return(NULL)
        #
        #     filter <- getFilter()
        #     if (type == "mean") {
        #         arg <- geneExprMean[filter]
        #         description <- "Mean per gene"
        #     } else if (type == "var") {
        #         arg <- geneExprVar[filter]
        #         description <- "Variance per gene"
        #     }
        #
        #     hc <- tryCatch({
        #         hcboxplot(arg)  %>%
        #             hc_chart(zoomType="y") %>%
        #             hc_tooltip(valueDecimals=2, followPointer=TRUE) %>%
        #             hc_xAxis(visible=FALSE) %>%
        #             hc_yAxis(title=list(text=description)) %>%
        #             hc_plotOptions(series=list(animation=FALSE))
        #     }, error=return, warning=return)
        #
        #     if (is(hc, "error") || is(hc, "warning"))
        #         return(NULL)
        #     else
        #         return(hc)
        # })
    })

    # Disable interface for gene filtering
    observeEvent(input$enableFiltering, {
        filter <- input$enableFiltering
        if (filter) {
            enable("filteringInterface")
            # show("assistantInterface", anim=TRUE)
        } else {
            disable("filteringInterface")
            # hide("assistantInterface", anim=TRUE)
        }
    })

    # Disable option to add counts to observations if not log2-transforming
    observe({
        # filter <- input$log2transformation
        # if (filter) {
        enable("priorCount")
        # show("assistantInterface", anim=TRUE)
        # } else {
        # disable("priorCount")
        # hide("assistantInterface", anim=TRUE)
        # }
    })

    # Filter and normalise gene expression
    observeEvent(input$processGeneExpr, {
        time <- startProcess("processGeneExpr")

        isolate({
            geneExpr      <- getGeneExpression(input$geneExpr)
            method        <- input$normalisation
            percentile    <- input$upperquartilePercentile
            sampleFilter  <- input$sampleFilter
            filter        <- input$enableFiltering
            priorCount    <- input$priorCount

            minMean        <- input$minMean
            maxMean        <- Inf # input$maxMean
            minVar         <- input$minVar
            maxVar         <- Inf # input$maxVar
            minCounts      <- input$minCounts
            minTotalCounts <- input$minTotalCounts

            voom <- input$voom

            convertToGeneSymbol <- input$convertToGeneSymbol
        })

        if (!is.null(sampleFilter) && sampleFilter != "") {
            samplesToKeep <- !colnames(geneExpr) %in% sampleFilter
            geneExpr <- geneExpr[ , samplesToKeep]
        }

        if (filter) {
            geneFilter <- getFilter()
        } else {
            geneFilter <- NULL
        }
        
        geneExprNorm <- normaliseGeneExpression(
            geneExpr, geneFilter, method, percentile, log2transform=TRUE,
            priorCount, performVoom=voom)

        if (convertToGeneSymbol) {
            rownames(geneExprNorm) <- convertGeneIdentifiers(
                org.Hs.eg.db, rownames(geneExprNorm))
        }

        attr(geneExprNorm, "filename") <- NULL
        if (!is.null(sampleFilter) && sampleFilter != "") {
            sampleFilterText <- paste(sampleFilter, collapse=", ")
        } else {
            sampleFilterText <- "None"
        }
        sampleFilterSettings <- c("Discarded samples"=sampleFilterText)

        if (filter) {
            geneFilterSettings <- c(
                "Gene filtering"="Enabled",
                "Mean >="=minMean, # "Mean <="=maxMean,
                "Variance >="=minVar, # "Variance <="=maxVar,
                "Counts for at least some samples >="=minCounts,
                "Total counts across samples >="=minTotalCounts)
        } else {
            geneFilterSettings <- c("Gene filtering"="Disabled")
        }

        if (!voom) {
            avgCountPerObservationText <- priorCount
            names(avgCountPerObservationText) <- c(
                "Average count added per observation")
        } else {
            avgCountPerObservationText <- NULL
        }

        convertToGeneSymbolText <- if (convertToGeneSymbol) "Yes" else "no"
        names(convertToGeneSymbolText) <- paste(
            "Replace unambiguous ENSEMBL gene identifiers with their gene",
            "symbols")

        settings <- c(list(
            "Original gene expression (file)"=attr(geneExpr, "filename"),
            "Original gene expression (label)"=isolate(input$geneExpr)
        ), sampleFilterSettings, geneFilterSettings, list(
            "Normalisation method"=method,
            "Mean-variance modelling (voom)"=if (voom) "Yes" else "No",
            "Log2-transformation"="Yes"),
        avgCountPerObservationText,
        convertToGeneSymbolText)
        attr(geneExprNorm, "settings") <- settings
        attr(geneExprNorm, "icon") <- list(symbol="cogs", colour="green")
        attr(geneExprNorm, "description") <- "Gene expression (normalised)"
        attr(geneExprNorm, "dataType") <- "Gene expression"

        if (is(geneExprNorm, "EList"))
            geneExprNorm$E <- inheritAttrs(geneExprNorm$E, geneExprNorm)
        setNormalisedGeneExpression(geneExprNorm)
        endProcess("processGeneExpr", time=time)
    })
    loadGeneExpressionSet(session, input, output)
}

attr(geNormalisationFilteringUI, "loader") <- "data"
attr(geNormalisationFilteringServer, "loader") <- "data"