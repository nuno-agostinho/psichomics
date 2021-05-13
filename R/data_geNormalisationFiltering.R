geneExprFilteringSetting <- function(ns, id, min=0, max=100, step=1, ...,
                                     label=id) {
    psiFilteringSetting(ns, id, min=min, max=max, step=step, ..., label=label)
}

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
        geneExprFilteringSetting(ns, "Mean", check=c(TRUE, FALSE)),
        geneExprFilteringSetting(ns, "Variance", check=c(TRUE, FALSE)),
        fluidRow(
            numericInputWithCheckbox(ns, "minCounts",
                                     div("Counts >=", icon("question-circle")),
                                     min=0, max=100, value=10, check=TRUE),
            numericInputWithCheckbox(ns, "minTotalCounts",
                                     div("Total counts >=",
                                         icon("question-circle")),
                                     min=0, max=100, value=15, check=TRUE),
            bsTooltip(ns("elementMinCounts"), placement="top",
                      options=list(container="body"),
                      paste("Minimum counts in a worthwhile number of samples:",
                            "for more information, check documentation for",
                            tags$code("edgeR::filterByExpr()"))),
            bsTooltip(ns("elementMinTotalCounts"), placement="top",
                      options=list(container="body"),
                      "Minimum total counts across all samples")))

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
                tagList(icon("vial"), "Sample filtering",
                        contextUI(ns("sampleFilterText"))),
                value="Sample filtering",
                selectizeInput(ns("sampleFilter"), "Samples to discard",
                               multiple=TRUE, width="100%",
                               choices=character(0))),
            bsCollapsePanel(
                tagList(icon("dna"), "Gene filtering",
                        contextUI(ns("filterText"))),
                value="Filtering", filters, filteringAssistant),
            bsCollapsePanel(
                tagList(icon("balance-scale"), "Normalisation",
                        contextUI(ns("normalisationText"))),
                value="Normalisation",
                helpText("Scale raw library sizes using",
                         tags$code("edgeR::calcNormFactors()"), ", unless",
                         tags$code("quantile"), "is selected."),
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
                         tags$code("limma::voom()"), "should be more powerful",
                         "and preferred.")),
            bsCollapsePanel(
                tagList(icon("retweet"), "CPM and log2",
                        contextUI(ns("logTransformText"))),
                value="Log-transformation",
                helpText("Compute log2-transformed counts per million",
                         "(log2CPM) using", tags$code("edgeR::cpm()"), "(or",
                         tags$code("limma::voom()"), ", if selected)."),
                numericInput(
                    ns("priorCount"), value=0.25, step=0.25, width="100%",
                    paste("Average count to add to each observation to avoid",
                          "zeroes after log-transformation"))),
            bsCollapsePanel(
                tagList(icon("address-book"), "Convert to gene symbol",
                        contextUI(ns("geneSymbolConversionText"))),
                value="Convert to gene symbol",
                checkboxInput(
                    ns("convertToGeneSymbol"), width="100%",
                    paste("Replace unambiguous ENSEMBL gene identifiers with",
                          "their gene symbols"), value=TRUE),
                selectizeInput(ns("orgDb"), "Gene database for conversion",
                               width="100%", choices=NULL,
                               options=list(highlight=FALSE)))))

    tagList(
        uiOutput(ns("modal")),
        errorDialog("Gene expression not loaded.",
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
#' @description
#' Gene expression is filtered and normalised in the following steps:
#'
#' \itemize{
#' \item{Filter gene expression;}
#' \item{Normalise gene expression with \code{\link[edgeR]{calcNormFactors}};}
#' \item{If \code{performVoom = FALSE}, compute counts per million (CPM) using
#' \code{\link[edgeR]{cpm}} and log2-transform values if
#' \code{log2transform = TRUE};}
#' \item{If \code{performVoom = TRUE}, use \code{\link[limma]{voom}} to compute
#' log2-CPM, quantile-normalise (if \code{method = "quantile"}) and estimate
#' mean-variance relationship to calculate observation-level weights.}
#' }
#'
#' @param geneExpr Matrix or data frame: gene expression
#' @param geneFilter Boolean: filtered genes (if \code{NULL}, skip filtering)
#' @param method Character: normalisation method, including \code{TMM},
#' \code{RLE}, \code{upperquartile}, \code{none} or \code{quantile} (see
#' Details)
#' @inheritParams edgeR::calcNormFactors
#' @param log2transform Boolean: perform log2-transformation?
#' @param priorCount Average count to add to each observation to avoid zeroes
#' after log-transformation
#' @param performVoom Boolean: perform mean-variance modelling
#' (using \code{\link[limma]{voom}})?
#'
#' @details \code{edgeR::calcNormFactors} will be used to normalise gene
#' expression if \code{method} is \code{TMM}, \code{RLE}, \code{upperquartile}
#' or \code{none}. If \code{performVoom = TRUE}, \code{\link[limma]{voom}} will
#' only normalise if \code{method = "quantile"}.
#'
#' Available normalisation methods:
#' \itemize{
#' \item{\code{TMM} is recommended for most RNA-seq data where more than half of
#' the genes are believed not differentially expressed between any pair of
#' samples;}
#' \item{\code{RLE} calculates the median library from the geometric mean of all
#' columns and the median ratio of each sample to the median library is taken as
#' the scale factor;}
#' \item{\code{upperquartile} calculates the scale factors from a given quantile
#' of the counts for each library, after removing genes with zero counts in all
#' libraries;}
#' \item{\code{quantile} forces the entire empirical distribution of each
#' column to be identical (only performed if \code{performVoom = TRUE}).}
#' }
#'
#' @importFrom edgeR DGEList [.DGEList calcNormFactors cpm
#' @importFrom limma voom
#'
#' @family functions for gene expression pre-processing
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
    geneFilterSettings <- attr(geneFilter, "settings")
    if (is.null(geneFilter)) geneFilter <- TRUE
    else if (!any(geneFilter)) return(NULL)

    originalGeneExpr <- geneExpr
    geneExpr <- DGEList(geneExpr)
    geneExprNorm <- geneExpr[geneFilter, , keep.lib.sizes=TRUE]

    updateProgress("Normalising gene expression")
    if (!performVoom && method == "quantile") method <- "none"
    if (method != "quantile") {
        geneExprNorm <- calcNormFactors(geneExprNorm, method=method, p=p)
    }

    if (!performVoom) {
        geneExprNorm <- cpm(geneExprNorm, log=log2transform,
                            prior.count=priorCount)

        avgCountPerObservationText <- priorCount
        names(avgCountPerObservationText) <- paste("Average count added per",
                                                   "observation")
    } else {
        norm <- if (method == "quantile") "quantile" else "none"
        log2transform <- TRUE
        geneExprNorm  <- voom(geneExprNorm, normalize.method=norm)

        avgCountPerObservationText <- NULL
    }

    updateProgress("Preparing gene expression data")
    if (!is(geneExprNorm, "EList")) geneExprNorm <- data.frame(geneExprNorm)
    colnames(geneExprNorm) <- colnames(geneExpr)

    geneExprNorm <- inheritAttrs(geneExprNorm, originalGeneExpr)
    if (is(geneExprNorm, "EList"))
        geneExprNorm$E <- inheritAttrs(geneExprNorm$E, originalGeneExpr)
    attr(geneExprNorm, "filename") <- NULL

    attr(geneExprNorm, "settings") <- c(
        "Original gene expression (file)"=attr(originalGeneExpr, "filename"),
        "Original gene expression (label)"=attr(originalGeneExpr, "label"),
        geneFilterSettings,
        "Normalisation method"=method,
        "Mean-variance modelling (voom)"=if (performVoom) "Yes" else "No",
        "Log2-transformation"=log2transform,
        avgCountPerObservationText)
    return(geneExprNorm)
}

#' @rdname normaliseGeneExpression
#' @export
normalizeGeneExpression <- normaliseGeneExpression

#' Set of functions to load splicing quantification
#'
#' @inherit geNormalisationFilteringServer
#'
#' @importFrom shiny tags
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

    # Load gene expression
    loadGeneExpression <- reactive({
        time <- startProcess("loadGeneExpr")

        updateProgress("Wait a moment", divisions=2)
        updateProgress("Loading gene expression")

        allFormats <- loadFileFormats()
        formats <- allFormats[sapply(allFormats, "[[",
                                     "dataType") == "Gene expression"]

        geneExpr <- tryCatch(
            loadFile(input$customGeneExpr, formats, multiple=TRUE),
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

            prepareGeneExpr <- function(geneExpr, set=FALSE) {
                tablename <- attr(geneExpr, "tablename")
                if (is.null(tablename)) tablename <- "Gene expression"
                if (set) {
                    name <- renameDuplicated(tablename,
                                             names(getCategoryData()))
                    setDataTable(name, geneExpr)
                } else {
                    res <- list(geneExpr)
                    names(res) <- tablename
                    return(res)
                }
            }

            if ( is.null(getData()) ) {
                name <- file_path_sans_ext( basename(input$customGeneExpr) )
                name <- gsub(" Gene expression.*$", "", name)
                if (name == "") name <- "Unnamed"

                if (is.data.frame(geneExpr)) {
                    geneExpr <- prepareGeneExpr(geneExpr)
                } else {
                    geneExpr <- unlist(lapply(geneExpr, prepareGeneExpr),
                                       recursive=FALSE)
                }
                data <- setNames(list(geneExpr), name)
                data <- processDatasetNames(data)
                setData(data)
                setCategory(name)

                samples <- colnames(geneExpr)
                parsed <- parseTCGAsampleInfo(samples)
                if ( !is.null(parsed) ) setSampleInfo(parsed)
            } else if (is.data.frame(geneExpr)) {
                prepareGeneExpr(geneExpr, set=TRUE)
            } else if (is.list(geneExpr)) {
                lapply(geneExpr, prepareGeneExpr, set=TRUE)
            }
            removeModal()
        }
        endProcess("loadGeneExpr", time)
    })
    observeEvent(input$loadCustomGE, loadGeneExpression())
}

#' Convert gene identifiers
#'
#' @param annotation \code{OrgDb} with genome wide annotation for an organism or
#'   \code{character} with species name to query \code{OrgDb}, e.g.
#'   \code{"Homo sapiens"}
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
#' @importFrom AnnotationHub AnnotationHub query
#'
#' @family functions for gene expression pre-processing
#' @return Character vector of the respective targets of gene identifiers. The
#' previous identifiers remain other identifiers have the same target (in case
#' \code{ignoreDuplicatedTargets = TRUE}) or if no target was found.
#' @export
#'
#' @examples
#' # Use species name to automatically look for a OrgDb database
#' sp <- "Homo sapiens"
#' genes <- c("ENSG00000012048", "ENSG00000083093", "ENSG00000141510",
#'            "ENSG00000051180")
#' convertGeneIdentifiers(sp, genes)
#' convertGeneIdentifiers(sp, genes, key="ENSEMBL", target="UNIPROT")
#'
#' # Alternatively, set the annotation database directly
#' ah <- AnnotationHub::AnnotationHub()
#' sp <- AnnotationHub::query(ah, c("OrgDb", "Homo sapiens"))[[1]]
#' columns(sp) # these attributes can be used to change the attributes
#'
#' convertGeneIdentifiers(sp, genes)
#' convertGeneIdentifiers(sp, genes, key="ENSEMBL", target="UNIPROT")
convertGeneIdentifiers <- function(annotation, genes, key="ENSEMBL",
                                   target="SYMBOL",
                                   ignoreDuplicatedTargets=TRUE) {
    if (is.character(annotation)) {
        ah <- AnnotationHub()
        annotation <- query(ah, c("OrgDb", annotation))[[1]]
        if (length(annotation) == 0) {
            stop(sprintf("No query found for species '%s'", annotation))
        }
    } else if (!is(annotation, "OrgDb")) {
        stop("Annotation needs to be a 'character' or 'OrgDb' object")
    }

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
#' Uses \code{\link[edgeR]{filterByExpr}} to determine genes with sufficiently
#' large counts to retain for statistical analysis.
#'
#' @param geneExpr Data frame or matrix: gene expression
#' @param minMean Numeric: minimum of read count mean per gene
#' @param maxMean Numeric: maximum of read count mean per gene
#' @param minVar Numeric: minimum of read count variance per gene
#' @param maxVar Numeric: maximum of read count variance per gene
#' @param minCounts Numeric: minimum number of read counts per gene for a
#' worthwhile number of samples (check \code{\link[edgeR]{filterByExpr}} for
#' more information)
#' @param minTotalCounts Numeric: minimum total number of read counts per gene
#'
#' @importFrom edgeR filterByExpr
#'
#' @family functions for gene expression pre-processing
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
    if (is.na(minMean))        minMean <- -Inf
    if (is.na(maxMean))        maxMean <- Inf
    if (is.na(minVar))         minVar <- -Inf
    if (is.na(maxVar))         maxVar <- Inf
    if (is.na(minCounts))      minCounts <- 0
    if (is.na(minTotalCounts)) minTotalCounts <- 0

    geneExprMean <- customRowMeans(geneExpr, fast=TRUE)
    geneExprVar  <- customRowVars(geneExpr, fast=TRUE)

    varMeanFilter <- geneExprMean >= minMean & geneExprMean <= maxMean &
        geneExprVar >= minVar & geneExprVar <= maxVar

    lowCountFilter <- suppressMessages(
        filterByExpr(geneExpr[varMeanFilter, ],
                     min.count=minCounts,
                     min.total.count=minTotalCounts))
    filteredGenes <- varMeanFilter
    filteredGenes[names(lowCountFilter[!lowCountFilter])] <- FALSE

    attr(filteredGenes, "settings") <- c(
        "Mean >="=minMean, "Mean <="=maxMean,
        "Variance >="=minVar, "Variance <="=maxVar,
        "Counts for at least some samples >="=minCounts,
        "Total counts across samples >="=minTotalCounts)
    return(filteredGenes)
}

#' Plot distribution of gene expression per sample
#'
#' @param geneExpr Data frame or matrix: gene expression
#' @inheritDotParams renderBoxplot
#'
#' @importFrom highcharter %>% hc_yAxis
#'
#' @family functions for gene expression pre-processing
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

#' Plot library size
#'
#' @param data Data frame or matrix: gene expression
#' @param log10 Boolean: log10-transform \code{data}?
#' @param title Character: plot title
#' @param subtitle Character: plot subtitle
#' @param colour Character: data colour
#'
#' @family functions for gene expression pre-processing
#' @return Library size distribution
#' @export
#'
#' @examples
#' df <- data.frame(geneA=c(2, 4, 5),
#'                  geneB=c(20, 3, 5),
#'                  geneC=c(5, 10, 21))
#' colnames(df) <- paste("Sample", 1:3)
#' plotLibrarySize(df)
plotLibrarySize <- function(
    data, log10=TRUE,
    title="Library size distribution across samples",
    subtitle="Library size: total number of mapped reads",
    colour="orange") {

    table <- colSums(data)
    if (log10) {
        table <- log10(table)
        xAxisLabel <- "log10(Library sizes)"
    } else {
        xAxisLabel <- "Library sizes"
    }
    yAxisLabel <- "Density"

    groups <- "All samples"
    attr(groups, "Colour") <- c("All samples"=colour)
    plot <- plotDistribution(table, groups,
                             rugLabels=TRUE, vLine=FALSE, legend=FALSE,
                             title=title, valueLabel="log10(library size)") %>%
        hc_xAxis(title=list(text=xAxisLabel)) %>%
        hc_yAxis(title=list(text=yAxisLabel)) %>%
        hc_subtitle(text=paste(subtitle))
    return(plot)
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
#' @importFrom highcharter hc_plotOptions hc_xAxis hc_chart
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

        checkSetting <- function(id) {
            if (startsWith(id, "min")) {
                res <- -Inf
            } else if (startsWith(id, "max")) {
                res <- Inf
            } else {
                res <- 0
            }

            enabled <- input[[paste0("enable", capitalize(id))]]
            if (isTRUE(enabled)) res <- input[[id]]
            return(res)
        }
        minMean        <- checkSetting("minMean")
        maxMean        <- checkSetting("maxMean")
        minVariance    <- checkSetting("minVariance")
        maxVariance    <- checkSetting("maxVariance")
        minCounts      <- checkSetting("minCounts")
        minTotalCounts <- checkSetting("minTotalCounts")

        sampleFilter   <- input$sampleFilter
        if (!is.null(sampleFilter) && sampleFilter != "") {
            samplesToKeep <- !colnames(geneExpr) %in% sampleFilter
            geneExpr <- geneExpr[ , samplesToKeep]
        }
        filtered <- filterGeneExpr(geneExpr, minMean, maxMean, minVariance,
                                   maxVariance, minCounts, minTotalCounts)
        return(filtered)
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
        geneExprMean <- customRowMeans(geneExpr, fast=TRUE)
        maxMean      <- ceiling(max(geneExprMean, na.rm=TRUE))
        updateNumericInput(session, "minMean", max=maxMean)
        updateNumericInput(session, "maxMean", max=maxMean, value=maxMean)

        # Update variance range
        geneExprVar <- customRowVars(geneExpr, fast=TRUE)
        maxVar      <- ceiling(max(geneExprVar, na.rm=TRUE))
        updateNumericInput(session, "minVar", max=maxVar)
        updateNumericInput(session, "maxVar", max=maxVar, value=maxVar)

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

    # Filter and normalise gene expression
    observeEvent(input$processGeneExpr, {
        time <- startProcess("processGeneExpr")

        isolate({
            geneExpr            <- getGeneExpression(input$geneExpr)
            method              <- input$normalisation
            percentile          <- input$upperquartilePercentile
            sampleFilter        <- input$sampleFilter
            priorCount          <- input$priorCount
            voom                <- input$voom
            convertToGeneSymbol <- input$convertToGeneSymbol
            orgDb               <- input$orgDb
        })

        # Filter samples
        if (!is.null(sampleFilter) && sampleFilter != "") {
            samplesToKeep <- !colnames(geneExpr) %in% sampleFilter
            geneExpr <- geneExpr[ , samplesToKeep]
            sampleFilterText <- paste(sampleFilter, collapse=", ")
        } else {
            sampleFilterText <- "None"
        }
        sampleFilterSettings <- c("Discarded samples"=sampleFilterText)

        # Filter and normalise
        geneFilter <- getFilter()
        attr(geneExpr, "label") <- isolate(input$geneExpr)
        geneExprNorm <- normaliseGeneExpression(
            geneExpr, geneFilter, method, percentile, log2transform=TRUE,
            priorCount, performVoom=voom)

        # Convert ENSEMBL gene id to gene symbols
        convertToGeneSymbolText <- "No"
        if (convertToGeneSymbol) {
            sp <- query(AnnotationHub(), c("OrgDb", orgDb))
            if (length(sp) >= 1) {
                db <- sp[1]
                genes <- rownames(geneExprNorm)
                rownames(geneExprNorm) <- convertGeneIdentifiers(db[[1]], genes)
                convertToGeneSymbolText <- sprintf("Yes, using %s (%s)",
                                                   db$title, db$species)
            } else {
                warning("Database for species not found")
            }
        }
        names(convertToGeneSymbolText) <- paste(
            "Replace unambiguous ENSEMBL gene identifiers with their gene",
            "symbols")

        # Prepare attributes
        description <- attr(geneExprNorm, "description")
        if (!is.null(description)) {
            description <- gsub(" \\(normalised\\)$", "", description)
            description <- paste(description, "(normalised)")
        } else {
            description <- "Gene expression (normalised)"
        }
        geneExprNorm <- addObjectAttrs(
            geneExprNorm,
            "settings"=c(attr(geneExprNorm, "settings"),
                         sampleFilterSettings, convertToGeneSymbolText),
            "icon"=list(symbol="cogs", colour="green"),
            "description"=description,
            "dataType"="Gene expression")

        if (is(geneExprNorm, "EList")) {
            geneExprNorm$E <- inheritAttrs(geneExprNorm$E, geneExprNorm)
        }
        setNormalisedGeneExpression(geneExprNorm)
        endProcess("processGeneExpr", time=time)
    })
    loadGeneExpressionSet(session, input, output)

    # Toggle filtering options
    toggleGEsetting <- function(id) {
        checkbox <- paste0("enable", capitalize(id))
        observe(toggleState(id, input[[checkbox]]))
    }
    toggleGEsetting("minMean")
    toggleGEsetting("maxMean")
    toggleGEsetting("minVariance")
    toggleGEsetting("maxVariance")
    toggleGEsetting("minCounts")
    toggleGEsetting("minTotalCounts")

    # Update context
    output$sampleFilterText <- renderText({
        sampleFilter <- input$sampleFilter
        if (is.null(sampleFilter) || sampleFilter == "") {
            text <- "No samples to discard"
        } else {
            len  <- length(sampleFilter)
            text <- sprintf("%s sample%s to discard",
                            len, ifelse(len == 1, "", "s"))
        }
        return(text)
    })

    output$filterText <- renderText({
        geneExpr <- input$geneExpr
        if (is.null(geneExpr) || geneExpr == "") return(NULL)
        geneExpr <- isolate(getGeneExpression(geneExpr))

        filter <- sum(getFilter())
        total  <- nrow(geneExpr)
        ratio  <- filter/total * 100

        msg <- sprintf("%s genes of %s (%s%%)", filter, total, round(ratio))
        return(msg)
    })

    output$normalisationText <- renderText({
        text <- input$normalisation
        if (text == "upperquartile") {
            text <- sprintf("%s (%s)", text, input$upperquartilePercentile)
        }
        if (input$voom) text <- paste(text, "+ voom")
        return(text)
    })

    observe({
        if (input$normalisation == "quantile") {
            updateCheckboxInput(session, "voom", value=TRUE)
            disable("voom")
        } else {
            enable("voom")
        }
    })

    output$logTransformText <- renderText({
        paste("Prior count:", input$priorCount)
    })

    output$geneSymbolConversionText <- renderText({
        ifelse(input$convertToGeneSymbol, "Enabled", "Disabled")
    })

    observe(toggleState("orgDb", input$convertToGeneSymbol))

    # Update OrgDb choices for gene symbol conversion
    observe({
        # Suppress snapshot date message
        orgDb <- suppressMessages(query(AnnotationHub(), "OrgDb")$species)
        orgDb <- unique(orgDb)
        if (is.null(orgDb)) orgDb <- "Homo sapiens"
        updateSelectizeInput(session, "orgDb", choices=orgDb,
                             selected="Homo sapiens")
    })
}

attr(geNormalisationFilteringUI, "loader") <- "data"
attr(geNormalisationFilteringServer, "loader") <- "data"
