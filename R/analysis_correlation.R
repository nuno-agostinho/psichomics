#' @rdname appUI
#'
#' @importFrom shiny NS tagList sidebarPanel mainPanel sliderInput actionButton
#' uiOutput
#' @importFrom highcharter highchartOutput
correlationUI <- function(id) {
    ns <- NS(id)

    corrParams <- bsCollapsePanel(
        tagList(icon("balance-scale-right"), "Correlation parameters"),
        value="corrParams", style="info",
        helpText(paste("Correlate gene expression against alternative",
                       "splicing quantification")),
        selectizeInput(ns("geneExpr"), "Gene expression", choices=NULL,
                       width="100%"),
        selectGroupsUI(ns("genes"), "Genes",
                       label="Genes from selected groups"),
        radioButtons(
            ns("ASeventsSelection"), "Alternative splicing event(s)",
            c("Selected event"="selectedEvent",
              "From selected groups"="byGroup")),
        conditionalPanel(
            sprintf("input[id='%s'] == '%s'",
                    ns("ASeventsSelection"), "byGroup"),
            selectGroupsUI(ns("ASevents"), "ASevents", label=NULL)),
        hr(), selectizeInput(
            ns("method"), "Correlation method", width="100%",
            c("Pearson's product-moment correlation"="pearson",
              "Kendall's rank correlation tau"="kendall",
              "Spearman's rank correlation rho"="spearman"),
            "spearman"),
        selectizeInput(ns("alternative"), "Alternative hypothesis",
                       c("Two-sided"="two.sided",
                         "Greater than zero"="greater",
                         "Less than zero"="less"), width="100%"),
        selectGroupsUI(ns("groupFilter"), "Samples",
                       label="Perform correlation analysis on...",
                       noGroupsLabel="All samples",
                       groupsLabel="Samples from selected groups"),
        processButton(ns("correlate"), label="Correlate"))

    generalPlotOptions <- bsCollapsePanel(
        "Plot options",
        radioButtons(
            ns("zoom"), "Range of axis for PSI values", width="100%",
            list("Fixed between 0 and 1"="fixed",
                 "Automatic based on available values"="auto")),
        numericInput(ns("height"), "Height of each plot (pixels)",
                     200, min=50, max=1000, step=50, width="100%"),
        selectizeInput(ns("cols"), "Plots per row", width="100%",
                       choices=c(1, 2, 3, 4, 6, 12), selected=3),
        numericInput(ns("fontSize"), "Font size", 12, min=1, max=50,
                     step=1, width="100%"), hr(),
        sliderInput(ns("size"), "Size of points", 0, 4, 1.5, 0.5,
                    width="100%"))

    loessOptions <- bsCollapsePanel(
        "Loess curve",
        checkboxInput(ns("loessSmooth"), "Plot Loess curve",
                      value=TRUE, width="100%"),
        div(id=ns("loessOptions"),
            selectizeInput(ns("loessFamily"), width="100%",
                           "Loess curve fitting family",
                           c(Gaussian="gaussian", Symmetric="symmetric")),
            sliderInput(ns("loessWidth"), "Width of Loess curve",
                        0, 2, 0.5, 0.25, width="100%"),
            colourInputMod(ns("loessColour"), "Colour for Loess curve",
                           "red", allowTransparent=TRUE, returnName=TRUE)))

    densityOptions <- bsCollapsePanel(
        "Density contours",
        checkboxInput(ns("plotDensity"), value=FALSE, width="100%",
                      "Plot contours based on a density estimate"),
        div(id=ns("densityOptions"),
            sliderInput(ns("densityWidth"), width="100%",
                        "Width of density contours", 0, 2, 0.5, 0.25),
            colourInputMod(
                ns("densityColour"), "Colour for density contours",
                "blue", allowTransparent=TRUE, returnName=TRUE)))

    scatterParams <- bsCollapsePanel(
        tagList(icon("sliders-h"), "Scatterplot options"),
        value="scatterplotOptions", style="info",
        radioButtons(
            ns("genesToPlot"), "Genes to plot", width="100%",
            c("All genes"="all", "Selected genes"="selected")),
        conditionalPanel(
            sprintf("input[id='%s'] == '%s'", ns("genesToPlot"), "selected"),
            selectizeInput(ns("selectedGenesToPlot"), label="Genes to plot",
                           choices=NULL, width="100%", multiple=TRUE,
                           options=list(placeholder="Select genes to plot"))),
        radioButtons(
            ns("ASeventsToPlot"), "Alternative splicing events to plot",
            c("All events"="all", "Selected events"="selected")),
        conditionalPanel(
            sprintf("input[id='%s'] == '%s'", ns("ASeventsToPlot"), "selected"),
            selectizeInput(
                ns("selectedASeventsToPlot"),
                label="Splicing events to plot",
                choices=NULL, width="100%", multiple=TRUE,
                options=list(placeholder="Select splicing events to plot"))),
        selectGroupsUI(
            ns("groupColour"), label="Sample colouring", "Samples",
            noGroupsLabel="Same colour for all samples",
            groupsLabel="Colour by selected groups", returnAllDataValue=TRUE,
            returnAllDataLabel="Display data outside selected groups"),
        colourInputMod(ns("colour"), "Base point colour", "black",
                       returnName=TRUE),
        sliderInput(ns("alpha"), "Point opacity", min=0, max=100, value=20,
                    step=1, post="%", width="100%"),
        bsCollapse(generalPlotOptions, loessOptions, densityOptions),
        processButton(ns("applyPlotStyle"), label="Plot"))

    options <- div(id=ns("options"), bsCollapse(open=c("corrParams"),
                                                corrParams, scatterParams))

    tagList(
        uiOutput(ns("modal")),
        sidebar(
            errorDialog(paste("No alternative splicing quantification or gene",
                              "expression data are available."),
                        id=ns("noData"), buttonLabel="Load data",
                        buttonIcon="plus-circle", buttonId=ns("loadData")),
            hidden(options)),
        mainPanel(
            hidden(dataTableOutput(ns("corTable"))),
            hidden(downloadButton(ns("saveTable"), "Save table", "btn-info")),
            uiOutput(ns("correlations"))))
}

#' Subset gene expression based on (full or partial) matching genes
#'
#' @param geneExpr Data frame or matrix: gene expression
#' @param gene Character: genes to look for
#'
#' @return Gene expression subset for the input genes
#' @keywords internal
subsetGeneExpressionFromMatchingGenes <- function(geneExpr, gene) {
    # Start by matching input genes with genes in gene expression
    exactMatch <- match(gene, rownames(geneExpr))
    unmatched  <- is.na(exactMatch)
    if (any(unmatched)) {
        # Parse input genes based on TCGA or non-TCGA style
        query <- gene[unmatched]
        tcgaStyledGenes <- grepl("|", query, fixed=TRUE)
        query[tcgaStyledGenes] <- gsub("\\|.*$", "", query[tcgaStyledGenes])

        # Parse genes in gene expression if TCGA-styled, i.e. "Gene|Probe"
        geneExprGenes <- rownames(geneExpr)
        isTcgaStyle <- all(grepl("|", head(rownames(geneExpr)), fixed=TRUE))
        if (isTcgaStyle) {
            geneExprGenes <- strsplit(geneExprGenes, "\\|")
            geneExprGenes <- sapply(geneExprGenes, "[[", 1)
        }

        # Retrieve gene based on first match of gene expression data
        bestMatch <- match(query, geneExprGenes)
    } else {
        bestMatch <- NULL
    }

    matched <- exactMatch
    matched[unmatched] <- bestMatch
    matched <- matched[!is.na(matched)]
    if (length(matched) == 0) stop("Gene expression not found for input genes.")
    if (is(geneExpr, "EList")) geneExpr <- geneExpr$E
    return(geneExpr[matched, , drop=FALSE])
}

#' Find splicing events based on given genes
#'
#' @param psi Data frame or matrix: alternative splicing quantification
#' @param gene Character: gene
#'
#' @return Character vector containing alternative splicing events
#' @keywords internal
findASeventsFromGene <- function(psi, gene) {
    # If no AS events are discriminated, find AS events for the given genes
    ASevents <- rownames(psi)
    query <- gene
    isTcgaStyle <- all(grepl("|", head(query), fixed=TRUE))
    if (isTcgaStyle) {
        query <- strsplit(query, "|", fixed=TRUE)
        query <- sapply(query, "[[", 1)
    }
    query <- sprintf("_%s|%s$|/%s/", query, query, query)
    ASevents <- unlist(lapply(query, grep, ASevents, value=TRUE))

    if (length(ASevents) == 0)
        stop("No alternative splicing events found based on the given gene.")
    else
        return(ASevents)
}

#' Correlate gene expression data against alternative splicing quantification
#'
#' Test for association between paired samples' gene expression (for any genes
#' of interest) and alternative splicing quantification.
#'
#' @param geneExpr Matrix or data frame: gene expression data
#' @param psi Matrix or data frame: alternative splicing quantification data
#' @param gene Character: gene symbol for genes of interest
#' @param ASevents Character: alternative splicing events to correlate with
#' gene expression of a gene (if \code{NULL}, the events will be automatically
#' retrieved from the given gene)
#' @param ... Extra parameters passed to \link[stats]{cor.test}
#'
#' @importFrom stats cor.test
#'
#' @family functions to correlate gene expression and alternative splicing
#' @return List of correlations where each element contains:
#' \item{\code{eventID}}{Alternative splicing event identifier}
#' \item{\code{cor}}{Correlation between gene expression and alternative
#' splicing quantification of one alternative splicing event}
#' \item{\code{geneExpr}}{Gene expression for the selected gene}
#' \item{\code{psi}}{Alternative splicing quantification for the alternative
#' splicing event}
#' @export
#'
#' @examples
#' annot <- readFile("ex_splicing_annotation.RDS")
#' junctionQuant <- readFile("ex_junctionQuant.RDS")
#' psi <- quantifySplicing(annot, junctionQuant, eventType=c("SE", "MXE"))
#'
#' geneExpr <- readFile("ex_gene_expression.RDS")
#' correlateGEandAS(geneExpr, psi, "ALDOA")
correlateGEandAS <- function(geneExpr, psi, gene, ASevents=NULL, ...) {
    # Filter and order samples available in both gene expression and alternative
    # splicing quantification data
    samples  <- intersect(colnames(geneExpr), colnames(psi))
    geneExpr <- geneExpr[ , samples]
    psi      <- psi[ , samples]

    geneExprSubset <- subsetGeneExpressionFromMatchingGenes(geneExpr, gene)

    if ( is.null(ASevents) || identical(ASevents, "") )
        ASevents <- findASeventsFromGene(psi, gene)

    # Calculate correlation betwenn GE and AS event(s)
    corrPerGene <- function(gene, ASevent, geneExpr, psi, ...) {
        updateProgress("Performing correlation analysis", console=FALSE)

        expr           <- geneExpr[gene, , drop=FALSE]
        exprNum        <- as.numeric(expr)
        names(exprNum) <- colnames(expr)

        eventPSI           <- psi[ASevent, , drop=FALSE]
        eventPSInum        <- as.numeric(eventPSI)
        names(eventPSInum) <- colnames(eventPSI)

        # Only compare samples available in both groups
        common <- intersect(names(exprNum), names(eventPSInum))
        if (length(common) == 0) {
            exprNum     <- NULL
            eventPSInum <- NULL
            cor         <- NULL

            msg <- paste("No common samples available to correlate %s",
                         "expression with %s splicing event.")
            warning(sprintf(msg, gene, ASevent))
        } else {
            exprNum     <- exprNum[common]
            eventPSInum <- eventPSInum[common]
            cor <- tryCatch(cor.test(exprNum, eventPSInum, ...), error=return)
        }
        res <- list(eventID=ASevent, gene=gene, cor=cor, geneExpr=exprNum,
                    psi=eventPSInum)
        return(res)
    }

    updateProgress("Performing correlation analyses",
                   divisions=length(gene) * length(ASevents))

    res <- lapply(ASevents, function(ASevent) {
        gene <- rownames(geneExprSubset)
        corr <- lapply(gene, corrPerGene, ASevent, geneExprSubset, psi, ...)
        names(corr) <- gene
        return(corr)
    })
    names(res) <- ASevents
    res        <- preserveAttributes(res)
    class(res) <- c("GEandAScorrelation", class(res))
    attr(res, "eventData") <- getSplicingEventData(psi)
    return(res)
}

#' @rdname plot.GEandAScorrelation
#' @param genes Character: genes
#' @param ASevents Character: AS events
#'
#' @importFrom stats na.omit
#'
#' @family functions to correlate gene expression and alternative splicing
#' @export
`[.GEandAScorrelation` <- function(x, genes=NULL, ASevents=NULL) {
    x <- unclass(x)

    if (!is.null(ASevents)) {
        if (is.numeric(ASevents)) {
            x <- x[ASevents]
        } else if (any(ASevents %in% names(x))) {
            ASevents <- ASevents[ASevents %in% names(x)]
            x <- x[ASevents]
        } else {
            x <- NULL
        }
    }
    # TODO: What if AS events do not match anything here?
    if (is.null(x)) return(NULL)

    if (!is.null(genes)) {
        for (eachASevent in seq(x)) {
            if (is.numeric(genes)) {
                x[[eachASevent]] <- x[[eachASevent]][genes]
            } else {
                ns <- names(x[[eachASevent]])
                matched <- c(match(genes, ns),
                             match(genes, gsub("\\|.*", "", ns)))
                matched <- na.omit(matched)
                geneNames <- ns[matched]
                if (any(geneNames %in% ns)) {
                    x[[eachASevent]] <- x[[eachASevent]][geneNames]
                } else {
                    x[[eachASevent]] <- NULL
                }
            }
        }
    }

    x <- lapply(x, function(item) Filter(Negate(is.null), item))
    x <- Filter(length, x)
    if (length(x) > 0 && !is.null(x)) {
        class(x) <- c("GEandAScorrelation", class(x))
    } else {
        x <- NULL
    }
    return(x)
}

#' Display results of correlation analyses
#'
#' Plot, print and display as table the results of gene expression and
#' alternative splicing
#'
#' @param x \code{GEandAScorrelation} object obtained after running
#'   \code{\link{correlateGEandAS}()}
#' @param loessSmooth Boolean: plot a smooth curve computed by
#'   \code{stats::loess.smooth}?
#' @param autoZoom Boolean: automatically set the range of PSI values based on
#'   available data? If \code{FALSE}, the axis relative to PSI values will range
#'   from 0 to 1
#' @param loessFamily Character: if \code{gaussian}, \code{loess} fitting is by
#'   least-squares, and if \code{symmetric}, a re-descending M estimator is used
#' @inheritDotParams stats::loess.smooth -x -y -family
#' @param colour Character: points' colour
#' @param alpha Numeric: points' alpha
#' @param size Numeric: points' size
#' @param loessColour Character: loess line's colour
#' @param loessAlpha Numeric: loess line's opacity
#' @param loessWidth Numeric: loess line's width
#' @param fontSize Numeric: plot font size
#' @param colourGroups List of characters: sample colouring by group
#' @param legend Boolean: show legend for sample colouring?
#' @param showAllData Boolean: show data outside selected groups as a single
#' group (coloured based on the \code{colour} argument)
#' @param density Boolean: contour plot of a density estimate
#' @param densityColour Character: line colour of contours
#' @param densityWidth Numeric: line width of contours
#'
#' @importFrom ggplot2 ggplot geom_point geom_line labs coord_cartesian ggtitle
#' aes theme_light scale_colour_manual geom_density_2d
#' @importFrom stats loess.smooth
#'
#' @method plot GEandAScorrelation
#' @family functions to correlate gene expression and alternative splicing
#' @return Plots, summary tables or results of correlation analyses
#' @export
#'
#' @aliases plotCorrelation
#'
#' @examples
#' annot <- readFile("ex_splicing_annotation.RDS")
#' junctionQuant <- readFile("ex_junctionQuant.RDS")
#' psi <- quantifySplicing(annot, junctionQuant, eventType=c("SE", "MXE"))
#'
#' geneExpr <- readFile("ex_gene_expression.RDS")
#' corr <- correlateGEandAS(geneExpr, psi, "ALDOA")
#'
#' # Quick display of the correlation results per splicing event and gene
#' print(corr)
#'
#' # Table summarising the correlation analysis results
#' as.table(corr)
#'
#' # Correlation analysis plots
#' colourGroups <- list(Normal=paste("Normal", 1:3),
#'                      Tumour=paste("Cancer", 1:3))
#' attr(colourGroups, "Colour") <- c(Normal="#00C65A", Tumour="#EEE273")
#' plot(corr, colourGroups=colourGroups, alpha=1)
plot.GEandAScorrelation <- function(
    x, autoZoom=FALSE, loessSmooth=TRUE, loessFamily=c("gaussian", "symmetric"),
    colour="black", alpha=0.2, size=1.5, loessColour="red", loessAlpha=1,
    loessWidth=0.5, fontSize=12, ..., colourGroups=NULL, legend=FALSE,
    showAllData=TRUE, density=FALSE, densityColour="blue", densityWidth=0.5) {
    loessFamily <- match.arg(loessFamily)

    plotCorrPerASevent <- function(single, eventData) {
        expr    <- single$geneExpr
        event   <- single$psi

        estimateMethod <- names(single$cor$estimate)
        if (!is.null(estimateMethod) && estimateMethod == "cor")
            estimateMethod <- "r"

        estimate <- trimWhitespace(signifDigits(unname(single$cor$estimate)))
        pvalue   <- trimWhitespace(signifDigits(unname(single$cor$p.value)))

        plot <- ggplot(mapping=aes(x=event, y=expr))
        if (!is.null(colourGroups)) {
            sampleColour <- rep(names(colourGroups),
                                sapply(colourGroups, length))
            groupNames   <- sampleColour[match(names(expr),
                                               unlist(colourGroups))]
            lookupColour <- attr(colourGroups, "Colour")

            if (!is.null(lookupColour)) {
                sampleColour <- lookupColour[groupNames]

                if (showAllData) {
                    names(sampleColour)[is.na(names(sampleColour))] <- "NA"
                    legendColours <- c(lookupColour, "NA"=colour)
                } else {
                    nonNAs        <- !is.na(sampleColour)
                    event         <- event[nonNAs]
                    expr          <- expr[nonNAs]
                    sampleColour  <- sampleColour[nonNAs]
                    legendColours <- lookupColour
                }

                plot <- plot +
                    geom_point(aes(colour=names(sampleColour)), na.rm=TRUE,
                               alpha=alpha, size=size) +
                    scale_colour_manual(name="", values=legendColours,
                                        guide=if(legend) "legend" else FALSE)
            } else {
                plot <- plot +
                    geom_point(aes(colour=sampleColour), na.rm=TRUE,
                               alpha=alpha, size=size)
            }
        } else {
            plot <- plot +
                geom_point(na.rm=TRUE, colour=colour, alpha=alpha, size=size)
        }

        if (!autoZoom)
            plot <- plot + coord_cartesian(xlim = c(0, 1))

        if (loessSmooth) {
            loess <- tryCatch(suppressWarnings(
                loess.smooth(event, expr, family=loessFamily, ...)),
                error=return)

            if (!is(loess, "error"))
                plot <- plot + geom_line(
                    aes(x=loess$x, y=loess$y),
                    colour=loessColour, alpha=loessAlpha, size=loessWidth)
        }

        if (density) {
            plot <- plot +
                geom_density_2d(show.legend=legend, colour=densityColour,
                                size=densityWidth, na.rm=TRUE)
        }

        eventId <- parseSplicingEvent(single$eventID, char=TRUE, pretty=TRUE,
                                      data=eventData)
        # Event ID example:
        # UHRF2 skipped exon (SE) (chr9: 6486925-6493826, positive strand)
        parenthesisRegex <- " \\((chr.*)\\)$"
        if (grepl(parenthesisRegex, eventId)) {
            eventDisplay <- gsub(parenthesisRegex, "", eventId)
            eventDetails <- paste0(gsub(
                paste0(".*", parenthesisRegex), "\\1", eventId), "\n")
        } else {
            eventDisplay <- eventId
            eventDetails <- ""
        }

        plot <- plot +
            ggtitle(eventDisplay,
                    sprintf("%s%s: %s (p-value: %s)", eventDetails,
                            estimateMethod, estimate, pvalue)) +
            labs(x="PSI", y=paste(single$gene, "gene expression"))
        return(plot + theme_light(fontSize))
    }

    lapply(x, lapply, plotCorrPerASevent, attr(x, "eventData"))
}

#' @export
plotCorrelation <- plot.GEandAScorrelation

#' @rdname plot.GEandAScorrelation
#' @export
print.GEandAScorrelation <- function(x, ...) {
    for (item in x) {
        for (elem in item) {
            consoleWidth <- options("width")
            cat(paste(rep("=", consoleWidth), collapse=""), fill=TRUE)
            cat(sprintf("%s splicing event\n%s gene expression",
                        elem$eventID, elem$gene), fill=TRUE)
            if (!is.null(elem$cor))
                print(elem$cor)
            else
                cat("\nNo correlation performed\n", fill=TRUE)
        }
    }
}

#' @rdname plot.GEandAScorrelation
#' @param pvalueAdjust Character: method used to adjust p-values (see Details)
#'
#' @details
#' The following methods for p-value adjustment are supported by using the
#' respective string in the \code{pvalueAdjust} argument:
#' \itemize{
#'      \item{\code{none}: do not adjust p-values}
#'      \item{\code{BH}: Benjamini-Hochberg's method (false discovery rate)}
#'      \item{\code{BY}: Benjamini-Yekutieli's method (false discovery rate)}
#'      \item{\code{bonferroni}: Bonferroni correction (family-wise error rate)}
#'      \item{\code{holm}: Holm's method (family-wise error rate)}
#'      \item{\code{hochberg}: Hochberg's method (family-wise error rate)}
#'      \item{\code{hommel}: Hommel's method (family-wise error rate)}
#' }
#' @export
as.table.GEandAScorrelation <- function (x, pvalueAdjust="BH", ...) {
    prepareCol <- function(object, FUN) unlist(lapply(object, lapply, FUN))

    gene     <- prepareCol(x, function(i) i[["gene"]])
    gene     <- prepareCol(x, function(i) i[["gene"]])
    eventID  <- prepareCol(x, function(i) i[["eventID"]])
    eventID  <- gsub("_", " ", eventID, fixed=TRUE) # Use prettifyEventID() ?

    estimate <- prepareCol(x, function(i) i[["cor"]][["estimate"]][[1]])
    pvalue   <- prepareCol(x, function(i) i[["cor"]][["p.value"]])
    method   <- prepareCol(x, function(i) i[["cor"]][["method"]])
    qvalue   <- p.adjust(pvalue, method=pvalueAdjust)
    qvalueLabel <- sprintf("p-value (%s adjusted)", pvalueAdjust)

    data <- data.frame(eventID, gene, method, estimate, pvalue, qvalue)
    if (length(unique(method)) > 1) {
        statCols <- c("Method", "Statistical value")
    } else {
        data$method <- NULL
        statCols    <- unique(method)
    }
    colnames(data) <- c("Alternative splicing event", "Gene", statCols,
                        "p-value", qvalueLabel)
    rownames(data) <- NULL
    return(data)
}

#' @rdname appServer
#'
#' @importFrom shiny renderUI observeEvent isolate tagList tags
#' @importFrom highcharter renderHighchart
#' @importFrom shinyjs show hide toggle
#' @importFrom graphics plot
correlationServer <- function(input, output, session) {
    selectGroupsServer(session, "ASevents", "ASevents")
    selectGroupsServer(session, "genes", "Genes")
    selectGroupsServer(session, "groupFilter", "Samples")
    selectGroupsServer(session, "groupColour", "Samples")

    observe({
        if (is.null( getInclusionLevels() ) || is.null( getGeneExpression() )) {
            show("noData")
            hide("options")
        } else {
            hide("noData")
            show("options")
        }
    })

    # Disable options for disabled features
    observe(toggle("loessOptions", condition=input$loessSmooth, anim=TRUE))
    observe(toggle("densityOptions", condition=input$plotDensity, anim=TRUE))

    observe({
        if (input$groupColourShowAllData ||
            input$groupColourSelection == "noGroups") {
            show("colour", anim=TRUE)
        } else {
            hide("colour", anim=TRUE)
        }
    })

    # Update gene expression data
    observe({
        geneExpr <- getGeneExpression()
        if ( !is.null(geneExpr) ) {
            updateSelectizeInput(session, "geneExpr",
                                 choices=rev(names(geneExpr)))
        }
    })

    displayCorrTable <- reactive({
        corr <- getCorrelation()
        if (is.null(corr)) return(NULL)
        data <- as.table(corr)

        show("corTable")
        output$corTable <- renderDataTable(
            data, style="bootstrap", server=TRUE, rownames=FALSE,
            selection="none", options=list(scrollX=TRUE))

        show("saveTable")
        output$saveTable <- downloadHandler(
            filename=function() paste(getCategory(), "Correlations"),
            content=function(con)
                write.table(data, con, quote=FALSE, sep="\t", row.names=FALSE))
    })

    performCorrelationAnalyses <- reactive({
        psi <- getInclusionLevels()

        if (input$ASeventsSelection == "selectedEvent") {
            ASevents <- getASevent()
        } else if (input$ASeventsSelection == "byGroup") {
            ASevents <- getSelectedGroups(input, "ASevents", "ASevents",
                                            filter=rownames(psi))
            ASevents <- unlist(ASevents)
        }
        geneExpr    <- getGeneExpression(input$geneExpr)
        gene        <- getSelectedGroups(input, "genes", "Genes",
                                            filter=rownames(geneExpr))
        gene        <- unlist(gene)

        method      <- input$method
        alternative <- input$alternative

        # Filter samples based on groups
        groupFilter <- getSelectedGroups(
            input, "groupFilter", "Samples",
            filter=intersect(colnames(geneExpr), colnames(psi)))
        groupFilter <- unname(unlist(groupFilter))
        if (is.null(groupFilter)) groupFilter <- TRUE
        geneExpr <- geneExpr[ , groupFilter, drop=FALSE]
        psi      <- psi[ , groupFilter, drop=FALSE]

        # Perform correlation analyses
        startProcess("correlate")
        corr <- suppressWarnings(
            tryCatch(correlateGEandAS(geneExpr, psi, gene, ASevents,
                                      method=method, alternative=alternative),
                     error=return))
        if ( !is(corr, "error") ) {
            setCorrelation(corr)
            displayCorrTable()
        }
        endProcess("correlate")
    })

    observeEvent(input$correlate, {
        ns <- session$ns
        isolate({
            psi         <- getInclusionLevels()

            if (input$ASeventsSelection == "selectedEvent") {
                ASevents <- getASevent()
            } else if (input$ASeventsSelection == "byGroup") {
                ASevents <- getSelectedGroups(input, "ASevents", "ASevents",
                                              filter=rownames(psi))
                ASevents <- unlist(ASevents)
            }

            geneExpr    <- getGeneExpression(input$geneExpr)
            gene        <- getSelectedGroups(input, "genes", "Genes",
                                             filter=rownames(geneExpr))
            gene        <- unlist(gene)

            cor <- getCorrelation()
        })

        if (is.null(psi)) {
            missingDataModal(session, "Inclusion levels",
                             ns("missingInclusionLevels"))
        } else if (is.null(geneExpr)) {
            errorModal(session, "No gene expression selected",
                       "Please selected gene expression data",
                       caller="Correlation analysis")
        } else if (is.null(gene) || identical(gene, "")) {
            errorModal(session, "No gene selected", "Please select genes",
                       caller="Correlation analysis")
        } else if (is.null(ASevents) || identical(ASevents, "")) {
            errorModal(session, "No alternative splicing event selected",
                       "Please select alternative splicing events",
                       caller="Correlation analysis")
        } else if (!is.null(cor)) {
            warningModal(session, "Correlation analyses already performed",
                         "Do you wish to discard the current results?",
                         footer=actionButton(ns("replace"), "Discard",
                                             class="btn-warning",
                                             "data-dismiss"="modal"),
                         caller="Correlation analyses")
        } else {
            geneSubset <- tryCatch(
                subsetGeneExpressionFromMatchingGenes(geneExpr, gene),
                error=return)
            if (is(geneSubset, "error")) {
                errorModal(session, "Selected genes not available",
                           "Gene expression dataset contains none of the",
                           "selected genes",
                           caller="Correlation analyses")
            } else {
                performCorrelationAnalyses()
            }
        }
    })

    # Replace previously performed differential analyses
    observeEvent(input$replace, performCorrelationAnalyses())

    # Plot correlation analyses
    plotShinyCorr <- reactive({
        ns <- session$ns
        corr <- getCorrelation()
        if (is.null(corr) || !is(corr, "GEandAScorrelation")) return(NULL)

        genes <- input$selectedGenesToPlot
        if (input$genesToPlot == "selected" && !is.null(genes))
            corr <- corr[genes=genes]

        ASevents <- input$selectedASeventsToPlot
        if (input$ASeventsToPlot == "selected" && !is.null(ASevents))
            corr <- corr[ASevents=ASevents]
        if (is.null(corr)) return(NULL)

        autoZoom    <- input$zoom
        autoZoom    <- identical(autoZoom, "auto")

        colour      <- input$colour
        alpha       <- input$alpha/100
        size        <- input$size
        fontSize    <- input$fontSize

        loessSmooth <- input$loessSmooth
        loessColour <- input$loessColour
        loessWidth  <- input$loessWidth
        loessFamily <- input$loessFamily

        plotDensity   <- input$plotDensity
        densityColour <- input$densityColour
        densityWidth  <- input$densityWidth
        showAllData   <- input$groupColourShowAllData

        # Colour samples based on groups
        groupColour <- getSelectedGroups(
            input, "groupColour", "Samples", filter=names(corr[[1]][[1]]$psi))

        plots <- plot(
            corr, colour=colour, alpha=alpha, size=size, fontSize=fontSize,
            autoZoom=autoZoom, loessFamily=loessFamily, loessSmooth=loessSmooth,
            loessColour=loessColour, loessWidth=loessWidth,
            colourGroups=groupColour, density=plotDensity,
            densityColour=densityColour, densityWidth=densityWidth,
            showAllData=showAllData)
        plots <- unlist(plots, recursive=FALSE)

        # Plot all groups
        output$correlations <- renderUI({
            distributeByCol <- function(id, len, cols, height) {
                ncols <- len
                nrows <- ceiling(len/cols)
                eachRow <- list()
                # Create rows
                for (i in seq(nrows)) {
                    eachCol <- list()
                    # Create columns
                    for (k in seq( min(cols, ncols) )) {
                        content <- plotOutput(paste0(id, k + cols * (i - 1)),
                                              height=height)
                        eachCol <- c(eachCol, list(column(12 / cols, content)))
                    }
                    eachRow <- c(eachRow, list(do.call(fluidRow, eachCol)))
                    ncols   <- ncols - cols
                }
                do.call(tagList, eachRow)
            }

            cols   <- as.numeric(input$cols)
            height <- input$height
            if ( length(cols) > 0 && height > 0 ) {
                height <- paste0(height, "px")
                tagList(
                    hr(),
                    distributeByCol(ns("plot"), length(plots), cols, height))
            }
        })

        lapply(seq(plots), function(i)
            output[[paste0("plot", i)]] <- renderPlot(plots[[i]]))
    })

    observeEvent(input$applyPlotStyle, {
        startProcess("applyPlotStyle")
        plotShinyCorr()
        endProcess("applyPlotStyle")
    })

    # Update choices to select genes and alternative splicing events to plot
    observe({
        corr <- getCorrelation()
        if (!is.null(corr)) {
            geneChoices    <- unique(unlist(sapply(corr, names)))
            ASeventChoices <- names(corr)
            names(ASeventChoices) <- parseSplicingEvent(names(corr), char=TRUE,
                                                        data=corr)
        } else {
            geneChoices    <- character(0)
            ASeventChoices <- character(0)
        }
        updateSelectizeInput(session, "selectedGenesToPlot",
                             choices=geneChoices)
        updateSelectizeInput(session, "selectedASeventsToPlot",
                             choices=ASeventChoices)
    })

    observeEvent(input$missingInclusionLevels,
                 missingDataGuide("Inclusion levels"))
    observeEvent(input$loadData, missingDataGuide("Inclusion levels"))
}

attr(correlationUI, "loader") <- "analysis"
attr(correlationUI, "name") <- "Correlation analysis"
attr(correlationServer, "loader") <- "analysis"
