#' @rdname appUI
#' 
#' @importFrom shiny NS tagList sidebarPanel mainPanel sliderInput actionButton
#' uiOutput
#' @importFrom highcharter highchartOutput
correlationUI <- function(id) {
    ns <- NS(id)
    
    options <- div(
        id=ns("options"),
        bsCollapse(
            open=c("corrParams", "corrOptions", "scatterplotOptions"),
            multiple=TRUE,
            bsCollapsePanel(
                tagList(icon("filter"), "Correlation parameters"), 
                value="corrParams", style="info",
                selectizeInput(ns("geneExpr"), "Gene expression", choices=NULL),
                selectizeGeneInput(ns("gene"), multiple=TRUE),
                # actionLink(ns("addRBPs"), 
                #            "Add RBPs from (Sebestyen et al., 2016)..."),
                selectizeInput(
                    ns("ASevents"), "Alternative splicing events", choices=NULL,
                    multiple=TRUE, 
                    options=list(plugins=list('remove_button', 'drag_drop'))), 
                selectGroupsUI(ns("groups"),
                               label="Perform correlation analysis on...",
                               noGroupsLabel="All samples",
                               groupsLabel="Samples from selected groups")),
            bsCollapsePanel(
                tagList(icon("wrench"), "Correlation options"), 
                value="corrOptions", style="info",
                selectizeInput(
                    ns("method"), "Correlation method", 
                    c("Pearson's product-moment correlation"="pearson", 
                      "Kendall's rank correlation tau"="kendall",
                      "Spearman's rank correlation rho"="spearman"),
                    "spearman"),
                selectizeInput(ns("alternative"), "Alternative hypothesis",
                               c("Two-sided"="two.sided",
                                 "Greater than zero"="greater", 
                                 "Less than zero"="less"))),
            bsCollapsePanel(
                tagList(icon("sliders"), "Scatterplot options"), 
                value="scatterplotOptions", style="info",
                radioButtons(
                    ns("zoom"), "Range of axis for PSI values",
                    list("Fixed between 0 and 1"="fixed",
                         "Automatic zoom based on available data"="auto")),
                selectizeInput(ns("cols"), "Number of plots per row",
                               choices=c(1:4, 6, 12), selected=3), hr(),
                sliderInput(ns("size"), "Size of points", 0, 4, 1.5, 0.5),
                colourInput(ns("colour"), "Colour for points", "#00000044",
                            allowTransparent=TRUE), hr(),
                checkboxInput(ns("loessSmooth"), "Plot Loess curve", 
                              value=TRUE),
                sliderInput(ns("loessWidth"), "Width of Loess curve", 
                            0, 2, 0.5, 0.25),
                colourInput(ns("loessColour"), "Colour for Loess curve", "red",
                            allowTransparent=TRUE))),
        processButton(ns("correlate"), label="Correlate"))
    
    tagList(
        sidebarPanel(
            errorDialog(paste("No alternative splicing quantification or gene",
                              "expression data are available."),
                        id=ns("noData"), buttonLabel="Load data",
                        buttonIcon="plus-circle", buttonId=ns("loadData")),
            hidden(options)),
        mainPanel(uiOutput(ns("correlations"))))
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
#' @inheritDotParams stats:::cor.test.default alternative:continuity
#' 
#' @importFrom stats cor.test
#' 
#' @export
#' @return List of correlations where each element is organised as such:
#' \item{\code{eventID}}{Alternative splicing event identifier}
#' \item{\code{cor}}{Correlation between gene expression and alternative 
#' splicing quantification of one alternative splicing event}
#' \item{\code{geneExpr}}{Gene expression for the selected gene}
#' \item{\code{psi}}{Alternative splicing quantification for the alternative 
#' splicing event}
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
    
    # Check if genes in gene expression are styled like in TCGA: "Gene|Probe"
    isTcgaStyle <- all(grepl("|", head(rownames(geneExpr)), fixed=TRUE))
    
    # Parse gene based on TCGA or non-TCGA style
    tcgaStyledGenes <- grepl("|", gene, fixed=TRUE)
    query <- character(length(tcgaStyledGenes))
    query[tcgaStyledGenes] <- gsub("|", "\\|", gene[tcgaStyledGenes],
                                   fixed=TRUE)
    pattern <- ifelse(isTcgaStyle, "^%s\\|", "^%s$")
    query[!tcgaStyledGenes] <- sprintf(pattern, gene[!tcgaStyledGenes])
    
    # Retrieve gene based on first match of gene expression data
    matched <- lapply(query, grep, rownames(geneExpr), value=TRUE)
    matched <- Filter(length, matched)
    if (length(matched) == 0) stop("Gene expression not found for given genes.")
    matched <- sapply(matched, "[[", 1)
    expr <- geneExpr[matched, ]
    
    if ( is.null(ASevents) || identical(ASevents, "") ) {
        # If no AS events are discriminated, find AS events for the given genes
        ASevents <- rownames(psi)
        query <- gene
        if (isTcgaStyle) {
            query <- strsplit(gene, "|", fixed=TRUE)
            query <- sapply(query, "[[", 1)
        }
        query <- sprintf("_%s|%s$|/%s/", query, query, query)
        ASevents <- unlist(lapply(query, grep, ASevents, value=TRUE))
    }
    
    if (length(ASevents) == 0)
        stop("No alternative splicing events found based on the given gene.")
    
    # Calculate correlation betwenn GE and AS event(s)
    corrPerGene <- function(gene, ASevent, geneExpr, psi, ...) {
        expr <- as.numeric(geneExpr[gene, ])
        eventPSI <- as.numeric(psi[ASevent, ])
        cor  <- tryCatch(cor.test(expr, eventPSI, ...), error=return)
        return(list(eventID=ASevent, gene=gene, cor=cor, geneExpr=expr, 
                    psi=eventPSI))
    }
    
    res <- lapply(ASevents, function(ASevent) {
        corr <- lapply(gene, corrPerGene, ASevent, geneExpr, psi, ...)
        names(corr) <- gene
        return(corr)
    })
    names(res) <- ASevents
    return(res)
}

#' Plot correlations
#'
#' Plot correlation results from \code{\link{correlateGEandAS}}
#'
#' @param corr List of correlations
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
#'
#' @importFrom ggplot2 ggplot geom_point geom_line labs coord_cartesian ggtitle
#'   aes theme_light
#' @importFrom stats loess.smooth
#'
#' @export
#' @return Renders plots for each correlation in \code{corr}
#'
#' @examples
#' annot <- readFile("ex_splicing_annotation.RDS")
#' junctionQuant <- readFile("ex_junctionQuant.RDS")
#' psi <- quantifySplicing(annot, junctionQuant, eventType=c("SE", "MXE"))
#'
#' geneExpr <- readFile("ex_gene_expression.RDS")
#' corr <- correlateGEandAS(geneExpr, psi, "ALDOA")
#' plotCorrelation(corr)
plotCorrelation <- function(corr, autoZoom=FALSE, loessSmooth=TRUE, 
                            loessFamily=c("gaussian", "symmetric"), 
                            colour="black", alpha=0.2, size=1.5,
                            loessColour="red", loessAlpha=1, loessWidth=0.5, 
                            ...) {
    plotCorrPerASevent <- function(single) {
        expr    <- single$geneExpr
        event   <- single$psi
        eventId <- parseSplicingEvent(single$eventID, char=TRUE, pretty=TRUE)
        eventDisplay <- gsub(" \\(.*\\)", "", eventId)
        eventDetails <- gsub(".*\\((.*)\\)", "\\1", eventId)
        
        estimateMethod <- names(single$cor$estimate)
        if (!is.null(estimateMethod) && estimateMethod == "cor") 
            estimateMethod <- "r"
        
        estimate <- trimWhitespace(signifDigits(unname(single$cor$estimate)))
        Pvalue   <- trimWhitespace(signifDigits(unname(single$cor$p.value)))
        
        plot <- ggplot(mapping=aes(x=event, y=expr)) + 
            geom_point(na.rm=TRUE, colour=colour, alpha=alpha, size=size) +
            ggtitle(eventDisplay, sprintf(
                "%s\n%s: %s (p-value: %s)", eventDetails, estimateMethod, 
                estimate, Pvalue)) +
            labs(x="PSI levels", y=paste(single$gene, "gene expression"))
        
        if (!autoZoom) plot <- plot + coord_cartesian(xlim = c(0, 1))
        if (loessSmooth) {
            loess <- tryCatch(suppressWarnings( 
                loess.smooth(event, expr, family="gaussian", ...)), 
                error=return)
            
            if (!is(loess, "error"))
                plot <- plot + geom_line(
                    aes(x=loess$x, y=loess$y),
                    colour=loessColour, alpha=loessAlpha, size=loessWidth)
        }
        return(plot + theme_light(12))
    }
    
    lapply(corr, lapply, plotCorrPerASevent)
}

#' @rdname appServer
#' 
#' @importFrom shiny renderUI observeEvent isolate tagList tags
#' @importFrom highcharter renderHighchart
#' @importFrom shinyjs show hide
correlationServer <- function(input, output, session) {
    selectGroupsServer(session, "groups", "Samples")
    
    observe({
        if (is.null( getInclusionLevels() ) || is.null( getGeneExpression() )) {
            show("noData")
            hide("options")
        } else {
            hide("noData")
            show("options")
        }
    })
    
    # Update available gene choices depending on gene expression data loaded
    # Reactive avoids updating if the input remains the same
    updateGeneChoices <- reactive({
        geneExpr <- getGeneExpression()[[input$geneExpr]]
        genes <- rownames(geneExpr)
        updateSelectizeInput(session, "gene", choices=genes, server=TRUE)
    })
    
    # Update gene expression data
    observe({
        geneExpr <- getGeneExpression()
        if ( !is.null(geneExpr) ) {
            updateSelectizeInput(session, "geneExpr",
                                 choices=rev(names(geneExpr)))
        }
    })
    
    # Update gene choices
    observe({
        geneExpr <- getGeneExpression()
        if ( !is.null(geneExpr) ) {
            updateGeneChoices()
            show("gene")
        } else {
            hide("gene")
        }
    })
    
    # Update alternative splicing events
    observe({
        psi  <- getInclusionLevels()
        if (!is.null(psi)) {
            updateSelectizeInput(session, "ASevents", choices=rownames(psi))
        }
    })
    
    # Update selected alternative splicing events based on selected gene
    observeEvent(input$gene, {
        geneExpr <- getGeneExpression()[[input$geneExpr]]
        gene <- input$gene
        ASevents <- isolate(input$ASevents)
        psi  <- getInclusionLevels()
        
        if (is.null(psi)) return(NULL)
        
        allEvents <- rownames(psi)
        names(allEvents) <- parseSplicingEvent(allEvents, char=TRUE)
        
        if (!is.null(ASevents)) return(NULL)
        
        # Automatically change AS events to those of selected genes if no AS
        # event has yet been selected
        if (!is.null(geneExpr) && !is.null(gene) && !identical(gene, "")) {
            
            isTcgaStyle <- all(grepl("|", head(rownames(geneExpr)), fixed=TRUE))
            
            query <- gene
            if (isTcgaStyle) query <- strsplit(gene, "|", fixed=TRUE)[[1]][[1]]
            
            if ( !identical(query, "?") ) {
                query <- sprintf("_%s|%s$|/%s/", query, query, query)
                ASevents <- grep(query, allEvents, value=TRUE)
            } else {
                ASevents <- character(0)
            }
            
            if (length(ASevents) == 0) {
                choices  <- c("No events found for the selected gene"="")
                selected <- NULL
            } else {
                choices  <- c("Select an alternative splicing event"="")
                selected <- ASevents
            }
            
            choices <- c(choices, allEvents)
            updateSelectizeInput(session, "ASevents", choices=choices, 
                                 selected=selected, server=TRUE)
        } else {
            choices <- c("Select an alternative splicing event"="", allEvents)
            updateSelectizeInput(session, "ASevents", choices=choices, 
                                 server=TRUE)
        }
    })
    
    observeEvent(input$correlate, {
        ns <- session$ns
        
        isolate({
            geneExpr    <- getGeneExpression()[[input$geneExpr]]
            psi         <- getInclusionLevels()
            gene        <- input$gene
            ASevents    <- input$ASevents
            method      <- input$method
            alternative <- input$alternative
            loessSmooth <- input$loessSmooth
            
            autoZoom    <- input$zoom
            autoZoom    <- identical(autoZoom, "auto")
            
            colour      <- input$colour
            size        <- input$size
            loessColour <- input$loessColour
            loessWidth  <- input$loessWidth
        })
        
        if (is.null(psi)) {
            missingDataModal(session, "Inclusion levels",
                             ns("missingInclusionLevels"))
            return(NULL)
        } else if (is.null(geneExpr)) {
            errorModal(session, "No gene expression selected", 
                       "Please selected gene expression data")
            return(NULL)
        } else if (is.null(gene) || identical(gene, "")) {
            errorModal(session, "No gene selected", "Please select a gene")
            return(NULL)
        } else if (is.null(ASevents) || identical(ASevents, "")) {
            errorModal(session, "No alternative splicing event selected",
                       "Please select one or more alternative splicing events")
            return(NULL)
        }
        
        # Filter samples based on groups
        startProcess("correlate")
        groups <- isolate(getSelectedGroups(
            input, "groups", "Samples", 
            filter=intersect(colnames(geneExpr), colnames(psi))))
        groups <- unname(unlist(groups))
        if (is.null(groups)) groups <- TRUE
        geneExpr <- geneExpr[ , groups, drop=FALSE]
        psi      <- psi[ , groups, drop=FALSE]
        
        corr <- suppressWarnings(
            correlateGEandAS(geneExpr, psi, gene, ASevents, method=method, 
                             alternative=alternative))
        plots <- plotCorrelation(corr, autoZoom, loessSmooth, colour=colour,
                                 size=size, loessColour=loessColour, 
                                 loessWidth=loessWidth)
        plots <- unlist(plots, recursive=FALSE)
        
        # Plot all groups
        output$correlations <- renderUI({
            distributeByCol <- function(id, len, cols) {
                ncols <- len
                nrows <- ceiling(len/cols)
                eachRow <- list()
                # Create rows
                for (i in seq(nrows)) {
                    eachCol <- list()
                    # Create columns
                    for (k in seq( min(cols, ncols) )) {
                        content <- plotOutput(paste0(id, k + cols * (i - 1)),
                                              height="200px")
                        eachCol <- c(eachCol, list(column(12 / cols, content)))
                    }
                    eachRow <- c(eachRow, list(do.call(fluidRow, eachCol)))
                    ncols   <- ncols - cols
                }
                do.call(tagList, eachRow)
            }
            
            cols <- as.numeric( isolate(input$cols) )
            if ( length(cols) > 0 )
                distributeByCol(ns("plot"), length(plots), cols)
        })
        
        lapply(seq(plots), function(i)
            output[[paste0("plot", i)]] <- renderPlot(plots[[i]]))
        
        endProcess("correlate")
    })
    
    observeEvent(input$missingInclusionLevels, 
                 missingDataGuide("Inclusion levels"))
    observeEvent(input$loadData, missingDataGuide("Inclusion levels"))
}

attr(correlationUI, "loader") <- "analysis"
attr(correlationUI, "name") <- c(
    "Correlation of gene expression and alternative splicing")
attr(correlationServer, "loader") <- "analysis"