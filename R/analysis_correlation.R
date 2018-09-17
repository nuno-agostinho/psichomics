#' @rdname appUI
#' 
#' @importFrom shiny NS tagList sidebarPanel mainPanel sliderInput actionButton
#' uiOutput
#' @importFrom highcharter highchartOutput
correlationUI <- function(id) {
    ns <- NS(id)
    
    corrParams <- bsCollapsePanel(
        tagList(icon("filter"), "Correlation parameters"), 
        value="corrParams", style="info",
        selectizeInput(ns("geneExpr"), "Gene expression", choices=NULL,
                       width="100%"),
        selectizeGeneInput(ns("gene"), multiple=TRUE),
        # actionLink(ns("addRBPs"), 
        #            "Add RBPs from (Sebestyen et al., 2016)..."),
        selectizeInput(
            ns("ASevents"), "Alternative splicing events", choices=NULL,
            multiple=TRUE, width="100%",
            options=list(plugins=list('remove_button', 'drag_drop'))), 
        hr(),
        selectizeInput(
            ns("method"), "Correlation method", width="100%",
            c("Pearson's product-moment correlation"="pearson", 
              "Kendall's rank correlation tau"="kendall",
              "Spearman's rank correlation rho"="spearman"),
            "spearman"),
        selectizeInput(ns("alternative"), "Alternative hypothesis",
                       c("Two-sided"="two.sided",
                         "Greater than zero"="greater", 
                         "Less than zero"="less"), width="100%"),
        selectGroupsUI(ns("groupFilter"),
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
                       choices=c(1:4, 6, 12), selected=3), 
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
        tagList(icon("sliders"), "Scatterplot options"), 
        value="scatterplotOptions", style="info",
        selectGroupsUI(
            ns("groupColour"), label="Sample colouring",
            noGroupsLabel="Same colour for all samples",
            groupsLabel="Colour by selected groups", returnAllDataValue=TRUE,
            returnAllDataLabel="Display data outside selected groups"),
        colourInputMod(ns("colour"), "Base point colour", "black",
                       returnName=TRUE),
        sliderInput(ns("alpha"), "Point opacity", min=0, max=100, value=20, 
                    step=1, post="%", width="100%"),
        bsCollapse(generalPlotOptions, loessOptions, densityOptions),
        processButton(ns("applyPlotStyle"), label="Apply"))
    
    options <- div(id=ns("options"), bsCollapse(open=c("corrParams"), 
                                                corrParams, scatterParams))
    
    tagList(
        sidebar(
            errorDialog(paste("No alternative splicing quantification or gene",
                              "expression data are available."),
                        id=ns("noData"), buttonLabel="Load data",
                        buttonIcon="plus-circle", buttonId=ns("loadData")),
            hidden(options)),
        mainPanel(
            uiOutput(ns("correlations")),
            hidden(dataTableOutput(ns("corTable"))),
            hidden(downloadButton(ns("saveTable"), "Save table", "btn-info"))))
}

#' Subset gene expression based on (full or partial) matching genes
#' 
#' @param geneExpr Data frame or matrix: gene expression
#' @param gene Character: genes to look for
#' 
#' @return Gene expression subset for the input genes
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
    return(geneExpr[matched, ])
}

#' Find splicing events based on given genes
#'
#' @param psi Data frame or matrix: alternative splicing quantification
#' @param gene Character: gene
#'
#' @return Character vector containing alternative splicing events
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
#' @inheritDotParams stats:::cor.test.default alternative:continuity
#' 
#' @importFrom stats cor.test
#' 
#' @export
#' @return List of correlations where each element contains:
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
    
    geneExprSubset <- subsetGeneExpressionFromMatchingGenes(geneExpr, gene)
    
    if ( is.null(ASevents) || identical(ASevents, "") )
        ASevents <- findASeventsFromGene(psi, gene)
    
    # Calculate correlation betwenn GE and AS event(s)
    corrPerGene <- function(gene, ASevent, geneExpr, psi, ...) {
        expr           <- geneExpr[gene, ]
        exprNum        <- as.numeric(expr)
        names(exprNum) <- colnames(expr)
        
        eventPSI           <- psi[ASevent, ]
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
    
    res <- lapply(ASevents, function(ASevent) {
        gene <- rownames(geneExprSubset)
        corr <- lapply(gene, corrPerGene, ASevent, geneExprSubset, psi, ...)
        names(corr) <- gene
        return(corr)
    })
    names(res) <- ASevents
    class(res) <- c("GEandAScorrelation", class(res))
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
#' 
#' colourGroups <- list(Normal=paste("Normal", 1:3), 
#'                      Tumour=paste("Cancer", 1:3))
#' attr(colourGroups, "Colour") <- c(Normal="#00C65A", Tumour="#EEE273")
#' plotCorrelation(corr, colourGroups=colourGroups, alpha=1)
plotCorrelation <- function(corr, autoZoom=FALSE, loessSmooth=TRUE,
                            loessFamily=c("gaussian", "symmetric"),
                            colour="black", alpha=0.2, size=1.5,
                            loessColour="red", loessAlpha=1, loessWidth=0.5,
                            fontSize=12, ..., colourGroups=NULL, legend=FALSE,
                            showAllData=TRUE, density=FALSE, 
                            densityColour="blue", densityWidth=0.5) {
    loessFamily <- match.arg(loessFamily)
    
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
        pvalue   <- trimWhitespace(signifDigits(unname(single$cor$p.value)))
        
        plot <- ggplot(mapping=aes(x=event, y=expr))
        if (!is.null(colourGroups)) {
            sampleColour <- rep(names(colourGroups), 
                                sapply(colourGroups, length))
            groupNames   <- sampleColour[match(names(expr),
                                               unlist(colourGroups))]
            lookupColour <- attr(colourGroups, "Colour")
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
        
        plot <- plot +
            ggtitle(eventDisplay, 
                    sprintf("%s\n%s: %s (p-value: %s)", eventDetails, 
                            estimateMethod, estimate, pvalue)) +
            labs(x="PSI", y=paste(single$gene, "gene expression"))
        return(plot + theme_light(fontSize))
    }
    
    lapply(corr, lapply, plotCorrPerASevent)
}

#' @rdname plotCorrelation
plot.GEandAScorrelation <- plotCorrelation

print.GEandAScorrelation <- function(object) {
    for (item in object) {
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

#' @rdname appServer
#' 
#' @importFrom shiny renderUI observeEvent isolate tagList tags
#' @importFrom highcharter renderHighchart
#' @importFrom shinyjs show hide toggle
correlationServer <- function(input, output, session) {
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
            updateSelectizeInput(
                session, "ASevents", 
                choices=c("Type to search for a splicing event..."="",
                          rownames(psi)))
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
    
    # Plot correlation analyses
    plotShinyCorr <- reactive({
        ns <- session$ns
        corr <- getCorrelation()
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
        
        plots <- plotCorrelation(
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
                    distributeByCol(ns("plot"), length(plots), cols, height), 
                    hr())
            }
        })
        
        lapply(seq(plots), function(i)
            output[[paste0("plot", i)]] <- renderPlot(plots[[i]]))
    })
    
    displayCorrTable <- reactive({
        corr <- getCorrelation()
        if (is.null(corr)) return(NULL)
        
        # Prepare table with correlation analyses
        eventID  <- unlist(lapply(corr, lapply, "[[", "eventID"))
        gene     <- unlist(lapply(corr, lapply, "[[", "gene"))
        estimate <- unlist(lapply(corr, lapply, 
                                  function(i) i[["cor"]][["estimate"]][[1]]))
        pvalue   <- unlist(lapply(corr, lapply, 
                                  function(i) i[["cor"]][["p.value"]]))
        method   <- unlist(lapply(corr, lapply, 
                                  function(i) i[["cor"]][["method"]]))
        qvalue   <- p.adjust(pvalue)
        
        method <- unique(method)
        if (length(method) != 1) 
            stop("Only one correlation method is currently supported.")
        
        data           <- data.frame(gsub("_", " ", eventID, fixed=TRUE), 
                                     gene, estimate, pvalue, qvalue)
        colnames(data) <- c("Alternative splicing event", "Protein",
                            method, "p-value", "p-value (BH adjusted)")
        
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
    
    observeEvent(input$correlate, {
        ns <- session$ns
        
        isolate({
            geneExpr    <- getGeneExpression()[[input$geneExpr]]
            psi         <- getInclusionLevels()
            gene        <- input$gene
            ASevents    <- input$ASevents
            method      <- input$method
            alternative <- input$alternative
        })
        
        if (is.null(psi)) {
            missingDataModal(session, "Inclusion levels",
                             ns("missingInclusionLevels"))
            return(NULL)
        } else if (is.null(geneExpr)) {
            errorModal(session, "No gene expression selected", 
                       "Please selected gene expression data",
                       caller="Correlation analysis")
            return(NULL)
        } else if (is.null(gene) || identical(gene, "")) {
            errorModal(session, "No gene selected", "Please select a gene",
                       caller="Correlation analysis")
            return(NULL)
        } else if (is.null(ASevents) || identical(ASevents, "")) {
            errorModal(session, "No alternative splicing event selected",
                       "Please select one or more alternative splicing events",
                       caller="Correlation analysis")
            return(NULL)
        }
        
        # Filter samples based on groups
        groupFilter <- isolate(getSelectedGroups(
            input, "groupFilter", "Samples", 
            filter=intersect(colnames(geneExpr), colnames(psi))))
        groupFilter <- unname(unlist(groupFilter))
        if (is.null(groupFilter)) groupFilter <- TRUE
        geneExpr <- geneExpr[ , groupFilter, drop=FALSE]
        psi      <- psi[ , groupFilter, drop=FALSE]
        
        # Perform correlation analyses
        startProcess("correlate")
        corr <- suppressWarnings(
            correlateGEandAS(geneExpr, psi, gene, ASevents, method=method, 
                             alternative=alternative))
        setCorrelation(corr)
        plotShinyCorr()
        displayCorrTable()
        endProcess("correlate")
    })
    
    observeEvent(input$applyPlotStyle, {
        startProcess("applyPlotStyle")
        plotShinyCorr()
        endProcess("applyPlotStyle")
    })
    
    observeEvent(input$missingInclusionLevels, 
                 missingDataGuide("Inclusion levels"))
    observeEvent(input$loadData, missingDataGuide("Inclusion levels"))
}

attr(correlationUI, "loader") <- "analysis"
attr(correlationUI, "name") <- "Correlation analysis"
attr(correlationServer, "loader") <- "analysis"