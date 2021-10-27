## TODO(NunoA): should default columns be a perfect match or just a partial
## match? A partial match would be better for certain situations

## TODO(NunoA): render UI for each data table instead of rendering UI for all
## so there's no refresh

contextUI <- function(id) {
    tags$span(class="pull-right", tags$small(textOutput(id, inline=TRUE)))
}

#' Set attributes to an object
#'
#' @param object Object
#' @param ... Named parameters to convert to attributes
#' @param replace Boolean: replace an attribute if already set?
#'
#' @return Object with attributes set
#' @keywords internal
#'
#' @examples
#' ll <- list(a="hey", b="there")
#' psichomics:::addObjectAttrs(ll, "words"=2, "language"="English")
addObjectAttrs <- function (object, ..., replace=TRUE) {
    args <- list(...)
    if (length(args) == 1 && is.list(args[[1]])) args <- args[[1]]
    if (length(args) > 0) {
        for (k in seq(args)) {
            attrName <- names(args[k])
            # Attribute is not set if it is not NULL and replace=TRUE
            if (is.null(attr(object, attrName)) || replace) {
                attr(object, attrName) <- args[[k]]
            }
        }
    }
    return(object)
}

calculateAxisStats <- function(data, x, y=NULL,
                               stats=c("range", "var", "median", "mean"),
                               cache=NULL, verbose=FALSE) {
    names(stats) <- stats
    input <- lapply(stats, grepl, c(x, y))

    x <- y <- NULL
    vars <- list()
    for (stat in stats) {
        if (any(input[[stat]])) {
            # Check if summary statistic was previously cached
            if (!is.null(cache[[stat]])) {
                vars[[stat]] <- cache[[stat]]
                if (verbose) {
                    message(sprintf("Loaded %s per row from cache", stat))
                }
            } else {
                if (verbose) message(sprintf("Calculating %s per row...", stat))
                FUN <- switch(stat,
                              "var"=customRowVars,
                              "mean"=customRowMeans,
                              "median"=customRowMedians,
                              "range"=customRowRanges)
                vars[[stat]]  <- FUN(data, na.rm=TRUE, fast=TRUE)
                cache[[stat]] <- vars[[stat]]
            }
        }
    }
    vars <- data.frame(vars, stringsAsFactors=FALSE)
    return(list(vars=vars, cache=cache))
}

#' Plot row-wise statistics
#'
#' Scatter plot to compare between the row-wise mean, median, variance or range
#' from a data frame or matrix. Also supports transformations of those
#' variables, such as \code{log10(mean)}. If \code{y = NULL}, a density plot is
#' rendered instead.
#'
#' @param data Data frame or matrix containing samples per column and, for
#' instance, gene or alternative splicing event per row
#' @param x,y Character: statistic to calculate and display in the plot per row;
#' choose between \code{mean}, \code{median}, \code{var} or \code{range}
#' (or transformations of those variables, e.g. \code{log10(var)}); if
#' \code{y = NULL}, the density of \code{x} will be plot instead
#' @param subset Boolean or integer: \code{data} points to highlight
#' @param xmin,xmax,ymin,ymax Numeric: minimum and maximum X and Y values to
#' draw in the plot
#' @param xlim,ylim Numeric: X and Y axis range
#' @param cache List of summary statistics for \code{data} previously calculated
#'   to avoid repeating calculations (output also returns cache in attribute
#'   named \code{cache} with appropriate data)
#' @param verbose Boolean: print messages of the steps performed
#' @param data2 Same as \code{data} argument but points in \code{data2} are
#'   highlighted (unless \code{data2 = NULL})
#' @param legend Boolean: show legend?
#' @param legendLabels Character: legend labels
#'
#' @importFrom ggplot2 geom_vline geom_hline xlim ylim ggtitle geom_density
#'   scale_fill_manual scale_colour_manual
#'
#' @family functions for gene expression pre-processing
#' @family functions for PSI quantification
#' @return Plot of \code{data}
#' @export
#'
#' @examples
#' library(ggplot2)
#'
#' # Plotting gene expression data
#' geneExpr <- readFile("ex_gene_expression.RDS")
#' plotRowStats(geneExpr, "mean", "var^(1/4)") +
#'     ggtitle("Mean-variance plot") +
#'     labs(y="Square Root of the Standard Deviation")
#'
#' # Plotting alternative splicing quantification
#' annot <- readFile("ex_splicing_annotation.RDS")
#' junctionQuant <- readFile("ex_junctionQuant.RDS")
#' psi <- quantifySplicing(annot, junctionQuant, eventType=c("SE", "MXE"))
#'
#' medianVar <- plotRowStats(psi, x="median", y="var", xlim=c(0, 1)) +
#'     labs(x="Median PSI", y="PSI variance")
#' medianVar
#'
#' rangeVar  <- plotRowStats(psi, x="range", y="log10(var)", xlim=c(0, 1)) +
#'     labs(x="PSI range", y="log10(PSI variance)")
#' rangeVar
plotRowStats <- function(data, x, y=NULL, subset=NULL, xmin=NULL, xmax=NULL,
                         ymin=NULL, ymax=NULL, xlim=NULL, ylim=NULL,
                         cache=NULL, verbose=FALSE, data2=NULL, legend=FALSE,
                         legendLabels=c("Original", "Highlighted")) {
    stats <- c("range", "var", "median", "mean")
    isValidX <- any(sapply(stats, grepl, x))
    isValidY <- !is.null(y) && any(sapply(stats, grepl, y))
    if (!isValidX && !isValidY) {
        stop("Arguments 'x' and 'y' must contain one of the strings: ",
             paste(stats, collapse=", "), " (alternatively, y may be NULL)")
    }

    subsetCol <- "orange"
    remainCol <- ifelse(!is.null(subset) || !is.null(data2),
                        "darkgrey", "orange")
    res   <- calculateAxisStats(data, x, y, stats, cache=cache, verbose=verbose)
    cache <- res$cache

    vars <- cbind(res$vars, colour=remainCol)
    if (!is.null(subset)) {
        vars <- rbind(vars, cbind(vars[subset], colour=subsetCol))
    }
    if (!is.null(data2)) {
        vars2 <- calculateAxisStats(data2, x, y, stats, verbose=verbose)$vars
        vars  <- rbind(vars, cbind(vars2, colour=subsetCol))
    }

    if (verbose) message("Preparing plot...")
    if (isValidY) {
        plot <- ggplot(vars, aes_string(x, y, colour="colour")) +
            geom_point(size=1, na.rm=TRUE, alpha=0.5, show.legend=legend) +
            labs(x=x, y=y)
    } else {
        plot <- ggplot(vars, aes_string(x, colour="colour", fill="colour")) +
            geom_density(na.rm=TRUE, adjust=0.5, alpha=0.1,
                         show.legend=legend) +
            labs(x=x)
    }

    values <- c("black"="black", "darkgrey"="darkgrey", "orange"="orange")
    legend.position <- ifelse(legend, "bottom", "none")
    plot <- plot +
        scale_fill_manual(name="", labels=legendLabels, values=values) +
        scale_colour_manual(name="", labels=legendLabels, values=values) +
        theme(legend.position=legend.position)

    if (!is.null(xlim)) plot <- plot + xlim(xlim)
    if (!is.null(ylim)) plot <- plot + ylim(ylim)

    # Intercept lines
    if (!is.null(xmin)) plot <- plot + geom_vline(xintercept=xmin, colour="red")
    if (!is.null(xmax)) plot <- plot + geom_vline(xintercept=xmax, colour="red")
    if (!is.null(ymin)) plot <- plot + geom_hline(yintercept=ymin, colour="red")
    if (!is.null(ymax)) plot <- plot + geom_hline(yintercept=ymax, colour="red")

    attr(plot, "cache") <- cache
    return(plot)
}

plotRowStatsShiny <- function(x, y, data, data2) {
    stats  <- getPSIsummaryStats()
    xLabel <- names(stats[stats == x])

    if (y == "none") {
        y <- NULL
        yLabel <- "Density"
    } else {
        yLabel <- names(stats[stats == y])
    }

    legendLabels <- c("Original", "Filtered in")
    cache <- isolate(getInclusionLevelsSummaryStatsCache())
    res <- plotRowStats(data, x, y, data2=data2, cache=cache,
                        legend=TRUE, legendLabels=legendLabels,
                        verbose=TRUE) +
        theme_light(14) +
        theme(legend.position="bottom") +
        labs(x=xLabel, y=yLabel)

    cache <- attr(res, "cache")
    if (!is.null(cache)) setInclusionLevelsSummaryStatsCache(cache)
    return(res)
}

#' Warn user about loaded data
#'
#' @param modalId Character: identifier of the modal
#' @param replaceButtonId Character: identifier of the button to replace data
#' @param keepButtonId Character: identifier of the button to append data
#' @param session Shiny session
#'
#' @return HTML elements for a warning modal reminding data is loaded
#' @keywords internal
loadedDataModal <- function(session, modalId, replaceButtonId, keepButtonId) {
    ns <- session$ns
    warningModal(session, "Data already loaded",
                 "Would you like to", tags$b("replace"), "the loaded data or",
                 tags$b("keep"), "both the previous and new data?",
                 footer = tagList(
                     actionButton(ns(keepButtonId), "data-dismiss"="modal",
                                  label="Keep both"),
                     actionButton(ns(replaceButtonId), class="btn-warning",
                                  "data-dismiss"="modal", label="Replace")),
                 modalId=modalId)
}

#' Process dataset names
#'
#' @details Avoid duplicated names and append the technology used for junction
#' quantification
#'
#' @param data List of lists of data frames
#'
#' @return Processed list of lists of data frames
#' @keywords internal
processDatasetNames <- function(data) {
    newData <- data
    # Avoid duplicate names in categories
    names(newData) <- renameDuplicated(names(data),
                                       names(data)[duplicated(names(data))])

    ns <- lapply(newData, names)
    for (each in names(ns)) {
        nse <- names(newData[[each]])

        # For read quantification, add the respective sequencing technology
        index <- nse %in% c("Junction quantification", "Gene expression")
        for (k in seq_along(nse)) {
            if (index[[k]]) {
                file <- attr(newData[[each]][[k]], "filename")
                if (is.null(file)) next

                if (grepl("illuminahiseq", file, fixed=TRUE))
                    names(newData[[each]])[[k]] <- paste(
                        names(newData[[each]])[[k]], "(Illumina HiSeq)")
                else if (grepl("illuminaga", file, fixed=TRUE))
                    names(newData[[each]])[[k]] <- paste(
                        names(newData[[each]])[[k]], "(Illumina GA)")
            }
        }

        # Avoid duplicate names in datasets from the same category
        nse <- names(newData[[each]])
        names(newData[[each]]) <- renameDuplicated(nse, nse[duplicated(nse)])
    }
    return(newData)
}

# Molecular data file input ----------------------------------------------------

#' @rdname fileBrowserInput
#' @importFrom shinyBS bsPopover
#' @importFrom shiny fileInput
fileBrowserInfoInput <- function(id, label, infoContent=NULL, clearable=FALSE) {
    if (!getOption("psichomics.shinyproxy", FALSE)) {
        input <- fileBrowserInput(
            id, label, placeholder="No file selected", clearable=clearable,
            info=TRUE, infoFUN=bsPopover, infoTitle=label,
            infoContent=infoContent)
    } else {
        input <- fileBrowserShinyproxyInput(
            id, label, clearable=FALSE, info=TRUE, infoFUN=bsPopover,
            infoTitle=label, infoContent=infoContent)
    }
    return(input)
}

#' File input for molecular data
#'
#' @param id Character: identifier for gene expression input
#' @inheritParams fileBrowserInput
#'
#' @return HTML elements
#' @keywords internal
geneExprFileInput <- function(id, clearable=FALSE) {
    info <- paste(
        tags$ul(
            class="popover-list",
            tags$li("Tab-separated values (TSV)"),
            tags$li("Read counts of genes (rows) across sample (columns)"),
            tags$li("The first column must contain gene symbols and be",
                    "named", tags$kbd("Gene ID"))),
        tags$hr(), helpText("Example:"), tags$table(
            class="table table-condensed",
            tags$thead(
                tableRow("Gene ID", "SMP-18", "SMP-03", "SMP-54",
                         th=TRUE)),
            tags$tbody(
                tableRow("AMP1", "24", "10", "43"),
                tableRow("BRCA1", "38", "46", "32"),
                tableRow("BRCA2", "43", "65", "21"))))
    input <- fileBrowserInfoInput(id, "Gene expression", clearable=clearable,
                                  infoContent=info)
    return(input)
}

#' @rdname geneExprFileInput
#' @importFrom shiny tags
ASquantFileInput <- function(id, clearable=FALSE){
    info <- paste(
        tags$ul(
            class="popover-list",
            tags$li("Tab-separated values (TSV)"),
            tags$li("PSI values of alternative splicing events (rows)",
                    "across samples (columns)"),
            tags$li(
                "The first column must contain alternative splicing event",
                "identifiers and be named", tags$kbd("AS Event ID")),
            tags$li(
                "PSI values must be between 0 and 1 or between 0 and 100;",
                "if the latter, values are scaled between 0 and 1")),
        tags$hr(), helpText("Example:"), tags$table(
            class="table table-condensed",
            tags$thead(
                tableRow("AS Event ID", "SMP-18", "SMP-03", th=TRUE)),
            tags$tbody(
                tableRow("someASevent001", "0.71", "0.30"),
                tableRow("anotherASevent653", "0.63", "0.37"),
                tableRow("yeatAnother097", "0.38", "0.62"))))
    input <- fileBrowserInfoInput(id, "Alternative splicing quantification",
                                  clearable=clearable, infoContent=info)
    return(input)
}

#' @rdname geneExprFileInput
#' @importFrom shinyBS bsPopover
#' @importFrom shiny tags
junctionQuantFileInput <- function(id, clearable=FALSE) {
    info <- paste(
        tags$ul(
            class="popover-list",
            tags$li("Tab-separated values (TSV)"),
            tags$li("Read counts of exon-exon junctions (rows) across",
                    "samples (columns)"),
            tags$li("The first column must contain junction identifiers",
                    "and be named", tags$kbd("Junction ID")),
            tags$li("Only chromosome number and capital letters X, Y, Z, W",
                    "and M, followed by the genomic regions are supported;",
                    "acceptable junction identifiers include:",
                    tags$kbd("10_18748_21822"), ",",
                    tags$kbd("chromosome 10 (18748 to 21822)"), "and",
                    tags$kbd("chr10:18748-21822")),
            tags$li("Optionally, indicate the strand with", tags$kbd("+"), "or",
                    tags$kbd("-"),
                    "at the end of the junction identifier; e.g.",
                    tags$kbd("10:3213:9402:+"), "and",
                    tags$kbd("chr10:3213-9402 -")),
            tags$li("Rows whose junction identifiers contain",
                    tags$kbd("alt"), ",", tags$kbd("random"), "or",
                    tags$kbd("Un"), "in chromosome names are discarded")),
        tags$hr(), helpText("Example:"), tags$table(
            class="table table-condensed",
            tags$thead(tableRow("Junction ID", "SMP-18", "SMP-03", th=TRUE)),
            tags$tbody(tableRow("10:6752-7393", "4", "0"),
                       tableRow("10:18748-21822", "8", "46"),
                       tableRow("10:24257-25325", "83", "65"))))
    input <- fileBrowserInfoInput(id, "Exon-exon junction read counts",
                                  clearable=clearable, infoContent=info)
    return(input)
}

#' @rdname geneExprFileInput
#' @importFrom shinyBS bsPopover
#' @importFrom shiny tags
sampleInfoFileInput <- function(id, clearable=FALSE) {
    info <- paste(
        tags$ul(
            class="popover-list",
            tags$li("Tab-separated values (TSV)"),
            tags$li(
                "Sample identifiers (rows) and their attributes (columns)"),
            tags$li("The first column must contain sample identifiers",
                    "and be named", tags$kbd("Sample ID")),
            tags$li("Optionally, indicate the subject associated to",
                    "each sample in a column named",
                    tags$kbd("Subject ID"))),
        tags$hr(), helpText("Example:"), tags$table(
            class="table table-condensed",
            tags$thead(
                tableRow("Sample ID", "Type", "Tissue", "Subject ID",
                         th=TRUE)),
            tags$tbody(
                tableRow("SMP-01", "Tumour", "Lung", "SUBJ-03"),
                tableRow("SMP-02", "Normal", "Blood", "SUBJ-12"),
                tableRow("SMP-03", "Normal", "Blood", "SUBJ-25"))))
    input <- fileBrowserInfoInput(id, "Sample information",
                                  clearable=clearable, infoContent=info)
    return(input)
}

#' @rdname geneExprFileInput
#' @importFrom shinyBS bsPopover
#' @importFrom shiny tags
subjectInfoFileInput <- function(id, clearable=FALSE) {
    info <- paste(
        tags$ul(
            class="popover-list",
            tags$li("Tab-separated values (TSV)"),
            tags$li("Subject identifiers (rows) and their attributes",
                    "(columns)"),
            tags$li("The first column must contain subject identifiers and",
                    "be named", tags$kbd("Subject ID"))),
        tags$hr(),
        helpText("Example:"), tags$table(
            class="table table-condensed",
            tags$thead(tableRow("Subject ID", "Age", "Gender", "Race",
                               th=TRUE)),
            tags$tbody(tableRow("SUBJ-01", "34", "Female", "Black"),
                       tableRow("SUBJ-02", "22", "Male", "Black"),
                       tableRow("SUBJ-03", "58", "Female", "Asian"))))
    input <- fileBrowserInfoInput(id, "Subject information",
                                  clearable=clearable, infoContent=info)
    return(input)
}

# Other ------------------------------------------------------------------------

#' @rdname appUI
#' @importFrom shinyjs hidden
dataUI <- function(id, tab) {
    ns <- NS(id)
    uiList <- getUiFunctions(
        ns, "data", bsCollapsePanel,
        priority=paste0(c("localData", "firebrowse", "gtexData", "recountData",
                          "inclusionLevels", "inclusionLevelsFilter",
                          "geNormalisationFiltering"), "UI"))

    tcga <- tags$abbr(title="The Cancer Genome Atlas", "TCGA")
    gtex <- tags$abbr(title="Genotype-Tissue Expression project", "GTEx")
    sra  <- tags$abbr(title="Sequence Read Archive", "SRA")

    analysesDescription <- tagList(
        fluidRow(
            column(3, style="padding: 5px !important;",
                   h4("Dimensionality reduction"),
                   helpText(tags$ul(
                       class="list-unstyled",
                       tags$li("Principal Component Analysis (PCA)"),
                       tags$li("Independent Component Analysis (ICA)")))),
            column(3, style="padding: 5px !important;",
                   h4("Differential splicing and expression analysis"),
                   helpText("Based on variance and median statistical tests")),
            column(3, style="padding: 5px !important;",
                   h4("Survival analysis"),
                   helpText(tags$ul(
                       class="list-unstyled",
                       tags$li("Analyse survival based on clinical attributes",
                               "(e.g. tumour stage, gender and race)"),
                       tags$li("Study the impact of a single alternative",
                               "splicing event or gene on subject survival")))),
            column(3, style="padding: 5px !important;",
                   h4("Gene, transcript and protein information"),
                   helpText("Check available annotation for splicing events",
                            "and genes including related research articles"))))

    customDataTutorial <- paste0("https://nuno-agostinho.github.io/psichomics/",
                                 "articles/custom_data.html")
    welcome <- div(
        id=ns("welcome"),
        linkToArticles(),
        h1("Welcome to psichomics"),
        "Integrative analyses of alternative splicing and gene expression",
        "based on transcriptomic and sample-associated data from multiple",
        "sources, including:",
        tags$ul(
            tags$li(tags$a(href=customDataTutorial, target="_blank",
                           "User-provided data")),
            tags$li("The Cancer Genome Atlas (TCGA) via Firebrowse"),
            tags$li("Genotype-Tissue Expression (GTEx) project"),
            tags$li("Sequence Read Archive (SRA) via recount2")),
        tags$hr(),
        tags$ol(
            id="list",
            tags$li(HTML(paste0(
                "Load gene expression values, alternative splicing ",
                "junction quantification and/or sample-associated data ",
                "from ", tcga, ", ", gtex, ", ", sra, " or user-provided data."
            ))),
            tags$li("Import or quantify alternative splicing. Alternative",
                    "splicing is calculated using the percent spliced-in (PSI)",
                    "metric.",
                    tags$br(), tags$small(
                        style="color: gray;",
                        "Note: retained intron (RI) events are not",
                        "calculated in psichomics.")),
            tags$li("Explore statistically significant and specific genes",
                    "and alternative splicing events using:")),
        analysesDescription,
        div(style="text-align:right",
            tags$a(href="http://imm.medicina.ulisboa.pt/group/distrans/",
                   target="_blank", "Disease Transcriptomics Lab, iMM"),
            tags$br(),
            tags$a(href="mailto:nunodanielagostinho@gmail.com",
                   "Nuno Saraiva-Agostinho", icon("envelope")),
            tags$br(),
            sprintf("psichomics %s, 2015-2021", packageVersion("psichomics"))))

    tab(title="Data", icon="table",
        sidebarLayout(
            sidebar( do.call(bsCollapse, c(id=ns("accordion"), uiList)) ),
            mainPanel( welcome, uiOutput(ns("tablesOrAbout")) ) ))
}

#' Creates a \code{tabPanel} template for a \code{datatable} with a title and
#' description
#'
#' @param ns Namespace function
#' @param title Character: tab title
#' @param tableId Character: id of the \code{datatable}
#' @param columns Character: column names of the \code{datatable}
#' @param visCols Boolean: visible columns
#' @param data Data frame: dataset of interest
#' @param description Character: description of the table (optional)
#' @param icon Character: list containing an item named \code{symbol}
#' (FontAwesome icon name) and another one named \code{colour} (background
#' colour)
#'
#' @importFrom shinyBS bsTooltip bsCollapse bsCollapsePanel
#' @importFrom DT dataTableOutput
#' @importFrom shiny hr br tabPanel selectizeInput column fluidRow p mainPanel
#' downloadButton
#'
#' @return HTML elements
#' @keywords internal
tabDataset <- function(ns, title, tableId, columns, visCols, data,
                       description=NULL, icon=NULL) {
    tablename <- ns(paste("table", tableId, sep="-"))

    downloadId <- paste(tablename, "download", sep="-")
    download <- downloadButton(downloadId, "Save table",
                               class="pull-right btn-info")

    if(!is.null(description)) {
        description <- p(tags$strong("Table description:"), description)
        download <- fluidRow(column(10, description), column(2, download))
    }

    # Get class of each column
    colType <- sapply(seq(ncol(data)), function(i) class(data[[i]]))
    colType[colType == "character"] <- "string"

    # Show class of each column
    choices <- columns
    names(choices) <- sprintf("%s (%s class)", columns, colType)

    visColsId <- paste(tablename, "columns", sep="-")
    visibleColumns <- suppressWarnings(selectizeInput(
        visColsId, label="Visible columns",  choices=choices, selected=visCols,
        multiple=TRUE, width="auto",
        options=list(plugins=list('remove_button', 'drag_drop'), render=I(
            "{ item: function(item, escape) {
            return '<div>' + escape(item.value) + '</div>'; } }"))))

    # Add a common HTML container to allow for multiple Highcharts plots
    multiPlotId        <- paste(tablename, "multiPlot", sep="-")
    multiHighchartsPlots <- fluidRow(column(12, uiOutput(multiPlotId)))

    if (is.null(icon)) {
        name <- title
    } else {
        colour <- switch(icon$colour,
                         "green"="progress-bar-success",
                         "blue"="progress-bar-info",
                         "orange"="progress-bar-warning",
                         "red"="progress-bar-danger")
        name <- tags$div(
            tags$span(class=paste("badge", colour), icon(icon$symbol)), title)
    }

    tabPanel(title=name, value=title, br(), download, br(), bsCollapse(
        open="Summary",
        bsCollapsePanel(tagList(icon("table"), "Data table"),
                        value="Data table", visibleColumns, hr(),
                        dataTableOutput(tablename)),
        bsCollapsePanel(tagList(icon("chart-pie"), "Summary"), value="Summary",
                        multiHighchartsPlots)))
}

prepareDatatableSettingsTable <- function(filename, settings) {
    if (!is.null(filename)) {
        filename <- prepareWordBreak(filename)
        filename <- tags$small(tags$b("Loaded based on file:"),
                               tags$var(filename))
    }

    if (!is.null(settings)) {
        settingsDf <- data.frame(names(settings), sapply(
            sapply(settings, paste, collapse=", "), prepareWordBreak))
        colnames(settingsDf) <- c("Attribute", "Item")
        settings <- table2html(settingsDf, rownames=FALSE, thead=TRUE,
                               class="table table-condensed table-striped")
        settings <- tags$small(tagList(tags$b("Dataset settings"), settings))
        settings <- gsub("&lt;", "<", settings, fixed=TRUE)
        settings <- gsub("&gt;", ">", settings, fixed=TRUE)
        settings <- HTML(settings)
        return(settings)
    }

    extra <- NULL
    if ( !is.null(filename) || !is.null(settings) ) {
        extra <- tagList(
            tags$hr(), filename,
            if (!is.null(filename) && !is.null(settings))
                tagList(tags$br(), tags$br()),
            settings)
    }
}

prepareDatatablePlots <- function(table, output, ns) {
    plots <- NULL
    if (is.null(attr(table, "plots"))) {
        isGeneExpr <- !is.null(attr(table, "dataType")) &&
            attr(table, "dataType") == "Gene expression"
        isPSI <- !is.null(attr(table, "dataType")) &&
            attr(table, "dataType") == "Inclusion levels"
        if (isGeneExpr) {
            if (is(table, "EList")) table <- table$E
            geneExprPerSamplePlot <- plotGeneExprPerSample(
                table, sortByMedian=TRUE,
                title="Gene expression distribution per sample")
            librarySizePlot <- plotLibrarySize(table)
            plots <- list(
                highchart=geneExprPerSamplePlot,
                highchart=librarySizePlot)
        } else if (isPSI) {
            cache     <- isolate(getInclusionLevelsSummaryStatsCache())
            medianVar <- plotRowStats(table, x="median", y="var",
                                      xlim=c(0,1 ), cache=cache, verbose=TRUE) +
                labs(x="Median PSI", y="PSI variance") +
                ggtitle(paste("Scatterplot of alternative splicing",
                              "quantification per event")) +
                theme_light(14)
            cache <- attr(medianVar, "cache")

            rangeVar  <- plotRowStats(table, x="range", y="log10(var)",
                                      xlim=c(0, 1), cache=cache, verbose=TRUE) +
                labs(x="PSI range", y="log10(PSI variance)") +
                ggtitle(paste("Scatterplot of alternative splicing",
                              "quantification per event")) +
                theme_light(14)
            cache <- attr(rangeVar, "cache")

            setInclusionLevelsSummaryStatsCache(cache)
            plots <- list(plot=medianVar, plot=rangeVar)
        }
        attr(table, "plots") <- plots
    }
    tablename <- attr(table, "tablenameID")
    plots     <- attr(table, "plots")

    renderedPlots <- lapply(seq(plots), function(i) {
        type <- names(plots)[[i]]
        FUN <- switch(type, highchart=renderHighchart, plot=renderPlot)
        res <- FUN(plots[[i]])
        attr(res, "type") <- type
        return(res)
    })

    plotList <- tagList(NULL)
    for (each in seq(renderedPlots)) {
        plot <- renderedPlots[[each]]
        type <- attr(plot, "type")
        id   <- paste0(gsub(" ", "_", tablename), "-", type, each)
        output[[id]] <- plot

        FUN  <- switch(type, highchart=highchartOutput, plot=plotOutput)
        item <- tagList(FUN(ns(id)))
        plotList <- tagAppendChild(plotList, item)
    }
    return(plotList)
}

#' Render a specific data tab (including data table and related interface)
#'
#' @param index Integer: index of the data to load
#' @param data Data frame: data with everything to load
#' @param name Character: name of the dataset
#' @param session Shiny session
#' @param input Shiny session input
#' @param output Shiny session output
#'
#' @importFrom shiny tags HTML
#' @importFrom DT renderDataTable
#' @importFrom shiny downloadHandler br
#' @importFrom utils write.table
#' @importFrom shinyjs show hide
#' @importFrom ggplot2 labs ggtitle theme_light
#'
#' @inherit psichomics return
#' @keywords internal
createDataTab <- function(index, data, name, session, input, output) {
    ns <- session$ns
    tablename <- paste("table", name, index, sep="-")

    table <- data[[index]]
    # Only show default columns if they are defined (don't cause problems)
    if (is(table, "EList")) table <- table$E
    subsetToShow <- table

    visCols <- input[[paste(tablename, "columns", sep="-")]]
    if (!is.null(visCols)) {
        match <- visCols %in% colnames(table)
        subsetToShow <- subset(table, select=visCols[match])
    }

    output[[tablename]] <- renderDataTable(
        subsetToShow, style="bootstrap", selection='none', filter="top",
        options=list(pageLength=10))

    downloadId <- paste(tablename, "download", sep="-")
    output[[downloadId]] <- downloadHandler(
        filename = paste(name, attr(table, "tablename")),
        content = function(file) {
            res <- cbind(rownames(table), table)
            names(res)[1] <- attr(table, "dataType")
            write.table(res, file, quote=FALSE, row.names=FALSE, sep="\t")
        }
    )

    multiPlotId        <- paste(tablename, "multiPlot", sep="-")
    createInfoInterface <- function(output, table) {
        rows     <- attr(table, "rows")
        rows     <- ifelse(!is.null(rows), rows, "rows")
        cols     <- attr(table, "columns")
        cols     <- ifelse(!is.null(cols), cols, "columns")
        settings <- prepareDatatableSettingsTable(attr(table, "filename"),
                                                  attr(table, "settings"))
        plotList <- prepareDatatablePlots(table, output, ns)

        tags$div(tags$h4(paste(ncol(table), cols)),
                 tags$h4(paste(nrow(table), rows)),
                 plotList,
                 settings)
    }

    attr(table, "tablenameID") <- tablename
    output[[multiPlotId]] <- renderUI(createInfoInterface(output, table))
}

prepareSubjectSampleMatch <- reactive({
    samples  <- getSampleId()
    subjects <- getSubjectId()
    match <- getSubjectFromSample(samples, subjects, sampleInfo=getSampleInfo())
    setClinicalMatchFrom("Inclusion levels", match)
})

#' @rdname appServer
#'
#' @importFrom shiny selectInput tabsetPanel tags h1 h2 HTML fluidRow column
#' tagList
#' @importFrom shinyjs show hide
dataServer <- function(input, output, session) {
    ns <- session$ns

    # Show welcome screen when there's no data loaded
    output$tablesOrAbout <- renderUI({
        if(is.null(getData())) {
            show("welcome", anim=TRUE, animType="fade")
        } else {
            hide("welcome", anim=TRUE, animType="fade")
            uiOutput(ns("datatabs"))
        }
    })

    # Render tables when data changes
    observe({
        data <- getData()
        if (!is.null(data)) {
            for (category in names(data)) {
                categoryData <- data[[category]]
                # Create data tab for each dataset in a data category
                lapply(seq_along(categoryData), createDataTab,
                       data=categoryData, category, session, input, output)
            }
        }
    })

    # Render tabs with data tables
    output$datatabs <- renderUI({
        categoryData <- getCategoryData()
        category <- getCategory()

        dataTablesUI <- lapply(
            seq_along(categoryData), function(i) {
                data <- categoryData[[i]]
                if (is(data, "EList")) data <- data$E

                # Display at most 100 columns if no visible columns are set
                visCols <- attr(data, "show")
                if (is.null(visCols) && ncol(data) > 100)
                    visCols <- colnames(data)[seq(100)]

                name <- names(categoryData)[i]
                if (grepl("Gene expression", name)) {
                    attr(data, "icon") <- list(symbol="dna", colour="orange")
                } else if (grepl("Junction quantification", name)) {
                    attr(data, "icon") <- list(symbol="dna", colour="orange")
                } else if (grepl("Sample metadata", name)) {
                    attr(data, "icon") <- list(symbol="vial", colour="blue")
                } else if (grepl("Clinical data", name)) {
                    attr(data, "icon") <- list(symbol="vial", colour="blue")
                }
                tabDataset(
                    ns, name, icon=attr(data, "icon"),
                    paste(category, i, sep="-"), colnames(data), visCols, data,
                    description=attr(data, "description"))
            }
        )
        do.call(tabsetPanel, c(id=ns("datasetTab"), dataTablesUI))
    })

    # Change the active dataset
    observe( setActiveDataset(input$datasetTab) )

    # Match clinical data with sample information
    observe({
        if ( !is.null(getSubjectId()) && !is.null(getSampleId()) ) {
            prepareSubjectSampleMatch()
        }
    })

    observe(setSelectedDataPanel(input$accordion))

    # Run server logic from the scripts
    getServerFunctions("data", priority=paste0(
        c("localData", "firebrowse", "gtexData",
          "inclusionLevels", "inclusionLevelsFilter",
          "geNormalisationFiltering"), "Server"))
}

attr(dataUI, "loader") <- "app"
attr(dataServer, "loader") <- "app"
