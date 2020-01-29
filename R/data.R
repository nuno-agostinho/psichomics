## TODO(NunoA): should default columns be a perfect match or just a partial
## match? A partial match would be better for certain situations

## TODO(NunoA): render UI for each data table instead of rendering UI for all
## so there's no refresh

#' Set attributes to an object
#' 
#' @param object Object
#' @param ... Named parameters to convert to attributes
#' 
#' @return Object with attributes set
#' @keywords internal
#' 
#' @examples 
#' ll <- list(a="hey", b="there")
#' psichomics:::addObjectAttrs(ll, "words"=2, "language"="English")
addObjectAttrs <- function (object, ...) {
    args <- as.list(match.call())[-c(1, 2)]
    for (k in seq(args)) attr(object, names(args[k])) <- args[[k]]
    return(object)
}

#' Plot sample statistics per row
#'
#' @param data Data frame or matrix
#' @param x,y Character: statistic to calculate and display in the plot per row;
#' choose between \code{mean}, \code{median}, \code{var} or \code{range}
#' (or transformations of those variables, e.g. \code{log10(var)})
#' @param subset Boolean or integer: \code{data} points to highlight (if
#' \code{NULL}, all points are highlighted)
#' @param xmin,xmax,ymin,ymax Numeric: minimum and maximum X and Y values to 
#' draw in the plot
#' @param xlim,ylim Numeric: X and Y axis range
#'
#' @importFrom ggplot2 geom_vline geom_hline xlim ylim ggtitle
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
plotRowStats <- function(data, x, y, subset=NULL, xmin=NULL, xmax=NULL, 
                         ymin=NULL, ymax=NULL, xlim=NULL, ylim=NULL) {
    stats <- c("range", "var", "median", "mean")
    if (!any(sapply(stats, grepl, x)) || !any(sapply(stats, grepl, y))) {
        stop("Arguments 'x' and 'y' must contain one of the strings:",
             paste(stats, collapse=", "))
    }
    
    calculateXandYvalues <- function(psi, stats) {
        names(stats) <- stats
        input <- lapply(stats, grepl, c(x, y))
        
        rowRanges <- function(mat, ...) {
            apply(mat, 1, max, ...) - apply(mat, 1, min, ...)
            # apply(mat, 1, function(k) max(k, ...) - min(k, ...))
        }
        
        x <- y <- NULL
        vars <- list()
        for (stat in stats) {
            if (any(input[[stat]])) {
                message(sprintf("Calculating %s per splicing event...", stat))
                FUN <- switch(stat,
                              "var"=rowVars,
                              "mean"=rowMeans,
                              "median"=rowMedians,
                              "range"=rowRanges)
                res <- FUN(psi, na.rm=TRUE)
                vars[[stat]] <- res
            }
        }
        vars <- data.frame(vars, stringsAsFactors=FALSE)
        return(vars)
    }
    vars <- calculateXandYvalues(data, stats)
    
    message("Preparing plot...")
    if (!is.null(subset)) {
        varsSubset <- vars[subset]
        plot <- ggplot(vars, aes_string(x, y)) +
            geom_point(size=1, na.rm=TRUE, alpha=0.5, colour="grey") +
            geom_point(data=varsSubset, size=1, na.rm=TRUE, alpha=0.5, 
                       colour="grey") +
            geom_density_2d(data=varsSubset, colour="orange", na.rm=TRUE) +
            labs(x=x, y=y)
    } else {
        plot <- ggplot(vars, aes_string(x, y)) +
            geom_point(size=1, na.rm=TRUE, alpha=0.5) +
            geom_point(size=1, na.rm=TRUE, alpha=0.5) +
            geom_density_2d(colour="orange", na.rm=TRUE) +
            labs(x=x, y=y)
    }
    
    if (!is.null(xlim)) plot <- plot + xlim(xlim)
    if (!is.null(ylim)) plot <- plot + ylim(ylim)
    
    # Intercept lines
    if (!is.null(xmin)) plot <- plot + geom_vline(xintercept=xmin, colour="red")
    if (!is.null(xmax)) plot <- plot + geom_vline(xintercept=xmax, colour="red")
    if (!is.null(ymin)) plot <- plot + geom_hline(yintercept=ymin, colour="red")
    if (!is.null(ymax)) plot <- plot + geom_hline(yintercept=ymax, colour="red")
    return(plot)
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

#' File input for gene expression
#' 
#' @param geneExprFileId Character: identifier for gene expression input
#' 
#' @return HTML elements
#' @keywords internal
geneExprFileInput <- function(geneExprFileId) {
    fileBrowserInput(
        geneExprFileId,
        "File with gene expression",
        placeholder="No file selected",
        info=TRUE, infoFUN=bsPopover, infoTitle=paste(
            "File containing the read counts of each gene (rows) per",
            "sample (columns)."),
        infoContent=paste(
            "The first column must contain gene symbols and be named", 
            tags$kbd("Gene ID"), tags$hr(), helpText("Example:"), tags$table(
                class="table table-condensed",
                tags$thead(
                    tableRow("Gene ID", "SMP-18", "SMP-03", "SMP-54", 
                             th=TRUE)),
                tags$tbody(
                    tableRow("AMP1", "24", "10", "43"),
                    tableRow("BRCA1", "38", "46", "32"),
                    tableRow("BRCA2", "43", "65", "21")))))
}

#' File input for alternative splicing quantification
#' 
#' @param ASquantFileId Character: identifier for alternative splicing 
#' quantification input
#' @param speciesId Character: identifier for species selection input
#' @param assemblyId Character: identifier for genome assembly selection input
#' 
#' @return HTML elements
#' @keywords internal
ASquantFileInput <- function(ASquantFileId, speciesId, assemblyId){
    tagList(
        fileBrowserInput(
            ASquantFileId, "File with alternative splicing quantification",
            placeholder="No file selected",
            info=TRUE, infoFUN=bsPopover, infoTitle=paste(
                "File containing the PSI value of each alternative splicing",
                "event (rows) per sample (columns)."),
            infoContent=paste(
                tags$ul(
                    class="popover-list",
                    tags$li(
                        "The first column must contain alternative splicing",
                        "event identifiers and should be named",
                        tags$kbd("AS Event ID")),
                    tags$li(
                        "An alternative splicing event must be represented by:",
                        tags$kbd(
                            paste0("EventType_Chromosome_Strand_Coordinate1_",
                                   "Coordinate2_..._Gene"))),
                    tags$li(
                        "PSI values may be handed between 0 and 1 or between 0",
                        "and 100. If the later, PSI values are then scaled",
                        "betwen 0 and 1.")))),
        selectizeInput(speciesId, "Species", choices="Human", width = "100%",
                       options=list(create=TRUE)),
        selectizeInput(assemblyId, "Assembly", choices=c("hg19", "hg38"),
                       width = "100%", options=list(create=TRUE)))
}

#' @rdname appUI
#' @importFrom shinyjs hidden
dataUI <- function(id, tab) {
    ns <- NS(id)
    uiList <- getUiFunctions(
        ns, "data", bsCollapsePanel,
        priority=paste0(c("localData", "firebrowse", "gtexData", "recountData",
                          "inclusionLevels", "geNormalisationFiltering"), "UI"))
    
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
                            "and genes including related research articles."))))
    
    welcome <- div(
        id=ns("welcome"),
        linkToArticles(),
        h1("Welcome to psichomics"), HTML(paste0(
            "Perform integrative analyses of alternative splicing and gene ",
            "expression based on transcriptomic and sample-associated data ",
            "from The Cancer Genome Atlas (", tcga, "), the Genotype-Tissue ",
            "Expression (", gtex, ") project, Sequence Read Archive (", sra, 
            ") or user-provided data.")),
        tags$br(), tags$br(), tags$ol(
            id="list",
            tags$li(HTML(paste0(
                "Load gene expression values, alternative splicing ",
                "junction quantification and/or sample-associated data ",
                "from ", tcga, ", ", gtex, ", ", sra, " or user-provided data.",
                tags$br(), tags$small(
                    style="color: gray;",
                    "More data types will soon be supported.")))),
            tags$li("Quantify alternative splicing events based on the",
                    "values from the percent spliced-in (PSI) metric.",
                    # "The following event types are available:",
                    # "skipped exon (SE), mutually exclusive exon",
                    # "(MXE), alternative 3' and 5' splice site (A3SS",
                    # "and A5SS) and alternative first and last exon",
                    # "(AFE and ALE).", tags$br(),
                    tags$br(), tags$small(
                        style="color: gray;",
                        "Note: as", tcga, "does not include exon-intron",
                        "junction quantification, intron retention",
                        "events are not measurable.")),
            tags$li("Explore statistically significant genes/events or",
                    "individual genes/events of interest using:")), 
        analysesDescription, br(), br(),
        p(style="text-align:right",
          tags$a(href="http://imm.medicina.ulisboa.pt/group/distrans/",
                 target="_blank", "Disease Transcriptomics Lab, iMM"), 
          "(", tags$a(href="mailto:nunodanielagostinho@gmail.com", 
                      "Nuno Saraiva-Agostinho", icon("envelope-o")),
          ", 2015-2020)", 
          br(), "Special thanks to my lab colleagues for their work-related",
          br(), "support and supporting chatter."))
    
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
    visibleColumns <- selectizeInput(
        visColsId, label="Visible columns",  choices=choices, selected=visCols, 
        multiple=TRUE, width="auto", 
        options=list(plugins=list('remove_button', 'drag_drop'), render=I(
            "{ item: function(item, escape) {
            return '<div>' + escape(item.value) + '</div>'; } }")))
    
    # Add a common HTML container to allow for multiple Highcharts plots
    multiPlotId        <- paste(tablename, "multiPlot", sep="-")
    loadingMultiPlotId <- paste(tablename, "loadingMultiPlot", sep="-")
    multiHighchartsPlots <- fluidRow(column(12, uiOutput(multiPlotId)))
    # div(id=loadingMultiPlotId, class="progress",
    #     div(class="progress-bar progress-bar-striped active",
    #         role="progressbar", style="width:100%",
    #         "Loading summary plots")))
    
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
        bsCollapsePanel(tagList(icon("pie-chart"), "Summary"), value="Summary",
                        multiHighchartsPlots)))
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
        })
    
    multiPlotId        <- paste(tablename, "multiPlot", sep="-")
    loadingMultiPlotId <- paste(tablename, "loadingMultiPlot", sep="-")
    
    createInfoInterface <- function(output, table) {
        rows <- attr(table, "rows")
        rows <- ifelse(!is.null(rows), rows, "rows")
        cols <- attr(table, "columns")
        cols <- ifelse(!is.null(cols), cols, "columns")
        
        filename <- attr(table, "filename")
        if (!is.null(filename)) {
            filename <- prepareWordBreak(filename)
            filename <- tags$small(tags$b("Loaded based on file:"),
                                   tags$var(filename))
        }
        
        settings <- attr(table, "settings")
        if (!is.null(settings)) {
            settingsDf <- data.frame(names(settings), sapply(
                sapply(settings, paste, collapse=", "), prepareWordBreak))
            colnames(settingsDf) <- c("Attribute", "Item")
            settings <- table2html(
                settingsDf, rownames=FALSE, thead=TRUE, 
                class="table table-condensed table-striped")
            settings <- tags$small(tagList(tags$b("Dataset settings"), 
                                           settings))
            settings <- gsub("&lt;", "<", settings, fixed=TRUE)
            settings <- gsub("&gt;", ">", settings, fixed=TRUE)
            settings <- HTML(settings)
        }
        
        extra <- NULL
        if ( !is.null(filename) || !is.null(settings) ) {
            extra <- tagList(
                tags$hr(), filename, 
                if (!is.null(filename) && !is.null(settings)) 
                    tagList(tags$br(), tags$br()), 
                settings)
        }
        
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
                
                librarySizePlot <- suppressWarnings(
                    plotDistribution(log10(colSums(table)),
                                     rugLabels=TRUE, vLine=FALSE) %>%
                        hc_xAxis(title=list(text="log10(Library sizes)")) %>%
                        hc_yAxis(title=list(text="Density")) %>%
                        hc_legend(enabled=FALSE) %>%
                        hc_title(
                            text="Library size distribution across samples") %>%
                        hc_subtitle(text=paste("Library size: number total",
                                               "mapped reads")))
                librarySizePlot$x$hc_opts$series[[1]]$color <- NULL
                librarySizePlot$x$hc_opts$series[[2]]$marker$fillColor <- NULL
                
                plots <- list(
                    highchart=geneExprPerSamplePlot,
                    highchart=librarySizePlot)
            } else if (isPSI) {
                medianVar <- plotRowStats(table, x="median", y="var", 
                                          xlim=c(0,1 )) +
                    labs(x="PSI median", y="PSI variance") +
                    ggtitle(paste("Scatterplot of alternative splicing",
                                  "quantification per event")) +
                    theme_light(14)
                rangeVar  <- plotRowStats(table, x="range", y="log10(var)", 
                                          xlim=c(0, 1)) +
                    labs(x="PSI range", y="log10(PSI variance)") +
                    ggtitle(paste("Scatterplot of alternative splicing",
                                  "quantification per event")) +
                    theme_light(14)
                plots <- list(plot=medianVar, plot=rangeVar)
            }
            attr(table, "plots") <- plots
        }
        tablename <- attr(table, "tablenameID")
        plots     <- attr(table, "plots")
        
        renderedPlots <- lapply(seq(plots), function(i) {
            FUN <- switch(names(plots)[[i]],
                          highchart=renderHighchart, plot=renderPlot)
            FUN(plots[[i]])
        })
        
        tags$div(tags$h4(paste(ncol(table), cols)),
                 tags$h4(paste(nrow(table), rows)),
                 # do.call(tagList, plotOutputs),
                 renderedPlots,
                 extra)
    }
    
    attr(table, "tablenameID") <- tablename
    output[[multiPlotId]] <- renderUI(createInfoInterface(output, table))
}

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
                
                tabDataset(
                    ns, names(categoryData)[i], icon=attr(data, "icon"),
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
        subjects   <- getSubjectId()
        samples    <- getSampleId()
        if ( !is.null(subjects) && !is.null(samples) ) {
            startProgress("Matching subjects to their samples...", 1)
            match <- getSubjectFromSample(samples, subjects, 
                                          sampleInfo=getSampleInfo())
            setClinicalMatchFrom("Inclusion levels", match)
            closeProgress("Matching process concluded")
        }
    })
    
    # Run server logic from the scripts
    getServerFunctions("data", priority=paste0(
        c("localData", "firebrowse", "gtexData",
          "inclusionLevels", "geNormalisationFiltering"), "Server"))
}

attr(dataUI, "loader") <- "app"
attr(dataServer, "loader") <- "app"