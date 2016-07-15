diffAnalysisTableUI <- function(id) {
    ns <- NS(id)
    tagList(
        uiOutput(ns("modal")),
        sidebarLayout(
            sidebarPanel(
                uiOutput(ns("groupsInfo")),
                checkboxGroupInput(
                    ns("statsChoices"),
                    "Choose statistical analyses to perform:",
                    # Basic stats is on and disabled by JavaScript
                    c("Variance and median"="basicStats",
                      "Wilcoxon Test (1 or 2 groups)"="wilcox",
                      "Kruskal-Wallis Rank Sum Test (2 or more groups)"="kruskal", 
                      "Levene's test (2 or more groups)"="levene"),
                    selected=c("basicStats", "wilcox", "kruskal", "levene")),
                # Disable checkbox of basic statistics
                HTML("<script>",
                     '$("[value=basicStats]").attr("disabled", true);', 
                     "</script>"),
                actionButton(ns("startAnalyses"), class="btn-primary", 
                             "Perform analyses")
            ), mainPanel(
                uiOutput(ns("showColumns")),
                dataTableOutput(ns("statsTable")),
                highchartOutput(ns("densitySparklines"))
            )
        )
    )
}

#' @importFrom highcharter highchart hc_chart hc_xAxis hc_plotOptions hc_tooltip
#' JS hc_add_series_scatter
densitySparkline <- function(hc, psi, type, bandwidth) {
    for (group in unique(type)) {
        row  <- psi[type == group]
        # Calculate the density of inclusion levels for each sample type with a
        # greatly reduced number of points for faster execution
        den <- density(row, n=10, bw=bandwidth, na.rm=TRUE)
        hc <- hc %>%
            hc_add_series_density(den, name=group, area=TRUE) %>%
            hc_new_sparkline()
    }
    return(hc)
}

#' Perform statistical analysis on a vector with elements from different groups
#' 
#' @param vector Numeric
#' @param group Character: group of each element in the vector
#' @param hc Highchart object used to get the density plots
#' @param threshold Integer: minimum number of data points to perform analysis
#' in a group (default is 1)
#' @param analyses Character: name of the analyses to perform (all by default)
#' 
#' @return A data frame row with the results
statsAnalysis <- function(vector, group, threshold=1, step=100,
                          analyses=c("wilcox", "kruskal", "levene")) {
    # Filter vector by a given threshold
    filterByThreshold <- function(thisType, allTypes, vector, threshold) {
        vector <- vector[thisType == allTypes]
        if ( sum(!is.na(vector)) >= threshold )
            return(vector)
    }
    names(vector) <- group
    vector <- lapply(unique(group), filterByThreshold, group, vector, 
                     threshold)
    
    vector  <- unlist(vector)
    group <- names(vector)
    vector  <- as.numeric(vector)
    len  <- length(unique(group))
    
    # Wilcoxon tests
    wilcox <- NULL
    if ("wilcox" %in% analyses) {
        if (len == 2) {
            typeOne <- group == unique(group)[1]
            wilcox  <- suppressWarnings(wilcox.test(vector[typeOne], 
                                                    vector[!typeOne]))
        } else if (len == 1) {
            wilcox <- suppressWarnings(wilcox.test(vector))
        }
    }
    
    # Kruskal-Wallis test
    kruskal <- NULL
    if ("kruskal" %in% analyses && len >= 2) {
        kruskal <- tryCatch(kruskal.test(vector, factor(group)),
                            error=return)
        if ("error" %in% class(kruskal)) kruskal <- NULL
    }
    
    # Levene's test
    levene <- NULL
    if ("levene" %in% analyses && len >= 2) {
        nas <- is.na(vector)
        levene <- tryCatch(levene.test(vector[!nas], factor(group[!nas])),
                           error=return)
        if ("error" %in% class(levene)) levene <- NULL
    }
    
    # Density sparklines
    hc <- getDensitySparklines()
    hc <- densitySparkline(hc, vector, group, 0.1)
    setDensitySparklines(hc)
    
    # Variance and median
    group <- split(vector, group)
    samples <- lapply(group, function(i) sum(!is.na(i))) # Number of samples
    med <- lapply(group, median, na.rm=TRUE) # Median
    var <- lapply(group, var, na.rm=TRUE) # Variance
    
    vector <- c(Samples=samples, Wilcox=wilcox, Kruskal=kruskal, 
                Levene=levene, Variance=var, Median=med)
    vector <- vector[!vapply(vector, is.null, logical(1))] # Remove NULL
    vector <- data.frame(vector, stringsAsFactors=FALSE)
    return(vector)
}

#' @importFrom lawstat levene.test
#' @importFrom stats kruskal.test median wilcox.test var
#' @importFrom DT renderDataTable
diffAnalysisTableServer <- function(input, output, session) {
    ns <- session$ns
    
    # Information on the data groups from TCGA
    output$groupsInfo <- renderUI({
        # Get event's inclusion levels
        psi <- getInclusionLevels()
        
        if (is.null(psi)) return(tagList(
            helpText(icon("exclamation-circle"), 
                     "No alternative splicing junction quantification loaded.")))
        
        # Separate samples by their type
        ids <- names(psi)
        type <- getSampleTypes(ids)
        
        bullet  <- "\u2022"
        groups <- NULL
        for (each in unique(type))
            groups <- tagList(groups, br(), bullet, each)
        
        return(tagList(
            helpText("The data contains the following groups:", groups),
            hr()))
    })
    
    
    observeEvent(input$startAnalyses, {
        isolate({
            # Get event's inclusion levels
            psi <- getInclusionLevels()
            statsChoices <- input$statsChoices
        })
        if (is.null(psi)) {
            errorModal(session, "No AS event quantification",
                       "Insert or quantify alternative splicing events by",
                       "going to the", icon("table"), tags$b("Data"),
                       "tab and opening",  icon("calculator"),
                       tags$b("Quantify alternative splicing events"), ".")
            return(NULL)
        }
        
        # Separate samples by their type
        ids <- names(psi)
        type <- getSampleTypes(ids)
        
        # cl <- parallel::makeCluster(getOption("cl.cores", getCores()))
        step <- 100 # Avoid updating after analysing each event
        startProgress("Performing statistical analysis", 
                      divisions=4+round(nrow(psi)/step))
        time <- Sys.time()
        
        hc <- highchart() %>%
            hc_xAxis(min=0, max=1) %>%
            hc_tooltip(
                headerFormat="<small>Inclusion levels: {point.x}</small></br>",
                pointFormat=paste(
                    span(style="color:{point.color}", "\u25CF "),
                    tags$b("{series.name}"), br()))
        setDensitySparklines(hc)
        
        count <- 0
        stats <- apply(psi, 1, function(...) {
            count <<- count + 1
            if (count %% step == 0)
                updateProgress("Performing statistical analysis", console=FALSE)
            return(statsAnalysis(...))
        }, factor(type), threshold=1, step=step, analyses=statsChoices)
        
        # Convert to data frame
        df <- do.call(rbind.fill, stats)
        df <- df[, !grepl("method|data.name", colnames(df))]
        updateProgress("Performing statistical analysis", console=FALSE)
        
        # Calculate delta variance and delta median if there are only 2 groups
        deltaVar <- df[, grepl("Variance", colnames(df)), drop=FALSE]
        if (ncol(deltaVar) == 2) {
            deltaVar <- deltaVar[, 2] - deltaVar[, 1]
            deltaMed <- df[, grepl("Median", colnames(df))]
            deltaMed <- deltaMed[, 2] - deltaMed[, 1]
            df <- cbind(df, deltaVar, deltaMed)
        }
        updateProgress("Performing statistical analysis", console=FALSE)
        
        hc <- getDensitySparklines()
        assign("hc2", hc, .GlobalEnv)
        sparklines <- hchart(hc)
        stats <- cbind(Density=sparklines, df)
        rownames(stats) <- rownames(psi)
        setDifferentialAnalyses(stats)
        updateProgress("Performing statistical analysis", console=FALSE)
        
        # parallel::stopCluster(cl)
        print(Sys.time() - time)
        closeProgress()
    })
    
    # Show table and respective interface when statistical table is available
    observe({
        stats <- getDifferentialAnalyses()
        if (is.null(stats)) return(NULL)
        
        # Columns to show in statistical table
        output$showColumns <- renderUI({
            tagList(
                downloadButton(ns("download"), "Download whole table"),
                selectizeInput(ns("columns"), "Show columns", multiple=TRUE,
                               choices=colnames(stats), width="auto",
                               selected=colnames(stats),
                               options=list(plugins=list('remove_button', 
                                                         'drag_drop'))),
                hr())
        })
        
        # Render statistical table with the selected columns
        output$statsTable <- renderDataTableSparklines({
            if (!is.null(input$columns) && all(input$columns %in% names(stats)))
                stats[, input$columns]
        }, style="bootstrap", selection="none", filter='top',
        options=list(pageLength=10))
        
        # Prepare table to be downloaded
        output$download <- downloadHandler(
            filename=paste(getCategories(), "Differential splicing analyses"),
            content=function(file)
                write.table(getDifferentialAnalyses()[-c("Density")], 
                            file, quote=FALSE, row.names=FALSE, sep="\t")
        )
    })
}

attr(diffAnalysisTableUI, "loader") <- "analysis"
attr(diffAnalysisTableUI, "name") <- "Differential analysis (exploratory)"
attr(diffAnalysisTableServer, "loader") <- "analysis"